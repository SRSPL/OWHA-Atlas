#' Quantify cell type positions along tissue axes for integrated ST objects
#' 
#' @description
#' Computes cell type positions along specified tissue axes (epidermis-dermis depth or tissue length)
#' for Seurat spatial transcriptomics objects. Supports both single and integrated multi-region datasets.
#' 
#' @param object Seurat object (can be integrated with multiple images)
#' @param axis Which tissue axis to quantify: "epidermis_dermis" (depth) or "length"
#' @param annotation_col Name of metadata column containing cell type annotations
#' @param celltypes_to_plot Optional vector of specific cell types to include in plot
#' @param image_names Named vector mapping region labels to image names (e.g., c("Region1" = "uw_Region1"))
#'                    If NULL, assumes object has "Region" metadata and single image per region
#' @param dataset Name of the dataset being analyzed. Used to load image scale factors for micron conversion
#' @param scalefactors_path Optional custom path to scalefactors JSON file. If NULL and dataset is provided,
#'                          will look for file at: <dataset>/<dataset>_Realignment/binned_outputs/square_008um/spatial/scalefactors_json.json
#'                          relative to the current working directory
#' @param combine_regions Logical, whether to combine all regions in a single plot (default: FALSE)
#' @param cols Named vector of colors for cell types
#' @param celltype_levels Optional vector specifying order of cell types
#' @param timepoint_label Optional label for timepoint (for metadata)
#' @param align_regions Logical, whether to align axes across regions using shared cell types (default: TRUE)
#' @param invert_axis Logical, whether to invert the selected axis (default: FALSE)
#' @param return_plot Logical, whether to generate and display plot (default: TRUE)
#' @param projection_data Optional pre-calculated projection data (for compatibility)
#' @param include_microns Logical, whether to calculate micron-based distances (requires dataset parameter, default: FALSE)
#' @param wound_center_microns Numeric value indicating the position of the wound center in microns (for length axis plots)
#' 
#' @return Data frame with the following columns:
#' \itemize{
#'   \item barcode: Cell barcode
#'   \item projection: Normalized position along selected axis (0-1)
#'   \item Region: Region identifier
#'   \item celltype: Cell type annotation
#'   \item weight: Cell weight (1 for singlets, 0.5 for doublets)
#'   \item metaclusters: Metacluster assignment (if available)
#'   \item timepoint: Timepoint label (if provided)
#'   \item weighted_id: Row identifier
#'   \item length_position, depth_position: Normalized positions (if include_microns = TRUE)
#'   \item length_microns, depth_microns: Micron-based distances (if include_microns = TRUE)
#' }
#' 
#' @details
#' This function performs PCA on spatial coordinates to identify primary tissue axes.
#' PC1 represents tissue length and PC2 represents epidermis-dermis depth.
#' When multiple regions are present, the function can optionally align axes using
#' shared cell type positions to ensure consistent orientation across regions.
#' 
#' @examples
#' \dontrun{
#' # Basic usage with single region
#' result <- quantify_celltype_axis_integrated(
#'   object = seurat_obj,
#'   axis = "epidermis_dermis",
#'   annotation_col = "cell_type"
#' )
#' 
#' # Multi-region with micron conversion
#' result <- quantify_celltype_axis_integrated(
#'   object = integrated_obj,
#'   axis = "length",
#'   annotation_col = "first_type",
#'   image_names = c("D4_R1" = "d4_Region1", "D4_R2" = "d4_Region2"),
#'   dataset = "YWV04_D4PW",
#'   include_microns = TRUE,
#'   wound_center_microns = 3500
#' )
#' }
#' 
#' @importFrom Seurat GetTissueCoordinates subset
#' @importFrom dplyr filter mutate group_by summarise bind_rows row_number case_when
#' @importFrom ggplot2 ggplot aes geom_hline geom_text theme_minimal labs facet_wrap scale_color_manual
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom stats prcomp cor
#' @importFrom jsonlite fromJSON
#' 
#' @export
quantify_celltype_axis_integrated <- function(
    object, 
    axis = c("epidermis_dermis", "length"),
    annotation_col = "first_type",
    celltypes_to_plot = NULL, 
    image_names = NULL,
    dataset = NULL,
    scalefactors_path = NULL,
    combine_regions = FALSE, 
    cols = NULL, 
    celltype_levels = NULL, 
    timepoint_label = NULL,
    align_regions = TRUE,
    invert_axis = FALSE,
    return_plot = TRUE,
    projection_data = NULL,
    include_microns = FALSE,
    wound_center_microns = NULL
) { 
  
  # Load required packages
  required_packages <- c("Seurat", "dplyr", "ggplot2", "ggbeeswarm", "jsonlite")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop(sprintf("Missing required packages: %s\nInstall with: install.packages(c(%s))",
                 paste(missing_packages, collapse = ", "),
                 paste(sprintf("'%s'", missing_packages), collapse = ", ")))
  }
  
  axis <- match.arg(axis)
  
  # Check if micron conversion is requested
  microns_per_pixel <- NULL
  if (include_microns) {
    if (is.null(dataset) && is.null(scalefactors_path)) {
      stop("Must provide either 'dataset' or 'scalefactors_path' parameter to calculate micron-based distances")
    }
    
    # Determine scalefactors path
    if (is.null(scalefactors_path)) {
      # Build default path relative to working directory
      scalefactors_path <- file.path(
        dataset,
        paste0(dataset, "_Realignment"),
        "binned_outputs",
        "square_008um",
        "spatial",
        "scalefactors_json.json"
      )
      message(sprintf("Looking for scale factors at: %s", scalefactors_path))
    }
    
    if (!file.exists(scalefactors_path)) {
      warning(sprintf("Scale factors file not found at %s. Proceeding without micron conversion.", scalefactors_path))
      include_microns <- FALSE
    } else {
      scalefactors <- jsonlite::fromJSON(txt = scalefactors_path)
      
      # Validate scalefactors structure
      if (!"microns_per_pixel" %in% names(scalefactors)) {
        warning("'microns_per_pixel' not found in scalefactors file. Proceeding without micron conversion.")
        include_microns <- FALSE
      } else {
        microns_per_pixel <- scalefactors$microns_per_pixel
        message(sprintf("Using scale factor: %.3f microns per pixel", microns_per_pixel))
      }
    }
  }
  
  # --- Handle pre-calculated data if provided ---
  if (!is.null(projection_data)) {
    required_cols <- c("region", "celltype")
    projection_col <- switch(axis,
                             "epidermis_dermis" = "depth_position",
                             "length" = "length_position",
                             stop("Invalid axis selection."))
    
    required_cols <- c(required_cols, projection_col)
    
    if (!all(required_cols %in% colnames(projection_data))) {
      stop(sprintf("Pre-calculated data must contain columns: %s", 
                   paste(required_cols, collapse = ", ")))
    }
    
    combined <- projection_data %>%
      dplyr::rename(projection = !!projection_col, Region = region) %>%
      dplyr::select(Region, projection, celltype, dplyr::everything())
    
    if (invert_axis) {
      combined$projection <- 1 - combined$projection
    }
    
    if (!is.null(celltypes_to_plot)) {
      combined <- combined %>% dplyr::filter(celltype %in% celltypes_to_plot)
    }
    
    combined$celltype <- if (!is.null(celltype_levels)) {
      factor(combined$celltype, levels = celltype_levels)
    } else {
      factor(combined$celltype)
    }
    
    combined$timepoint <- timepoint_label
    combined$weighted_id <- 1
    combined$weight <- 1
    
    message("Using pre-calculated projections. Skipping internal calculation.")
    
  } else {
    
    # --- Determine regions and images ---
    if (is.null(image_names)) {
      # Assume Region metadata exists and each region has one image
      if (!"Region" %in% colnames(object@meta.data)) {
        stop("Either 'image_names' must be provided or object must have 'Region' in metadata")
      }
      regions <- unique(as.character(object$Region))
      
      # Try to automatically detect image names
      available_images <- names(object@images)
      if (length(available_images) == 0) {
        stop("No images found in Seurat object")
      }
      
      if (length(available_images) < length(regions)) {
        stop(sprintf("Found %d regions but only %d images. Please provide 'image_names' parameter.",
                     length(regions), length(available_images)))
      }
      
      # Create simple mapping
      image_names <- setNames(available_images[seq_along(regions)], regions)
      message("Auto-detected image mapping:")
      for (r in names(image_names)) {
        message(sprintf("  %s -> %s", r, image_names[r]))
      }
    } else {
      regions <- names(image_names)
      message("Using provided image mapping:")
      for (r in regions) {
        message(sprintf("  %s -> %s", r, image_names[r]))
      }
    }
    
    # Validate annotation column
    if (!annotation_col %in% colnames(object@meta.data)) {
      stop(sprintf("Annotation column '%s' not found in object metadata. Available columns: %s",
                   annotation_col, paste(colnames(object@meta.data), collapse = ", ")))
    }
    
    # --- Step 1: Compute projections for all regions/images ---
    proj_list <- list()
    meta_list <- list()
    scale_info <- list()
    
    for (region in regions) {
      message("Processing ", region, " (image: ", image_names[region], ")")
      
      image_name <- image_names[region]
      
      # Validate image exists
      if (!image_name %in% names(object@images)) {
        warning(sprintf("Image '%s' not found in object. Available images: %s. Skipping region %s.",
                        image_name, paste(names(object@images), collapse = ", "), region))
        next
      }
      
      # Get coordinates from the specific image
      coords <- tryCatch({
        Seurat::GetTissueCoordinates(object, image = image_name)
      }, error = function(e) {
        warning(sprintf("Could not get coordinates for image %s: %s", image_name, e$message))
        return(NULL)
      })
      
      if (is.null(coords) || nrow(coords) == 0) {
        warning(sprintf("No coordinates found for image %s. Skipping.", image_name))
        next
      }
      
      cells_with_coords <- rownames(coords)
      
      # Subset to cells with valid coordinates
      obj_subset <- subset(object, cells = cells_with_coords)
      
      message(sprintf("  %d cells with coordinates", length(cells_with_coords)))
      
      # Get coordinates and compute PCA
      coords_mat <- as.matrix(coords[, c("x", "y")])
      mode(coords_mat) <- "numeric"
      coords_centered <- scale(coords_mat, center = TRUE, scale = FALSE)
      
      pca <- prcomp(coords_centered)
      
      # Project onto both axes (for scale info)
      length_vec <- pca$rotation[, 1]
      length_vec <- length_vec / sqrt(sum(length_vec^2))
      depth_vec <- pca$rotation[, 2]
      depth_vec <- depth_vec / sqrt(sum(depth_vec^2))
      
      length_proj <- as.numeric(coords_centered %*% length_vec)
      depth_proj <- as.numeric(coords_centered %*% depth_vec)
      
      # Choose projection based on selected axis
      proj_values <- if (axis == "epidermis_dermis") depth_proj else length_proj
      
      proj_list[[region]] <- list(
        selected = proj_values,
        length = length_proj,
        depth = depth_proj,
        cells = cells_with_coords
      )
      
      # Store scale information if microns requested
      if (include_microns) {
        length_range_pixels <- max(length_proj) - min(length_proj)
        depth_range_pixels <- max(depth_proj) - min(depth_proj)
        
        scale_info[[region]] <- list(
          length_min = min(length_proj),
          length_max = max(length_proj),
          length_range_pixels = length_range_pixels,
          length_range_microns = length_range_pixels * microns_per_pixel,
          depth_min = min(depth_proj),
          depth_max = max(depth_proj),
          depth_range_pixels = depth_range_pixels,
          depth_range_microns = depth_range_pixels * microns_per_pixel,
          microns_per_pixel = microns_per_pixel
        )
        
        message(sprintf("  Length axis spans: %.1f pixels = %.1f µm", 
                        length_range_pixels, length_range_pixels * microns_per_pixel))
        message(sprintf("  Depth axis spans: %.1f pixels = %.1f µm", 
                        depth_range_pixels, depth_range_pixels * microns_per_pixel))
      }
      
      meta_list[[region]] <- obj_subset@meta.data[cells_with_coords, , drop = FALSE]
      
      # Verify dimensions
      if (length(proj_values) != nrow(meta_list[[region]])) {
        stop(sprintf("Dimension mismatch in %s: %d projections vs %d metadata rows",
                     region, length(proj_values), nrow(meta_list[[region]])))
      }
    }
    
    if (length(proj_list) == 0) {
      stop("No valid regions processed. Check image names and data.")
    }
    
    # Update regions to only processed ones
    regions <- names(proj_list)
    
    # --- Step 2: Align axes across regions using shared cell types ---
    if (align_regions && length(regions) > 1) {
      message("\nAligning axes across regions using shared cell types...")
      
      mean_proj <- lapply(regions, function(region) {
        celltype_vals <- meta_list[[region]][[annotation_col]]
        valid_idx <- !is.na(celltype_vals) & celltype_vals != "" & celltype_vals != "Unknown"
        
        if (sum(valid_idx) == 0) {
          warning(sprintf("No valid cell types in %s for alignment", region))
          return(NULL)
        }
        
        df <- data.frame(
          celltype = celltype_vals[valid_idx],
          proj = proj_list[[region]]$selected[valid_idx]
        )
        
        df %>% 
          dplyr::group_by(celltype) %>% 
          dplyr::summarise(mean_proj = mean(proj, na.rm = TRUE), .groups = "drop")
      })
      names(mean_proj) <- regions
      
      mean_proj <- mean_proj[!sapply(mean_proj, is.null)]
      
      if (length(mean_proj) > 1) {
        ref_region <- names(mean_proj)[1]
        ref_means <- mean_proj[[ref_region]]
        
        for (region in names(mean_proj)[-1]) {
          curr_means <- mean_proj[[region]]
          shared <- intersect(ref_means$celltype, curr_means$celltype)
          
          if (length(shared) > 1) {
            ref_vals <- ref_means$mean_proj[match(shared, ref_means$celltype)]
            curr_vals <- curr_means$mean_proj[match(shared, curr_means$celltype)]
            
            if (cor(ref_vals, curr_vals, use = "complete.obs") < 0) {
              message("  Flipping axis for ", region, " to align with ", ref_region)
              proj_list[[region]]$selected <- -proj_list[[region]]$selected
              
              # Also flip in the length/depth specific projections
              if (axis == "epidermis_dermis") {
                proj_list[[region]]$depth <- -proj_list[[region]]$depth
                if (include_microns) {
                  old_min <- scale_info[[region]]$depth_min
                  old_max <- scale_info[[region]]$depth_max
                  scale_info[[region]]$depth_min <- -old_max
                  scale_info[[region]]$depth_max <- -old_min
                }
              } else {
                proj_list[[region]]$length <- -proj_list[[region]]$length
                if (include_microns) {
                  old_min <- scale_info[[region]]$length_min
                  old_max <- scale_info[[region]]$length_max
                  scale_info[[region]]$length_min <- -old_max
                  scale_info[[region]]$length_max <- -old_min
                }
              }
            }
          } else {
            message(sprintf("  Warning: Only %d shared cell type(s) between %s and %s", 
                            length(shared), ref_region, region))
          }
        }
      }
    }
    
    # --- Step 3: Normalize projections and combine ---
    combined_list <- lapply(regions, function(region) {
      proj_values <- proj_list[[region]]$selected
      length_proj <- proj_list[[region]]$length
      depth_proj <- proj_list[[region]]$depth
      cells <- proj_list[[region]]$cells
      
      # Apply manual inversion if requested
      if (invert_axis) {
        proj_values <- -proj_values
      }
      
      # Normalize to 0-1
      proj_norm <- (proj_values - min(proj_values)) / (max(proj_values) - min(proj_values))
      
      meta_proj <- meta_list[[region]]
      
      # Build base data frame
      result_df <- data.frame(
        barcode = cells,
        projection = proj_norm,
        Region = region,
        row.names = cells
      )
      
      # Add micron-based distances if requested
      if (include_microns) {
        length_norm <- (length_proj - min(length_proj)) / (max(length_proj) - min(length_proj))
        depth_norm <- (depth_proj - min(depth_proj)) / (max(depth_proj) - min(depth_proj))
        
        # For depth: invert so 0 = dermis, 1 = epidermis (matching gene function)
        if (axis == "epidermis_dermis" && invert_axis) {
          depth_norm <- 1 - depth_norm
        }
        
        result_df$length_position <- length_norm
        result_df$depth_position <- depth_norm
        result_df$length_raw <- length_proj
        result_df$depth_raw <- depth_proj
        result_df$length_microns <- length_norm * scale_info[[region]]$length_range_microns
        
        # Depth microns from epidermis (0 µm = epidermis)
        result_df$depth_microns <- (1 - depth_norm) * scale_info[[region]]$depth_range_microns
      }
      
      # Add cell type annotation
      if (annotation_col %in% colnames(meta_proj)) {
        result_df$celltype <- meta_proj[cells, annotation_col]
      } else {
        warning(sprintf("Annotation column '%s' not found in %s metadata", annotation_col, region))
        result_df$celltype <- NA
      }
      
      # Filter and weight cells
      result_df <- result_df %>%
        dplyr::mutate(
          weight = dplyr::case_when(
            !is.na(celltype) & 
              (is.na(meta_proj[cells, "second_type"]) | meta_proj[cells, "singlet_score"] > 0.5) ~ 1,
            !is.na(celltype) & 
              !is.na(meta_proj[cells, "second_type"]) & meta_proj[cells, "singlet_score"] <= 0.5 ~ 0.5,
            TRUE ~ NA_real_
          ),
          metaclusters = if("metaclusters" %in% colnames(meta_proj)) meta_proj[cells, "metaclusters"] else NA,
          timepoint = timepoint_label,
          weighted_id = dplyr::row_number()
        ) %>%
        dplyr::filter(
          !is.na(celltype), 
          !is.na(weight), 
          celltype != "Unknown"
        )
      
      return(result_df)
    })
    
    combined <- dplyr::bind_rows(combined_list)
    
    if (nrow(combined) == 0) {
      stop("No valid cells remaining after filtering. Check annotation column and cell type labels.")
    }
    
    # Filter to specific cell types if requested
    if (!is.null(celltypes_to_plot)) {
      before_filter <- nrow(combined)
      combined <- combined %>% dplyr::filter(celltype %in% celltypes_to_plot)
      after_filter <- nrow(combined)
      
      if (after_filter == 0) {
        stop(sprintf("No cells match celltypes_to_plot. Available cell types: %s",
                     paste(unique(combined$celltype), collapse = ", ")))
      }
      message(sprintf("Filtered to %d cells (from %d) matching specified cell types", 
                      after_filter, before_filter))
    }
    
    # Set cell type levels
    combined$celltype <- if (!is.null(celltype_levels)) {
      factor(combined$celltype, levels = celltype_levels)
    } else {
      factor(combined$celltype)
    }
    
    # Print summary
    message("\n=== Summary ===")
    for (region in regions) {
      n_cells <- sum(combined$Region == region)
      message(sprintf("%s: %d cells", region, n_cells))
    }
    message(sprintf("Total: %d cells across %d cell types", 
                    nrow(combined), length(unique(combined$celltype))))
    
    # Attach scale information as attributes
    if (include_microns) {
      attr(combined, "scale_info") <- scale_info
      attr(combined, "microns_per_pixel") <- microns_per_pixel
    }
  }
  
  # --- Step 4: Generate plot if requested ---
  if (return_plot) {
    # Determine which variable to plot
    if (include_microns) {
      plot_var <- if (axis == "epidermis_dermis") "depth_microns" else "length_microns"
      ylab <- if (axis == "epidermis_dermis") {
        "Distance from epidermis (µm)"
      } else {
        "Position along tissue length (µm)"
      }
    } else {
      plot_var <- "projection"
      ylab <- if (axis == "epidermis_dermis") {
        "Position along epidermis–dermis axis"
      } else {
        "Position along tissue length"
      }
    }
    
    # Create base plot
    p <- if (combine_regions) {
      ggplot2::ggplot(combined, ggplot2::aes(x = celltype, y = .data[[plot_var]], color = celltype)) +
        ggbeeswarm::geom_quasirandom(groupOnX = TRUE, alpha = 0.8, size = 1.2) +
        ggplot2::labs(x = annotation_col, y = ylab) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          legend.position = "none", 
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    } else {
      ggplot2::ggplot(combined, ggplot2::aes(x = celltype, y = .data[[plot_var]], color = celltype)) +
        ggbeeswarm::geom_quasirandom(groupOnX = TRUE, alpha = 0.8, size = 1.2) +
        ggplot2::facet_wrap(~Region, scales = "free_x") +
        ggplot2::labs(x = annotation_col, y = ylab) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          legend.position = "none", 
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    }
    
    # Add colors if provided
    if (!is.null(cols)) {
      p <- p + ggplot2::scale_color_manual(values = cols)
    }
    
    # Add wound center line if applicable
    if (!is.null(wound_center_microns) && axis == "length" && include_microns) {
      p <- p + 
        ggplot2::geom_hline(yintercept = wound_center_microns, 
                            linetype = "dashed", 
                            color = "red", 
                            linewidth = 0.8) +
        ggplot2::annotate("text", 
                          x = Inf, 
                          y = wound_center_microns, 
                          label = "Wound Center", 
                          hjust = 1.1, 
                          vjust = -0.5, 
                          color = "red", 
                          size = 3.5)
    }
    
    print(p)
  }
  
  invisible(combined)
}