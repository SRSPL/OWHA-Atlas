#' Test for Cell Type Spatial Localization Along Tissue Axes
#' 
#' @description
#' Tests whether cell types show significant spatial localization (non-uniform distribution)
#' along a tissue axis by comparing observed variance in position to a uniform null distribution.
#' Supports multi-region replicates with statistical combination using Stouffer's method.
#' 
#' @param proj_values Numeric vector of normalized positions along tissue axis (0-1 scale).
#'   Must be same length as cell_labels
#' @param cell_labels Character or factor vector of cell type assignments. Must be same length as proj_values
#' @param regions Optional character or factor vector of replicate/region identifiers. If NULL, treats all
#'   cells as single replicate. Must be same length as proj_values if provided
#' @param n_iter Integer, number of bootstrap iterations for null distribution simulation (default: 10000).
#'   Used for small sample sizes (< small_n_threshold)
#' @param small_n_threshold Integer, threshold for switching from bootstrap to analytical chi-square test
#'   (default: 500). Cell types with fewer cells use bootstrap; larger use chi-square approximation
#' @param seed Integer, random seed for reproducibility (default: 123)
#' @param z_threshold Numeric, minimum absolute z-score for significance. If NULL (default), 
#'   automatically set to 10th percentile of observed |z| values. Cell types below this threshold
#'   are labeled as "ns" regardless of p-value
#' 
#' @return Data frame with one row per cell type containing:
#' \describe{
#'   \item{celltype}{Cell type label}
#'   \item{var_proj}{Observed variance in position (pooled across regions)}
#'   \item{z_var}{Combined z-score from Stouffer's method (positive = localized, negative = dispersed)}
#'   \item{p_var}{Two-tailed p-value for localization test}
#'   \item{p_var_adj}{Benjamini-Hochberg adjusted p-value for multiple testing}
#'   \item{sig_label}{Significance indicator: "***" (p<0.001), "**" (p<0.01), "*" (p<0.05), "ns" (not significant)}
#' }
#' 
#' @details
#' The function tests the null hypothesis that cells are uniformly distributed along the tissue axis
#' (variance = 1/12 for uniform [0,1] distribution). For each cell type:
#' 
#' **Small samples (n < small_n_threshold):**
#' - Simulates null distribution by bootstrap sampling from uniform(0,1)
#' - Calculates z-score and p-value empirically
#' 
#' **Large samples (n >= small_n_threshold):**
#' - Uses chi-square approximation for faster computation
#' - Z-score calculated as (obs_var - 1/12) / SE(var)
#' 
#' **Multiple regions:**
#' When regions are provided, the test:
#' 1. Calculates z-scores independently for each region
#' 2. Combines using Stouffer's method: Z_combined = sum(Z_i) / sqrt(k)
#' 3. Accounts for variation across biological replicates
#' 
#' **Interpretation:**
#' - Positive z-scores indicate spatial localization (variance < uniform)
#' - Negative z-scores indicate spatial dispersion (variance > uniform)
#' - Magnitude indicates strength of spatial pattern
#' 
#' @examples
#' \dontrun{
#' # Single region analysis
#' results <- cluster_localization_test_reps(
#'   proj_values = data$position,
#'   cell_labels = data$celltype,
#'   n_iter = 10000
#' )
#' 
#' # Multi-region replicate analysis
#' results <- cluster_localization_test_reps(
#'   proj_values = data$projection,
#'   cell_labels = data$celltype,
#'   regions = data$Region,
#'   n_iter = 10000,
#'   seed = 456
#' )
#' 
#' # View most localized cell types
#' head(results[order(-results$z_var), ])
#' }
#' 
#' @references
#' Stouffer's method for combining independent tests:
#' Stouffer, S.A., et al. (1949). The American Soldier, Vol. 1: Adjustment during Army Life.
#' 
#' @importFrom stats var pchisq pnorm p.adjust quantile runif
#' 
#' @export
cluster_localization_test_reps <- function(
    proj_values, 
    cell_labels, 
    regions = NULL,
    n_iter = 10000, 
    small_n_threshold = 500,
    seed = 123, 
    z_threshold = NULL
) {
  
  # Input validation
  if (!is.numeric(proj_values)) {
    stop("'proj_values' must be a numeric vector")
  }
  
  if (length(proj_values) != length(cell_labels)) {
    stop(sprintf("Length mismatch: proj_values (%d) and cell_labels (%d) must have same length",
                 length(proj_values), length(cell_labels)))
  }
  
  if (!is.null(regions) && length(regions) != length(proj_values)) {
    stop(sprintf("Length mismatch: regions (%d) must match proj_values (%d)",
                 length(regions), length(proj_values)))
  }
  
  # Check for valid projection range
  if (any(!is.na(proj_values) & (proj_values < 0 | proj_values > 1))) {
    warning("Some proj_values are outside [0,1] range. Results may be unreliable. Consider normalizing data.")
  }
  
  # Check for missing data
  if (any(is.na(proj_values))) {
    n_missing <- sum(is.na(proj_values))
    warning(sprintf("Found %d NA values in proj_values (%.1f%%). These will be excluded.",
                    n_missing, 100 * n_missing / length(proj_values)))
  }
  
  if (any(is.na(cell_labels))) {
    n_missing <- sum(is.na(cell_labels))
    warning(sprintf("Found %d NA values in cell_labels (%.1f%%). These will be excluded.",
                    n_missing, 100 * n_missing / length(cell_labels)))
  }
  
  # Parameter validation
  if (n_iter < 100) {
    warning("n_iter < 100 may produce unreliable p-values. Consider using n_iter >= 1000")
  }
  
  if (small_n_threshold < 30) {
    warning("small_n_threshold < 30 may lead to unreliable chi-square approximations")
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Handle regions
  if (is.null(regions)) {
    regions <- rep("R0", length(proj_values))
    message("No regions provided. Treating all cells as single replicate.")
  } else {
    n_regions <- length(unique(regions[!is.na(regions)]))
    message(sprintf("Analyzing %d region(s)", n_regions))
  }
  regions <- as.character(regions)
  
  # Get unique cell types
  unique_cells <- sort(unique(as.character(cell_labels[!is.na(cell_labels)])))
  n_celltypes <- length(unique_cells)
  
  if (n_celltypes == 0) {
    stop("No valid cell types found after removing NAs")
  }
  
  message(sprintf("Testing %d cell type(s)", n_celltypes))
  
  # Initialize results data frame
  results <- data.frame(
    celltype = character(),
    n_cells = integer(),
    n_regions = integer(),
    var_proj = numeric(),
    z_var = numeric(),
    p_var = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop over each cell type
  for (cell in unique_cells) {
    inds <- which(cell_labels == cell & !is.na(proj_values))
    n_cell <- length(inds)
    
    if (n_cell == 0) {
      warning(sprintf("Cell type '%s' has no valid data points. Skipping.", cell))
      next
    }
    
    if (n_cell < 5) {
      warning(sprintf("Cell type '%s' has only %d cells. Results may be unreliable.", cell, n_cell))
    }
    
    # Observed variance pooled across all regions
    obs_var <- var(proj_values[inds], na.rm = TRUE)
    
    if (is.na(obs_var)) {
      warning(sprintf("Could not compute variance for cell type '%s'. Skipping.", cell))
      next
    }
    
    # ---- Region-specific analysis ----
    region_results <- list()
    regions_tested <- character()
    
    for (r in unique(regions[inds])) {
      inds_r <- inds[regions[inds] == r & !is.na(regions[inds])]
      n_r <- length(inds_r)
      
      if (n_r < 2) {
        # Need at least 2 cells to compute variance
        next
      }
      
      obs_var_r <- var(proj_values[inds_r], na.rm = TRUE)
      
      if (is.na(obs_var_r)) {
        next
      }
      
      # Simulate null variance distribution for this region
      if (n_r < small_n_threshold) {
        # Bootstrap method for small samples
        null_vars <- replicate(n_iter, var(runif(n_r)))
        null_center <- mean(null_vars)
        null_sd <- sd(null_vars)
        
        if (null_sd == 0) {
          warning(sprintf("Zero standard deviation in null for %s, region %s. Skipping.", cell, r))
          next
        }
        
        z_r <- (obs_var_r - null_center) / null_sd
        # Two-tailed p-value
        p_r <- (sum(abs(null_vars - null_center) >= abs(obs_var_r - null_center)) + 1) / (n_iter + 1)
        
      } else {
        # Analytical chi-square approximation for large samples
        null_var <- 1/12  # Variance of uniform(0,1)
        chi_stat <- (n_r - 1) * obs_var_r / null_var
        
        # Two-tailed test
        p_r <- 2 * min(pchisq(chi_stat, df = n_r - 1),
                       1 - pchisq(chi_stat, df = n_r - 1))
        
        # Standard error of variance for uniform distribution
        se_var <- sqrt((2 * (null_var^2)) / (n_r - 1))
        z_r <- (obs_var_r - null_var) / se_var
      }
      
      region_results[[r]] <- list(z = z_r, p = p_r, n = n_r)
      regions_tested <- c(regions_tested, r)
    }
    
    if (length(region_results) == 0) {
      warning(sprintf("Cell type '%s' has insufficient data in all regions. Skipping.", cell))
      next
    }
    
    # ---- Combine regions using Stouffer's method ----
    zs <- sapply(region_results, `[[`, "z")
    n_regions_tested <- length(zs)
    
    # Stouffer's Z: sum of z-scores divided by sqrt(number of tests)
    combined_z <- sum(zs) / sqrt(n_regions_tested)
    
    # Two-tailed p-value from standard normal
    combined_p <- 2 * pnorm(-abs(combined_z))
    
    # Store results
    results <- rbind(results,
                     data.frame(
                       celltype = cell,
                       n_cells = n_cell,
                       n_regions = n_regions_tested,
                       var_proj = obs_var,
                       z_var = combined_z,
                       p_var = combined_p,
                       stringsAsFactors = FALSE
                     ))
  }
  
  if (nrow(results) == 0) {
    stop("No cell types had sufficient data for testing")
  }
  
  # Multiple testing correction (Benjamini-Hochberg FDR)
  results$p_var_adj <- p.adjust(results$p_var, method = "BH")
  
  # Generate significance labels based on adjusted p-values
  stars <- cut(results$p_var_adj,
               breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
               labels = c("***", "**", "*", "ns"),
               right = FALSE)
  results$sig_label <- as.character(stars)
  
  # Optional z-threshold filter to avoid spurious low-variance calls
  if (is.null(z_threshold)) {
    # Use 10th percentile of absolute z-scores as threshold
    z_threshold <- quantile(abs(results$z_var), 0.10, na.rm = TRUE)
    message(sprintf("Auto-detected z-threshold: %.2f", z_threshold))
  } else {
    message(sprintf("Using provided z-threshold: %.2f", z_threshold))
  }
  
  # Mark as non-significant if below threshold
  n_filtered <- sum(abs(results$z_var) < z_threshold & results$sig_label != "ns")
  results$sig_label[abs(results$z_var) < z_threshold] <- "ns"
  
  if (n_filtered > 0) {
    message(sprintf("Filtered %d cell type(s) below z-threshold", n_filtered))
  }
  
  # Sort by significance and effect size
  results <- results[order(results$p_var_adj, -abs(results$z_var)), ]
  rownames(results) <- NULL
  
  # Summary statistics
  n_sig <- sum(results$sig_label != "ns")
  message("\n=== Summary ===")
  message(sprintf("Tested: %d cell types", nrow(results)))
  message(sprintf("Significant (adj p < 0.05): %d (%.1f%%)", 
                  n_sig, 100 * n_sig / nrow(results)))
  message(sprintf("  *** (p < 0.001): %d", sum(results$sig_label == "***")))
  message(sprintf("  **  (p < 0.01):  %d", sum(results$sig_label == "**")))
  message(sprintf("  *   (p < 0.05):  %d", sum(results$sig_label == "*")))
  
  # Add effect size interpretation column
  results$pattern <- ifelse(results$z_var > 0, "localized", "dispersed")
  results$pattern[results$sig_label == "ns"] <- "uniform"
  
  return(results)
}


#' Plot Localization Test Results
#' 
#' @description
#' Visualizes results from cluster_localization_test_reps as a volcano plot
#' showing z-scores vs adjusted p-values
#' 
#' @param results Data frame from cluster_localization_test_reps
#' @param label_top Integer, number of top cell types to label (default: 10)
#' @param z_threshold Numeric, z-score threshold line to draw (default: NULL)
#' 
#' @return ggplot2 object
#' 
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_hline geom_vline scale_color_manual labs theme_minimal
#' @importFrom ggrepel geom_text_repel
#' 
#' @export
plot_localization_results <- function(results, label_top = 10, z_threshold = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting. Install with: install.packages('ggplot2')")
  }
  
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    message("Package 'ggrepel' recommended for better label placement. Install with: install.packages('ggrepel')")
  }
  
  # Add -log10 p-value for plotting
  results$neg_log10_p <- -log10(results$p_var_adj)
  
  # Identify top cell types to label
  results$label <- ""
  if (label_top > 0) {
    top_idx <- head(order(-abs(results$z_var)), label_top)
    results$label[top_idx] <- results$celltype[top_idx]
  }
  
  # Color by significance
  results$color_group <- factor(results$sig_label, 
                                levels = c("***", "**", "*", "ns"))
  
  p <- ggplot2::ggplot(results, ggplot2::aes(x = z_var, y = neg_log10_p, 
                                             color = color_group)) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::scale_color_manual(
      values = c("***" = "#D62728", "**" = "#FF7F0E", "*" = "#2CA02C", "ns" = "gray60"),
      name = "Significance"
    ) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "gray40") +
    ggplot2::labs(
      x = "Z-score (positive = localized, negative = dispersed)",
      y = "-log10(adjusted p-value)",
      title = "Cell Type Spatial Localization"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  
  # Add z-threshold lines if provided
  if (!is.null(z_threshold)) {
    p <- p + 
      ggplot2::geom_vline(xintercept = c(-z_threshold, z_threshold), 
                          linetype = "dotted", color = "blue")
  }
  
  # Add labels
  if (label_top > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5
    )
  } else if (label_top > 0) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = label),
      size = 3,
      vjust = -0.5
    )
  }
  
  return(p)
}


# ============================================================================
# Example Usage
# ============================================================================

# Single region analysis
# results <- cluster_localization_test_reps(
#   proj_values = data$projection,
#   cell_labels = data$celltype,
#   n_iter = 10000,
#   seed = 123
# )

# Multi-region replicate analysis (recommended for spatial transcriptomics)
# results <- cluster_localization_test_reps(
#   proj_values = df_uw$projection,
#   cell_labels = df_uw$celltype,
#   regions = df_uw$Region,
#   n_iter = 10000,
#   small_n_threshold = 500,
#   seed = 123
# )

# View results
# print(head(results))

# Plot results
# plot_localization_results(results, label_top = 15)