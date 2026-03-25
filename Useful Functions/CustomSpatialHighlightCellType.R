library(Seurat)
library(dplyr)
library(ggplot2)

CustomSpatialHighlightPlot <- function(
    seurat_obj,
    dataset,
    region,
    annotation_col,
    cells_highlight,
    pt.size.highlight = 1.2,
    pt.size.unselected = 0.3,
    stroke.size = 0.2,
    pal,
    custom_theme,
    scalebar_length_um = 100,   # length of bar in microns
    scalebar_position = c("bottomright", "bottomleft", "topright", "topleft") # location on image
) {
  # Load scalefactors JSON for this dataset
  scalefactors_path <- file.path(
    paste0("~/SRSP Laboratory Dropbox/SRSP Lab/Team/Immunodynamics/Simon Van Deursen/DATA/Spatial Transcriptomics/",
           dataset, "/", dataset, "_Realignment/binned_outputs/square_008um/spatial/scalefactors_json.json")
  )
  scalefactors <- jsonlite::fromJSON(txt = scalefactors_path)
  
  Idents(seurat_obj) <- annotation_col
  annotation_sym <- sym(annotation_col)
  highlight_cells <- unlist(CellsByIdentities(seurat_obj, idents = cells_highlight))
  
  coords <- GetTissueCoordinates(seurat_obj, image = region)
  
  # Flip x/y to match SpatialDimPlot orientation
  tmp <- coords$x
  coords$x <- coords$y
  coords$y <- tmp
  
  barcodes <- rownames(coords)
  coords$highlight <- barcodes %in% highlight_cells
  coords$cell_type <- as.character(seurat_obj@meta.data[barcodes, annotation_col])
  coords$cell_type[!coords$cell_type %in% cells_highlight] <- "Unselected"
  
  cat("Number of highlighted cells: ", sum(coords$highlight), "\n")
  cat("Number of unselected cells: ", sum(!coords$highlight), "\n")
  
  # --- Scale bar calculation ---
  microns_per_pixel <- scalefactors$microns_per_pixel
  scalebar_length_px <- scalebar_length_um / microns_per_pixel
  
  # Position for the bar
  scalebar_position <- match.arg(scalebar_position)
  x_min <- min(coords$x)
  x_max <- max(coords$x)
  y_min <- min(coords$y)
  y_max <- max(coords$y)
  
  margin <- 0.02 * (x_max - x_min)
  bar_height <- 0.01 * (y_max - y_min)
  
  if (scalebar_position == "bottomright") {
    x_start <- x_max - scalebar_length_px - margin
    x_end <- x_max - margin
    y_bar <- y_max - margin
    text_vjust <- -0.8
  } else if (scalebar_position == "bottomleft") {
    x_start <- x_min + margin
    x_end <- x_min + scalebar_length_px + margin
    y_bar <- y_max - margin
    text_vjust <- -0.8
  } else if (scalebar_position == "topright") {
    x_start <- x_max - scalebar_length_px - margin
    x_end <- x_max - margin
    y_bar <- y_min + margin
    text_vjust <- 1.8
  } else {
    x_start <- x_min + margin
    x_end <- x_min + scalebar_length_px + margin
    y_bar <- y_min + margin
    text_vjust <- 1.8
  }
  
  # --- Plot ---
  p <- ggplot() +
    geom_point(
      data = subset(coords, !highlight),
      aes(x = x, y = y, fill = cell_type),
      shape = 21,
      size = pt.size.unselected,
      alpha = 1,
      color = "black",
      stroke = stroke.size
    ) +
    geom_point(
      data = subset(coords, highlight),
      aes(x = x, y = y, fill = cell_type),
      shape = 21,
      size = pt.size.highlight,
      alpha = 1,
      color = "black",
      stroke = stroke.size
    ) +
    # --- Fixed scale bar (no warning, proper label position) ---
    annotate(
      "segment",
      x = x_start, xend = x_end,
      y = y_bar, yend = y_bar,
      color = "black", linewidth = 1.2, lineend = "round"
    ) +
    annotate(
      "text",
      x = mean(c(x_start, x_end)),
      y = y_bar - 0.03 * (y_max - y_min),  # 10% of y-range above the bar
      label = paste0(scalebar_length_um, " µm"),
      vjust = 0.5, hjust = 0.5, size = 3.5
  ) +
    scale_fill_manual(values = pal, drop = FALSE) +
    scale_y_reverse() +
    coord_equal() +
    custom_theme +
    theme_void()
  
  return(p)
}
