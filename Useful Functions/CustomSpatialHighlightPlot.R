library(Seurat)
library(dplyr)
library(ggplot2)

CustomSpatialHighlightPlot <- function(
    seurat_obj,
    region,
    annotation_col,
    cells_highlight,
    pt.size.highlight = 1.2,
    pt.size.unselected = 0.3,
    stroke.size = 0.2,
    pal,
    custom_theme
) {
  Idents(seurat_obj) <- annotation_col
  annotation_sym <- sym(annotation_col)
  highlight_cells <- unlist(CellsByIdentities(seurat_obj, idents = cells_highlight))
  
  coords <- GetTissueCoordinates(seurat_obj, image = region)
  # Flip y-axis to match SpatialDimPlot orientation
  tmp <- coords$x
  coords$x <- coords$y
  coords$y <- tmp
  
  barcodes <- rownames(coords)
  
  coords$highlight <- barcodes %in% highlight_cells
  coords$cell_type <- as.character(seurat_obj@meta.data[barcodes, annotation_col])
  coords$cell_type[!coords$cell_type %in% cells_highlight] <- "Unselected"
  
  cat("Number of highlighted cells: ", sum(coords$highlight), "\n")
  cat("Number of unselected cells: ", sum(!coords$highlight), "\n")
  
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
    scale_fill_manual(values = pal, drop = FALSE) +
    scale_y_reverse() +
    custom_theme +
    theme_void()
  
  return(p)
}