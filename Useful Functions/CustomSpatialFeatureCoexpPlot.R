library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

CustomSpatialFeatureCoexpPlot <- function(
    object,
    image,
    genes,
    gene_colors,
    assay         = NULL,
    layer         = "data",
    image.scale   = "hires",
    pt.size.expr  = 1.5,     # size for expressing spots
    pt.size.zero  = 0.5,     # size for zero-expression spots
    alpha.expr = NA,
    alpha.zero = NA,
    stroke        = 0.1,
    bg.col        = "white"  # color for zero-expression spots
) {
  
  if (!is.null(assay)) DefaultAssay(object) <- assay
  
  # get image + raster
  if (!(image %in% names(object@images)))
    stop("Image '", image, "' not found")
  spat_img <- object@images[[image]]
  img_rast <- GetImage(spat_img, mode = "raster")
  
  # get spot coords
  coords <- as.data.frame(GetTissueCoordinates(spat_img, scale = image.scale))
  coords$barcode <- rownames(coords)
  
  # flip coordinates so the image is in the same orientation as default plots
  tmp <- coords$x
  coords$x <- coords$y
  coords$y <- tmp
  
  # fetch expression
  expr <- FetchData(object, vars = genes, cells = coords$barcode,
                    layer = layer, clean = FALSE)
  
  # 4) blend colors
  expr_scaled <- as.data.frame(lapply(expr, scales::rescale))
  rgb_mat     <- t(grDevices::col2rgb(gene_colors)/255)
  wts         <- expr_scaled / rowSums(expr_scaled + 1e-6)
  blend       <- as.matrix(wts) %*% rgb_mat
  spot_cols   <- grDevices::rgb(blend[,1], blend[,2], blend[,3])
  names(spot_cols) <- rownames(expr)
  
  # mark zero-expression spots
  zero_bcs <- rownames(expr)[rowSums(expr)==0]
  spot_cols[zero_bcs] <- bg.col
  coords$col <- spot_cols[coords$barcode]
  
  # split coords
  coords_zero <- coords[coords$barcode %in% zero_bcs, ]
  coords_expr <- coords[!coords$barcode %in% zero_bcs, ]
  
  # plot limits
  xlim <- range(coords$x); ylim <- range(coords$y)
  
  # --- Main spatial plot
  spatial_plot <- ggplot() +
    # annotation_raster(
    #   img_rast,
    #   xmin = xlim[1], xmax = xlim[2],
    #   ymin = ylim[1], ymax = ylim[2],
    #   interpolate = TRUE
    # ) +
    geom_point(
      data = coords_zero,
      aes(x = x, y = y),
      color = "grey75",
      alpha = alpha.zero,
      shape = 21,
      fill  = bg.col,
      size  = pt.size.zero,
      stroke= stroke
    ) +
    geom_point(
      data = coords_expr,
      aes(x = x, y = y, fill = col),
      shape = 21,
      alpha = alpha.expr,
      size  = pt.size.expr,
      stroke= stroke
    ) +
    scale_fill_identity() +
    coord_fixed() +
    scale_y_reverse() +
    theme_void()
  
  # --- Make gene legend strips
  make_gene_legend <- function(gene, color) {
    df <- data.frame(x = 1:100, y = 1, val = seq(0, 1, length.out = 100))
    ggplot(df, aes(x = x, y = y, fill = val)) +
      geom_tile() +
      scale_fill_gradientn(colors = c("white", color), name = gene) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.width = unit(1.2, "cm"),
        plot.margin = margin(2, 2, 2, 2)
      ) +
      guides(fill = guide_colorbar(title.position = "top", barwidth = 5))
  }
  
  legend_plots <- lapply(seq_along(genes), function(i) {
    make_gene_legend(genes[i], gene_colors[i]) +
      theme(legend.title = element_text(face = "italic"))
  })
  
  legends_combined <- plot_grid(plotlist = legend_plots, nrow = 1)
  
  # --- Combine and return
  final_plot <- plot_grid(spatial_plot, legends_combined,
                          ncol = 1, rel_heights = c(1, 0.15))
  return(final_plot)
}
