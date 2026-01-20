# This function will save png and svg outputs for non-ggplot plots
# Updated save_plot function with circlize canvas parameters
save_plot <- function(filename, plot_expr, 
                      png_width = 2000, png_height = 2000, png_res = 300,
                      svg_width = 9, svg_height = 9) {
  
  # Remove file extension if provided
  filename <- sub("\\.(png|svg)$", "", filename)
  
  # Check if it's already an object (like ComplexHeatmap or pheatmap) vs an expression to evaluate
  if (inherits(plot_expr, c("Heatmap", "HeatmapList", "AdditiveUnit"))) {
    # It's a ComplexHeatmap object - save directly
    png(paste0(filename, ".png"), width = png_width, height = png_height, res = png_res)
    draw(plot_expr)
    dev.off()
    
    svg(paste0(filename, ".svg"), width = svg_width, height = svg_height)
    draw(plot_expr)
    dev.off()
    
  } else if (inherits(plot_expr, "pheatmap")) {
    # It's a pheatmap object - use grid graphics
    png(paste0(filename, ".png"), width = png_width, height = png_height, res = png_res)
    grid::grid.newpage()
    grid::grid.draw(plot_expr$gtable)
    dev.off()
    
    svg(paste0(filename, ".svg"), width = svg_width, height = svg_height)
    grid::grid.newpage()
    grid::grid.draw(plot_expr$gtable)
    dev.off()
    
  } else {
    # It's a base R graphics expression - capture and evaluate
    if (is.call(substitute(plot_expr)) || is.name(substitute(plot_expr))) {
      plot_expr <- substitute(plot_expr)
    }
    
    # Save PNG
    circlize::circos.clear()
    png(paste0(filename, ".png"), width = png_width, height = png_height, res = png_res)
    par(mfrow = c(1, 1), xpd = TRUE)
    # Set circlize parameters for more space
    circlize::circos.par(canvas.xlim = c(-1.5, 1.5), canvas.ylim = c(-1.5, 1.5),
                         track.margin = c(0.005, 0.005))
    eval(plot_expr, envir = parent.frame())
    dev.off()
    circlize::circos.clear()
    
    # Save SVG
    circlize::circos.clear()
    svg(paste0(filename, ".svg"), width = svg_width, height = svg_height)
    par(mfrow = c(1, 1), xpd = TRUE)
    # Set circlize parameters for more space
    circlize::circos.par(canvas.xlim = c(-1.5, 1.5), canvas.ylim = c(-1.5, 1.5),
                         track.margin = c(0.005, 0.005))
    eval(plot_expr, envir = parent.frame())
    dev.off()
    circlize::circos.clear()
  }
  
  message("Saved: ", filename, ".png and ", filename, ".svg")
}