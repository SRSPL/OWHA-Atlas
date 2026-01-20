# Define a function to save plots in multiple formats
ggsave_multidevice <- function(plot = ggplot2::last_plot(), filename_base, 
                               devices = c("svg", "png"), width = 7, height = 7,
                               units = "in", dpi = 300) {
  for (device in devices) {
    ggsave(
      filename = paste0(filename_base, ".", device),
      plot = plot,
      width = width,
      height = height,
      units = units,
      dpi = dpi,
      bg = "transparent"
    )
  }
}
