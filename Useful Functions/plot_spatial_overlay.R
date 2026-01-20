library(tidyverse)
library(dplyr)
library(Seurat)
library(ggplot2)

plot_spatial_overlay <- function(
    obj,
    gene,
    highlight_group,
    annotation_col = NULL,
    image             = "slice1",
    low_color         = "white",
    high_color        = "red",
    highlight_color   = "darkblue",
    pt.size.highlight = 1.2,
    pt.size.unselected= 0.3,
    stroke.size       = 0.05,
    min_expr          = 0.1
) {
  coords_df <- GetTissueCoordinates(obj, image = image) %>%
    as.data.frame() %>%
    rownames_to_column("barcode") %>%
    dplyr::rename(x = 2, y = 3)
  
  # Invert x and y to match expected orientation
  tmp <- coords_df$x
  coords_df$x <- coords_df$y
  coords_df$y <- tmp
  
  DefaultAssay(obj) <- "Spatial"
  
  expr_df <- FetchData(obj, vars = gene) %>%
    rownames_to_column("barcode") %>%
    dplyr::rename(expr = !!gene)
  
  label_df <- tibble(barcode = colnames(obj), label = obj@meta.data[[annotation_col]])
  
  df <- coords_df %>%
    left_join(expr_df, by = "barcode") %>%
    left_join(label_df, by = "barcode") %>%
    mutate(
      is_na     = is.na(expr),
      highlight = label %in% highlight_group,
      low_expr  = expr < min_expr,
      expressed = !is_na & !low_expr
    )
  
  ggplot() +
    # NA cells as background
    geom_point(
      data = df %>% filter(is_na),
      aes(x = x, y = y),
      shape = 21, fill = "grey90", color = "black", stroke = 0.05,
      size = pt.size.unselected, alpha = 0.5
    ) +
    # Low-expression, unhighlighted cells → white fill
    geom_point(
      data = df %>% filter(!highlight & !expressed),
      aes(x = x, y = y),
      shape = 21, fill = "white", color = "black", stroke = stroke.size,
      size = pt.size.unselected, alpha = 0.75
    ) +
    # Highlighted cells with low/no expression → solid highlight fill
    geom_point(
      data = df %>% filter(highlight & !expressed),
      aes(x = x, y = y),
      shape = 21, fill = highlight_color, color = "black", stroke = stroke.size,
      size = pt.size.highlight, alpha = 1
    ) +
    # Non-highlighted expressed cells → gradient fill
    geom_point(
      data = df %>% filter(expressed & !highlight),
      aes(x = x, y = y, fill = expr),
      shape = 21, color = "black", stroke = stroke.size,
      size = pt.size.highlight, alpha = 0.9
    ) +
    # Highlighted + expressed → expression fill + highlight outline
    geom_point(
      data = df %>% filter(expressed & highlight),
      aes(x = x, y = y, fill = expr),
      shape = 21, color = highlight_color, stroke = 0.5,
      size = pt.size.highlight, alpha = 0.9
    ) +
    scale_fill_gradient(low = low_color, high = high_color) +
    scale_y_reverse() +
    coord_fixed() +
    ggtitle(paste0(gene, " — ", paste(highlight_group, collapse = ", "))) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}