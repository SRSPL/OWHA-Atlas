```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
# Sets root directory for chunks to the directory where the file is located
knitr::opts_knit$set(root.dir = "/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Jonathan Current Analysis/Spatial")

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(dplyr)
library(tidyr)
library(purrr)
library(writexl)
library(SeuratWrappers)

library(hdf5r)

library(sfarrow)
library(Rfast2)
library(rlang)
#install.packages("ggsignif")  # only once
library(ggsignif)
library(colorway)
library(khroma)
library(ggbreak)


#font_import() # Needed to import fonts from computer 
#Load needed functions
source("~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/scRNA SEQ of wounds/Skin Wound Healing Atlas/Useful Functions/GGsave_multidevice.R")

# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e9)

#Save a custom theme for ggplots, including setting Arial as default font and
#making background transparent
custom_theme <- theme(
  text = element_text(family = "Arial", face = "bold"),
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  legend.background = element_rect(fill = "transparent")
)
```
Load in the datasets
```{r Load data}
object_uw <- readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Shared Data/Spatial Data/RDS Files/YWV01_UW_FineAnnotation_2025-08-05.rds")

object_d4 <- readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Shared Data/Spatial Data/RDS Files/YWV04_D4PW_FineAnnotation_2025-08-05.rds")

object_d7 <- readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Shared Data/Spatial Data/RDS Files/YWV03_D7PW_FineAnnotation_2025-08-05.rds")

object_d30 <- readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Shared Data/Spatial Data/RDS Files/YWV02_D30PW_FineAnnotation_2025-08-05.rds")
```
For ease of plotting, generate a list of objects
```{r Generate dataset list for plotting}
# List of all conditions with metadata
plot_info <- list(
  "object_uw"   = list(object = object_uw,  regions = c("Region1", "Region2"), pt.size.factor = 1.2),
  "object_d4" = list(object = object_d4, regions = c("Region1", "Region2"),   pt.size.factor = 1.25),
  "object_d7" = list(object = object_d7, regions = c("Region1", "Region2"),   pt.size.factor = 1.25),
  "object_d30"= list(object = object_d30, regions = c("Region1", "Region2"), pt.size.factor = 1)
)
```

```{r Module Scoring To Identify Populations of Intrest Basal IV/Prolif Endo}

library(Seurat)
library(patchwork)

# === Parameters to update ===
celltype_name <- "ProlifEndo"  # e.g. "Basal_IV", "ProlifEndo", etc.'


celltype_name <- "Basal IV"



# Define your gene module
genes_of_interest <- c(
  "Cdh13", "Antxr1", "Col4a6", "P3h2", "Sox6", "Sema5a","Frmd4b", "Cblb", "Sema3c","Trp63")

#Markes for Basal IV
genes_of_interest <- c("Sox6", "Cdh13", "Nfia", "Mef2a", "Trp63", "Fgfr2", "Slc24a3")

genes_of_interest <- c(
  # Proliferation / mitosis-related markers
  "Kif4", "Kif15", "Kif11", "Cenpe", "Tpx2", "Neil3", 
  "Anln", "Knl1", "Kif20b", "Ect2", "Top2a", "Mki67", "Igfbp3")






score_name <- paste0(celltype_name, "_Score")
max_value <- 0.7
image_region <- "slice1"
pt_size <- 1.5

# === Add module score ===
object_uw  <- AddModuleScore(object_uw, features = list(genes_of_interest), name = score_name, assay = "Spatial")
object_d4  <- AddModuleScore(object_d4, features = list(genes_of_interest), name = score_name, assay = "Spatial")
object_d7  <- AddModuleScore(object_d7, features = list(genes_of_interest), name = score_name, assay = "Spatial")
object_d30 <- AddModuleScore(object_d30, features = list(genes_of_interest), name = score_name, assay = "Spatial")

score_column <- paste0(score_name, "1")

# === Generate spatial plots ===
p_uw <- SpatialFeaturePlot(
  object_uw,
  features = score_column,
  images = image_region,
  pt.size.factor = pt_size,
  min.cutoff = 0,
  max.cutoff = max_value
) + ggtitle(paste(celltype_name, "- UW")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

p_d4 <- SpatialFeaturePlot(
  object_d4,
  features = score_column,
  images = image_region,
  pt.size.factor = pt_size,
  min.cutoff = 0,
  max.cutoff = max_value
) + ggtitle(paste(celltype_name, "- D4")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

p_d7 <- SpatialFeaturePlot(
  object_d7,
  features = score_column,
  images = image_region,
  pt.size.factor = pt_size,
  min.cutoff = 0,
  max.cutoff = max_value
) + ggtitle(paste(celltype_name, "- D7")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

p_d30 <- SpatialFeaturePlot(
  object_d30,
  features = score_column,
  images = image_region,
  pt.size.factor = pt_size,
  min.cutoff = 0,
  max.cutoff = max_value
) + ggtitle(paste(celltype_name, "- D30")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))

# === Combine and save plots ===
combined_plot <- wrap_plots(p_uw, p_d4, p_d7, p_d30, ncol = 2)

dir.create("spatial_module_score_plots", showWarnings = FALSE)
ggsave(
  filename = paste0("spatial_module_score_plots/", score_name, "_combined.png"),
  plot = combined_plot,
  width = 14, height = 10, dpi = 320
)

# === Optional violin plot of module score by celltype ===
VlnPlot(object_d7, features = score_column, group.by = "first_type") + NoLegend()
```

```{r Marking Cell Lables Using Module Scores}
library(Seurat)
library(patchwork)

# Parameters
basal_score_col <- "Basal IV_Score1"
prolif_score_col <- "ProlifEndo_Score1"

basal_score_threshold <- 0.2
prolif_score_threshold <- 0.2

basal_marker_genes <- c("Itga6", "Krt14")
prolif_marker_gene <- c("Pecam1")

# Function to label cells for one Seurat object
label_cells <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "Spatial"  # use RNA assay or adjust as needed
  
  module_scores <- seurat_obj@meta.data[, c(basal_score_col, prolif_score_col), drop = FALSE]
  
  # Average basal marker expression per cell
  basal_marker_expr <- Matrix::colMeans(GetAssayData(seurat_obj, slot = "data")[basal_marker_genes, , drop = FALSE])
  prolif_marker_expr <- GetAssayData(seurat_obj, slot = "data")[prolif_marker_gene, ]
  
  cell_labels <- rep("Other", ncol(seurat_obj))
  
  basal_cells <- which(
    module_scores[, basal_score_col] >= basal_score_threshold &
      basal_marker_expr >= 1
  )
  prolif_cells <- which(
    module_scores[, prolif_score_col] >= prolif_score_threshold &
      prolif_marker_expr >= 1
  )
  
  cell_labels[basal_cells] <- "BasalIV"
  cell_labels[prolif_cells] <- "ProlifEndo"
  
  seurat_obj$cell_label <- cell_labels
  return(seurat_obj)
}


object_uw  <- label_cells(object_uw)
object_d4  <- label_cells(object_d4)
object_d7  <- label_cells(object_d7)
object_d30 <- label_cells(object_d30)

# Define your regions (update this based on your object)
regions <- c("Region1","Region2")  # replace or extend as needed

#Make Spatial Plots of Cell Locations
# Define your color map
label_colors <- c(
  "Other"       = "white",
  "BasalIV"    = "#364B9A",  # blue
  "ProlifEndo" = "#F6C800"   # gold
)

# Set plot parameters
pt.size <- 1.5  # adjust as needed
info_name <- "cell_label"  # or another identifier if you'd like

# Loop over each region and each object
for (region in regions) {
  plots <- list()
  
  for (obj_name in c("object_uw", "object_d4", "object_d7", "object_d30")) {
    obj <- get(obj_name)
    
    p <- SpatialDimPlot(
      obj,
      group.by = "cell_label",
      images = region,
      pt.size.factor = pt.size,
      cols = label_colors,
      label = FALSE,
      stroke = 0.05
    ) +
      ggtitle(paste(obj_name, region, info_name)) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) +
      NoLegend()
    
    plots[[obj_name]] <- p
  }
  
  
  
  # Combine and save
  combined <- wrap_plots(plots, ncol = 2)
  ggsave(
    filename = paste0("spatial_module_score_plots/", region, "BasalIVProlif_labelled_combined.png"),
    plot = combined,
    width = 12,
    height = 10,
    dpi = 320
  )
}



# Define your regions (update this based on your object)
regions <- c("Region1","Region2")  # replace or extend as needed

#Make Spatial Plots of Cell Locations
# Define your color map
label_colors <- c(
  "Other"       = "white",
  "BasalIV"    = "#364B9A",  # blue
  "ProlifEndo" = "#F6C800"   # gold
)

# Set plot parameters
pt.size <- 1.5  # adjust as needed
info_name <- "cell_label"  # or another identifier if you'd like

immune_types <- c("BasalIV","ProlifEndo")

# # Define palette for immune types
# pal_fun <- color("smoothrainbow")
# pal <- pal_fun(length(immune_types), range = c(0.1, 0.9))
# names(pal) <- immune_types
# pal <- c(pal, Unselected = "white")


CustomSpatialHighlightPlot <- function(
    seurat_obj,
    region,
    cells_highlight,
    pt.size.highlight = 1.2,
    pt.size.unselected = 0.3,
    stroke.size = 0.1,
    label_colors,
    custom_theme
) {
  Idents(seurat_obj) <- "cell_label"
  highlight_cells <- unlist(CellsByIdentities(seurat_obj, idents = cells_highlight))
  coords <- GetTissueCoordinates(seurat_obj, image = region)
  # Flip y-axis to match SpatialDimPlot orientation
  tmp <- coords$x
  coords$x <- coords$y
  coords$y <- tmp
  barcodes <- rownames(coords)
  coords$highlight <- barcodes %in% highlight_cells
  coords$cell_type <- as.character(seurat_obj$cell_label[barcodes])
  coords$cell_type[!coords$cell_type %in% cells_highlight] <- "Other"
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
    scale_fill_manual(values = label_colors, drop = FALSE) +
    scale_y_reverse() +
    custom_theme +
    theme_void()
  return(p)
}

for (region in regions) {
  plots <- list()
  
  for (obj_name in c("object_uw", "object_d4", "object_d7", "object_d30")) {
    obj <- get(obj_name)
    
    # Set identity
    Idents(obj) <- "cell_label"
    
    # Check which labels are present in this object
    valid_immune_types <- intersect(c("BasalIV", "ProlifEndo"), levels(obj))
    
    if (length(valid_immune_types) == 0) {
      cat("Skipping", obj_name, "for", region, "- no valid immune types.\n")
      next
    }
    
    # Use your custom plot function
    p <- CustomSpatialHighlightPlot(
      seurat_obj = obj,
      region = region,
      cells_highlight = valid_immune_types,
      pt.size.highlight = 0.7,
      pt.size.unselected = 0.5,
      stroke.size = 0.01,
      label_colors,
      custom_theme = custom_theme
    ) +
      ggtitle(paste(obj_name, region, "BasalIV & Prolf.Endo")) +
      custom_theme +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) +
      NoLegend()
    
    plots[[obj_name]] <- p
  }
  
  # Only combine and save if any plots were generated
  if (length(plots) > 0) {
    combined <- wrap_plots(plots, ncol = 2)
    
    ggsave(
      filename = paste0("spatial_module_score_plots/", region, "_BasalIVProlifEndo_labelled_combined.png"),
      plot = combined,
      width = 12,
      height = 10,
      dpi = 320
    )
  } else {
    cat("No valid plots for region", region, "\n")
  }
}


# Loop over each region and each object
for (region in regions) {
  plots <- list()
  
  for (obj_name in c("object_uw", "object_d4", "object_d7", "object_d30")) {
    obj <- get(obj_name)
    
    p <- SpatialDimPlot(
      obj,
      group.by = "cell_label",
      images = region,
      pt.size.factor = pt.size ,
      cols = label_colors,
      label = FALSE,
      stroke = 0.05
    ) +
      ggtitle(paste(obj_name, region, info_name)) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) 
    plots[[obj_name]] <- p
  }
  
  
  
  # Combine and save
  combined <- wrap_plots(plots, ncol = 2)
  ggsave(
    filename = paste0("spatial_module_score_plots/", region, "BasalIVProlifEndo_labelled_combined.png"),
    plot = combined,
    width = 12,
    height = 10,
    dpi = 320
  )
}

```

#font_import() # Needed to import fonts from computer 
#Load needed functions
source("~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/scRNA SEQ of wounds/Skin Wound Healing Atlas/Useful Functions/GGsave_multidevice.R")

# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e9)

#Save a custom theme for ggplots, including setting Arial as default font and
#making background transparent
custom_theme <- theme(
  text = element_text(family = "Arial", face = "bold"),
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  legend.background = element_rect(fill = "transparent")
)
```
Load in the datasets
```{r Load data}
object_uw <- readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Shared Data/Spatial Data/RDS Files/YWV01_UW_FineAnnotation_2025-08-05.rds")

object_d4 <- readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Shared Data/Spatial Data/RDS Files/YWV04_D4PW_FineAnnotation_2025-08-05.rds")

object_d7 <- readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Shared Data/Spatial Data/RDS Files/YWV03_D7PW_FineAnnotation_2025-08-05.rds")

object_d30 <- readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Shared Data/Spatial Data/RDS Files/YWV02_D30PW_FineAnnotation_2025-08-05.rds")
```


```{r}
# List of all conditions with metadata
plot_info <- list(
  "object_uw"   = list(object = object_uw,  regions = c("Region1", "Region2"), pt.size.factor = 1.2),
  "object_d4" = list(object = object_d4, regions = c("Region1", "Region2"),   pt.size.factor = 1.25),
  "object_d7" = list(object = object_d7, regions = c("Region1", "Region2"),   pt.size.factor = 1.25),
  "object_d30"= list(object = object_d30, regions = c("Region1", "Region2"), pt.size.factor = 1)
)
```

```{r}
# ============================================================
# CustomSpatialFeatureCoexpPlot (fixed gradient version)
# ============================================================

CustomSpatialFeatureCoexpPlot <- function(
    object,
    image,
    genes,
    gene_colors,
    assay         = NULL,
    layer         = "data",
    image.scale   = "lowres",
    pt.size.expr  = 1.5,   # size for expressing spots
    pt.size.zero  = 0.5,   # size for zero-expression spots
    stroke        = 0.1,
    bg.col        = "white" # color for zero-expression spots
) {
  if (!is.null(assay)) DefaultAssay(object) <- assay
  
  # get image + raster
  if (!(image %in% names(object@images)))
    stop("Image '", image, "' not found in Seurat object")
  spat_img <- object@images[[image]]
  img_rast <- GetImage(spat_img, mode = "raster")
  
  # get spot coordinates
  coords <- as.data.frame(GetTissueCoordinates(spat_img, scale = image.scale))
  coords$barcode <- rownames(coords)
  
  # flip coordinates for correct orientation
  tmp <- coords$x
  coords$x <- coords$y
  coords$y <- tmp
  
  # fetch expression data
  expr <- FetchData(object, vars = genes, cells = coords$barcode,
                    layer = layer, clean = FALSE)
  
  # scale expression between 0–1
  expr_scaled <- as.data.frame(lapply(expr, scales::rescale))
  
  # ===========================================
  # Handle color assignment
  # ===========================================
  if (length(genes) == 1) {
    # ---- Single gene: gradient from white → color ----
    expr_vals <- expr[[1]]
    expr_scaled_single <- scales::rescale(expr_vals)
    
    spot_cols <- scales::col_numeric(
      palette = c("white", gene_colors[1]),
      domain = range(expr_scaled_single, na.rm = TRUE)
    )(expr_scaled_single)
    
    names(spot_cols) <- rownames(expr)
    
  } else {
    # ---- Multi-gene: blend colors ----
    rgb_mat <- t(grDevices::col2rgb(gene_colors)/255)
    wts <- expr_scaled / rowSums(expr_scaled + 1e-6)
    blend <- as.matrix(wts) %*% rgb_mat
    spot_cols <- grDevices::rgb(blend[,1], blend[,2], blend[,3])
    names(spot_cols) <- rownames(expr)
  }
  
  # mark zero-expression spots
  zero_bcs <- rownames(expr)[rowSums(expr) == 0]
  spot_cols[zero_bcs] <- bg.col
  coords$col <- spot_cols[coords$barcode]
  
  # split coords
  coords_zero <- coords[coords$barcode %in% zero_bcs, ]
  coords_expr <- coords[!coords$barcode %in% zero_bcs, ]
  
  # plot limits
  xlim <- range(coords$x); ylim <- range(coords$y)
  
  # --- Main spatial plot ---
  spatial_plot <- ggplot() +
    annotation_raster(
      img_rast,
      xmin = xlim[1], xmax = xlim[2],
      ymin = ylim[1], ymax = ylim[2],
      interpolate = TRUE
    ) +
    geom_point(
      data = coords_zero,
      aes(x = x, y = y),
      color = "grey75",
      shape = 21,
      fill  = bg.col,
      size  = pt.size.zero,
      stroke= stroke
    ) +
    geom_point(
      data = coords_expr,
      aes(x = x, y = y, fill = col),
      shape = 21,
      size  = pt.size.expr,
      stroke= stroke
    ) +
    scale_fill_identity() +
    coord_fixed() +
    scale_y_reverse() +
    theme_void()
  
  # --- Legend strips per gene ---
  make_gene_legend <- function(gene, color) {
    df <- data.frame(x = 1:100, y = 1, val = seq(0, 1, length.out = 100))
    ggplot(df, aes(x = x, y = y, fill = val)) +
      geom_tile() +
      scale_fill_gradientn(colors = c("white", color), name = gene) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "italic"),
        legend.text  = element_text(size = 8),
        legend.key.width = unit(1.2, "cm"),
        plot.margin = margin(2, 2, 2, 2)
      ) +
      guides(fill = guide_colorbar(title.position = "top", barwidth = 5))
  }
  
  legend_plots <- lapply(seq_along(genes), function(i) {
    make_gene_legend(genes[i], gene_colors[i])
  })
  legends_combined <- plot_grid(plotlist = legend_plots, nrow = 1)
  
  # --- Combine main plot + legends ---
  final_plot <- plot_grid(spatial_plot, legends_combined,
                          ncol = 1, rel_heights = c(1, 0.15))
  
  return(final_plot)
}

# ============================================================
# Example usage (loop for single-gene gradient plots)
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

feats <- c("Tenm4", "Adgrl2")
feat_colors <- c("seagreen")  # single color; will apply gradient white→seagreen

for (gene in feats) {
  message("Processing: ", gene)
  
  p <- CustomSpatialFeatureCoexpPlot(
    object = object_d7,
    image = "Region1",
    genes = gene,
    gene_colors = feat_colors,
    assay = "Spatial",
    pt.size.expr = 1.5,
    pt.size.zero = 0,
    bg.col = "white",
    stroke = 0
  )
  
  ggsave(
    filename = paste0("D7PWRegion1_", gene, "_Expression.png"),
    plot = p, width = 6, height = 6, dpi = 300
  )
  
  ggsave(
    filename = paste0("D7PWRegion1_", gene, "_Expression.svg"),
    plot = p, width = 6, height = 6
  )
}

```

```{r}
# ===============================
# Load required libraries
# ===============================
library(Seurat)
library(patchwork)
library(dplyr)

# ===============================
# Define gene sets per cell type
# ===============================
gene_sets <- list(
  "BasalIV" = c("Sox6", "Cdh13", "Nfia", "Mef2a", "Trp63", "Fgfr2", "Slc24a3"),
  "PericyteI" = c("Dlk1", "Abcc9", "Plac8", "Steap4", "Agtr1a", "Adra2a", "Ifitm1", "Pla1a", "Ccl11", "Gpx3", "Rgs5"),
  "RepairSCII" = c("Morrbid", "Cfap299", "Col8a1", "Bche", "Apod", "Gpc3", "Col20a1", "Etv1", "Tmem158", "Mpz", "Sox10", "Plp1"),
  "ProliferatingECs" = c("Kif4", "Kif15", "Kif11", "Cenpe", "Tpx2", "Neil3", "Anln", "Knl1", "Kif20b", "Ect2", "Top2a", "Mki67", "Igfbp3")
)

# ===============================
# Parameters
# ===============================
image_region <- "slice1"
pt_size <- 1.2
max_value <- 0.7

# ===============================
# Add module scores for all objects
# ===============================
seurat_list <- list(
  UW = object_uw,
  D4 = object_d4,
  D7 = object_d7,
  D30 = object_d30
)

for (celltype in names(gene_sets)) {
  score_name <- paste0(celltype, "_Score")
  for (cond in names(seurat_list)) {
    seurat_list[[cond]] <- AddModuleScore(
      seurat_list[[cond]],
      features = list(gene_sets[[celltype]]),
      name = score_name,
      assay = "Spatial"
    )
  }
}

# ===============================
# Define labeling rules
# ===============================
label_params <- list(
  "BasalIV" = list(score_col = "BasalIV_Score1", markers = c("Itga6", "Krt14"), threshold = 0.25),
  "PericyteI" = list(score_col = "PericyteI_Score1", markers = c("Rgs5"), threshold = 0.35),
  "RepairSCII" = list(score_col = "RepairSCII_Score1", markers = c("Mpz"), threshold = 0.30),
  "ProliferatingECs" = list(score_col = "ProliferatingECs_Score1", markers = c("Pecam1"), threshold = 0.25)
)

# ===============================
# Labeling function
# ===============================
label_cells <- function(obj, label_params) {
  obj$cell_label <- "Other"
  
  for (ct in names(label_params)) {
    p <- label_params[[ct]]
    score_col <- p$score_col
    threshold <- p$threshold
    
    # Check if all marker genes are in data
    valid_markers <- intersect(p$markers, rownames(obj))
    if (length(valid_markers) == 0) next
    
    expr_df <- FetchData(obj, vars = c(score_col, valid_markers))
    high_cells <- rownames(expr_df)[
      expr_df[[score_col]] > threshold &
        rowMeans(expr_df[, valid_markers, drop = FALSE]) > 1
    ]
    
    obj$cell_label[high_cells] <- ct
  }
  
  return(obj)
}

# Apply labeling to all conditions
for (cond in names(seurat_list)) {
  seurat_list[[cond]] <- label_cells(seurat_list[[cond]], label_params)
}

# Unpack labeled objects
object_uw  <- seurat_list$UW
object_d4  <- seurat_list$D4
object_d7  <- seurat_list$D7
object_d30 <- seurat_list$D30

# ===============================
# Plot labeled cells and module scores
# ===============================

# Define plot function for a given condition
plot_spatial_labeling <- function(obj, label_col = "cell_label", region = image_region, pt_size = 1.2) {
  SpatialPlot(obj, group.by = label_col, images = region, pt.size.factor = pt_size) +
    ggtitle(paste0(unique(obj$orig.ident), " - Cell Labels"))
}

plot_module_score <- function(obj, feature, title, region = image_region, pt_size = 1.2, max_value = 0.7) {
  SpatialFeaturePlot(obj, features = feature, images = region, pt.size.factor = pt_size) +
    scale_fill_viridis_c(limits = c(0, max_value)) +
    ggtitle(title)
}

# === Example: plotting RepairSCII module across timepoints ===
plots_repair_sc <- wrap_plots(
  plot_module_score(object_uw,  "RepairSCII_Score1", "UW - RepairSCII Module"),
  plot_module_score(object_d4,  "RepairSCII_Score1", "D4 - RepairSCII Module"),
  plot_module_score(object_d7,  "RepairSCII_Score1", "D7 - RepairSCII Module"),
  plot_module_score(object_d30, "RepairSCII_Score1", "D30 - RepairSCII Module"),
  ncol = 4
)

# === Example: plotting labeled cells ===
plots_labels <- wrap_plots(
  plot_spatial_labeling(object_uw),
  plot_spatial_labeling(object_d4),
  plot_spatial_labeling(object_d7),
  plot_spatial_labeling(object_d30),
  ncol = 4
)

# View plots
plots_repair_sc
plots_labels

```

```{r}
# ===============================
# Load plotting libraries
# ===============================
library(tidyverse)
library(ggbeeswarm)
library(scCustomize)

# ===============================
# Define file paths and objects
# ===============================
files <- c(
  "YWV01_UW_BanksyDistributionFromWoundSite.csv",
  "YWV04_D4PW_BanksyDistributionFromWoundSite.csv",
  "YWV03_D7PW_BanksyDistributionFromWoundSite.csv",
  "YWV02_D30PW_BanksyDistributionFromWoundSite.csv"
)

objects <- list(
  "YWV01_UW_BanksyDistributionFromWoundSite.csv" = object_uw,
  "YWV04_D4PW_BanksyDistributionFromWoundSite.csv" = object_d4,
  "YWV03_D7PW_BanksyDistributionFromWoundSite.csv" = object_d7,
  "YWV02_D30PW_BanksyDistributionFromWoundSite.csv" = object_d30
)

base_path <- "/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Immunome Paper/Figure 3 Spatial/SpaceFold/"

# ===============================
# Combine data across timepoints
# ===============================
df_list <- list()
for (file_name in files) {
  df <- read.csv(file.path(base_path, file_name), row.names = "X")
  df$barcode <- rownames(df)
  seurat_obj <- objects[[file_name]]
  
  labels <- data.frame(barcode = colnames(seurat_obj), cell_label = seurat_obj$cell_label)
  df <- left_join(df, labels, by = "barcode")
  
  df <- df %>%
    filter(cell_label %in% c("BasalIV", "PericyteI", "RepairSCII", "ProliferatingECs")) %>%
    mutate(
      timepoint = case_when(
        grepl("UW", file_name) ~ "UW",
        grepl("D4PW", file_name) ~ "D4PW",
        grepl("D7PW", file_name) ~ "D7PW",
        grepl("D30PW", file_name) ~ "D30PW"
      )
    )
  df_list[[file_name]] <- df
}

combined_df <- bind_rows(df_list)
combined_df$cell_label <- factor(combined_df$cell_label,
                                 levels = c("BasalIV", "PericyteI", "RepairSCII", "ProliferatingECs")
)
combined_df$timepoint <- factor(combined_df$timepoint,
                                levels = c("UW", "D4PW", "D7PW", "D30PW")
)

# ===============================
# Define colors
# ===============================
cell_colors <- c(
  "BasalIV" = "#364B9A",
  "PericyteI" = "#D12E82",
  "RepairSCII" = "#8DD3C7",
  "ProliferatingECs" = "#F6C800"
)

# ===============================
# BeeSwarm Plot
# ===============================
p_bee <- ggplot(combined_df, aes(x = projection, y = cell_label, fill = cell_label)) +
  geom_quasirandom(
    groupOnY = TRUE, shape = 21, size = 1.8, stroke = 0.25,
    color = "black", alpha = 0.9
  ) +
  scale_fill_manual(values = cell_colors) +
  facet_wrap(~timepoint, scales = "free_x") +
  labs(x = "Position along tissue axis", y = "Cell population") +
  theme_ggprism_mod(base_size = 14) +
  theme(legend.position = "right")

ggsave("FourCellTypes_BeeSwarm.png", plot = p_bee, width = 10, height = 5, dpi = 300)

# ===============================
# Expression Histograms per Cell Type
# ===============================
genes_by_celltype <- list(
  "BasalIV" = c("Sema3c"),
  "PericyteI" = c("Nrp1","Nrp2"),
  "RepairSCII" = c("Ngfr","Nrp1"),
  "ProliferatingECs" = c("Mki67","Pecam1")
)

plots_hist <- list()
for (ct in names(genes_by_celltype)) {
  df_ct <- combined_df %>% filter(cell_label == ct)
  p <- ggplot(df_ct, aes(x = projection, fill = cell_label)) +
    geom_histogram(binwidth = 0.05, color = "black", alpha = 0.9) +
    scale_fill_manual(values = cell_colors) +
    labs(title = ct, x = "Position along tissue axis", y = "Cell count") +
    theme_ggprism_mod(base_size = 14)
  plots_hist[[ct]] <- p
  ggsave(paste0(ct, "_Projection_Histogram.png"), plot = p, width = 10, height = 6, dpi = 300)
}

# ===============================
# Spatial Plots per Timepoint
# ===============================
for (cond in names(seurat_list)) {
  p <- SpatialDimPlot(
    seurat_list[[cond]],
    group.by = "cell_label",
    pt.size.factor = 1.2,
    alpha = c(0.2, 1),
    image.alpha = 0.6,
    cols = cell_colors
  ) + ggtitle(paste0(cond, " — Labeled Cell Types"))
  ggsave(paste0(cond, "_Spatial_Labeling.png"), plot = p, width = 20, height = 20, dpi = 300)
}

```
```{r}
# ===============================
# Load required libraries
# ===============================
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(tidyverse)
library(scCustomize)

# ===============================
# Define gene sets per cell type
# ===============================
gene_sets <- list(
  "BasalIV" = c("Sox6", "Cdh13", "Nfia", "Mef2a", "Trp63", "Fgfr2", "Slc24a3"),
  "PericyteI" = c("Dlk1", "Abcc9", "Plac8", "Steap4", "Agtr1a", "Adra2a", "Ifitm1", "Pla1a", "Ccl11", "Gpx3", "Rgs5"),
  "RepairSCII" = c("Morrbid", "Cfap299", "Col8a1", "Bche", "Apod", "Gpc3", "Col20a1", "Etv1", "Tmem158", "Mpz", "Sox10", "Plp1"),
  "ProliferatingECs" = c("Kif4", "Kif15", "Kif11", "Cenpe", "Tpx2", "Neil3", "Anln", "Knl1", "Kif20b", "Ect2", "Top2a", "Mki67", "Igfbp3")
)

# ===============================
# Parameters
# ===============================
image_region <- "slice1"
pt_size <- 1.2
max_value <- 0.7

# ===============================
# Add module scores for all objects
# ===============================
seurat_list <- list(
  UW = object_uw,
  D4 = object_d4,
  D7 = object_d7,
  D30 = object_d30
)

for (celltype in names(gene_sets)) {
  score_name <- paste0(celltype, "_Score")
  for (cond in names(seurat_list)) {
    seurat_list[[cond]] <- AddModuleScore(
      seurat_list[[cond]],
      features = list(gene_sets[[celltype]]),
      name = score_name,
      assay = "Spatial"
    )
  }
}

# ===============================
# Define labeling rules
# ===============================
label_params <- list(
  "BasalIV" = list(score_col = "BasalIV_Score1", markers = c("Itga6", "Krt14"), threshold = 0.25),
  "PericyteI" = list(score_col = "PericyteI_Score1", markers = c("Rgs5"), threshold = 0.30),
  "RepairSCII" = list(score_col = "RepairSCII_Score1", markers = c("Mpz"), threshold = 0.40),
  "ProliferatingECs" = list(score_col = "ProliferatingECs_Score1", markers = c("Pecam1"), threshold = 0.25)
)

# ===============================
# Labeling function
# ===============================
label_cells <- function(obj, label_params) {
  obj$cell_label <- "Other"
  
  for (ct in names(label_params)) {
    p <- label_params[[ct]]
    score_col <- p$score_col
    threshold <- p$threshold
    valid_markers <- intersect(p$markers, rownames(obj))
    if (length(valid_markers) == 0) next
    
    expr_df <- FetchData(obj, vars = c(score_col, valid_markers))
    high_cells <- rownames(expr_df)[
      expr_df[[score_col]] > threshold &
        rowMeans(expr_df[, valid_markers, drop = FALSE]) > 1
    ]
    obj$cell_label[high_cells] <- ct
  }
  return(obj)
}

# Apply labeling
for (cond in names(seurat_list)) {
  seurat_list[[cond]] <- label_cells(seurat_list[[cond]], label_params)
}

# Unpack labeled objects
object_uw  <- seurat_list$UW
object_d4  <- seurat_list$D4
object_d7  <- seurat_list$D7
object_d30 <- seurat_list$D30

# ===============================
# Spatial Plots
# ===============================
plot_spatial_labeling <- function(obj, label_col = "cell_label") {
  SpatialPlot(obj, group.by = label_col, images = image_region, pt.size.factor = pt_size) +
    ggtitle(paste0(unique(obj$orig.ident), " - Cell Labels"))
}

plot_module_score <- function(obj, feature, title) {
  SpatialFeaturePlot(obj, features = feature, images = image_region, pt.size.factor = pt_size) +
    scale_fill_viridis_c(limits = c(0, max_value)) +
    ggtitle(title)
}

plots_labels <- wrap_plots(
  plot_spatial_labeling(object_uw),
  plot_spatial_labeling(object_d4),
  plot_spatial_labeling(object_d7),
  plot_spatial_labeling(object_d30),
  ncol = 4
)

plots_repair_sc <- wrap_plots(
  plot_module_score(object_uw, "RepairSCII_Score1", "UW - RepairSCII Module"),
  plot_module_score(object_d4, "RepairSCII_Score1", "D4 - RepairSCII Module"),
  plot_module_score(object_d7, "RepairSCII_Score1", "D7 - RepairSCII Module"),
  plot_module_score(object_d30, "RepairSCII_Score1", "D30 - RepairSCII Module"),
  ncol = 4
)

# ===============================
# Beeswarm (Projection) Plots
# ===============================
# combined_df must contain columns: barcode, projection, timepoint, cell_label

```
```{r}


# =========================
# Load packages
# =========================
library(tidyverse)
library(ggbeeswarm)
library(scCustomize)

# =========================
# Define files and Seurat objects
# =========================
files <- c(
  "YWV01_UW_BanksyDistributionFromWoundSite.csv",
  "YWV04_D4PW_BanksyDistributionFromWoundSite.csv",
  "YWV03_D7PW_BanksyDistributionFromWoundSite.csv",
  "YWV02_D30PW_BanksyDistributionFromWoundSite.csv"
)

objects <- list(
  "YWV01_UW_BanksyDistributionFromWoundSite.csv" = object_uw,
  "YWV04_D4PW_BanksyDistributionFromWoundSite.csv" = object_d4,
  "YWV03_D7PW_BanksyDistributionFromWoundSite.csv" = object_d7,
  "YWV02_D30PW_BanksyDistributionFromWoundSite.csv" = object_d30
)

base_path <- "/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Immunome Paper/Figure 3 Spatial/SpaceFold/"

# =========================
# Collect data across files
# =========================
df_list <- list()
for (file_name in files) {
  
  # Read CSV
  df <- read.csv(file.path(base_path, file_name), row.names = "X")
  df$barcode <- rownames(df)
  
  # Get Seurat object
  seurat_obj <- objects[[file_name]]
  
  # Pull cell labels
  labels <- data.frame(
    barcode = colnames(seurat_obj),
    cell_label = seurat_obj$cell_label
  )
  
  # Merge
  df <- left_join(df, labels, by = "barcode")
  
  # Filter for the four cell types
  df <- df %>% filter(cell_label %in% c("BasalIV", "PericyteI", "RepairSCII", "ProliferatingECs"))
  
  # Add timepoint
  df$timepoint <- case_when(
    grepl("UW", file_name) ~ "UW",
    grepl("D4PW", file_name) ~ "D4PW",
    grepl("D7PW", file_name) ~ "D7PW",
    grepl("D30PW", file_name) ~ "D30PW",
    TRUE ~ "Other"
  )
  
  # Store
  df_list[[file_name]] <- df
}

# Combine all timepoints
combined_df <- bind_rows(df_list)

# =========================
# Factor levels
# =========================
combined_df$cell_label <- factor(
  combined_df$cell_label,
  levels = c("BasalIV", "PericyteI", "RepairSCII", "ProliferatingECs")
)
combined_df$timepoint <- factor(
  combined_df$timepoint,
  levels = c("UW", "D4PW", "D7PW", "D30PW")
)

# =========================
# Define color palette
# =========================
cell_colors <- c(
  "BasalIV" = "#364B9A",
  "PericyteI" = "#D12E82",
  "RepairSCII" = "#8DD3C7",
  "ProliferatingECs" = "#F6C800"
)

# =========================
# Faceted plot by timepoint
# =========================
p_faceted <- ggplot(combined_df, aes(x = projection, y = cell_label, fill = cell_label)) +
  geom_quasirandom(
    groupOnY = TRUE,
    shape = 21,
    size = 1.8,
    stroke = 0.25,
    color = "black",
    alpha = 0.9
  ) +
  scale_fill_manual(values = cell_colors) +
  facet_wrap(~timepoint, scales = "free_x") +
  labs(x = "Position along tissue axis", y = "Cell population") +
  theme_ggprism_mod(base_size = 14) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 12)
  )

# Save
ggsave("FourCellTypes_Faceted.svg", plot = p_faceted, width = 10, height = 5, dpi = 300)
ggsave("FourCellTypes_Faceted.png", plot = p_faceted, width = 10, height = 5, dpi = 300)

# =========================
# Facet in one column
# =========================
p_faceted_col <- ggplot(combined_df, aes(x = projection, y = cell_label, fill = cell_label)) +
  geom_quasirandom(
    groupOnY = TRUE,
    shape = 21,
    size = 1.8,
    stroke = 0.25,
    color = "black",
    alpha = 0.9
  ) +
  scale_fill_manual(values = cell_colors) +
  facet_wrap(~timepoint, ncol = 1, scales = "free_x") +
  labs(x = "Position along tissue axis", y = "Cell population") +
  theme_ggprism_mod(base_size = 14) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 12)
  )

# Save
ggsave("FourCellTypes_FacetedColumn.svg", plot = p_faceted_col, width = 8, height = 10, dpi = 300)
ggsave("FourCellTypes_FacetedColumn.png", plot = p_faceted_col, width = 8, height = 10, dpi = 300)


```


