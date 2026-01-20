

---
  title: "Spatial Figure Making"
author: "Jonathan Chin Cheong"
date: "`r Sys.Date()`"
output: html_document
---
  
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
```{r Quantifying Spatail Data By Metaclusters}
# Define Metacluster Colors list
metacluster_colors <- c(
  "IFE Keratinocytes" = "#364B9A",
  "HF Keratinocytes"  = "#56B4E9",
  "Sebocytes"         = "#8DD3C7",
  "Fibroblasts"       = "indianred2",
  "Immune cells"      = "#198F1B",
  "Endothelial cells" = "#F6C800",
  "Pericytes"         = "#D12E82",
  "Muscle cells"      = "orange2",
  "Schwann cells"     = "purple",
  "Melanocytes"       = "pink1",
  "Adipocytes"        = "deeppink",
  "Red Blood Cells"   = "#A7222B",
  "Unknown" = "grey50"
)

# Loop through all Seurat objects
for (cond in names(plot_info)) {
  
  seurat_obj <- plot_info[[cond]]$object
  
  # Make sure Idents is set to alignment_clusters
  Idents(seurat_obj) <- "banksy_cluster_0.5"
  
  # Pull metadata
  meta_df <- seurat_obj@meta.data %>%
    dplyr::count(banksy_cluster_0.5, metaclusters) %>%
    group_by(banksy_cluster_0.5) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
  # Stacked bar plot: alignment_clusters (y-axis), stacked by metaclusters
  p <- ggplot(meta_df, aes(y = factor(banksy_cluster_0.5), x = proportion, fill = metaclusters)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = metacluster_colors) +
    theme_ggprism_mod() +
    labs(
      y = "Alignment Cluster",
      x = "Proportion of Metaclusters",
      fill = "Metacluster"
    ) +
    theme(
      strip.text = element_text(size = 18, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),  # Cluster labels
      axis.text.x = element_text(size = 14, color = "black"),
      axis.title.x = element_text(size = 16, color = "black"),
      axis.title.y = element_text(size = 16, color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    ) +
    custom_theme
  
  print(p)
  
  # Save plot as SVG
  ggsave(
    filename = paste0(cond, "_AlignmentCluster_MetaclusterComposition.svg"),
    plot = p,
    width = 8, height = 10, dpi = 320
  )
}


```


```{r}

```{r Quantifying Spatail Data By Metaclusters}
# Define Metacluster Colors list
metacluster_colors <- c(
  "IFE Keratinocytes" = "#364B9A",
  "HF Keratinocytes"  = "#56B4E9",
  "Sebocytes"         = "#8DD3C7",
  "Fibroblasts"       = "indianred2",
  "Immune cells"      = "#198F1B",
  "Endothelial cells" = "#F6C800",
  "Pericytes"         = "#D12E82",
  "Muscle cells"      = "orange2",
  "Schwann cells"     = "purple",
  "Melanocytes"       = "pink1",
  "Adipocytes"        = "deeppink",
  "Red Blood Cells"   = "#A7222B",
  "Unknown" = "grey50"
)

# Loop through all Seurat objects
for (cond in names(plot_info)) {
  
  seurat_obj <- plot_info[[cond]]$object
  
  # Make sure Idents is set to alignment_clusters
  Idents(seurat_obj) <- "banksy_cluster_0.5"
  
  # Pull metadata
  meta_df <- seurat_obj@meta.data %>%
    dplyr::count(banksy_cluster_0.5, metaclusters) %>%
    group_by(banksy_cluster_0.5) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
  # Stacked bar plot: alignment_clusters (y-axis), stacked by metaclusters
  p <- ggplot(meta_df, aes(y = factor(banksy_cluster_0.5), x = proportion, fill = metaclusters)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = metacluster_colors) +
    theme_ggprism_mod() +
    labs(
      y = "Alignment Cluster",
      x = "Proportion of Metaclusters",
      fill = "Metacluster"
    ) +
    theme(
      strip.text = element_text(size = 18, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),  # Cluster labels
      axis.text.x = element_text(size = 14, color = "black"),
      axis.title.x = element_text(size = 16, color = "black"),
      axis.title.y = element_text(size = 16, color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    ) +
    custom_theme
  
  print(p)
  
  # Save plot as SVG
  ggsave(
    filename = paste0(cond, "_AlignmentCluster_MetaclusterComposition.svg"),
    plot = p,
    width = 8, height = 10, dpi = 320
  )
}


```

```{r Bar Plots shwoing Spatial Makeup of Cell Types by Banksy Clusters}
# Define the cell types you want to plot
celltypes_of_interest <- c("Fascia", "Myofibroblasts II", "Reticular III"
)

celltypes_of_interest <- c("Spinous III", "Pericyte I", "Proliferating ECs"
)


# Loop through all Seurat objects
for (cond in names(plot_info)) {
  
  seurat_obj <- plot_info[[cond]]$object
  
  # Make sure Idents are set to banksy clusters
  Idents(seurat_obj) <- "banksy_cluster_0.5"
  
  # Build palette for this object based on the clusters it actually has
  cluster_levels <- levels(seurat_obj$banksy_cluster_0.5)
  full_pal <- setNames(
    viridis(length(cluster_levels), option = "H"),
    cluster_levels
  )
  
  # Pull metadata, focusing only on your selected cell types
  meta_df <- seurat_obj@meta.data %>%
    filter(first_type %in% celltypes_of_interest) %>%
    dplyr::count(first_type, banksy_cluster_0.5) %>%
    group_by(first_type) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
  # Force ordering of y-axis based on celltypes_of_interest
  meta_df$first_type <- factor(meta_df$first_type, levels = rev(celltypes_of_interest))
  
  # Stacked proportional bar plot (keeps order, bold labels)
  p <- ggplot(meta_df, aes(
    y = first_type,
    x = proportion,
    fill = factor(banksy_cluster_0.5)
  )) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = full_pal) +
    theme_ggprism_mod() +
    labs(
      y = "Cell Type",
      x = "Proportion of Banksy Clusters",
      fill = "Banksy Cluster"
    ) +
    theme(
      strip.text = element_text(size = 18, color = "black"),
      axis.text.y = element_text(size = 20, color = "black", face = "bold"), # bold groups
      axis.text.x = element_text(size = 14, color = "black"),
      axis.title.x = element_text(size = 16, color = "black"),
      axis.title.y = element_text(size = 16, color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    ) +
    custom_theme
  
  print(p)
  
  # Save as SVG
  ggsave(
    filename = paste0(cond, "_SelectedCellTypes_BanksyClusterCompositionSpatialMajorWithLegendsD30PW.svg"),
    plot = p,
    width = 10, height = 2, dpi = 320
  )
}

```
```{r}
# Define the cell types you want to plot
celltypes_of_interest <- c("Adipocyte II","Fascia", "Muscle Progenitors","Pericyte I",
                           "Proliferating ECs", "Repair SC II"
                           
)

# Loop through all Seurat objects
for (cond in names(plot_info)) {
  
  seurat_obj <- plot_info[[cond]]$object
  
  # Make sure Idents are set to banksy clusters
  Idents(seurat_obj) <- "banksy_cluster_0.5"
  
  # Build palette for this object based on the clusters it actually has
  cluster_levels <- levels(seurat_obj$banksy_cluster_0.5)
  full_pal <- setNames(
    viridis(length(cluster_levels), option = "H"),
    cluster_levels
  )
  
  # Pull metadata, focusing only on your selected cell types
  meta_df <- seurat_obj@meta.data %>%
    filter(first_type %in% celltypes_of_interest) %>%
    dplyr::count(first_type, banksy_cluster_0.5) %>%
    group_by(first_type) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
  # Force ordering of y-axis based on celltypes_of_interest
  meta_df$first_type <- factor(meta_df$first_type, levels = rev(celltypes_of_interest))
  
  # Stacked proportional bar plot (keeps order, bold labels)
  p <- ggplot(meta_df, aes(
    y = first_type,
    x = proportion,
    fill = factor(banksy_cluster_0.5)
  )) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = full_pal) +
    theme_ggprism_mod() +
    labs(
      y = "Cell Type",
      x = "Proportion of Banksy Clusters",
      fill = "Banksy Cluster"
    ) +
    theme(
      strip.text = element_text(size = 18, color = "black"),
      axis.text.y = element_text(size = 14, color = "black", face = "bold"), # bold groups
      axis.text.x = element_text(size = 14, color = "black"),
      axis.title.x = element_text(size = 16, color = "black"),
      axis.title.y = element_text(size = 16, color = "black"),
      panel.border = element_rect(fill = NA, color = "black"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    ) +
    custom_theme
  
  print(p)
  
  # Save as SVG
  ggsave(
    filename = paste0(cond, "_SelectedCellTypes_BanksyClusterCompositionSpatialMinorControllersLegends.svg"),
    plot = p,
    width = 10, height = 8, dpi = 320
  )
}
```

```{r}
# Load packages
library(tidyverse)
library(ggbeeswarm)   # for geom_quasirandom
library(viridis)      # for viridis palette
library(scCustomize)
# Import the csv

read
df <- read.csv("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Immunome Paper/Figure 3 Spatial/SpaceFold/YWV03_D7PW_BanksyDistributionFromWoundSite.csv", row.names = "X")

# Make celltype a factor so colors are consistent
df$celltype <- factor(df$celltype)

# Number of unique cell types
n_celltypes <- length(unique(df$celltype))

# Viridis color palette
celltype_colors <- viridis(n_celltypes, option = "H")

# Plot
p <- ggplot(df, aes(x = projection, y = celltype, fill = celltype)) +
  geom_quasirandom(
    groupOnY = TRUE, 
    alpha = 1, 
    size = 1.8, 
    shape = 21,        # filled circle with outline
    color = "black",   # black outline
    stroke = 0.3       # thin outline
  ) +
  scale_fill_manual(values = celltype_colors) +
  labs(y = "Cluster", x = "Position along tissue length") +
  theme_ggprism_mod(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 10)
  )

print(p)

# Save as SVG
ggsave("D4PWBanksyCluster_axis_plot.svg", plot = p,
       width = 3, height = 5, units = "in", dpi = 300)

# Save as PNG
ggsave("D4PWBanksyClusters_axis_plot.png", plot = p,
       width = 3, height = 5, units = "in", dpi = 300)
```

```{r}
# Load packages
library(tidyverse)
library(ggbeeswarm)   # for geom_quasirandom
library(viridis)      # for viridis palette
library(scCustomize)

# Import the csv
df <- read.csv(
  "/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Immunome Paper/Figure 3 Spatial/SpaceFold/YWV03_D7PW_BanksyDistributionFromWoundSite.csv",
  row.names = "X"
)

# Make celltype a factor with all levels (for consistent color mapping)
df$celltype <- factor(df$celltype)

# Number of unique cell types
n_celltypes <- length(unique(df$celltype))

# Viridis color palette (same across all annotations)
celltype_colors <- viridis(n_celltypes, option = "H")
names(celltype_colors) <- levels(df$celltype)

# ---- Filter only for the selected clusters ----
df_sub <- df %>% filter(celltype %in% c("28","23","22","9"))

# ---- Normalize projection to 0-1 ----
df_sub <- df_sub %>%
  mutate(projection_scaled = (projection - min(projection)) / (max(projection) - min(projection)))

# ---- Plot ----
p <- ggplot(df_sub, aes(x = projection_scaled, y = celltype, fill = celltype)) +
  geom_quasirandom(
    groupOnY = TRUE,
    alpha = 1,
    size = 1.8,
    shape = 21,        # filled circle with outline
    color = "black",   # black outline
    stroke = 0.3       # thin outline
  ) +
  scale_fill_manual(values = celltype_colors) +  # keep full color scheme
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) + # show 0–1
  labs(y = "Cluster", x = "Position along tissue length (scaled 0–1)") +
  theme_ggprism_mod(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 10)
  )

print(p)

# Save as SVG
ggsave("D7PW_BanksyClusters_axis_plot_subsetLong.svg", plot = p,
       width = 7, height = 2, units = "in", dpi = 300)

# Save as PNG
ggsave("D7PW_BanksyClusters_axis_plot_subsetLong.png", plot = p,
       width = 7, height = 2, units = "in", dpi = 300)

```






library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

results_list <- list()

for (cond in names(plot_info)) {
  
  seurat_obj <- plot_info[[cond]]$object
  Idents(seurat_obj) <- "banksy_cluster_0.5"
  
  meta_df <- object_d4@meta.data %>%
    select(banksy_cluster_0.5, metaclusters, Region) %>%
    filter(!is.na(Region),
           !is.na(metaclusters),
           !is.na(banksy_cluster_0.5))
  
  # Count: alignment cluster × metacluster × Region
  counts_df <- meta_df %>%
    group_by(banksy_cluster_0.5, metaclusters, Region) %>%
    summarize(n = n(), .groups = "drop")
  
  # Pivot wide: rows = Region, columns = metaclusters
  wide_df <- counts_df %>%
    pivot_wider(
      names_from = metaclusters,
      values_from = n,
      values_fill = list(n = 0)
    )
  
  # Split by alignment cluster
  test_results <- wide_df %>%
    group_split(banksy_cluster_0.5) %>%
    map_df(function(df_cluster) {
      
      cluster_id <- unique(df_cluster$banksy_cluster_0.5)
      
      # Build contingency table
      cont_table <- df_cluster %>%
        ungroup() %>%
        select(-banksy_cluster_0.5) %>%
        column_to_rownames("Region") %>%
        as.matrix()
      
      # ---- SKIP INVALID CASES ----
      # Skip if fewer than 2 Regions
      if (nrow(cont_table) < 2) {
        return(NULL)
      }
      
      # Skip if a Region has zero total cells in that alignment cluster
      # (meaning the cluster only appears in one replicate)
      region_totals <- rowSums(cont_table)
      if (any(region_totals == 0)) {
        return(NULL)
      }
      
      # Force numeric + finite + >=0
      cont_table[is.na(cont_table)] <- 0
      cont_table[!is.finite(cont_table)] <- 0
      
      # ---- Fisher test only ----
      test <- fisher.test(cont_table, simulate.p.value = TRUE, B = 1e5)
      
      tibble(
        alignment_cluster = cluster_id,
        p_value = test$p.value,
        method = "Fisher’s exact"
      )
    })
  
  results_list[[cond]] <- test_results
}

# Combine all
all_stats <- bind_rows(results_list, .id = "condition")

# Optional: FDR correction
all_stats$p_adj <- p.adjust(all_stats$p_value, method = "BH")

all_stats

# ---- Combine results ----
all_stats <- bind_rows(results_list, .id = "condition")

# ---- Add FDR (BH) correction ----
all_stats$p_adj <- p.adjust(all_stats$p_value, method = "BH")

# ---- Add significance stars ----
all_stats$signif <- cut(
  all_stats$p_adj,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "n.s")
)

# ---- Order columns cleanly ----
stats_table <- all_stats %>%
  select(
    condition,
    alignment_cluster,
    p_value,
    p_adj,
    signif,
    method
  ) %>%
  arrange(condition, alignment_cluster)

# ---- Save as CSV ----
write.csv(
  stats_table,
  file = "AlignmentCluster_Metacluster_FisherStats.csv",
  row.names = FALSE
)

# If you want also an Excel version:
# install.packages("openxlsx")
# library(openxlsx)
# write.xlsx(stats_table, "AlignmentCluster_Metacluster_FisherStats.xlsx")
```


```{r}
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# BANKSY clusters you want to test
clusters_of_interest <- c("28", "23", "22", "9")



# Pull relevant metadata
meta_df <- object_d7@meta.data %>%
  filter(
    banksy_cluster_0.5 %in% clusters_of_interest,
    !is.na(Region),
    !is.na(first_type)
  )

results_list <- list()

# ---- run per cluster ----
for (cl in clusters_of_interest) {
  
  df_cluster <- meta_df %>% filter(banksy_cluster_0.5 == cl)
  
  # Build Region × first_type contingency table
  cont_table <- df_cluster %>%
    count(Region, first_type) %>%
    tidyr::pivot_wider(
      names_from = first_type,
      values_from = n,
      values_fill = 0
    ) %>%
    column_to_rownames("Region") %>%
    as.matrix()
  
  # Must have 2 regions (otherwise invalid)
  if (nrow(cont_table) < 2) next
  
  # Skip clusters that only appear in 1 region
  if (any(rowSums(cont_table) == 0)) next
  
  # Fisher’s exact test
  fisher_res <- fisher.test(cont_table)
  
  results_list[[cl]] <- tibble(
    cluster = cl,
    p_value = fisher_res$p.value,
    method = "Fisher’s exact (Region × CellType)"
  )
}

# ---- Combine all clusters ----
all_stats <- bind_rows(results_list)

# Add BH correction
all_stats$p_adj <- p.adjust(all_stats$p_value, method = "BH")

# Add significance stars
all_stats$significance <- cut(
  all_stats$p_adj,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "n.s.")
)

all_stats

```
```{r}
library(Seurat)
library(dplyr)
library(ggplot2)
library(khroma)

# -----------------------------
# 1. Define stress TF genes
# -----------------------------
stress_tf.genes <- intersect(genes_use, tf.genes) # use your previously defined lists

genes_use <- stress.genes

genes_use <- stress_tf.genes
# -----------------------------
# 2. Loop over spatial objects to compute StressTFScore
# -----------------------------
for(name in names(plot_info)){
  obj <- plot_info[[name]]$object
  
  # Filter genes present in object
  genes_present <- intersect(genes_use, rownames(obj))
  
  obj <- AddModuleScore(
    obj,
    features = list(genes_present),
    name = "StressScore"
  )
  
  # Save back
  plot_info[[name]]$object <- obj
}

library(Seurat)
library(ggplot2)
library(khroma)

# Define color palette for modality (update with your actual modalities)
modality_pal <- color("smooth rainbow")(37) 

# Loop over plot_info objects
for(name in names(plot_info)){
  obj_info <- plot_info[[name]]
  obj <- obj_info$object
  
  # Add timepoint and modality metadata
  obj$timepoint <- gsub("object_", "", name)  # e.g., "uw", "d4", etc.
  
  
  
  # Normalize StressTFScore1 to UW mean (calculate across all cells of UW)
  if("UW" %in% obj$timepoint){
    mean_UW <- mean(obj$StressScore1, na.rm = TRUE)
  } else {
    # fallback: use previously computed UW mean from the UW object
    mean_UW <- mean(plot_info$object_uw$object$StressScore1, na.rm = TRUE)
  }
  obj$StressTFScore_FC <- obj$StressScore1 / mean_UW
  
  # Violin plot by modality
  p_mod <- VlnPlot(
    obj,
    features = "StressScore1",
    pt.size = 0,
    cols = modality_pal
  ) + theme(
    text = element_text(family = "Arial", face = "bold", size = 14),
    axis.text = element_text(family = "Arial", face = "bold", size = 12),
    axis.title = element_text(family = "Arial", face = "bold", size = 12),
    axis.title.x = element_blank(),
    legend.text = element_text(family = "Arial", face = "bold", size = 12),
    legend.title = element_text(family = "Arial", face = "bold", size = 12)
  )
  
  ggsave(paste0(name, "_StressTFScore_FC_Vln_Modality.png"), p_mod, width = 6, height = 5, dpi = 300)
  
  # Save the object back
  plot_info[[name]]$object <- obj
}



library(dplyr)
library(ggplot2)
library(khroma)

# Prepare a data frame of averages per orig.ident across timepoints
avg_df <- lapply(names(plot_info), function(name) {
  obj <- plot_info[[name]]$object
  
  # Assign timepoint in meta.data
  obj@meta.data$timepoint <- gsub("object_", "", name)
  
  # Compute average per orig.ident
  df <- obj@meta.data %>%
    group_by(orig.ident) %>%
    summarise(
      StressTFScore1 = mean(StressTFScore1, na.rm = TRUE),
      .groups = "drop"
    )
  df$timepoint <- gsub("object_", "", name)
  return(df)
}) %>% bind_rows()

# Optional: order timepoints
avg_df$timepoint <- factor(avg_df$timepoint, levels = c("uw","d4","d7","d15","d30"))

# Color palette
pal <- color("smooth rainbow")(length(unique(avg_df$timepoint)), range = c(0.6,0.9))

# Bar plot with jittered points per orig.ident
p <- ggplot(avg_df, aes(x = timepoint, y = StressTFScore1, fill = timepoint)) +
  stat_summary(fun = mean, geom = "bar", position = "dodge") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_jitter(aes(color = timepoint), width = 0.15, size = 2, alpha = 0.7) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  labs(x = "Timepoint", y = "Average Stress TF Score") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", face = "bold", size = 14),
    axis.text = element_text(family = "Arial", face = "bold", size = 12),
    axis.title = element_text(family = "Arial", face = "bold", size = 12),
    legend.position = "none"
  )

# Save plot
ggsave("StressTFScore1_Barplot_Timepoints.png", p, width = 6, height = 5, dpi = 300)


library(dplyr)
library(ggplot2)
library(khroma)

# Combine all plot_info objects into a single dataframe
violin_df <- lapply(names(plot_info), function(name) {
  obj <- plot_info[[name]]$object
  df <- obj@meta.data
  df$timepoint <- gsub("object_", "", name)
  return(df)
}) %>% bind_rows()

# Optional: order timepoints
violin_df$timepoint <- factor(violin_df$timepoint, levels = c("uw","d4","d7","d15","d30"))

# Color palette
pal <- color("smooth rainbow")(length(unique(violin_df$timepoint)), range = c(0.6,0.9))
library(ggpubr)

comparisons <- list(
  c("uw", "d4"),
  c("uw", "d7"),
  c("uw", "d15"),
  c("uw", "d30")
)


# Violin plot
p <- ggplot(
  violin_df,
  aes(x = timepoint, y = StressScore1, fill = timepoint)
) +
  geom_violin(trim = FALSE, alpha = 1) +
  
  # very faint single-cell points
  geom_jitter(
    aes(color = timepoint),
    width = 0.15,
    size = 0.1,
    alpha = 0.01
  ) +
  
  # mean per timepoint (red dot)
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2,
    color = "red"
  ) +
  
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  
  labs(
    x = "Timepoint",
    y = "Stress TF Score"
  ) +
  
  stat_compare_means(
    ref.group = "uw",
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  
  theme_ggprism_mod() +
  theme(
    text = element_text(family = "Arial", face = "bold", size = 14),
    axis.text = element_text(family = "Arial", face = "bold", size = 12),
    axis.title = element_text(family = "Arial", face = "bold", size = 12),
    legend.position = "none"
  )
# Save plot
ggsave("StressTFScoreScore_MouseSpatialViolin_TimepointsSmall.png", p, width = 4, height = 4, dpi = 300)
ggsave("StressTFScoreScore_MouseSpatialViolin_Timepoints.svg", p, width = 6, height = 5, dpi = 300)


ggsave("StressTFScore1_Violin_Timepoints.png", p, width = 6, height = 5, dpi = 300)

SpatialFeaturePlot(plot_info$object_d30$object, "StressTFScore1",images = "slice1")


library(Seurat)
library(ggplot2)

# Loop through each object in plot_info
for(name in names(plot_info)) {
  obj <- plot_info[[name]]$object
  
  # Optionally pick the first image if multiple exist
  img_name <- if("slice1" %in% names(obj@images)) "slice1" else names(obj@images)[1]
  
  # Create SpatialFeaturePlot
  p <- SpatialFeaturePlot(obj, features = "StressTFScore1", images = img_name) +
    ggtitle(paste0("StressTFScore1 - ", gsub("object_", "", name)))
  
  # Print to viewer
  print(p)
  
  # Save each plot
  ggsave(filename = paste0("StressTFScore1_Spatial_", gsub("object_", "", name), ".png"),
         plot = p, width = 6, height = 5, dpi = 300)
}
```

```{r}
library(Seurat)
library(dplyr)
library(ggplot2)
library(khroma)

# ------------------------------
# 1. Loop through plot_info to make modality violin plots
# ------------------------------
modality_pal <- color("smooth rainbow")(2)  # adjust number of modalities

for(name in names(plot_info)){
  obj_info <- plot_info[[name]]
  obj <- obj_info$object
  
  # Add timepoint metadata
  obj$timepoint <- gsub("object_", "", name)
  
  # Violin plot by modality
  p_mod <- VlnPlot(
    obj,
    features = "StressScore1",
    pt.size = 0
  ) + theme(
    text = element_text(family = "Arial", face = "bold", size = 14),
    axis.text = element_text(family = "Arial", face = "bold", size = 12),
    axis.title = element_text(family = "Arial", face = "bold", size = 12),
    axis.title.x = element_blank(),
    legend.text = element_text(family = "Arial", face = "bold", size = 12),
    legend.title = element_text(family = "Arial", face = "bold", size = 12)
  )
  
  ggsave(paste0(name, "_StressScore1_Vln_Modality.png"), p_mod, width = 6, height = 5, dpi = 300)
  
  # Save updated object back
  plot_info[[name]]$object <- obj
}

# ------------------------------
# 2. Average per orig.ident across timepoints (bar plot)
# ------------------------------
avg_df <- lapply(names(plot_info), function(name) {
  obj <- plot_info[[name]]$object
  obj@meta.data$timepoint <- gsub("object_", "", name)
  
  df <- obj@meta.data %>%
    group_by(orig.ident) %>%
    summarise(
      StressScore1 = mean(StressScore1, na.rm = TRUE),
      .groups = "drop"
    )
  df$timepoint <- gsub("object_", "", name)
  return(df)
}) %>% bind_rows()

avg_df$timepoint <- factor(avg_df$timepoint, levels = c("uw","d4","d7","d15","d30"))

pal <- color("smooth rainbow")(length(unique(avg_df$timepoint)), range = c(0.6,0.9))

p_bar <- ggplot(avg_df, aes(x = timepoint, y = StressScore1, fill = timepoint)) +
  stat_summary(fun = mean, geom = "bar", position = "dodge") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_jitter(aes(color = timepoint), width = 0.15, size = 2, alpha = 0.7) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  labs(x = "Timepoint", y = "Average StressScore1") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", face = "bold", size = 14),
    axis.text = element_text(family = "Arial", face = "bold", size = 12),
    axis.title = element_text(family = "Arial", face = "bold", size = 12),
    legend.position = "none"
  )

ggsave("StressScore1_Barplot_Timepoints.png", p_bar, width = 6, height = 5, dpi = 300)

# ------------------------------
# 3. Combined violin plot across all timepoints
# ------------------------------
violin_df <- lapply(names(plot_info), function(name) {
  obj <- plot_info[[name]]$object
  df <- obj@meta.data
  df$timepoint <- gsub("object_", "", name)
  return(df)
}) %>% bind_rows()

# Rename timepoints to desired labels
violin_df$timepoint <- recode(
  violin_df$timepoint,
  "uw"  = "UW",
  "d4"  = "D4PW",
  "d7"  = "D7PW",
  "d30" = "D30PW"
)

# Set factor order
violin_df$timepoint <- factor(
  violin_df$timepoint,
  levels = c("UW", "D4PW", "D7PW", "D30PW")
)

# Color palette
pal <- color("smooth rainbow")(length(levels(violin_df$timepoint)), range = c(0.6, 0.9))

# Violin plot (NO POINTS)
p_violin <- ggplot(violin_df, aes(x = timepoint, y = StressScore1, fill = timepoint)) +
  geom_violin(trim = FALSE, alpha = 0.9) +
  geom_jitter(
    aes(color = timepoint),
    width = 0.12,
    size = 0.11,
    alpha = 0.15
  ) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  labs(x = "Timepoint", y = "StressTFScore_FC") +
  theme_ggprism_mod() +
  theme(
    text = element_text(family = "Arial", face = "bold", size = 14),
    axis.text = element_text(family = "Arial", face = "bold", size = 12),
    axis.title = element_text(family = "Arial", face = "bold", size = 12),
    legend.position = "none"
  )

ggsave(
  "StressScore1_Violin_Timepoints.png",
  p_violin,
  width = 6,
  height = 5,
  dpi = 300
)


# ------------------------------
# 4. SpatialFeaturePlots for all objects
# ------------------------------
for(name in names(plot_info)) {
  obj <- plot_info[[name]]$object
  img_name <- if("slice1" %in% names(obj@images)) "slice1" else names(obj@images)[1]
  
  p <- SpatialFeaturePlot(obj, features = "StressScore1", images = img_name) +
    ggtitle(paste0("StressScore1 - ", gsub("object_", "", name)))
  
  print(p)
  
  ggsave(filename = paste0("StressScore1_Spatial_", gsub("object_", "", name), ".png"),
         plot = p, width = 6, height = 5, dpi = 300)
}

```

```{r}
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(khroma)

# ------------------------------
# 1. Loop through plot_info to make modality violin plots
# ------------------------------
modality_pal <- color("smooth rainbow")(2)  # adjust number of modalities

for(name in names(plot_info)){
  obj_info <- plot_info[[name]]
  obj <- obj_info$object
  
  # Add timepoint metadata
  obj$timepoint <- gsub("object_", "", name)
  
  # Violin plot by modality
  p_mod <- VlnPlot(
    obj,
    features = "StressScore1",
    pt.size = 0
  ) + theme(
    text = element_text(family = "Arial", face = "bold", size = 14),
    axis.text = element_text(family = "Arial", face = "bold", size = 12),
    axis.title = element_text(family = "Arial", face = "bold", size = 12),
    axis.title.x = element_blank(),
    legend.text = element_text(family = "Arial", face = "bold", size = 12),
    legend.title = element_text(family = "Arial", face = "bold", size = 12)
  )
  
  ggsave(paste0(name, "_StressScore1_Vln_Modality.png"), p_mod, width = 6, height = 5, dpi = 300)
  
  # Save updated object back
  plot_info[[name]]$object <- obj
}

# ------------------------------
# 2. Average per orig.ident across timepoints (bar plot)
# ------------------------------
avg_df <- lapply(names(plot_info), function(name) {
  obj <- plot_info[[name]]$object
  obj@meta.data$timepoint <- gsub("object_", "", name)
  
  df <- obj@meta.data %>%
    group_by(orig.ident) %>%
    summarise(
      StressScore1 = mean(StressScore1, na.rm = TRUE),
      .groups = "drop"
    )
  df$timepoint <- gsub("object_", "", name)
  return(df)
}) %>% bind_rows()

avg_df$timepoint <- factor(avg_df$timepoint, levels = c("uw","d4","d7","d15","d30"))

pal <- color("smooth rainbow")(length(unique(avg_df$timepoint)), range = c(0.6,0.9))

p_bar <- ggplot(avg_df, aes(x = timepoint, y = StressScore1, fill = timepoint)) +
  stat_summary(fun = mean, geom = "bar", position = "dodge") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_jitter(aes(color = timepoint), width = 0.15, size = 2, alpha = 0.7) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  labs(x = "Timepoint", y = "Average StressScore1") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", face = "bold", size = 14),
    axis.text = element_text(family = "Arial", face = "bold", size = 12),
    axis.title = element_text(family = "Arial", face = "bold", size = 12),
    legend.position = "none"
  )

ggsave("StressScore1_Barplot_Timepoints.png", p_bar, width = 6, height = 5, dpi = 300)

# ------------------------------
# 3. Combined violin plot across all timepoints
# ------------------------------
violin_df <- lapply(names(plot_info), function(name) {
  obj <- plot_info[[name]]$object
  df <- obj@meta.data
  df$timepoint <- gsub("object_", "", name)
  return(df)
}) %>% bind_rows()

violin_df$timepoint <- factor(violin_df$timepoint, levels = c("uw","d4","d7","d15","d30"))

pal <- color("smooth rainbow")(length(unique(violin_df$timepoint)), range = c(0.6,0.9))

p_violin <- ggplot(violin_df, aes(x = timepoint, y = StressScore1, fill = timepoint)) +
  geom_violin(trim = FALSE, alpha = 0.3) +
  geom_jitter(aes(color = timepoint), width = 0.15, size = 0.5, alpha = 0.7) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  labs(x = "Timepoint", y = "StressScore1") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", face = "bold", size = 14),
    axis.text = element_text(family = "Arial", face = "bold", size = 12),
    axis.title = element_text(family = "Arial", face = "bold", size = 12),
    legend.position = "none"
  )

ggsave("StressScore1_Violin_Timepoints.png", p_violin, width = 6, height = 5, dpi = 300)


p_metaclusters <- ggplot(violin_df, aes(x = metaclusters, y = StressScore1, fill = metaclusters)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(aes(color = metaclusters), width = 0.15, size = 0.1, alpha = 0.2) +
  scale_fill_manual(values = metacluster_colors) +
  scale_color_manual(values = metacluster_colors) +
  labs(x = "Metacluster", y = "StressScore1") +
  theme_ggprism_mod() +
  theme(
    text = element_text(family = "Arial", face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 0.5, family = "Arial", face = "bold", size = 12),
    axis.text.y = element_text(family = "Arial", face = "bold", size = 12),
    axis.title = element_text(family = "Arial", face = "bold", size = 12),
    legend.position = "none"
  )

ggsave("StressScore1_Violin_Metaclusters.png", p_metaclusters, width = 8, height = 5, dpi = 300)




# ------------------------------
# 4. SpatialFeaturePlots and UMAP Feature Plots for all objects
# ------------------------------
for(name in names(plot_info)) {
  obj <- plot_info[[name]]$object
  timepoint_name <- gsub("object_", "", name)
  
  # --- Spatial Feature Plot ---
  img_name <- if("slice1" %in% names(obj@images)) "slice1" else names(obj@images)[1]
  
  p_spatial <- SpatialFeaturePlot(obj, features = "StressScore1", images = img_name) +
    ggtitle(paste0("StressScore1 - Spatial - ", timepoint_name))
  
  print(p_spatial)
  
  ggsave(filename = paste0("StressScore1_Spatial_", timepoint_name, ".png"),
         plot = p_spatial, width = 6, height = 5, dpi = 300)
  
  # --- UMAP Feature Plot (NEW) ---
  # Generates a plot showing StressScore1 expression on the UMAP embedding.
  # Assumes UMAP reduction coordinates are available in the Seurat object.
  p_umap <- FeaturePlot_scCustom(obj, features = "StressScore1", reduction = "full.umap.sketch", na_cutoff = 0.05) +
    ggtitle(paste0("StressScore1 - UMAP - ", timepoint_name)) +
    theme(
      text = element_text(family = "Arial", face = "bold", size = 14),
      plot.title = element_text(hjust = 0.5)
    )
  
  print(p_umap)
  
  ggsave(filename = paste0("StressScore1_UMAP_", timepoint_name, ".png"),
         plot = p_umap, width = 6, height = 5, dpi = 300)
}
```


```{r}
# Define the clusters you want to plot
library(ggplot2)
library(dplyr)
library(cowplot)   # for get_legend()

# Define clusters
clusters_of_interest <- c("20","1","0")

# Load your predefined color palette
color_pals <- readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 1/color_palette.rds")

# Loop through all Seurat objects
for (cond in names(plot_info)) {
  
  seurat_obj <- plot_info[[cond]]$object
  
  # Make sure Idents are set to banksy clusters
  Idents(seurat_obj) <- "banksy_cluster_0.5"
  
  # Pull metadata only for the clusters of interest
  meta_df <- seurat_obj@meta.data %>%
    filter(banksy_cluster_0.5 %in% clusters_of_interest) %>%
    dplyr::count(banksy_cluster_0.5, first_type) %>%
    group_by(banksy_cluster_0.5) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
  # Ensure consistent ordering of clusters on y-axis (since flipped)
  meta_df$banksy_cluster_0.5 <- factor(
    meta_df$banksy_cluster_0.5,
    levels = rev(clusters_of_interest)  # reversed so "28" at top
  )
  
  # Ensure consistent cell type ordering in legend
  meta_df$first_type <- factor(meta_df$first_type, levels = names(color_pals))
  
  # Restrict palette to only the cell types present in this subset
  type_pal <- color_pals[levels(meta_df$first_type)]
  
  # Base plot
  p <- ggplot(meta_df, aes(
    x = banksy_cluster_0.5,
    y = proportion,
    fill = first_type
  )) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = type_pal, drop = FALSE) +
    theme_ggprism_mod() +
    labs(
      x = "Banksy Cluster",
      y = "Proportion of Cell Types",
      fill = "Cell Type"
    ) +
    theme(
      strip.text = element_text(size = 20, color = "black"),
      axis.text.x = element_text(size = 20, color = "black"),
      axis.text.y = element_text(size = 20, color = "black", face = "bold"),
      axis.title.x = element_text(size = 20, color = "black"),
      axis.title.y = element_text(size = 20, color = "black"),
      panel.border = element_rect(fill = NA, color = "white"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    ) +
    coord_flip()
  
  # Save barplot without legend
  p_nolegend <- p + NoLegend()
  ggsave(
    filename = paste0(cond, "_SelectedBanksyClusters_CellTypeCompositionD30PW.svg"),
    plot = p_nolegend,
    width = 15, height = 4, dpi = 320
  )
  
  # Extract legend and save separately
  legend <- get_legend(p + theme(legend.position = "right"))
  legend_plot <- ggdraw(legend)
  
  ggsave(
    filename = paste0(cond, "_Legend_CellTypeCompositionD30PW.svg"),
    plot = legend_plot,
    width = 20, height = 20, dpi = 320
  )
}
```


```{r}
library(Seurat) # Assumed for Idents() and NoLegend(), ensure this is loaded in your environment
library(ggplot2)
library(dplyr)
library(cowplot) # for get_legend()
library(readr)   # You might need this if using the full filepath/readRDS

# --- USER DEFINITIONS (UPDATE THESE) ---


# 2. Your 'plot_info' object should be defined here or already loaded in your environment.
# Assuming 'plot_info' is a named list where keys are the condition names (e.g., "UW", "D4PW")
# and each element contains a Seurat object under the key 'object'.
# For example, if you want to explicitly run for UW, D4PW, D7PW, D30PW:
conditions_to_run <- c("UW", "D4PW", "D7PW", "D30PW")
# Check if your plot_info keys match these conditions.


color_pals = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 1/color_palette.rds")
# --- SCRIPT START ---

# Loop through all Seurat objects in the plot_info list
for (cond in names(plot_info)) {
  
  # Check if the current condition is one of the desired ones (if you need to filter)
  # If you only want to run UW, D4PW, D7PW, D30PW:
  # if (!(cond %in% conditions_to_run)) { next }
  
  message(paste("Processing condition:", cond))
  
  seurat_obj <- plot_info[[cond]]$object
  
  # Make sure Idents are set to banksy clusters
  Idents(seurat_obj) <- "banksy_cluster_0.5"
  
  # Pull metadata for *ALL* clusters (no filter applied)
  meta_df <- seurat_obj@meta.data %>%
    # The filter step is removed to include all clusters
    dplyr::count(banksy_cluster_0.5, first_type) %>%
    group_by(banksy_cluster_0.5) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
  # Get the names of all unique clusters for dynamic ordering
  # FIX: Use a more robust check to ensure only valid numeric cluster labels are included,
  # and filter out any non-numeric strings or explicit NAs.
  cluster_labels_to_sort <- as.character(unique(meta_df$banksy_cluster_0.5))
  
  # Filter out labels that result in NA when converted to numeric, then sort
  all_clusters <- cluster_labels_to_sort[!is.na(suppressWarnings(as.numeric(cluster_labels_to_sort)))] %>%
    as.numeric() %>%
    sort(decreasing = FALSE) %>%
    as.character()
  
  # Ensure consistent ordering of clusters on y-axis (since flipped)
  # We use rev(all_clusters) so the lowest numbered cluster is at the bottom
  meta_df$banksy_cluster_0.5 <- factor(
    meta_df$banksy_cluster_0.5,
    levels = rev(all_clusters)
  )
  
  # Ensure consistent cell type ordering in legend
  meta_df$first_type <- factor(meta_df$first_type, levels = names(color_pals))
  
  # Restrict palette to only the cell types present in this subset
  type_pal <- color_pals[levels(meta_df$first_type)]
  
  # Base plot
  p <- ggplot(meta_df, aes(
    x = banksy_cluster_0.5,
    y = proportion,
    fill = first_type
  )) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = type_pal, drop = FALSE) +
    # Assuming theme_ggprism_mod() is defined in your environment
    theme_ggprism_mod() +
    labs(
      x = "Banksy Cluster",
      y = "Proportion of Cell Types",
      fill = "Cell Type"
    ) +
    theme(
      strip.text = element_text(size = 20, color = "black"),
      axis.text.x = element_text(size = 20, color = "black"),
      axis.text.y = element_text(size = 20, color = "black", face = "bold"),
      axis.title.x = element_text(size = 20, color = "black"),
      axis.title.y = element_text(size = 20, color = "black"),
      panel.border = element_rect(fill = NA, color = "white"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    ) +
    coord_flip()
  
  # Calculate height dynamically: 1.5 inches per cluster for good spacing
  # REMOVED DYNAMIC HEIGHT CALCULATION AS REQUESTED.
  fixed_height <- 20 # Using a fixed height of 20 inches
  
  # Save barplot without legend
  # File name updated to reflect ALL clusters are included, and condition is dynamic
  p_nolegend <- p + theme(legend.position = "none") # Replaced NoLegend() for robustness
  
  ggsave(
    filename = paste0(cond, "_AllBanksyClusters_CellTypeComposition.svg"),
    plot = p_nolegend,
    width = 15, height = 25, dpi = 320 # Use fixed height
  )
  
  message(paste("Saved barplot for", cond, "to", paste0(cond, "_AllBanksyClusters_CellTypeComposition.svg")))
  
  # Extract legend and save separately
  legend <- get_legend(p + theme(legend.position = "right"))
  legend_plot <- ggdraw(legend)
  
  # Save legend
  ggsave(
    filename = paste0(cond, "_Legend_CellTypeComposition.svg"),
    plot = legend_plot,
    width = 20, height = 25, dpi = 320 # Keep high dimensions for detailed legend
  )
  
  message(paste("Saved legend for", cond, "to", paste0(cond, "_Legend_CellTypeComposition.svg")))
}
```