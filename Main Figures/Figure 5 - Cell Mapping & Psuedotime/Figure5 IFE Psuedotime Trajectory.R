```{r Keraitnocyte IFE Focus Psuedotime}


alldata.1 = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Sequencing Repository/Data/OWHA RDS files/Subclusters/KeratinocytesAnnotated_2026-01-06.rds")

alldata.2 = subset(alldata.1, ident = c("Basal IV", "Basal III", "Basal II","Spinous I", "Spinous II", "Spinous III", "Spinous IV","Ker.Cycling III", "Ker.Cycling II", "Ker.Cycling I"))
#alldata.2 <- subset(alldata.1, subset = Sox6 > 0)


Idents(alldata.2) <- alldata.2$celltype


# Ensure your condition labels are consistent
alldata.2@meta.data <- alldata.2@meta.data %>%
  mutate(condition = case_when(
    timepoint == "UW" ~ "Unwounded",
    timepoint %in% c("D4PW","D7PW") ~ "Wounded",
    TRUE ~ as.character(timepoint)
  ))

# Set identity to condition for subsetting
Idents(alldata.2) <- alldata.2$condition
# Subset Unwounded and Wounded datasets
alldata_unwounded <- subset(alldata.2, idents = "Unwounded")
alldata_wounded   <- subset(alldata.2, idents = "Wounded")

#Remiove Low Count Layers They Cause an Error To Be Thown
subset_by_min_cells <- function(seurat_obj, min_cells = 100) {
  # Count cells per orig.ident
  idents_to_keep <- seurat_obj$orig.ident %>% 
    table() %>%
    .[. >= min_cells] %>%
    names()
  
  # Subset Seurat object
  seurat_obj_subset <- subset(seurat_obj, subset = orig.ident %in% idents_to_keep)
  return(seurat_obj_subset)
}

alldata.2 <- subset_by_min_cells(alldata.2, min_cells = 20)
alldata_wounded <- subset_by_min_cells(alldata_wounded, min_cells = 20)
alldata_unwounded <- subset_by_min_cells(alldata_unwounded, min_cells = 20)

#Seperate By Modality and Run Analysis
Idents(alldata_unwounded) <- alldata_unwounded$modality
Idents(alldata_wounded) <- alldata_wounded$modality
Idents(alldata.2) <- alldata.2$modality

alldata_unwounded <- subset(alldata_unwounded, idents = "snRNAseq")
alldata_wounded <- subset(alldata_wounded, idents = "snRNAseq")

alldata.2 <- subset(alldata.2, idents = "snRNAseq")

process_seurat_subset <- function(seurat_obj) {
  seurat_obj <- Split_Layers(seurat_obj, split.by = "orig.ident")
  seurat_obj <- SCTransform(seurat_obj, assay = "RNA", new.assay.name = "SCT", vars.to.regress = "percent_mito", verbose = TRUE)
  DefaultAssay(seurat_obj) <- "SCT"
  seurat_obj <- RunPCA(seurat_obj, verbose = TRUE)
  seurat_obj <- RunHarmony(seurat_obj, group.by.vars = c("orig.ident"), assay.use = "SCT", theta = c(3), verbose = TRUE)
  
  # Calculate number of PCs explaining 90% of variance
  
  
  # Clustering and UMAP
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = 0.8, group.singletons = FALSE, cluster.name = "harmony_clusters") %>%
    RunUMAP(reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")
  
  return(seurat_obj)
}

process_seurat_subset <- function(seurat_obj) {
  #Optional: split layers by sample (if this is required for your pipeline)
  seurat_obj <- Split_Layers(seurat_obj, split.by = "orig.ident")
  
  # SCTransform normalization
  seurat_obj <- SCTransform(
    seurat_obj,
    assay = "RNA",
    new.assay.name = "SCT",
    vars.to.regress = "percent_mito",
    verbose = TRUE
  )
  
  DefaultAssay(seurat_obj) <- "SCT"
  
  # Dimensionality reduction
  seurat_obj <- RunPCA(seurat_obj, verbose = TRUE)
  
  # Calculate number of PCs explaining 90% of variance
  
  # Clustering and UMAP based on PCA
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20) %>%
    FindClusters(resolution = 0.8, group.singletons = T, cluster.name = "pca_clusters") %>%
    RunUMAP(dims = 1:20, reduction.name = "umap.pca")
  
  return(seurat_obj)
}



#Run Function and Harmony Intergration on
alldata_unwounded <- process_seurat_subset(alldata_unwounded)
alldata_wounded   <- process_seurat_subset(alldata_wounded)

alldata.2 <- process_seurat_subset(alldata.2)


Idents(alldata_unwounded) <- "harmony_clusters"
Idents(alldata_wounded) <- "harmony_clusters"

DimPlot(alldata_unwounded, reduction = "umap.harmony", cols = colors_pals, label = TRUE) + ggtitle("Unwounded")
DimPlot(alldata_wounded, reduction = "umap.harmony", cols = colors_pals, label = TRUE) + ggtitle("Wounded")

#Remove Clusters Including Cortex Markers and Doublets as well
# Remove cluster 8 from alldata_unwounded
#alldata_unwounded <- subset(alldata_unwounded, idents = c(13), invert = TRUE)

# Remove clusters 7, 9, 15 from alldata_wounded
#alldata_wounded <- subset(alldata_wounded, idents = c(17,18), invert = TRUE)
#alldata.2 <- subset(alldata.2, idents = c(), invert = T)



alldata_unwounded <- process_seurat_subset(alldata_unwounded)
alldata_wounded   <- process_seurat_subset(alldata_wounded)


p_unwounded <- DimPlot_scCustom(
  alldata_unwounded,
  reduction = "umap.harmony",
  group.by = "celltype",
  colors_use = colors_pals,
  label = TRUE
) + ggtitle("Unwounded")

p_wounded <- DimPlot_scCustom(
  alldata_wounded,
  reduction = "umap.harmony",
  group.by = "celltype",
  colors_use = colors_pals,
  label = TRUE
) + ggtitle("Wounded")

ggsave(
  filename = "UMAP_Unwounded_Celltypes.png",
  plot = p_unwounded,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "UMAP_Wounded_Celltypes.png",
  plot = p_wounded,
  width = 7,
  height = 5,
  dpi = 300
)

DimPlot_scCustom(alldata.2, reduction = "umap.harmony", group.by = "celltype", colors_use =  colors_pals, label = TRUE) + ggtitle("Basal IV All Timepoints")

Idents(alldata_unwounded) <- "celltype"
Idents(alldata_wounded) <- "celltype"


#Run Psuedtotime
library(monocle3)
library(SeuratWrappers)



run_pseudotime_max_Krt14 <- function(seurat_obj, umap_name = "umap.pca") {
  library(monocle3)
  library(SeuratWrappers)
  
  # 1. Identify the root cell: highest Krt14 expression (SCT assay)
  expr_vec <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")["Sox5", ]
  root_cell <- names(which.max(expr_vec))
  
  # 2. Convert to CDS and set UMAP
  cds <- as.cell_data_set(seurat_obj)
  reducedDims(cds)$UMAP <- Embeddings(seurat_obj, umap_name)
  colData(cds)$celltype <- Idents(seurat_obj)
  
  # 3. Monocle3 pseudotime steps
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- order_cells(cds, root_cells = root_cell)
  
  # 4. Extract pseudotime
  pseudotime_vals <- cds@principal_graph_aux[["UMAP"]]$pseudotime
  pseudotime_df <- data.frame(
    pseudotime = pseudotime_vals,
    celltype = seurat_obj$celltype[names(pseudotime_vals)],
    timepoint = seurat_obj$timepoint[names(pseudotime_vals)]
  )
  
  return(list(cds = cds, pseudotime_df = pseudotime_df, root_cell = root_cell))
}



# Run on Unwounded
pt_unwounded <- run_pseudotime_max_Krt14(alldata_unwounded)

# Run on Wounded
pt_wounded <- run_pseudotime_max_Krt14(alldata_wounded)

pt_wounded <- run_pseudotime_max_Krt14(alldata.2)


# Combine for joint comparison if needed
pseudotime_combined <- rbind(
  cbind(pt_unwounded$pseudotime_df, group = "Unwounded"),
  cbind(pt_wounded$pseudotime_df, group = "Wounded")
)


pt_unwounded$cds <- order_cells(pt_unwounded$cds)
# Plot Unwounded pseudotime
p_unwounded <- plot_cells(
  pt_unwounded$cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
) + ggtitle("Unwounded: Pseudotime") + NoAxes()


pt_wounded$cds <- order_cells(pt_wounded$cds)

# Plot Wounded pseudotime
p_wounded <- plot_cells(
  pt_wounded$cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
) + ggtitle("Wounded: Pseudotime")+ NoAxes()
ggsave("Monocle3_YvonUpdatesUWOwtihLegend.svg", plot = p_unwounded, width = 5, height = 5, dpi = 300)
ggsave("Monocle3_YvonUpdatesWoundedwithLegend.svg", plot = p_wounded, width = 5, height = 5, dpi = 300)






# Save
ggsave("Monocle3_UnwoundedNuceliOnly_Pseudotime_ManualRoot.png", plot = p_unwounded, width = 5, height = 5, dpi = 300)
ggsave("Monocle3_Sox6+OnlyWoundedD4PWD7PWNucleiOnly_Pseudotime_ManualRoot.png", plot = p_wounded, width = 5, height = 5, dpi = 300)


timepoint_colors <- c(
  "UW"    = "#D41159",   # bold red for unwounded
  "D1PW"  = "#F8766D",   # red-orange
  "D2PW"  = "#E7861B",   # orange
  "D4PW"  = "#FFD700",   # golden yellow
  "D7PW"  = "#7CAE00",   # yellow-green
  "D15PW" = "#00BFC4",   # bright teal
  "D30PW" = "purple3"    # deep violet-blue
)


# Timepoint grouping with legend
p1 <- DimPlot(alldata_unwounded, group.by = "timepoint", cols = timepoint_colors) +
  NoAxes() + customtheme
p2 <- DimPlot(alldata_wounded, group.by = "timepoint",
              cols = timepoint_colors, shuffle = TRUE) +
  NoAxes() + customtheme
p3 <- DimPlot(alldata.2, group.by = "timepoint",
              cols = timepoint_colors, shuffle = TRUE) +
  NoAxes() + customtheme


ggsave("IFEPseudotime_UnwoundedNucleiOnlyTimepointRootCells_umap.png", plot = p1, width = 5, height = 5, dpi = 300)
ggsave("IFEPseudotime_D4PWD7PW_NucleiOnlyTimepointRoot_umap.png", plot = p2, width = 5, height = 5, dpi = 300)
ggsave("IFEPseudotime_AllTimepoints_NucleiOnlyTimepointSox6+Cells_umap.png", plot = p3, width = 5, height = 5, dpi = 300)

# Timepoint grouping without legend
p1 <- DimPlot(alldata_unwounded, group.by = "timepoint", cols = timepoint_colors) +
  NoAxes() + customtheme + NoLegend()
p2 <- DimPlot(alldata_wounded, group.by = "timepoint",
              cols = timepoint_colors, shuffle = TRUE) +
  NoAxes() + customtheme + NoLegend()

ggsave("IFEPseudotime_UnwoundedNucleiOnlySox6+CellsBlank_umap.png", plot = p1, width = 5, height = 5, dpi = 300)
ggsave("IFEPseudotime_D4PWD7PW_NucleiOnlyTimepointSox6+CellsBlank_umap.png", plot = p2, width = 5, height = 5, dpi = 300)


Idents(alldata_unwounded) <- "celltype"
Idents(alldata_wounded) <- "celltype"
Idents(alldata.2) <- "celltype"


# Color palette grouping with legend
s1 <- DimPlot(alldata_wounded) +
  NoAxes() + customtheme
s2 <- DimPlot(alldata_unwounded) +
  NoAxes() + customtheme
p3 <- DimPlot(alldata.2) +
  NoAxes() + customtheme

ggsave("IFEPseudotime_WoundedD4PWD7PWNucleiOnlyLegendRoot_umap.png", plot = s1, width = 6, height = 5, dpi = 300)
ggsave("IFEPseudotime_UnWoundedBNucleiOnlyWithLegendRoot_umap.png", s2, width = 6, height = 5, dpi = 300)
ggsave("IFEPseudotime_AllTimepointsNucleiOnlyWithLegendRoot_umap.png", p3, width = 6, height = 5, dpi = 300)

# Color palette grouping without legend
p1 <- DimPlot(alldata_unwounded, reduction = "umap.harmony") +
  NoAxes() + customtheme + NoLegend()
p2 <- DimPlot(alldata.2, reduction = "umap.harmony") +
  NoAxes() + customtheme + NoLegend()
p3 <- DimPlot(alldata.2, reduction = "umap.harmony") +
  NoAxes() + customtheme + NoLegend()

ggsave("IFEPseudotime_UnwoundedNucleiOnlyLegendRootSox6_umapBlank.png", plot = p1, width = 6, height = 5, dpi = 300)
ggsave("IFEPseudotime_WoundedNucleiOnlyWithLegendRootSox6_umapBlank.png", plot = p2, width = 6, height = 5, dpi = 300)
ggsave("IFEPseudotime_AllTimepointsNucleiOnlyLegendSox6+Cells_umapBlank.png", plot = p3, width = 6, height = 5, dpi = 300)



p_unwounded + p1 + s2
ggsave("IFEPseudotime_UnwoundedNucleiOnlyWithCornfiedLegendRootManual_ALlPlots.png", width = 15, height = 5, dpi = 300)

p_wounded + p2 + s1
ggsave("IFEPseudotime_UWD4PWD7PWBasalIVOnlyLegendRootManual_ALlPlots.png", width = 15, height = 5, dpi = 300)

gg_violin <- ggplot(pseudotime_combined, aes(x = group, y = pseudotime, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +  # Add p-value
  scale_fill_manual(values = c("Unwounded" = "#1f77b4", "Wounded" = "#d62728")) +
  theme_classic() +
  ylab("Pseudotime") +
  ggtitle("Pseudotime Comparison Between Conditions") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 13, face = "bold")
  )

# Print the plot
print(gg_violin)

saveRDS(alldata_unwounded,"IFEPseudotimeUnwoundedNucleiOnly.rds")
saveRDS(alldata_wounded,"IFEPseudotimeWoundedNucleiOnlyUWOD4PW.rds")
saveRDS(alldata.2,"IFEPseudotimeWoundedAllTimePoints.rds")


# Make sure you're selecting valid gene IDs present in the CDS
valid_genes <- rownames(pt_unwounded$cds)  # gene IDs
sig_gene_ids <- rownames(sig_genes_unwounded)[rownames(sig_genes_unwounded) %in% valid_genes]

# Now subset and run find_gene_modules

get_gene_modules <- function(cds, gene_ids, num_dim = 30, resolution = 1e-2) {
  # Subset to genes of interest
  gene_ids <- intersect(gene_ids, rownames(cds))
  cds_sub <- cds[gene_ids, ]
  
  # Preprocess & reduce
  cds_sub <- preprocess_cds(cds_sub, num_dim = num_dim, method = "PCA")
  cds_sub <- reduce_dimension(cds_sub, reduction_method = "UMAP", preprocess_method = "PCA")
  
  # Module finding
  mods <- find_gene_modules(cds_sub, resolution = resolution)
  return(mods)
}

gene_modules_unwounded <- get_gene_modules(
  pt_unwounded$cds,
  sig_gene_ids,
  num_dim = 30,
  resolution = 1e-2
)

gene_modules_wounded <- get_gene_modules(
  pt_wounded$cds,
  sig_gene_ids,
  num_dim = 30,
  resolution = 1e-2
)

plot_module_trend <- function(cds, gene_modules) {
  gene_module_expr <- aggregate_gene_expression(cds, gene_modules)
  
  ggplot(gene_module_expr, aes(x = pseudotime, y = expression, color = module)) +
    geom_smooth(se = FALSE, method = "loess", span = 0.3) +
    theme_minimal() +
    labs(y = "Module expression (mean z-score)", x = "Pseudotime")
}

ggplot(gene_modules_unwounded, aes(x = dim_1, y = dim_2, color = factor(module))) +
  geom_point(size = 1.5, alpha = 0.8) +
  theme_minimal() +
  labs(color = "Module", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "right")

saveRDS(alldata_unwounded, "IFEPseudtoimeUnwounded.rds")
saveRDS(alldata_wounded, "IFEPseudtoimeWounded.rds")

alldata_unwounded = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Jonathan Current Analysis/Pseudotime/IFEPseudtoimeUnwounded.rds")
alldata_wounded = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Jonathan Current Analysis/Pseudotime/IFEPseudtoimeWounded.rds")

library(ggplot2)
library(patchwork)

# Unwounded
g1_unw <- FeaturePlot_scCustom(alldata_unwounded, "Krt10")
g2_unw <- FeaturePlot_scCustom(alldata_unwounded, "Krt1")
g3_unw <- FeaturePlot_scCustom(alldata_unwounded, "Krt14")
g4_unw <- FeaturePlot_scCustom(alldata_unwounded, "Krt16")


plot_unw <- g1_unw + g2_unw + g3_unw + g4_unw + plot_layout(ncol = 4)

ggsave("Krt_ExpressionRootIFE_NucleiOnlyCornfiedUnwounded.png",
       plot = plot_unw,
       width = 15, height = 4, dpi = 300)

# Wounded
g1_w <- FeaturePlot_scCustom(alldata_wounded, "Krt10")
g2_w <- FeaturePlot_scCustom(alldata_wounded, "Krt1")
g3_w <- FeaturePlot_scCustom(alldata_wounded, "Krt14")
g4_w <- FeaturePlot_scCustom(alldata_wounded, "Krt16")

plot_w <- g1_w + g2_w + g3_w + g4_w + plot_layout(ncol = 4)

ggsave("Krt_ExpressionRootBasalIV_NucleiOnlyUWOD4PWD7PW.png",
       plot = plot_w,
       width = 15, height = 4, dpi = 300)


# Unwounded
g1_unw <- FeaturePlot_scCustom(alldata.2, "Krt10")
g2_unw <- FeaturePlot_scCustom(alldata.2, "Krt1")
g3_unw <- FeaturePlot_scCustom(alldata.2, "Krt14")
g4_unw <- FeaturePlot_scCustom(alldata.2, "Krt16")

plot_unw <- g1_unw + g2_unw + g3_unw + g4_unw + plot_layout(ncol = 4)

ggsave("Krt_ExpressionBasalIVCells_NucleiOnlyAllTimepoints.png",
       plot = plot_unw,
       width = 15, height = 4, dpi = 300)


Idents(alldata_unwounded) <- alldata_unwounded$celltype
DimPlot_scCustom(alldata_unwounded)


Idents(alldata_wounded) <- alldata_wounded$celltype
DimPlot_scCustom(alldata_unwounded)

levels(alldata_unwounded)

# Define custom colors
celltype_colors <- c(
  "Ker.Cycling III" = "darkblue",
  "Basal II"   = "#364B9A",  # deep blue
  "Basal III"  = "#4A6FE3",  # medium blue
  "Basal IV"   = "#8DB1FF",  # light blue
  "Spinous I"  = "#D1495B",  # deep red-orange
  "Spinous II" = "#E68653",  # medium orange
  "Spinous III"= "#F6C667",  # golden orange
  "Spinous IV" = "#F9E3A3"   # pale orange
)


# Make sure Idents are set
Idents(alldata_unwounded) <- alldata_unwounded$celltype
Idents(alldata_wounded)   <- alldata_wounded$celltype

# Define colors
celltype_colors <- c(
  "Ker.Cycling III" = "darkblue",
  "Basal II"   = "#364B9A",  # deep blue
  "Basal III"  = "#4A6FE3",  # medium blue
  "Basal IV"   = "#8DB1FF",  # light blue
  "Spinous I"  = "#D1495B",  # deep red-orange
  "Spinous II" = "#E68653",  # medium orange
  "Spinous III"= "#F6C667",  # golden orange
  "Spinous IV" = "#F9E3A3"   # pale orange
)

# Generate plots
p_unwounded <- DimPlot_scCustom(
  alldata_unwounded,
  
  group.by = "celltype",
  colors_use = celltype_colors,
  label = F,
  pt.size = 1.2
) + ggtitle("Unwounded")

p_wounded <- DimPlot_scCustom(
  alldata_wounded,
  group.by = "celltype",
  colors_use = celltype_colors,
  label = F,pt.size = 1.2
  
) + ggtitle("Wounded")

# Save as PDF
ggsave("UnwoundedIFEPsuedotime_DimPlotYvon.pdf", p_unwounded, width = 6, height = 5)
ggsave("WoundedIFEPsuedotime_DimPlotYvon.pdf", p_wounded, width = 6, height = 5)

# Save as SVG
ggsave("UnwoundedIFEPsuedotime_DimPlotYvon.svg", p_unwounded, width = 6, height = 5)
ggsave("WoundedIFEPsuedotime_DimPlotYvon.svg", p_wounded, width = 6, height = 5)

saveRDS(alldata_unwounded, "IFEwithKerCyclingIIIUnwounded.rds")
saveRDS(alldata_wounded, "IFEwithKerCyclingIIIUnwoundedD4PWD7PW.rds")

alldata_unwounded = readRDS("IFEwithKerCyclingIIIUnwounded.rds")
alldata_wounded = readRDS("IFEwithKerCyclingIIIUnwoundedD4PWD7PW.rds")

```


```{r Pseudtome Chunk For Module Scoring and Gene Ontology}
library(dplyr)
library(tibble)

plot_cells(pt_wounded$cds,
           color_cells_by = "pseudotime",
           label_branch_points = TRUE,
           label_groups_by_cluster = FALSE)

selected_cells_wo <-choose_graph_segments(pt_wounded$cds,return_list = T)

#selected_cells_wo <- choose_cells(pt_wounded$cds, return_list = T)

branch_genes <- c("Krt1")

rowData(pt_wounded$cds)$gene_short_name <- rownames(pt_wounded$cds)

# Subset CDS to genes and selected cells
branch_lineage_cds <- pt_wounded$cds[
  rowData(pt_wounded$cds)$gene_short_name %in% branch_genes,
  colnames(pt_wounded$cds) %in% selected_cells_wo
]


plot_genes_in_pseudotime(branch_lineage_cds,
                         color_cells_by = "timepoint",
                         min_expr = 0.5,
                         label_by_short_name = TRUE)

seurat_branch <- subset(alldata_wounded, cells = selected_cells_wo)

library(Seurat)
library(dplyr)

# Example: get gene lists per module
module_gene_lists <- gene_modules_wounded %>%
  group_by(module) %>%
  summarise(genes = list(id)) %>%
  deframe()

for (mod_name in names(module_gene_lists)) {
  genes <- module_gene_lists[[mod_name]]
  
  # Ensure genes are present in your data
  genes_present <- intersect(genes, rownames(seurat_branch))
  
  # Add module score
  seurat_branch <- AddModuleScore(seurat_branch,
                                  features = list(genes_present),
                                  name = paste0("ModuleScore_", mod_name))
}

pseudotime_vals <- pt_wounded$cds@principal_graph_aux[["UMAP"]]$pseudotime

# Convert to a data frame with cell names
pseudotime_df <- data.frame(
  cell = names(pseudotime_vals),
  pseudotime = pseudotime_vals,
  stringsAsFactors = FALSE
)

pseudotime_df_branch <- pseudotime_df %>%
  filter(cell %in% Cells(seurat_branch))

mod_score_cols <- grep("^ModuleScore_", colnames(seurat_branch@meta.data), value = TRUE)

library(stringr)  # for str_remove

pseudotime_df_branch <- pseudotime_df_branch %>%
  left_join(
    seurat_branch@meta.data %>%
      rownames_to_column("cell") %>%
      select(cell, timepoint),
    by = "cell"
  )


library(stringr)

module_scores_df <- seurat_branch@meta.data %>%
  rownames_to_column("cell") %>%
  select(cell, all_of(mod_score_cols)) %>%
  inner_join(pseudotime_df_branch, by = "cell") %>%
  pivot_longer(
    cols = all_of(mod_score_cols),
    names_to = "module",
    values_to = "score"
  ) %>%
  # Remove only the "ModuleScore_" prefix but keep the numeric suffix
  mutate(module = str_remove(module, "^ModuleScore_"))



p <- ggplot(module_scores_df, aes(x = pseudotime, y = score, color = timepoint)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  facet_wrap(~ module, scales = "fixed") +   # fixed y-axis scales for all facets
  theme_ggprism_mod(base_size = 14) +
  labs(
    title = "Module Scores Across Pseudotime in Selected Branch",
    x = "Pseudotime",
    y = "Module Score",
    color = "Timepoint"
  ) +
  theme(
    legend.position = "bottom"
  )

# Save the plot
ggsave("Module_scores_WoundedUnwoundedBranchPseudotime.png", plot = p, width = 8, height = 8, dpi = 300)

#Per Module



library(dplyr)
library(tibble)

# Group genes by module
module_gene_lists <- gene_modules_wounded %>%
  group_by(module) %>%
  summarise(gene_ids = list(id)) %>%
  deframe()


library(org.Mm.eg.db)  # For mouse; use org.Hs.eg.db for human
library(clusterProfiler)

# Function to map symbols to Entrez
map_to_entrez <- function(symbols) {
  bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) %>%
    distinct(SYMBOL, .keep_all = TRUE)
}

go_results <- lapply(names(module_gene_lists), function(mod) {
  genes <- module_gene_lists[[mod]]
  entrez <- map_to_entrez(genes)
  
  ego <- enrichGO(
    gene         = entrez$ENTREZID,
    OrgDb        = org.Mm.eg.db,
    ont          = "BP",  # Biological Process
    keyType      = "ENTREZID",
    pAdjustMethod= "BH",
    pvalueCutoff = 0.05,
    readable     = TRUE
  )
  
  ego@result$module <- mod  # Tag results with module
  return(ego@result)
})

# Combine all modules
go_df <- do.call(rbind, go_results)
go_df_filtered <- go_df %>%
  filter(Count >= 5)
write.csv(go_df_filtered, file = "GO_enrichment_all_modules.csv", row.names = FALSE)



library(dplyr)
library(forcats)
library(ggplot2)
library(stringr)

# Step 1: Filter top N GO terms per module (e.g., top 5)
go_top_df <- go_df %>%
  group_by(module) %>%
  slice_min(order_by = p.adjust, n = 5) %>%
  ungroup()

# Step 2: Wrap long GO terms and create sorted factor levels
go_top_df$Description <- str_wrap(go_top_df$Description, width = 50)

# Reorder GO terms by module and significance (so terms are grouped and sorted within module)
go_top_df <- go_top_df %>%
  mutate(Description = fct_reorder2(Description, module, -p.adjust))

# Step 3: Create dotplot
ggplot(go_top_df, aes(x = -log10(p.adjust), y = fct_reorder(Description, -p.adjust))) +
  geom_point(aes(size = Count), color = "steelblue") +
  facet_wrap(~ module, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Top GO Terms per Module",
    x = "-log10(adj. p-value)",
    y = NULL,
    size = "Gene Count"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )



```

```{r Pseudotime Monocle Differential Expression, message=FALSE, warning=FALSE}
library(monocle3)
library(dplyr)
library(ggplot2)
library(purrr)
library(Seurat)
library(tibble)

# -----------------------
# 1. Select a branch manually
# -----------------------
selected_cells <- choose_graph_segments(pt_wounded$cds, return_list = TRUE)

# Subset to branch
branch_cds <- pt_wounded$cds[, colnames(pt_wounded$cds) %in% selected_cells$cells]

# Ensure gene_short_name exists
rowData(branch_cds)$gene_short_name <- rownames(branch_cds)

# -----------------------
# 2. Run graph_test for DE along pseudotime
# -----------------------
subset_test_res <- graph_test(branch_cds, neighbor_graph = "principal_graph", cores = 8)

deg_ids <- subset_test_res %>%
  filter(q_value < 0.05) %>%
  arrange(desc(morans_I)) %>%
  rownames()

# Remove ribosomal genes
filtered_deg_ids <- deg_ids[!grepl("^Rpl|^Rps", deg_ids)]

# Top genes (or set manually)
pseudo_genes <- head(filtered_deg_ids[filtered_deg_ids != "Gm42418"], 20)
# pseudo_genes <- c("Plac8","Ifitm6","Pf4","Clec4e","C1qa","Arg1","Pcolce2")

# -----------------------
# 3. Build pseudotime df
# -----------------------
pseudotime_vals <- pseudotime(pt_wounded$cds) # safer than slotting directly

pt_df <- tibble(
  cell = names(pseudotime_vals),
  pseudotime = pseudotime_vals
) %>%
  filter(cell %in% selected_cells$cells)

# -----------------------
# 4. Join pseudotime with Seurat expression
# -----------------------
genes_to_plot <- c("Sema3c","Krt6a","Sox6")  # define your genes

seurat_branch <- subset(alldata_wounded, cells = selected_cells$cells)

plot_df_all <- map_dfr(genes_to_plot, function(g) {
  FetchData(
    seurat_branch,
    vars = c(g, "timepoint"),
    assay = "SCT",
    slot = "data"
  ) %>%
    rownames_to_column("cell") %>%
    # rename the single gene column to "expression"
    rename_with(~ "expression", all_of(g)) %>%
    inner_join(pt_df, by = "cell") %>%
    mutate(gene = g)
})


# -----------------------
# 5. Plot expression by pseudotime and timepoint
# -----------------------
p1 <- ggplot(plot_df_all, aes(x = pseudotime, y = expression, color = timepoint)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_ggprism_mod() +
  scale_color_manual(
    values = c(
      UW    = "#00BFC4",
      D1PW  = "#F8766D",
      D2PW  = "#7CAE00",   # replaced "" with green
      D4PW  = "#E7861B",
      D7PW  = "#009F75",
      D15PW = "#619CFF",
      D30PW = "grey40"
    )
  ) +
  labs(
    title = "Pseudotime Genes",
    x = "Pseudotime",
    y = "SCT-normalized expression",
    color = "Timepoint"
  ) +
  theme(legend.position = "bottom")

ggsave("PsuedotimeWoundedBranch_TimepointColoredRedBranchYvon.svg",
       p1, width = 8, height = 4, dpi = 300)
ggsave("PsuedotimeWoundedBranch_TimepointColoredRedBranchYvon.png",
       p1, width = 8, height = 4, dpi = 300)

# -----------------------
# 6. Plot expression by pseudotime with pseudotime color scale
# -----------------------
p2 <- ggplot(plot_df_all, aes(x = pseudotime, y = expression, color = pseudotime)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "black") +
  facet_wrap(~ gene, scales = "free_y") +
  theme_ggprism_mod(base_size = 14) +
  scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
  labs(
    title = "Gene Expression Across Pseudotime",
    x = "Pseudotime",
    y = "SCT-normalized expression"
  ) +
  theme(legend.position = "bottom")

ggsave("PseudotimeGenes_BranchColoredByPseudotime.svg",
       p2, width = 8, height = 4, dpi = 300)
```

```{r Psuedotime Unwounded Time}
library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(purrr)
library(tidyr)
library(viridis)

# =========================
# 1. Select branch cells
# =========================
selected_cells <- choose_graph_segments(pt_unwounded$cds, return_list = TRUE)

branch_cds <- pt_unwounded$cds[, colnames(pt_unwounded$cds) %in% selected_cells$cells]
rowData(branch_cds)$gene_short_name <- rownames(branch_cds)

# =========================
# 2. Graph test along pseudotime
# =========================
subset_test_res <- graph_test(branch_cds, neighbor_graph = "principal_graph", cores = 8)

deg_ids <- subset_test_res %>%
  filter(q_value < 0.05) %>%
  arrange(desc(morans_I)) %>%
  rownames()

filtered_deg_ids <- deg_ids[!grepl("^Rpl|^Rps", deg_ids)]

pseudo_genes <- head(filtered_deg_ids[filtered_deg_ids != "Gm42418"], 20)
# Or define manually:
# pseudo_genes <- c("Plac8","Ifitm6","Pf4","Clec4e","C1qa","Arg1","Pcolce2")

genes_to_plot <- c("Sema3c","Krt6a","Sox6")

# =========================
# 3. Build pseudotime dataframe
# =========================
pseudotime_vals <- pseudotime(pt_unwounded$cds)

pt_df <- tibble(
  cell = names(pseudotime_vals),
  pseudotime = pseudotime_vals
) %>%
  filter(cell %in% selected_cells$cells)

# =========================
# 4. Fetch gene expression + metadata
# =========================
seurat_branch <- subset(alldata_unwounded, cells = selected_cells$cells)

plot_df_all <- FetchData(
  seurat_branch,
  vars = c(genes_to_plot, "timepoint"),
  assay = "SCT",
  slot = "data"
) %>%
  rownames_to_column("cell") %>%
  inner_join(pt_df, by = "cell") %>%
  pivot_longer(cols = all_of(genes_to_plot),
               names_to = "gene",
               values_to = "expression")

# =========================
# 5. Plot: expression by pseudotime + timepoint
# =========================
p1 <- ggplot(plot_df_all, aes(x = pseudotime, y = expression, color = timepoint)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_ggprism_mod() +
  scale_color_manual(
    values = c(
      UW    = "#00BFC4",
      D1PW  = "#F8766D",
      D2PW  = "#7CAE00",   # fixed (was "")
      D4PW  = "#E7861B",
      D7PW  = "#009F75",
      D15PW = "#619CFF",
      D30PW = "grey40"
    )
  ) +
  labs(
    title = "Pseudotime Genes (Unwounded branch)",
    x = "Pseudotime",
    y = "SCT-normalized expression",
    color = "Timepoint"
  ) +
  theme(legend.position = "bottom")

ggsave("Pseudotime_UnwoundedBranch_TimepointColored.svg",
       p1, width = 8, height = 4, dpi = 300)
ggsave("Pseudotime_UnwoundedBranch_TimepointColored.png",
       p1, width = 8, height = 4, dpi = 300)


# =========================
# 6. Plot: expression by pseudotime gradient
# =========================
p2 <- ggplot(plot_df_all, aes(x = pseudotime, y = expression, color = pseudotime)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "black") +
  facet_wrap(~ gene, scales = "free_y") +
  theme_ggprism_mod(base_size = 14) +
  scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
  labs(
    title = "Gene Expression Across Pseudotime (Unwounded branch)",
    x = "Pseudotime",
    y = "SCT-normalized expression"
  ) +
  theme(legend.position = "bottom")

ggsave("Pseudotime_UnwoundedBranch_PseudotimeColored.svg",
       p2, width = 8, height = 4, dpi = 300)
```