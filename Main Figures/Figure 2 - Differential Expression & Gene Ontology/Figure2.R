---
  title: "Figure 2 UMAP Making"
author: "Jonathan Chin Cheong"
date: "2024-10-21"
output: html_document
---
  
  ```{r}
set.seed(123)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 1-2")
```

```{r pressure, echo=FALSE}
#Load Libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Rfast2)
library(viridis)
library(scCustomize)
library(harmony)
library(khroma)
```
```{r Find Markers Heatmap}

# Load libraries
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)  # ⬅️ Mouse annotation database
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(khroma)

# Set identity and assay
Idents(alldata) <- "timepoint"
DefaultAssay(alldata) <- "RNA"
alldata<- JoinLayers(alldata)

# Recode and factor timepoints
alldata@meta.data$timepoint <- recode(alldata@meta.data$timepoint,
                                      "Unwounded" = "UW",
                                      "Wounded_D1PW" = "D1PW",
                                      "Wounded_D2PW" = "D2PW",
                                      "Wounded_D4PW" = "D4PW",
                                      "Wounded_D7PW" = "D7PW",
                                      "Wounded_D15PW" = "D15PW",
                                      "Wounded_D30PW" = "D30PW"
)
alldata@meta.data$timepoint <- factor(alldata@meta.data$timepoint,
                                      levels = c("UW", "D1PW", "D2PW", "D4PW", "D7PW", "D15PW", "D30PW"))

# Find markers
cellmarkers_mouse <- FindAllMarkers(alldata, min.pct = 0.30, recorrect_umi = FALSE, only.pos = TRUE)

# Subset DE genes per timepoint
timepoints <- levels(alldata@meta.data$timepoint)
de_genes_list <- lapply(timepoints, function(tp) {
  subset(cellmarkers_mouse, cluster == tp & avg_log2FC > 1)
})
names(de_genes_list) <- timepoints

# Convert to ENTREZ (mouse)
gene_lists_entrez <- lapply(de_genes_list, function(df) {
  bitr(df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
})

# GO enrichment
ck <- compareCluster(geneCluster = gene_lists_entrez, fun = enrichGO, OrgDb = org.Mm.eg.db,
                     pvalueCutoff = 0.05, ont = "BP")
ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
ck <- simplify(ck)

# Unique pathways
ck_df <- as.data.frame(ck)
unique_df <- ck_df %>% group_by(ID) %>% filter(n_distinct(Cluster) == 1) %>% ungroup()
ck@compareClusterResult <- unique_df


# Dotplot
s1 <- dotplot(ck, showCategory = 8, font.size = 12, label_format = 100) +
  scale_color_gradientn(colors = c("blue", "white", "red")) +
  theme(plot.margin = grid::unit(c(0, 0, 0, 0), "mm")) +
  theme(axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, face = "bold", color = "black"))

s1 + customtheme
ggsave("Mouse_GeneOntology_GlobalPathways_AllTimepoints.png", width = 10, height = 7)

# Save GO results
write.csv(as.data.frame(ck), file = "Mouse_GlobalGOcompareCluster_results.csv", row.names = FALSE)
ck <- read.csv("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 2/Mouse_GlobalGOcompareCluster_results.csv")

library(ggplot2)
library(dplyr)
library(khroma)

# ---- Define Nightfall colors for timepoints ----
timepoint_colors <- color("nightfall")(7)
names(timepoint_colors) <- c("UW", "D1PW", "D2PW", "D4PW", "D7PW", "D15PW", "D30PW")

# ---- Ensure Cluster is a factor in correct order ----
timepoint_order <- c("UW", "D1PW", "D2PW", "D4PW", "D7PW", "D15PW", "D30PW")
ck_df$Cluster <- factor(ck_df$Cluster, levels = timepoint_order)

# ---- Optional: Keep only unique pathways (appear in 1 cluster) ----
unique_df <- ck_df %>%
  group_by(ID) %>%
  filter(n_distinct(Cluster) == 1) %>%
  ungroup()

# ---- Select top 5 GO terms per cluster by adjusted p-value ----
top_terms <- unique_df %>%
  group_by(Cluster) %>%
  slice_min(order_by = p.adjust, n = 5) %>%
  ungroup()

# ---- Create bar plot ----
ggplot(top_terms, aes(x = zScore, y = reorder(Description, zScore), fill = Cluster)) +
  geom_col() +
  scale_fill_manual(values = timepoint_colors) +
  facet_wrap(~ Cluster, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12)
  ) +
  labs(
    x = "-log10(adj. p-value)",
    y = "GO Term",
    fill = "Timepoint",
    title = "Top Enriched GO Terms per Timepoint"
  ) + theme_ggprism_mod()

# ---- Save figure ----
ggsave("Mouse_GO_TopTerms_Barplot_Nightfall.svg", width = 12, height = 12)



# Heatmap section
# Top 30 genes per timepoint
top_genes_list <- lapply(names(de_genes_list), function(tp) {
  top <- head(de_genes_list[[tp]][order(-de_genes_list[[tp]]$avg_log2FC), ], 50)
  transform(top, condition = tp)
})
top.genes <- do.call(rbind, top_genes_list)
mouse_genes <- unique(top.genes$gene)
mouse_genes <- setdiff(mouse_genes, "Gm42418")

# Scale and extract matrix
alldata <- ScaleData(alldata, features = mouse_genes)
mouse_genes <- intersect(mouse_genes, rownames(alldata))
matDE1 <- GetAssayData(alldata, slot = "scale.data")[mouse_genes, ]

# Order meta
meta <- alldata@meta.data[colnames(matDE1), ]
meta$order <- paste(meta$timepoint, meta$metaclusters, sep = "_")
ordered_cells <- rownames(meta[order(meta$timepoint, meta$metaclusters), ])
matDE1 <- matDE1[, ordered_cells, drop = FALSE]

# Annotations
metacluster_colors <- c(
  "IFE Keratinocytes" = "#364B9A",
  "HF Keratinocytes"  = "#56B4E9",
  "Sebocytes"         = "#8DD3C7",
  "Fibroblasts"       = "indianred2",
  "Immune cells"      = "#198F1B",
  "Endothelial cells" = "#F6C800",
  "Pericytes"         = "#D12E82",
  "Muscle cells"            = "orange2",
  "Schwann cells"     = "purple",
  "Melanocytes"       = "pink1",
  "Adipocytes"        = "deeppink",
  "Red Blood Cells"   = "#A7222B"
)
timepoint_colors <- color("nightfall")(7)
names(timepoint_colors) <- c("UW", "D1PW", "D2PW", "D4PW", "D7PW", "D15PW", "D30PW")

col_annot <- HeatmapAnnotation(
  df = meta[ordered_cells, c("metaclusters", "timepoint")],
  col = list(
    metaclusters = metacluster_colors,
    timepoint = timepoint_colors
  )
)

col_fun <- colorRamp2(c(-1, 0, 2), c("blue", "white", "red"))

# Save heatmap
png(file = "Mouse_HeatmapByTimepoint-Metacluster_TopGenes.png",
    width = 12, height = 15, units = 'in', bg = "transparent", res = 300)
hm <- Heatmap(
  matDE1,
  name = "logFC\n(Z-score)",
  top_annotation = col_annot,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  row_names_max_width = unit(10, "cm"),
  show_column_names = FALSE,
  show_heatmap_legend = TRUE,
  col = col_fun,
  column_names_rot = 45,
  column_gap = unit(5, 'mm')
)
draw(hm, merge_legend = TRUE)
dev.off()


```




```{r Global Gene Ontology Across Timepoints}

alldata = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/scRNA SEQ of wounds/Skin Wound Healing Atlas/IntegratedMultimodal_FineCellTypes_061225.rds")
colors_pals = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 1/color_palette.rds")

alldata$celltype <- Idents(alldata)
# Update timepoint names
alldata@meta.data$timepoint <- recode(alldata@meta.data$timepoint,
                                      "Unwounded" = "UW",
                                      "Wounded_D1PW" = "D1PW",
                                      "Wounded_D2PW" = "D2PW",
                                      "Wounded_D4PW" = "D4PW",
                                      "Wounded_D7PW" = "D7PW",
                                      "Wounded_D15PW" = "D15PW",
                                      "Wounded_D30PW" = "D30PW"
)

# ---- Libraries ----
library(Seurat)
library(DElegate)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(khroma)
library(CellChat)

# ---- DE + GO on full alldata ----
run_DE_GO_analysis <- function(alldata) {
  
  # ---- Run DE for pairwise comparisons ----
  de_results.04 <- DElegate::findDE(
    object = alldata,
    compare = c("D4PW", "UW"),
    method = "deseq",
    group_column = "timepoint",
    replicate_column = "orig.ident"
  )
  
  de_results.47 <- DElegate::findDE(
    object = alldata,
    compare = c("D7PW", "D4PW"),
    method = "deseq",
    group_column = "timepoint"
  )
  
  de_results.715 <- DElegate::findDE(
    object = alldata,
    compare = c("D15PW", "D7PW"),
    method = "deseq",
    group_column = "timepoint"
  )
  
  de_results.1530 <- DElegate::findDE(
    object = alldata,
    compare = c("D30PW", "D15PW"),
    method = "deseq",
    group_column = "timepoint"
  )
  
  # ---- Filter DE results ----
  x <- 0     # logFC threshold
  y <- 0.01   # padj threshold
  
  res.04   <- subset(de_results.04, padj < y & log_fc > x)
  res.47   <- subset(de_results.47, padj < y & log_fc > x)
  res.715  <- subset(de_results.715, padj < y & log_fc > x)
  res.1530 <- subset(de_results.1530, padj < y & log_fc > x)
  
  # ---- Map SYMBOL to ENTREZ ----
  ids.04   <- bitr(res.04$feature, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  ids.47   <- bitr(res.47$feature, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  ids.715  <- bitr(res.715$feature, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  ids.1530 <- bitr(res.1530$feature, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  genelist <- list(
    D4PW_vs_UW      = ids.04$ENTREZID,
    D7PW_vs_D4PW    = ids.47$ENTREZID,
    D15PW_vs_D7PW   = ids.715$ENTREZID,
    D30PW_vs_D15PW  = ids.1530$ENTREZID
  )
  
  # ---- GO enrichment ----
  ck <- compareCluster(
    geneCluster = genelist,
    fun = enrichGO,
    OrgDb = org.Mm.eg.db,
    pvalueCutoff = 0.05,
    ont = "BP"
  )
  ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  
  ck_df <- as.data.frame(ck)
  
  # ---- Keep unique pathways ----
  unique_df <- ck_df %>%
    group_by(ID) %>%
    filter(n_distinct(Cluster) == 2) %>%
    ungroup()
  
  # ---- Select top 5 per comparison ----
  top_terms <- unique_df %>%
    group_by(Cluster) %>%
    slice_min(order_by = p.adjust, n = 8) %>%
    ungroup()
  
  # ---- Order comparisons ----
  comp_order <- c("D4PW_vs_UW", "D7PW_vs_D4PW", "D15PW_vs_D7PW", "D30PW_vs_D15PW")
  top_terms$Cluster <- factor(top_terms$Cluster, levels = comp_order)
  
  # ---- Colors (nightfall palette) ----
  comp_colors <- color("nightfall")(length(comp_order))
  names(comp_colors) <- comp_order
  
  # ---- Plot by Z-score ----
  p <- ggplot(top_terms,
              aes(x = zScore,
                  y = reorder(Description, zScore),
                  fill = Cluster)) +
    geom_col() +
    scale_fill_manual(values = comp_colors) +
    facet_wrap(~ Cluster, scales = "free_y", ncol = 1) +
    theme_ggprism_mod(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 12),
      legend.position = "none"
    ) +
    labs(
      x = "Z-score",
      y = "GO Term",
      title = ""
    )
  
  # ---- Save outputs ----
  ggsave("GO_TopTerms_Gloabal_Barplot_NightfallNew.svg",
         p, width = 10, height = 10, dpi = 300)
  write.csv(ck_df, "GO_Global_resultsNew.csv", row.names = FALSE)
  
  print(p)
  
}

# ---- Run ----
options(future.globals.maxSize = 8000 * 1024^2)
run_DE_GO_analysis(alldata)




p <- ggplot(top_terms,
            aes(x = zScore,
                y = reorder(Description, zScore),
                fill = Cluster)) +
  geom_col() +
  scale_fill_manual(values = comp_colors) +
  facet_wrap(~ Cluster, scales = "free_y", ncol = 1) +
  theme_ggprism_mod(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(
    x = "Z-score",
    y = "GO Term",
    title = "Top Enriched GO Terms across Timepoints (Keratinocytes)"
  )
```
```{r}
# ---- Libraries ----
# ========================================
# Differential Expression Heatmap (Grouped by Peak Timepoint)
# ========================================

# ---- Libraries ----
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(khroma)
library(ggplot2)
library(grid)

# ---- Create a named list of DE results ----
de_results_list <- list(
  D4_vs_UW   = res.04,
  D7_vs_D4   = res.47,
  D15_vs_D7  = res.715,
  D30_vs_D15 = res.1530
)

# ---- Define Nightfall colors for timepoints ----
timepoint_colors <- color("nightfall")(7)
names(timepoint_colors) <- c("UW", "D1PW", "D2PW", "D4PW", "D7PW", "D15PW", "D30PW")

# ---- Define metacluster colors ----
metacluster_colors <- c(
  "IFE Keratinocytes" = "#364B9A",
  "HF Keratinocytes"  = "#56B4E9",
  "Sebocytes"         = "#8DD3C7",
  "Fibroblasts"       = "#F95B3C",
  "Immune cells"      = "#198F1B",
  "Endothelial cells" = "#F6C800",
  "Pericytes"         = "#D12E82",
  "Muscle cells"      = "#EC6F2D",
  "Schwann cells"     = "#7A4FCF",
  "Melanocytes"       = "#FA5B77",
  "Adipocytes"        = "#19A1D1",
  "Red Blood Cells"   = "#A7222B"
)

# ---- Prepare top genes list, filter by logFC and rate thresholds ----
logFC_min <- -10
logFC_max <- 10

top_genes_list <- lapply(names(de_results_list), function(tp) {
  df <- de_results_list[[tp]] %>%
    filter(
      log_fc > logFC_min & log_fc <= logFC_max,
      (rate1 > 0.15 | rate2 > 0.15)   # rate filter
    ) %>%
    arrange(-log_fc) %>%
    head(200)
  df$timepoint <- tp
  return(df)
})

top_genes_df <- do.call(rbind, top_genes_list)

# Remove unwanted genes
genes_to_plot <- setdiff(unique(top_genes_df$feature), "Gm42418")

# ---- Scale data in Seurat ----
alldata <- ScaleData(alldata, features = genes_to_plot)
genes_to_plot <- intersect(genes_to_plot, rownames(alldata))
mat <- GetAssayData(alldata, slot = "scale.data")[genes_to_plot, ]

# ---- Order cells by timepoint and metacluster ----
meta <- alldata@meta.data[colnames(mat), ]
meta$order <- paste(meta$timepoint, meta$metaclusters, sep = "_")
ordered_cells <- rownames(meta[order(meta$timepoint, meta$metaclusters), ])
mat <- mat[, ordered_cells, drop = FALSE]

# ---- Compute average expression per timepoint ----
time_levels <- c("UW", "D1PW", "D2PW", "D4PW", "D7PW", "D15PW", "D30PW")

avg_by_time <- sapply(time_levels, function(tp) {
  cells_tp <- rownames(meta[meta$timepoint == tp, ])
  if (length(cells_tp) > 1) {
    rowMeans(mat[, cells_tp, drop = FALSE])
  } else if (length(cells_tp) == 1) {
    mat[, cells_tp, drop = FALSE]
  } else {
    # handle missing timepoints gracefully
    rep(NA, nrow(mat))
  }
})

# ---- Align row names and remove incomplete rows ----
rownames(avg_by_time) <- rownames(mat)
avg_by_time <- avg_by_time[complete.cases(avg_by_time), , drop = FALSE]

# ---- Determine peak timepoint for each gene ----
gene_peak_tp <- apply(avg_by_time, 1, function(x) names(which.max(x)))
gene_peak_tp <- factor(gene_peak_tp, levels = time_levels)

# ---- Ensure matching order ----
gene_peak_tp <- gene_peak_tp[rownames(avg_by_time)]
mat <- mat[rownames(avg_by_time), , drop = FALSE]

# ---- Reorder genes by peak timepoint and mean expression ----
gene_order <- order(gene_peak_tp, rowMeans(avg_by_time))
mat <- mat[gene_order, , drop = FALSE]
# ---- Column annotations ----
col_annot <- HeatmapAnnotation(
  df = meta[ordered_cells, c("metaclusters", "timepoint")],
  col = list(
    metaclusters = metacluster_colors,
    timepoint = timepoint_colors
  )
)

# ---- Heatmap color scale ----
col_fun <- colorRamp2(c(-1, 0, 2), c("blue", "white", "red"))

# ---- Plot heatmap ----
png(
  filename = "Mouse_Heatmap_TopGenes_GroupedByTimepoint.png",
  width = 12, height = 15, units = "in", res = 300, bg = "transparent"
)

hm <- Heatmap(
  mat,
  name = "Expression\n(Z-score)",
  top_annotation = col_annot,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  row_names_max_width = unit(10, "cm"),
  show_column_names = FALSE,
  show_heatmap_legend = TRUE,
  col = col_fun,
  column_names_rot = 45,
  column_gap = unit(5, "mm")
)

draw(hm, merge_legend = TRUE)
dev.off()

```
```{r}
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ---- Combine DE results into one matrix ----
de_results_list <- list(
  D4_vs_UW   = de_results.04,
  D7_vs_D4   = de_results.47,
  D15_vs_D7  = de_results.715,
  D30_vs_D15 = de_results.1530
)

# ---- Create gene x timepoint matrix and cap logFC ----
logFC_cap <- 5

all_genes <- unique(unlist(lapply(de_results_list, function(x) x$feature)))

logFC_mat <- sapply(names(de_results_list), function(tp) {
  df <- de_results_list[[tp]]
  vals <- setNames(df$log_fc, df$feature)
  vals[all_genes]
})
logFC_mat[is.na(logFC_mat)] <- 0
logFC_mat <- pmin(pmax(logFC_mat, -logFC_cap), logFC_cap)  # cap logFC at ±5
rownames(logFC_mat) <- all_genes

# ---- Define colors for timepoints ----
timepoint_colors <- color("nightfall")(length(names(de_results_list)))
names(timepoint_colors) <- names(de_results_list)

col_fun <- colorRamp2(c(-logFC_cap, 0, logFC_cap), c("blue", "white", "red"))

# ---- Save heatmap as PNG ----
png(filename = "Mouse_LogFC_Heatmap_Capped.png",
    width = 12, height = 15, units = "in", res = 300, bg = "transparent")

Heatmap(
  logFC_mat,
  name = "logFC",
  cluster_rows = TRUE,          # cluster genes to show temporal patterns
  cluster_columns = FALSE,      # keep timepoints in chronological order
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_side = "top",
  row_names_gp = gpar(fontsize = 6),
  clustering_method_rows = "complete"
)

dev.off()
```

