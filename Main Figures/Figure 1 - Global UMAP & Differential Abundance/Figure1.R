---
  
  title: "Figure 2 – UMAPs, QC, and Cell Composition"
author: "Jonathan Chin Cheong"
date: "2024-10-21"
output: html_document
---------------------
  
  ## 0. Global Setup
  
  set.seed(123)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(
  root.dir = "/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 1-2"
)

## 1. Load Libraries

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Rfast2)
library(viridis)
library(scCustomize)
library(harmony)
library(nichenetr)
library(khroma)
library(colorspace)
library(ggpattern)
library(clusterProfiler)
library(org.Mm.eg.db)

## 2. Load Integrated Dataset

alldata <- readRDS(
  "/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Sequencing Repository/Data/OWHA RDS files/IntegratedMultimodal_FineCellTypes_061225.rds"
)

## 3. Define Cell Type Hierarchy

### 3.1 Metacluster → Fine Cell Types

metacluster_groups <- list(
  "IFE Keratinocytes" = c("Basal II","Basal III","Basal IV","uHF I","uHF II",
                          "Spinous I","Spinous II","Spinous III","Spinous IV","Cornified",
                          "Ker.Cycling I","Ker.Cycling II","Ker.Cycling III"),
  "HF Keratinocytes"  = c("IRS I","IRS II","Medulla","Cortex I","Cortex II",
                          "Outer Bulge I","Outer Bulge II","Outer Bulge III","HF I","HF II"),
  "Sebocytes"         = c("Transitional Basal Cells","Sebocyte II","Sebocyte I","Sebocyte III"),
  "Fibroblasts"       = c("Papillary I","Papillary II","Papillary III","Papillary IV",
                          "Reticular I","Reticular II","Reticular III","Reticular IV","Reticular V",
                          "Fascia","Myofibroblasts I","Myofibroblasts II","Myofibroblasts III",
                          "Lef1+ Fibroblasts","Tnc+ Fibroblasts","Myoc+ Fibroblasts",
                          "Lox+ Fibroblasts","Nfkb+ Fibroblasts"),
  "Immune cells"      = c("Monocytes I","Monocytes II","Monocytes III",
                          "Macrophages I","Macrophages II","Macrophages III","Macrophages IV",
                          "Cycling Macrophages","OC-like cells","moDCs","cDC1","cDC2",
                          "Migratory DCs","Langerhans cells","pDCs","Neutrophils I",
                          "Neutrophils II","Neutrophils III","Basophils","Mast cells",
                          "NK cells","γδT cells","T cells","Tregs","ILCs","B cells"),
  "Endothelial cells" = c("Tacr1+ ECs","Capillary ECs","Mature ECs","LEC I","LEC II",
                          "Artery ECs","Proliferating ECs","Progenitor ECs","EC I","EC II",
                          "Vein ECs","PVMs","Proliferating LECs","LEC III","EC III","VWF-hi ECs"),
  "Pericytes"         = c("Pericytes","VSMCs","Pericyte I","Pericyte II","Pericyte III","Pericyte IV"),
  "Muscle cells"      = c("Muscle Progenitors","Myonuclei","Skeletal Muscle"),
  "Schwann cells"     = c("Myelinating SC I","Myelinating SC II","Repair SC I","Repair SC II","Repair SC III"),
  "Melanocytes"       = c("Oca2- Melanocytes","Oca2-lo Melanocytes","Oca2+ Melanocytes"),
  "Adipocytes"        = c("Adipocyte I","Adipocyte II"),
  "Red Blood Cells"   = c("Red Blood Cells")
)

### 3.2 Metacluster Base Colors

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
  "Red Blood Cells"   = "#A7222B"
)

## 4. Generate Fine-Grained Color Palette

fine_colors <- list()
for (meta in names(metacluster_groups)) {
  pal <- colorway::build_palette(
    central_color = metacluster_colors[[meta]],
    n_colors      = length(metacluster_groups[[meta]]),
    hue_range     = 0.2,
    val_range     = 0.4,
    sat_range     = -0.3
  )
  names(pal$palette) <- metacluster_groups[[meta]]
  fine_colors[[meta]] <- pal$palette
}

color_palette <- unlist(fine_colors)
names(color_palette) <- sub("^[^.]+\.", "", names(color_palette))

# Manual refinements (fibroblasts, pericytes, adipocytes, Schwann cells)

# -- kept explicit for figure reproducibility

color_palette[c(
  "Skeletal Muscle"=" #FF6600",
  "Papillary I"="#A66C8F","Papillary II"="#B178A8","Papillary III"="#A04B72","Papillary IV"="#8C2B5D",
  "Reticular I"="#BD5B99","Reticular II"="#D1647A","Reticular III"="#E07682",
  "Reticular IV"="#9B3040","Reticular V"="#6B1E28",
  "Fascia"="#E37471",
  "Myofibroblasts I"="#A74337","Myofibroblasts II"="#B85D4A","Myofibroblasts III"="#FFAC90",
  "Lef1+ Fibroblasts"="#BC6C45","Tnc+ Fibroblasts"="#6A4F2A",
  "Myoc+ Fibroblasts"="#E3A571","Lox+ Fibroblasts"="#D49F63","Nfkb+ Fibroblasts"="#A97132",
  "Pericytes"="#A7DD0B","VSMCs"="#A2E813","Pericyte I"="#8669B8","Pericyte II"="#5E3490",
  "Pericyte III"="#3A1A65","Pericyte IV"="#B8A3DA",
  "Adipocyte I"="#FF1493","Adipocyte II"="#FF69B4",
  "Myelinating SC I"="#9C27B0","Myelinating SC II"="#BE4FD9",
  "Repair SC I"="#D18EFF","Repair SC II"="#E6B3FF","Repair SC III"="#A63FCC"
)] <- color_palette[names(color_palette)]

saveRDS(color_palette, "color_palette.rds")

## 5. Global UMAPs (Figure 1–2)

Idents(alldata) <- "celltype_fine"
DimPlot(alldata, cols = color_palette, reduction = "umap.rpca") + custom_theme

ggsave("UMAPwithLegend.png", width = 15, height = 11)

ggsave(
  "UMAPNoLegend.png",
  DimPlot(alldata, cols = color_palette, reduction = "umap.rpca", raster = FALSE) &
    NoLegend() & custom_theme,
  width = 11, height = 11
)

## 6. Condition, Modality, and QC Annotations

alldata$condition <- ifelse(alldata$timepoint == "UW", "Unwounded", "Wounded")
alldata$modality_group <- ifelse(alldata$modality %in% c("scRNAseq","CITEseq"), "whole cell", "nuclei")
levels(alldata$timepoint) <- c("UWO","D1PW","D2PW","D4PW","D7PW","D15PW","D30PW")

## 7. Split & Downsampled UMAPs

target_n <- 43686
set.seed(1234)

downsampled_cells <- alldata@meta.data %>%
  mutate(cell_id = rownames(.)) %>%
  group_by(condition) %>%
  slice_sample(n = target_n) %>%
  pull(cell_id)

alldata_ds <- subset(alldata, cells = downsampled_cells)

p <- DimPlot(
  alldata_ds,
  cols = color_palette,
  reduction = "umap.rpca",
  split.by = "condition",
  raster = FALSE
) & NoLegend() + NoAxes() + custom_theme

ggsave("UMAP_SplitCondition_Downsampled.png", p, width = 12, height = 6)

## 8. Cell Composition Over Time and Modality

*(Bar plots and stacked proportions; see original logic preserved below)*
  
  cols <- metacluster_colors
cells <- as.data.frame(table(alldata$metaclusters, alldata$timepoint))

ggplot(cells, aes(Var2, Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual(values = cols) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## 9. Longitudinal Line Plots (Counts & Proportions)

metadata <- alldata@meta.data

result_df <- metadata %>%
  group_by(orig.ident, timepoint, metaclusters) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(TotalCount = sum(Count), Proportion = Count / TotalCount) %>%
  ungroup() %>%
  group_by(timepoint, metaclusters) %>%
  summarise(
    AverageCount = mean(Count),
    Count_SD = sd(Count),
    Percentage = mean(Proportion) * 100,
    SD_Percentage = sd(Proportion) * 100,
    .groups = "drop"
  )

## 10. Differential Expression & GO Enrichment

cell_types <- c(
  "IFE Keratinocytes","HF Keratinocytes","Immune cells","Sebocytes",
  "Fibroblasts","Schwann cells","Endothelial cells","Pericytes",
  "Muscle cells","Melanocytes","Adipocytes","Red Blood Cells"
)

for (ct in cell_types) {
  Idents(alldata) <- "metaclusters"
  alldata_ct <- subset(alldata, ident = ct)
  
  de <- DElegate::findDE(
    alldata_ct,
    compare = c("Wounded_D4PW","Unwounded"),
    method = "deseq",
    group_column = "timepoint",
    replicate_column = "orig.ident"
  )
  
  sig <- subset(de, padj < 0.01 & log_fc > 1)
  ids <- bitr(sig$feature, "SYMBOL", "ENTREZID", org.Mm.eg.db)
  
  ego <- enrichGO(ids$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP")
  print(dotplot(ego, showCategory = 5) + ggtitle(ct))
}


