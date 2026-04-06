# title: "RNA Decontamination with DecontX"
# author: "Simon Van Deursen"
# date: "04-06-2026"
# ---- Set up ----
setwd("~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/OWHA Files/Manuscript_Multilineage dissection of wound healing crosstalks/OWHA/Cell submission/Cell Stem Cell/Revisions/Files for Revision/Data Analysis/DecontX")
# Install the decontX package if not already installed
if(!require(decontX)){
  BiocManager::install("decontX")
}
# Load necessary packages
library(Seurat)
library(celda)
library(SingleCellExperiment)
library(ggplot2)
library(decontX)
# We can check the analysis vignette here to determine how to run decontX
# vignette('decontX', package = 'decontX')

# The goal of this script is to run DecontX on our integrated multimodal atlas,
# to remove any unwanted RNA contamination from the count matrix

# ---- Load data -----
# Load dataset
alldata <- readRDS("~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Sequencing Repository/Data/OWHA Raw Files/Processed RDS Files/IntegratedMultimodal_FineCellTypes_061225.rds")

# ---- Convert to SCE object -----
DefaultAssay(alldata) <- "RNA"
alldata <- JoinLayers(alldata)
alldata.sce <- as.SingleCellExperiment(alldata, assay = "RNA")

# ---- Run DecontX -----
#run decontX on the object, specifying that clustering has already been done and the batches
alldata.sce <- decontX(alldata.sce, z = alldata.sce$rpca_clusters, batch = alldata.sce$orig.ident)
#Add the decontaminated counts back to the dataset
alldata[["decontXcounts"]] <- CreateAssayObject(counts = decontXcounts(alldata.sce))
rm(alldata.sce)

# ---- Normalize DecontX results-----
#as decontX takes the raw counts for decontamination, the new decontXcounts should be normalized and scaled
DefaultAssay(alldata) <- 'decontXcounts'
alldata[['decontXcounts']] <- split(alldata[['decontXcounts']], f = alldata$orig.ident)
alldata <- NormalizeData(alldata)
alldata <- FindVariableFeatures(alldata)
#alldata <- ScaleData(alldata)
alldata <- JoinLayers(alldata)

# ---- Plot comparison of results-----
# First we can group by the integrated clusters
grouping_var <- "rpca_clusters"
# Plot VlnPlots to compare expression of markers before and after DecontX
v1 <- VlnPlot(alldata, c("Krt14","Col1a1","Col1a2","Ptprc"), assay = "RNA", stack = T,
              flip = T, group.by = grouping_var) + labs(title = "RNA Assay")

v2 <- VlnPlot(alldata, c("Krt14","Col1a1","Col1a2","Ptprc"), assay = "SCT", stack = T,
              flip = T, group.by = grouping_var) + labs(title = "SCT Assay")

v3 <- VlnPlot(alldata, c("Krt14","Col1a1","Col1a2","Ptprc"), assay = "decontXcounts",
              stack = T, flip = T, group.by = grouping_var) + labs(title = "DecontX Assay")

v1 + v2 + v3 & NoLegend() & theme(axis.title.x = element_blank(),
                                  axis.text.x = element_text(size = 8))

# Save outputs grouped by cluster
ggsave(paste0("OWHA_DecontX_", grouping_var,"_VlnPlot.png"), width = 21, height = 7)

# Next, group by the metaclusters and fill wiht
# Save color palette
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

grouping_var <- "metaclusters"
# Plot VlnPlots to compare expression of markers before and after DecontX
v1 <- VlnPlot(alldata, c("Krt14","Col1a1","Col1a2","Ptprc"), assay = "RNA", stack = T,
              flip = T, group.by = grouping_var, fill.by = "ident",
              cols = metacluster_colors) + labs(title = "RNA Assay")

v2 <- VlnPlot(alldata, c("Krt14","Col1a1","Col1a2","Ptprc"), assay = "SCT", stack = T,
              flip = T, group.by = grouping_var, fill.by = "ident",
              cols = metacluster_colors) + labs(title = "SCT Assay")

v3 <- VlnPlot(alldata, c("Krt14","Col1a1","Col1a2","Ptprc"), assay = "decontXcounts",
              stack = T, flip = T, group.by = grouping_var, fill.by = "ident",
              cols = metacluster_colors) + labs(title = "DecontX Assay")

v1 + v2 + v3 & NoLegend() & theme(axis.title.x = element_blank())

# Save outputs grouped by celltype metacluster
ggsave(paste0("OWHA_DecontX_", grouping_var,"_VlnPlot.png"), width = 21, height = 7)

## ----- snRNA-seq only analysis-----
# Because snRNA-seq associated cell lysis can introduce additional ambient RNAs,
# check the nuclei data specifically for decontX
Idents(alldata) <- "modality"
grouping_var <- "rpca_clusters"

v1 <- VlnPlot(alldata, c("Krt14","Col1a1","Col1a2","Ptprc"),
              idents = "snRNAseq",
              assay = "RNA", 
              stack = T, flip = T, 
              group.by = grouping_var) + labs(title = "RNA Assay")

v2 <- VlnPlot(alldata, c("Krt14","Col1a1","Col1a2","Ptprc"),
              idents = "snRNAseq",
              assay = "SCT", 
              stack = T, flip = T, 
              group.by = grouping_var) + labs(title = "SCT Assay")

v3 <- VlnPlot(alldata, c("Krt14","Col1a1","Col1a2","Ptprc"),
              idents = "snRNAseq",
              assay = "decontXcounts", 
              stack = T, flip = T, 
              group.by = grouping_var) + labs(title = "DecontX Assay")

v1 + v2 + v3 & NoLegend() & theme(axis.title.x = element_blank(),
                                  axis.text.x = element_text(size = 8))
# Save outputs grouped by cluster
ggsave(paste0("OWHA_DecontX_", grouping_var,"snRNAseqOnly_VlnPlot.png"), width = 21, height = 7)

# Repeat again for metaclusters
grouping_var <- "metaclusters"

v1 <- VlnPlot(alldata, c("Krt14","Col1a1","Col1a2","Ptprc"),
              idents = "snRNAseq",
              assay = "RNA", 
              stack = T, flip = T, 
              group.by = grouping_var,
              fill.by = "ident",
              cols = metacluster_colors) + labs(title = "RNA Assay")

v2 <- VlnPlot(alldata, c("Krt14","Col1a1","Col1a2","Ptprc"),
              idents = "snRNAseq",
              assay = "SCT", 
              stack = T, flip = T, 
              group.by = grouping_var,
              fill.by = "ident",
              cols = metacluster_colors) + labs(title = "SCT Assay")

v3 <- VlnPlot(alldata, c("Krt14","Col1a1","Col1a2","Ptprc"),
              idents = "snRNAseq",
              assay = "decontXcounts", 
              stack = T, flip = T, 
              group.by = grouping_var,
              fill.by = "ident",
              cols = metacluster_colors) + labs(title = "DecontX Assay")

v1 + v2 + v3 & NoLegend() & theme(axis.title.x = element_blank())
# Save outputs grouped by cluster
ggsave(paste0("OWHA_DecontX_", grouping_var,"snRNAseqOnly_VlnPlot.png"), width = 21, height = 7)