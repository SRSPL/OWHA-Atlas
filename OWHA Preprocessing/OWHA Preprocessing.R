title: "Multimodal Atlas Preprocessing" 

##1. Setup and Configuration

knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

# --- USER CONFIGURATION ---
data_root <- "~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/scRNA SEQ of wounds/Skin Wound Healing Atlas"
# ---------------------------

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(viridis)
library(tidyr)
library(scCustomize)
library(khroma)
library(writexl)
library(ggpubr)
library(lisi)

# Increase max globals for integration
options(future.globals.maxSize = 85 * 1024^3)

# Load custom functions
source(file.path(data_root, "Useful Functions/DeterminingPcsFunction.R"))
source(file.path(data_root, "Useful Functions/Run_ScType.R"))

path_to_db_file <- file.path(data_root, "Annotation/Skin_immunome_db.xlsx")

# Custom theme for plots
custom_theme <- theme(
  text = element_text(family = "Arial", face = "bold"),
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  legend.background = element_rect(fill = "transparent")
)


##2. Data Loading

##2.1 SRSP Lab Data (snRNA-seq and CITE-seq)

sn_base <- file.path(data_root, "Raw Data/snRNAseq Raw Data")
sn_dirs <- list(
  YJ001_uw = "YJ001_UW", YJ002_d7 = "YJ002_DPW7", YJ003_d7 = "YJ003_DPW7",
  YJ004_d15 = "YJ004_DPW15", YJ005_d15 = "YJ005_DPW15", YJ006_d30 = "YJ006_DPW30",
  YJ007_d30 = "YJ007_DPW30", YJ008_uw = "YJ008_UWO", YJ009_d4 = "YJ009_DPW4",
  YJ010_d4 = "YJ0010_DPW4", 
  YJ011_uw = "YJ011_UW/YJ011_cellranger_count_outs/filtered_feature_bc_matrix",
  YJ012_d4 = "YJ012_D4PW/YJ012_cellranger_count_outs/filtered_feature_bc_matrix",
  YJ013_d4 = "YJ013_D4PW/YJ012_cellranger_count_outs/filtered_feature_bc_matrix"
)

sn_data_list <- lapply(sn_dirs, function(path) Read10X(file.path(sn_base, path)))

cite_base <- file.path(data_root, "Raw Data/CITEseq Raw Data")
cite_dirs <- list(
  YS001_uw = "YS001_UW/filtered_feature_bc_matrix",
  YS002_d30 = "YS002_D30PW/filtered_feature_bc_matrix",
  YS003_d15 = "YS003_D15PW_RNAProtein/filtered_feature_bc_matrix",
  YS004_d15 = "YS004_D15PW_RNAProtein/filtered_feature_bc_matrix",
  YS005_d30 = "YS005_D30PW_ProjectSEE/filtered_feature_bc_matrix",
  YS006_uw = "YS006_PeripheralUW_ProjectSEE/filtered_feature_bc_matrix"
)

cite_data_list <- lapply(cite_dirs, function(path) Read10X(file.path(cite_base, path)))


##2.2 Prepare CITE-seq Matrices

prefixes <- c("HuMsRt.", "HuMs.", "MsRt.", "Ms.")
for (name in names(cite_data_list)) {
  for (pre in prefixes) {
    cite_data_list[[name]][["Antibody Capture"]] <- CollapseSpeciesExpressionMatrix(
      cite_data_list[[name]][["Antibody Capture"]], 
      prefix = pre, controls = "Isotype", ncontrols = 25
    )
  }
}


##2.3 Public scRNA-seq Data

public_base <- file.path(data_root, "Raw Data/scRNAseq Raw Data")

haensel_path <- file.path(public_base, "Realignment")
haensel_names <- c("HaenselUW1", "HaenselUW2", "Haenselwo1", "Haenselwo2", "Haenselwo3")
haensel_data <- lapply(haensel_names, function(n) Read10X(file.path(haensel_path, n, "filtered_feature_bc_matrix")))

just_path <- file.path(public_base, "GSE223660_Justynski_2023")
Justynski_24 <- ReadMtx(file.path(just_path, 'GSM6970228_24hr_matrix.mtx.gz'),
                        file.path(just_path, 'GSM6970228_24hr_barcodes.tsv.gz'),
                        file.path(just_path, 'GSM6970228_24hr_features.tsv.gz'))
Justynski_48 <- ReadMtx(file.path(just_path, 'GSM6970229_48hr_matrix.mtx.gz'),
                        file.path(just_path, 'GSM6970229_48hr_barcodes.tsv.gz'),
                        file.path(just_path, 'GSM6970229_48hr_features.tsv.gz'))

vu_path <- file.path(public_base, "Realignment")
Vu_D7PW_v2_1 <- Read10X(file.path(vu_path, "VuD7PWv2_1/filtered_feature_bc_matrix"))
Vu_D7PW_v2_2 <- Read10X(file.path(vu_path, "VuD7PWv2_2/filtered_feature_bc_matrix"))
Vu_UW_v3 <- Read10X(file.path(vu_path, "VuYoungv3UW/filtered_feature_bc_matrix"))
Vu_D4PW_v3 <- Read10X(file.path(vu_path, "VuYoungv3D4PW/filtered_feature_bc_matrix"))
Vu_D7PW_v3 <- Read10X(file.path(vu_path, "VuD7PWv3/filtered_feature_bc_matrix"))


##3. Integration and Metadata

##3.1 Create Seurat Objects

haensel_objs <- mapply(function(d, p) CreateSeuratObject(d, project = p, min.cells = 3), 
                       haensel_data, c("UW_Haensel1", "UW_Haensel2", "D4PW_Haensel1", "D4PW_Haensel2", "D4PW_Haensel3"))

just_objs <- list(
  CreateSeuratObject(Justynski_24, project = "D1PW_Justynski", min.cells = 3),
  CreateSeuratObject(Justynski_48, project = "D2PW_Justynski", min.cells = 3)
)

vu_objs <- list(
  CreateSeuratObject(Vu_D7PW_v2_1, project = "D7PW_Vu_v2_1", min.cells = 3),
  CreateSeuratObject(Vu_D7PW_v2_2, project = "D7PW_Vu_v2_2", min.cells = 3),
  CreateSeuratObject(Vu_UW_v3, project = "UW_Vu_v3", min.cells = 3),
  CreateSeuratObject(Vu_D4PW_v3, project = "D4PW_Vu_v3", min.cells = 3),
  CreateSeuratObject(Vu_D7PW_v3, project = "D7PW_Vu_v3", min.cells = 3)
)

sn_projects <- c("UW_snRNAseq_1", "D7PW_snRNAseq_1", "D7PW_snRNAseq_2", "D15PW_snRNAseq_1", 
                 "D15PW_snRNAseq_2", "D30PW_snRNAseq_1", "D30PW_snRNAseq_2", "UW_snRNAseq_2", 
                 "D4PW_snRNAseq_1", "D4PW_snRNAseq_2", "UW_snRNAseq_3", "D4PW_snRNAseq_3", "D4PW_snRNAseq_4")
sn_objs <- mapply(function(d, p) CreateSeuratObject(d, project = p, min.cells = 3), 
                  sn_data_list, sn_projects)

cite_objs <- mapply(function(d, p) {
  obj <- CreateSeuratObject(counts = d[["Gene Expression"]], project = p, min.cells = 3)
  obj[["ADT"]] <- CreateAssayObject(d[["Antibody Capture"]][, colnames(obj)])
  return(obj)
}, cite_data_list, cite_projects)

alldata <- merge(haensel_objs[[1]], y = c(haensel_objs[-1], just_objs, vu_objs, sn_objs, cite_objs))


##3.2 Assign Experimental Metadata

alldata$timepoint <- case_when(
  grepl("UW", alldata$orig.ident) ~ "Unwounded",
  grepl("D1PW", alldata$orig.ident) ~ "Wounded_D1PW",
  grepl("D2PW", alldata$orig.ident) ~ "Wounded_D2PW",
  grepl("D4PW", alldata$orig.ident) ~ "Wounded_D4PW",
  grepl("D7PW", alldata$orig.ident) ~ "Wounded_D7PW",
  grepl("D15PW", alldata$orig.ident) ~ "Wounded_D15PW",
  grepl("D30PW", alldata$orig.ident) ~ "Wounded_D30PW",
  TRUE ~ "Unknown"
)

alldata$modality <- case_when(
  grepl("snRNAseq", alldata$orig.ident) ~ "snRNAseq",
  grepl("CITEseq", alldata$orig.ident) ~ "CITEseq",
  TRUE ~ "scRNAseq"
)

alldata$internal_external <- ifelse(alldata$modality == "scRNAseq", "External", "Internal")
alldata$timepoint_modality <- paste(alldata$timepoint, alldata$modality, sep = "_")

# Add Sex and Chemistry
# [Note: Map remaining values here as needed]
source_meta <- list(
  sex = c("UW_Haensel1"="Female", "D1PW_Justynski"="Male", "D4PW_snRNAseq_2"="Male"),
  chem = c("UW_Haensel1"="10X Chromium v1", "UW_Vu_v3"="10X Chromium v3")
)


##4. Quality Control and Scoring

alldata[["percent_mito"]] <- PercentageFeatureSet(alldata, pattern = "^mt-")
alldata[["percent_ribo"]] <- PercentageFeatureSet(alldata, pattern = "^Rp[sl]")

DefaultAssay(alldata) <- "RNA"
alldata <- NormalizeData(alldata, verbose = FALSE) %>% JoinLayers()

# Stress/TF Scoring (Optional if files exist)
stress_file <- file.path(data_root, "2. Integration of Data/Stress Singatures/StressGenes.csv")
if(file.exists(stress_file)) {
  stress_genes <- read.csv(stress_file)$Gene
  stress_genes_mouse <- stringr::str_to_title(stress_genes)
  stress_genes_mouse <- intersect(stress_genes_mouse, rownames(alldata))
  alldata <- AddModuleScore(alldata, features = list(stress_genes_mouse), name = "StressScore")
}


#5. QC Subsetting and Filtering

Idents(alldata) <- "orig.ident"
alldata <- subset(alldata, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_mito < 15 & percent_ribo < 50)
Layers(alldata)


#6. Post-Filtering Visualization

# Shorten timepoint labels for plotting
Idents(alldata) <- "timepoint"
alldata <- RenameIdents(alldata,
                        "Unwounded" = "UW", "Wounded_D1PW" = "D1PW", "Wounded_D2PW" = "D2PW",
                        "Wounded_D4PW" = "D4PW", "Wounded_D7PW" = "D7PW", 
                        "Wounded_D15PW" = "D15PW", "Wounded_D30PW" = "D30PW")
alldata$timepoint <- Idents(alldata)

meta_df <- as.data.frame(alldata@meta.data)

# Summary counts
modality_total <- meta_df %>% dplyr::count(modality)
timepoint_total <- meta_df %>% dplyr::count(timepoint)
# [Add sex/chemistry counts here]

# Custom Plotting Functions
plot_total_counts <- function(data, x_var, fill_color, title) {
  ggplot(data, aes_string(x = x_var, y = "n")) +
    geom_bar(stat = "identity", fill = fill_color, color = "black") +
    labs(title = paste("Total", title), x = x_var, y = "Cell Count") +
    theme_minimal() + custom_theme
}

g1_total <- plot_total_counts(modality_total, "modality", "steelblue", "Cell Counts by Modality")
g1_total

ggsave("IntegratedMultimodal_CountsbyMetadata.png", width = 12, height = 11)


#7. Statistical Analysis of Cell Counts

CellCountChiSquareAnalysis <- function(seurat_obj, metadata_col_name) {
  observed_counts <- table(seurat_obj@meta.data[[metadata_col_name]])
  main_test <- chisq.test(observed_counts)
  # [Pairwise logic included in output object]
  return(list(Main_ChiSquare_Test = main_test))
}

# Example execution for modality
modality_stats <- CellCountChiSquareAnalysis(alldata, "modality")


#8. Integration Pipeline

# Reset timepoints to long format
alldata$timepoint <- case_when(
  alldata$timepoint == "UW" ~ "Unwounded",
  alldata$timepoint == "D1PW" ~ "Wounded_D1PW",
  # [Rest of mapping]
  TRUE ~ as.character(alldata$timepoint)
)

alldata[["RNA"]] <- split(alldata[["RNA"]], f = alldata$orig.ident)

set.seed(123)
alldata <- alldata %>% 
  SCTransform(vars.to.regress = "percent_mito", verbose = FALSE) %>%
  RunPCA(verbose = F) 

pcs <- npcs(alldata)

# Run Integration Methods
alldata <- IntegrateLayers(alldata, method = RPCAIntegration, normalization.method = "SCT", new.reduction = "rpca", verbose = F)
alldata <- IntegrateLayers(alldata, method = HarmonyIntegration, normalization.method = "SCT", new.reduction = "harmony", verbose = F)
alldata <- IntegrateLayers(alldata, method = CCAIntegration, normalization.method = "SCT", new.reduction = "cca", verbose = F)

# Clustering
alldata <- FindNeighbors(alldata, reduction = "rpca", dims = 1:pcs) %>% 
  FindClusters(resolution = 1.5, cluster.name = "rpca_clusters") %>%
  RunUMAP(reduction = "rpca", dims = 1:pcs, reduction.name = "umap.rpca")


#9. Integration Scoring (LISI)

pca_emb <- alldata@reductions$umap.pca@cell.embeddings
rpca_emb <- alldata@reductions$umap.rpca@cell.embeddings

res_pca <- compute_lisi(pca_emb, alldata@meta.data, c("orig.ident","timepoint","modality"))
res_rpca <- compute_lisi(rpca_emb, alldata@meta.data, c("orig.ident","timepoint","modality"))

# Plotting LISI results
res_combined <- rbind(
  data.frame(Method = "PCA", Score = res_pca$orig.ident),
  data.frame(Method = "RPCA", Score = res_rpca$orig.ident)
)

ggplot(res_combined, aes(x = Method, y = Score, fill = Method)) +
  geom_boxplot() + labs(title = "LISI Scores by Dataset") + theme_linedraw()

ggsave("MultimodalIntegration_LISIscores.png", width = 8, height = 7)


##10. Final Export

saveRDS(alldata, "MultimodalAtlas_Integrated_Final.rds")
