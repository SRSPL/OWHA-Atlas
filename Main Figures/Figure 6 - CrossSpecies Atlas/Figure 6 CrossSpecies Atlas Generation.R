
  ```{r}
set.seed(123)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 4")
```


```{r setup, include=FALSE}
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Rfast2)
library(viridis)
library(scCustomize)
library(harmony)
```

```{r}
alldata.1 = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/scRNA SEQ of wounds/Jonathan/Single-Cell RNA seq/Liu 2024/LiuHumanData_Labeled2025-01-24.rds")

alldata = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/scRNA SEQ of wounds/Skin Wound Healing Atlas/IntegratedMultimodal_FineCellTypes_061225.rds")

alldata.1$metaclusters <- Idents(alldata.1)
alldata.1$timepoint <- factor(alldata.1$timepoint,
                              levels = c("Unwounded",
                                         "D1PW",
                                         "D7PW",
                                         "D30PW"))
alldata.1 <- RenameIdents(alldata.1,
                          "Keratinocytes" = "Keratinocytes",           # or "HF Keratinocytes", adjust as needed
                          "Fibroblasts" = "Fibroblasts",
                          "Immune Cells" = "Immune cells",
                          "Mast Cells" = "Immune cells",
                          "Endothelial" = "Endothelial cells",
                          "Lymphatic Endothelial" = "Endothelial cells",
                          "Pericytes" = "Pericytes",
                          "Melanocytes" = "Melanocytes"
)

```

```{r}
# load the human and mouse converted gene symbols (made by myself)
hs_ms_genes <- read.csv("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 4/HumantoMouseOrthologs.csv")

# Day 3 after wounding (Dongqing Snhg26)
ms_seu <- readRDS("../s1-scData/01-our ms data/s2_miceD3data_onlyWT.rds")
## convert the mouse gene symbol to human gene symbol (orthology)
exp_mtx <- as.matrix(GetAssayData(alldata.1, assay = "RNA", slot = "counts"))
v2genes <- data.frame(ms_gene = rownames(exp_mtx)) %>% left_join(., hs_ms_genes)

table(is.na(v2genes$hs_gene))
##lots of NAs for which there are no human genes matching
## Remove NAs
v2genes <- v2genes[!is.na(v2genes$hs_gene),,F]
## Filter the expression matrix for genes which a mouse counterpart is available
exp_mtx <- exp_mtx[v2genes$ms_gene,]
## Now change the rownames of the matrix to the mouse gene names
rownames(exp_mtx) <- v2genes$hs_gene
#check the duplicated gene names
duplicated(rownames(exp_mtx)) %>% table()
exp_mtx <- exp_mtx[!duplicated(rownames(exp_mtx)),]
dim(exp_mtx);identical(colnames(exp_mtx), colnames(ms_seu))
## Create the seurat object with mouse genes.
ms_seu <- CreateSeuratObject(counts = exp_mtx, meta.data = ms_seu@meta.data)

saveRDS(ms_seu, file = "../s1-scData/01-our ms data/s2_miceD3data_onlyWT_humanGenes.rds")
```


```{r}

DefaultAssay(alldata.1) <- "RNA"
DefaultAssay(alldata) <- "RNA"
# Extract count matrix from human Seurat object
human_counts <- GetAssayData(alldata.1, assay = "RNA", slot = "counts")

# Create a named vector for mapping: human -> mouse
gene_map <- setNames(hs_ms_genes$mouse_Symbol, hs_ms_genes$human_Symbol)

# Rename rows (genes) using the mapping
rownames(human_counts) <- gene_map[rownames(human_counts)]

# Optional: remove rows where mapping was NA (no mouse ortholog)
human_counts <- human_counts[!is.na(rownames(human_counts)), ]

# Clone the original Seurat object
alldata.1.mouse <- alldata.1

# Replace the assay with the renamed count matrix
alldata.1.mouse[["RNA"]] <- CreateAssayObject(counts = human_counts)

# Optional: confirm dimensions
dim(alldata.1.mouse)
VlnPlot(alldata.1.mouse, "Krt14")

# Find intersecting genes
shared_genes <- intersect(rownames(alldata), rownames(alldata.1.mouse))

# Subset both objects to shared genes
alldata <- subset(alldata, features = shared_genes)
alldata.1.mouse <- subset(alldata.1.mouse, features = shared_genes)

saveRDS(alldata.1.mouse, "LiuAnnotatedConvertedtoMouseGenes04232025.rds")
```

```{r}
alldata.1 = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 4/LiuAnnotatedConvertedtoMouseGenes04232025.rds")


alldata$species <- "Mouse"
alldata.1$species <- "Human"

# Remove SCT assay from both objects to simplify
DefaultAssay(alldata) <- "RNA"

alldata[["SCT"]] <- NULL
alldata.1[["SCT"]] <- NULL

# Now merge
cross_seu <- merge(alldata, alldata.1)

```

```{r}
Layers(cross_seu)
cross_seu <- JoinLayers(cross_seu)
cross_seu[["RNA"]] <- split(cross_seu[["RNA"]], f = cross_seu$orig.ident)
```


```{r}
set.seed(123)
cross_seu <- cross_seu %>% 
  SCTransform(vars.to.regress = "percent_mito", conserve.memory = T) %>%
  RunPCA(verbose = F) 
pcs <- npcs(cross_seu)
cross_seu <- cross_seu %>%
  FindNeighbors(reduction = "pca", dims = 1:pcs, verbose = FALSE) %>%
  FindClusters(resolution = 1.5, cluster.name = "pca_clusters", verbose = FALSE) %>%
  RunUMAP(reduction = "pca", reduction.name = "umap.pca", dims = 1:pcs, verbose = FALSE)

```

```{r}

cross_seu <- IntegrateLayers(
  object = cross_seu, 
  method = HarmonyIntegration, 
  normalization.method = "SCT", 
  orig.reduction = "pca", 
  new.reduction = "harmony",
  theta = 3,
  verbose = F)

cross_seu <- FindNeighbors(cross_seu, reduction = "harmony", dims = 1:pcs) %>% 
  FindClusters(resolution = 0.5, cluster.name = "harmony_clusters") %>%
  RunUMAP(reduction = "harmony", dims = 1:pcs, reduction.name = "umap.harmony")
```

```{r}
cross_seu = readRDS("CrossSpeciesFineCellTypeLabels_2025-08-05.rds")



DimPlot(cross_seu, split.by = "timepoint", reduction = "umap.harmony")

Idents(cross_seu) <- cross_seu$metaclusters
DimPlot(cross_seu, reduction = "umap.harmony")

```
```{r}
library(dplyr)
library(ggplot2)

cross_seu$species[is.na(cross_seu$species)] <- "Human"

# 1. Extract metadata
meta <- cross_seu@meta.data %>%
  mutate(cell = rownames(.)) # keep cell ID for subsampling

# 2. Downsample 40,000 cells from each species
set.seed(42)
meta_sub <- meta %>%
  group_by(species) %>%
  sample_n(size = min(80000, n()), replace = FALSE) %>%
  ungroup()

# 3. Count cells per cluster and species
cluster_counts <- meta_sub %>%
  count(harmony_clusters, species)

# 4. Calculate proportions within each cluster
cluster_props <- cluster_counts %>%
  group_by(harmony_clusters) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# 5. Plot pie charts (faceted by harmony cluster)
ggplot(cluster_props, aes(x = "", y = prop, fill = species)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ harmony_clusters, ncol = 6) +
  theme_void() +
  scale_fill_manual(values = c("Mouse" = "#56B4E9", "Human" = "#E69F00")) +
  labs(title = "Species Proportions per Harmony Cluster (80k subsampled)",
       fill = "Species")




```
```{r}
library(dplyr)
library(ggplot2)

# 1. Extract metadata
meta <- cross_seu@meta.data %>%
  mutate(cell = rownames(.)) # keep cell ID for subsampling

# 2. Downsample 40,000 cells from each species
set.seed(42)
meta_sub <- meta %>%
  group_by(species) %>%
  sample_n(size = min(80000, n()), replace = FALSE) %>%
  ungroup()

# 3. Count cells per cluster and species
cluster_counts <- meta_sub %>%
  count(harmony_clusters, species)

# 4. Calculate proportions within each cluster
cluster_props <- cluster_counts %>%
  group_by(harmony_clusters) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# 5. Plot pie charts (faceted by harmony cluster)
ggplot(cluster_props, aes(x = "", y = prop, fill = species)) +
  geom_col(width = 1, color = "black") +
  coord_polar(theta = "y") +
  facet_wrap(~ harmony_clusters, ncol = 6) +
  theme_void() +
  scale_fill_manual(values = c("Mouse" = "#56B4E9", "Human" = "#E69F00")) +
  labs(title = "Species Proportions per Harmony Cluster (80k subsampled)",
       fill = "Species") + customtheme


ggsave("PieChartsOfSpeciesPropotion.svg", height = 7, width = 7)
```
```{r}
library(dplyr)
library(ggplot2)
library(ggalluvial)

# 1. Subsample 40,000 cells from each species
meta_sub <- cross_seu@meta.data %>%
  mutate(cell = rownames(.)) %>%
  group_by(species) %>%
  sample_n(size = min(80000, n()), replace = FALSE) %>%
  ungroup()

# 2. Separate and label mouse/human cells
mouse_meta <- meta_sub %>%
  filter(species == "Mouse") %>%
  dplyr::select(cell, mouse_label = metaclusters, cluster = harmony_clusters)

human_meta <- meta_sub %>%
  filter(species == "Human") %>%
  dplyr::select(cell, human_label = metaclusters, cluster = harmony_clusters)

# 3. Count flows between labels and clusters
sankey_df <- full_join(mouse_meta, human_meta, by = "cluster") %>%
  count(mouse_label, cluster, human_label) %>%
  filter(!is.na(mouse_label), !is.na(human_label))  # clean

# 4. Convert to factors to control order (optional, useful for display)

# 5. Plot with fixed stratum height using fake constant values in an extra plot layer
ggplot(sankey_df,
       aes(axis1 = mouse_label,
           axis2 = paste0("Cl_", cluster),
           axis3 = human_label,
           y     = n)) +
  geom_alluvium(aes(fill = mouse_label), width = 1/12, alpha = 0.8) +
  
  # Fixed-height stratum text and bars
  geom_stratum(width = 1/12, fill = "grey90", color = "black", na.rm = TRUE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5)+
  
  scale_x_discrete(
    limits = c("Mouse metacluster", "Integrated cluster", "Human metacluster"),
    expand = c(0.1, 0.05)
  ) +
  theme_minimal() +
  labs(
    title = "Mouse → Integrated Cluster → Human Cell Flow (Fixed Metacluster Width)",
    y = "Cell Count",
    x = ""
  ) +
  theme(
    axis.text.x = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.title = element_text(size = 14, face = "bold")
  ) & customtheme

ggsave("ThreeAxisSankeyPlots.svg", width = 10, height = 10)
```

```{r}
library(dplyr)
library(ggplot2)
library(ggalluvial)

# 1. Extract just human cells’ metadata, renaming “Mast cells” → “Immune Cells”
human_df <- cross_seu@meta.data %>%
  filter(tolower(species) == "human") %>%
  mutate(
    # recode metaclusters: rename “Mast cells” to “Immune Cells”
    human_label = recode(metaclusters,
                         "Mast Cells" = "Immune Cells",
                         .default    = metaclusters),
    cluster_lbl = paste0("Cl_", harmony_clusters)
  )

# 2. Count per human_label → cluster
sankey_human <- human_df %>%
  count(human_label, cluster_lbl, name = "n") %>%
  arrange(human_label, cluster_lbl)

# 3. Compute proportions within each human_label
sankey_human <- sankey_human %>%
  group_by(human_label) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# 4. Re-order the cluster_lbl factor numerically
sankey_human <- sankey_human %>%
  mutate(
    cluster_num = as.numeric(sub("^Cl_", "", cluster_lbl))
  ) %>%
  mutate(
    cluster_lbl = factor(
      cluster_lbl,
      levels = paste0("Cl_", sort(unique(cluster_num)))
    )
  ) %>%
  select(-cluster_num)

# 5. Plot: fixed bar height (y = 1) and ribbon width = prop
ggplot(sankey_human,
       aes(axis1  = human_label,
           axis2  = cluster_lbl,
           y      = n,         # all bars same height
           weight = prop)) +   # ribbon thickness ∝ proportion
  geom_alluvium(aes(fill = human_label), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/12, fill = "grey85", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("human_label", "cluster_lbl"),
    labels = c("Human metacluster", "Integrated cluster"),
    expand = c(0.05, 0.05)
  ) +
  theme_minimal() +
  labs(
    title = "Human metacluster → Integrated Harmony Cluster",
    x     = NULL,
    y     = NULL,
    fill  = "Human metacluster"
  ) +
  theme(
    axis.text.x = element_text(face = "bold"),
    plot.title  = element_text(face = "bold", size = 14)
  ) & customtheme

ggsave("HumanAxisSankeyPlots.png", width = 15, height = 15)

```

```{r}
library(dplyr)
library(ggplot2)
library(ggalluvial)

# 1. Extract just mouse cells’ metadata
mouse_df <- cross_seu@meta.data %>%
  filter(tolower(species) == "mouse") %>%
  mutate(
    mouse_label = metaclusters,
    cluster_lbl = paste0("Cl_", harmony_clusters)
  )

# 2. Count per mouse_label → cluster
sankey_mouse <- mouse_df %>%
  count(mouse_label, cluster_lbl, name = "n") %>%
  arrange(mouse_label, cluster_lbl)

# 3. Compute proportions within each mouse_label
sankey_mouse <- sankey_mouse %>%
  group_by(mouse_label) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


sankey_mouse <- sankey_mouse %>%
  mutate(
    # extract the number after "Cl_", convert to numeric, sort, then re-prefix
    cluster_num = as.numeric(sub("^Cl_", "", cluster_lbl))
  ) %>%
  mutate(
    cluster_lbl = factor(
      cluster_lbl,
      levels = paste0("Cl_", sort(unique(cluster_num)))
    )
  ) %>%
  select(-cluster_num)

# 4. Plot: fixed bar height (y = 1), ribbon width = prop
ggplot(sankey_mouse,
       aes(axis1  = mouse_label,
           axis2  = cluster_lbl,
           y      = n,         # all bars same height
           weight = prop)) +   # ribbon thickness ∝ proportion
  geom_alluvium(aes(fill = mouse_label), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/12, fill = "grey85", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("mouse_label", "cluster_lbl"),
    labels = c("Mouse metacluster", "Integrated cluster"),
    expand = c(0.05, 0.05)
  ) +
  theme_minimal() +
  labs(
    title = "Mouse metacluster → Integrated Harmony Cluster",
    x     = NULL,
    y     = NULL,
    fill  = "Mouse metacluster"
  ) +
  theme(
    axis.text.x = element_text(face = "bold"),
    plot.title  = element_text(face = "bold", size = 14)
  ) &customtheme
ggsave("MouseAxisSankeyPlots.png", width = 15, height = 15)
```
```{r}
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(Seurat)

Idents(cross_seu) <- cross_seu$metaclusters

# ---- Downsample 80,000 cells ----
set.seed(123)  # for reproducibility
cross_seu_down <- subset(cross_seu, cells = sample(colnames(cross_seu), 80000))

# ---- Store renamed identities in metadata ----
cross_seu_down$celltype <- Idents(cross_seu_down)

# ---- Create combined species + celltype label ----
cross_seu_down$sp_celltype <- paste(cross_seu_down$species, cross_seu_down$celltype, sep = "_")

# Use cluster identity
cross_seu_down$cluster <- cross_seu_down$seurat_clusters

# ---- Generate pie_df with proportion data ----
pie_df <- cross_seu_down@meta.data %>%
  group_by(cluster, sp_celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(prop = count / sum(count))

# ---- Define color palette ----
n_slices <- n_distinct(pie_df$sp_celltype)

# Split into human and mouse types
human_types <- sort(unique(grep("^Human_", pie_df$sp_celltype, value = TRUE)))
mouse_types <- sort(unique(grep("^Mouse_", pie_df$sp_celltype, value = TRUE)))

# Custom palettes
human_pal <- c("gold4","red3","springgreen4","royalblue4","orchid4","mediumpurple4")
mouse_pal <- c("olivedrab2","gold2","coral2","skyblue2","royalblue2","springgreen3",
               "orchid2","orange2","mediumpurple2","tomato4","darkseagreen2", "turquoise2")

if (length(human_types) > 8) {
  human_pal <- colorRampPalette(human_pal)(length(human_types))
}
if (length(mouse_types) > 8) {
  mouse_pal <- colorRampPalette(mouse_pal)(length(mouse_types))
}

pal <- c(human_pal, mouse_pal)
names(pal) <- c(human_types, mouse_types)

# ---- Plot ----
ggplot(pie_df, aes(x = "", y = prop, fill = sp_celltype)) +
  geom_col(width = 1, color = "black") +  # black outlines
  coord_polar(theta = "y") +
  facet_wrap(~ cluster, ncol = 6) +
  scale_fill_manual(values = pal) +
  theme_void() +
  labs(
    title = "Integrated Clusters Composition by Species+Metacluster (Downsampled 80k)",
    fill  = "Species_Metacluster"
  ) +
  theme(
    strip.text   = element_text(face = "bold"),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    legend.text  = element_text(size = 8),
    legend.title = element_text(face = "bold")
  ) +
  customtheme  # Apply your custom theme if defined elsewhere


ggsave("PieChartsofCellTpyeSpeciesContributionMetaclusters.svg", width = 7, height = 7 )


library(dplyr)
library(colorspace)  # for darken()

# Updated baseline palette
metacluster_colors <- c(
  "IFE Keratinocytes"     = "#364B9A",
  "HF Keratinocytes"      = "#56B4E9",
  "Sebocytes"             = "#8DD3C7",
  "Fibroblasts"           = "#F95B3C",
  "Immune Cells"          = "#198F1B",
  "Mast Cells"            = "chartreuse",
  "Endothelial cells"     = "#F6C800",
  "Endothelial"           = "#F6C800",   # keep for compatibility
  "Pericytes"             = "#D12E82",
  "Muscle cells"          = "#EC6F2D",
  "Schwann cells"         = "#7A4FCF",
  "Melanocytes"           = "#FA5B77",
  "Lymphatic Endothelial" = "#EC6F2D",
  "Red Blood Cells"       = "#A7222B",
  "Keratinocytes"         = "#364B9A",
  "Immune cells"          = "#198F1B",
  "Adipocytes"            = "deeppink"
)

# Extract unique labels from your Seurat object
celltypes <- unique(cross_seu_down$sp_celltype)

# Strip prefixes for base matching
clean_names <- gsub("^(Human_|Mouse_)", "", celltypes)

# Build palette
pal <- setNames(character(length(celltypes)), celltypes)

for (ct in celltypes) {
  ct_clean <- gsub("^(Human_|Mouse_)", "", ct)
  
  if (!ct_clean %in% names(metacluster_colors)) {
    warning(paste("No base color for:", ct_clean))
    next
  }
  
  base_col <- metacluster_colors[ct_clean]
  
  if (startsWith(ct, "Mouse_")) {
    pal[ct] <- base_col
  } else if (startsWith(ct, "Human_")) {
    pal[ct] <- darken(base_col, amount = 0.4)
  }
}

# ---- Check palette ----
print(pal)

# ---- Use in plot ----
ggplot(pie_df, aes(x = "", y = prop, fill = sp_celltype)) +
  geom_col(width = 1, color = "black") +
  coord_polar(theta = "y") +
  facet_wrap(~ cluster, ncol = 6) +
  scale_fill_manual(values = pal) +
  theme_void() +
  labs(
    title = "Integrated Clusters Composition by Species+Metacluster (Downsampled 80k)",
    fill  = "Species_Metacluster"
  ) +
  theme(
    strip.text   = element_text(face = "bold"),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    legend.text  = element_text(size = 8),
    legend.title = element_text(face = "bold")
  ) +
  customtheme

# ---- Debug missing colors ----
bad_colors <- pal[pal == "" | is.na(pal)]
bad_colors

ggsave("PieChartsofCellTpyeSpeciesContributionMetaclusters.svg", width = 7, height = 7 )
```
```{r}
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(Seurat)

# ---- Rename identities in Seurat object ----
cross_seu <- RenameIdents(cross_seu,
                          "Mast Cells" = "Immune",
                          "Immune Cells" = "Immune",
                          "Lymphatic Endothelial" = "Endothelial")

# Store renamed identities in metadata
cross_seu$celltype <- Idents(cross_seu)

# ---- Create combined species + celltype label ----
cross_seu$sp_celltype <- paste(cross_seu$species, cross_seu$celltype, sep = "_")

# ---- Generate pie_df grouped by species and timepoint ----
pie_df <- cross_seu@meta.data %>%
  group_by(species, timepoint, sp_celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(species, timepoint) %>%
  mutate(prop = count / sum(count),
         group = paste(species, timepoint, sep = "_"))

# ---- Define color palette ----
n_slices <- n_distinct(pie_df$sp_celltype)

# Split into human and mouse types
human_types <- sort(unique(grep("^Human_", pie_df$sp_celltype, value = TRUE)))
mouse_types <- sort(unique(grep("^Mouse_", pie_df$sp_celltype, value = TRUE)))

# Custom palettes
human_pal <- c("gold4", "red3", "springgreen4", "royalblue4", "orchid4", "mediumpurple4")
mouse_pal <- c("olivedrab2", "gold2", "coral2", "skyblue2", "royalblue2", "springgreen3",
               "orchid2", "orange2", "mediumpurple2", "tomato4", "darkseagreen2", "turquoise2")

if (length(human_types) > length(human_pal)) {
  human_pal <- colorRampPalette(human_pal)(length(human_types))
}
if (length(mouse_types) > length(mouse_pal)) {
  mouse_pal <- colorRampPalette(mouse_pal)(length(mouse_types))
}

pal <- c(human_pal, mouse_pal)
names(pal) <- c(human_types, mouse_types)

# ---- Plot ----
ggplot(pie_df, aes(x = "", y = prop, fill = sp_celltype)) +
  geom_col(width = 1, color = "black") +  # black outlines on wedges
  coord_polar(theta = "y") +
  facet_wrap(~ group, ncol = 4) +  # facet by species and timepoint combo
  scale_fill_manual(values = pal) +
  theme_void() +
  labs(
    title = "Cell Type Composition by Species and Timepoint",
    fill  = "Species_Celltype"
  ) +
  theme(
    strip.text   = element_text(face = "bold"),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    legend.text  = element_text(size = 8),
    legend.title = element_text(face = "bold")
  ) +
  customtheme  # apply your theme if defined
ggsave("PieChartsCellCompostionofTimepoints.png", width = 7, height = 7)
```
```{r}
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(Seurat)


# ---- Create combined species + celltype label ----
cross_seu$sp_celltype <- paste(cross_seu$species, cross_seu$timepoint, sep = "_")

# Use cluster identity (you can change this if you're using a custom cluster column)
cross_seu$cluster <- cross_seu$seurat_clusters

# ---- Generate pie_df with proportion data ----
pie_df <- cross_seu@meta.data %>%
  group_by(cluster, sp_celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(prop = count / sum(count))

# ---- Define color palette ----
n_slices <- n_distinct(pie_df$sp_celltype)

# Split into human and mouse types
human_types <- sort(unique(grep("^Human_", pie_df$sp_celltype, value = TRUE)))
mouse_types <- sort(unique(grep("^Mouse_", pie_df$sp_celltype, value = TRUE)))

human_pal <- c("gold4","red3","springgreen4","royalblue4")

mouse_pal <- c("royalblue2","mediumpurple2","gold2","skyblue2","red1","orange2","springgreen3")

if (length(human_types) > 8) {
  human_pal <- colorRampPalette(human_pal)(length(human_types))
}
if (length(mouse_types) > 8) {
  mouse_pal <- colorRampPalette(mouse_pal)(length(mouse_types))
}

pal <- c(human_pal, mouse_pal)
names(pal) <- c(human_types, mouse_types)

# ---- Plot ----
# ---- Plot ----
ggplot(pie_df, aes(x = "", y = prop, fill = sp_celltype)) +
  geom_col(width = 1, color = "black") +  # black outlines on each wedge
  coord_polar(theta = "y") +
  facet_wrap(~ cluster, ncol = 6) +
  scale_fill_manual(values = pal) +
  theme_void() +
  labs(
    title = "Integrated Clusters Composition by Species+Metacluster",
    fill  = "Species_Metacluster"
  ) +
  theme(
    strip.text   = element_text(face = "bold"),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    legend.text  = element_text(size = 8),
    legend.title = element_text(face = "bold")
  ) +
  customtheme  # Apply your custom theme if defined elsewhere

ggsave()
```

```{r UMAP Plots for Datasets}
set.seed(123)  # for reproducibility

# Get metadata
meta <- cross_seu@meta.data
meta$cell <- rownames(meta)

# Make sure species column exists and has "Human" / "Mouse"
table(meta$species)

# Sample 40,000 cells from each species
human_cells <- meta$cell[meta$species == "Human"]
mouse_cells <- meta$cell[meta$species == "Mouse"]

human_sample <- sample(human_cells, size = 40000)
mouse_sample <- sample(mouse_cells, size = 40000)

# Combine sampled cells
sampled_cells <- c(human_sample, mouse_sample)

# Subset Seurat object
cross_seu_down<- subset(cross_seu, cells = sampled_cells)

levels(cross_seu_down)
Idents(cross_seu_down) <- cross_seu_down$metaclusters



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
  "Red Blood Cells"   = "#A7222B",
  "Keratinocytes" = "#364B9A"
)


levels(cross_seu)
Idents(cross_seu_down) <- cross_seu_down$pred

# Define correct level order based on color map
correct_levels <- names(metacluster_colors)

# Set factor levels for both Seurat objects
cross_seu <- SetIdent(cross_seu, value = factor(Idents(cross_seu), levels = correct_levels))
cross_seu_down <- SetIdent(cross_seu_down, value = factor(Idents(cross_seu_down), levels = correct_levels))


Idents(cross_seu) <- cross_seu$fine_names

Idents(cross)


levels(cross_seu$celltype_fine)

cross_seu$fine_names <- factor(cross_seu$fine_names, levels = names(colors_pals))


DimPlot_scCustom(cross_seu, label = F,colors_use = colors_pals, reduction = "umap.harmony", order = F, raster = F, split.by = "species", group.by = "fine_names") & NoAxes() + NoLegend() & customtheme
ggsave("UMAPCrossSpeciesIntergratedDatasetFineCellTypeSplitBySpecies.png", width = 14, height = 7)




DimPlot_scCustom(cross_seu, label = F, colors_use = colors_pals, reduction = "umap.harmony") & NoAxes() & customtheme
ggsave("UMAPCrossSpeciesIntergratedDatasetLegendOnlyFineCellType.png", width = 7, height = 7)

cols = c("springgreen4", "gold2", "royalblue2", "skyblue2", "turquoise3", "coral2", "darkseagreen2", "mediumpurple2", "orange2", "orchid2", "tomato3","olivedrab2", "royalblue2")


cols = c("royalblue2","coral2","skyblue2","springgreen3","orange2","turquoise3","mediumpurple2","gold2","darkseagreen2","orchid2",
         "tomato3","olivedrab2","royalblue4","springgreen4","red3","gold4","mediumpurple4","orchid4")

DimPlot_scCustom(cross_seu_down, label = F, metacluster_colors, reduction = "umap.harmony", split.by = "species", order = F) & NoAxes()
ggsave("UMAPCrossSpeciesIntergratedDatasetSplitbyHumanMouse.png", width = 15, height = 7)


DimPlot_scCustom(cross_seu_down, label = F, metacluster_colors, reduction = "umap.harmony", split.by = "species") & NoAxes() & customtheme
ggsave("UMAPCrossSpeciesIntergratedDatasetSplitbyHumanMouseLegendOnly.png", width = 15, height = 7)


DimPlot_scCustom(cross_seu, label = T, group.by = "harmony_clusters", reduction = "umap.harmony", order = F, raster = F, ggplot_default_colors = T,label.size = 5) & NoLegend() & NoAxes() & customtheme
ggsave("UMAPCrossSpeciesIntergratedDatasetClustersLabeledFineCellTypes.png", width = 7, height = 7)

DimPlot_scCustom(cross_seu, label = T, group.by = "harmony_clusters", reduction = "umap.harmony", order = F, raster = F, ggplot_default_colors = T,label.size = 5) & NoAxes() & customtheme
ggsave("UMAPCrossSpeciesIntergratedDatasetClustersLabeledLegendOnly.png", width = 7, height = 7)
```
```{r}
features = c("Krt14","Itga6", #Keratincoytes Labels
             "Krt15","Sox9",
             "Pparg", #Sebocytes
             "Col1a2","Pdgfra", #Fibroblast Markers
             "Ptprc", # Immune Cell Markers
             "Pecam1", #Endothelial Cell Markers
             "Rgs5", #Pericytes
             "Mylpf", "Des", # Muscle Markers
             "Mpz", #Schwann Cell Markers
             "Mlana", #Melanocyte Marker
             "Adipoq", #Adipocyte Markers
             "Hba-a1","Hbb-bt" #RBC marker
)


p <- DotPlot(cross_seu, features, dot.scale = 7, cols = "RdYlBu", split.by = "species") + 
  RotatedAxis() +
  theme(
    text = element_text(size = 15, face = "bold", family = "Arial")
  ) +
  scale_size(range = c(0, 7)) + # Ensures size consistency
  geom_point(aes(size = pct.exp), shape = 21, color = "black") # Outlines circles in black

# Display the plot
p

ggsave("DotPlotCrossSpeciesAcrossCellTypesNew.png",dpi = 300,height = 6, width =9 )

```

```{r}
library(Seurat)
library(patchwork)

# List of genes
genes <- c("KRT19","FOXA1","BEST2","AQP5","ATP1A1", "PPARG")

# Generate individual FeaturePlots
plots <- lapply(genes, function(gene) {
  FeaturePlot(alldata.1, features = gene, reduction = "umap.harmony", label = TRUE)
})

# Combine plots into a single row
combined_plot <- wrap_plots(plots, nrow = 1)

# Save the combined plot
ggsave("HumanDataSweatGlandFeaturePlots.png", combined_plot, width = 25, height = 7)

```

```{r}
saveRDS(cross_seu, file = paste0("CrossSpeciesAtlas_MouseHuman_Healing_", Sys.Date(), ".rds"))


Idents(cross_seu) <- cross_seu$species
# Remove "Human_" and "Mouse_" from timepoint values
cross_seu$timepoint <- gsub("^Human_|^Mouse_", "", cross_seu$timepoint)
cross_seu$timepoint <- gsub("^Wounded_", "", cross_seu$timepoint)
# View the updated table of timepoints
table(cross_seu$timepoint)

# Clean and recode timepoints
cross_seu_down$timepoint <- recode(
  cross_seu_down$timepoint,
  "Wounded_D1PW" = "D1PW",
  "Wounded_D2PW" = "D2PW",
  "Wounded_D4PW" = "D4PW",
  "Wounded_D7PW" = "D7PW",
  "Wounded_D15PW" = "D15PW",
  "Wounded_D30PW" = "D30PW"
)


# Define the desired order
desired_order <- c("Unwounded","D1PW", "D2PW", "D4PW", "D7PW", "D15PW", "D30PW")

# Convert to factor with this order
cross_seu_down$timepoint <- factor(cross_seu_down$timepoint, levels = desired_order)

# Check ordering with table
table(cross_seu$timepoint)



# Subset the Seurat object
cross_seu_down <- subset(cross_seu, cells = cross_seu_down)

DimPlot(cross_seu_down, split.by = "timepoint", group.by = "species",reduction = "umap.harmony", shuffle = T, cols = c("#E69F00","#56B4E9"),
        pt.size = 0.4) & customtheme & NoAxes()
ggsave("UMAPCrossSpeciesIntergratedDatasetSplitbyHumanMouseOverTime.png", width = 20, height = 4)

#56B4E9", "Human" = "#E69F00"

```


```{r}
alldata.1 = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 4/CrossSpeciesAtlasMouseHumanHealing05142025.rds")


p <- VlnPlot(
  alldata.1, 
  features = c("Lgals1","Col5a1","Col6a3",
               "Mmp23b","Thbs2"),
  split.by = "timepoint",
  flip = TRUE,
  stack = TRUE
) + 
  theme(plot.margin = unit(c(1, 1, 1.5, 1), "cm"))  # top, right, bottom, left

# Save the plot
ggsave("testfigure.png", plot = p, width = 15, height = 5)
```

```{r}

cross_seu$species_timepoint <- paste(cross_seu$species, cross_seu$timepoint, sep = "_")

# Desired order
desired_order <- c(
  "Mouse_UWO", "Mouse_D1PW", "Mouse_D2PW", "Mouse_D7PW", "Mouse_D15PW", "Mouse_D30PW",
  "Human_UWO", "Human_D1PW", "Human_D7PW", "Human_D30PW"
)

# Apply factor with the desired levels to the metadata column
cross_seu$species_timepoint <- factor(cross_seu$species_timepoint, levels = desired_order)


FeaturePlot_scCustom(cross_seu, "Sema3c", split.by = "species_timepoint") & NoAxes()
ggsave("FeaturePlotSema3cSplitBySpecies_TimepointsCrossSpecies.png", height =5, width = 35)



# Split the Seurat object by species
species_list <- SplitObject(alldata.1, split.by = "species")


species_list <- lapply(species_list, function(x) {
  Idents(x) <- "celltype"
  return(x)
})

species_list$Mouse$timepoint <- factor(species_list$Mouse$timepoint, levels = c("UWO", "D1PW", "D2PW", "D7PW", "D15PW", "D30PW"))
species_list$Human$timepoint <- factor(species_list$Human$timepoint, levels = c("UWO", "D1PW", "D7PW", "D30PW"))


# Mouse plot
p_mouse <- VlnPlot_scCustom(species_list$Mouse, features = "Lgals1", split.by = "timepoint", pt.size = 0) + ggtitle("Mouse")

# Human plot
p_human <- VlnPlot_scCustom(species_list$Human, features = "Lgals1", split.by = "timepoint", pt.size = 0) + ggtitle("Human")

# Optionally display side by side
p_mouse + p_human

ggsave("VlnPlotLgals1SplitBySpeciestimepoint.png")

```
```{r Cross Species GSEA for Human Data Alone.}

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(forcats)


alldata.1 = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/scRNA SEQ of wounds/Jonathan/Single-Cell RNA seq/Liu 2024/LiuHumanData_Labeled2025-01-24.rds")

# Get all markers

Idents(alldata.1) <- "timepoint"
cellmarkers <- FindAllMarkers(alldata.1, min.pct = 0.2, recorrect_umi = FALSE, only.pos = TRUE)

# Define a function to prepare ranked gene list for GSEA
prep_gsea <- function(df) {
  df <- df[order(df$avg_log2FC, decreasing = TRUE), ]
  gene.df <- bitr(df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb =org.Hs.eg.db)
  df <- merge(df, gene.df, by.x = "gene", by.y = "SYMBOL")
  gene_list <- df$avg_log2FC
  names(gene_list) <- df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  return(gene_list)
}

# Prepare gene lists
glist.uwo <- prep_gsea(subset(cellmarkers, cluster == "Unwounded"))
glist.wo1 <- prep_gsea(subset(cellmarkers, cluster == "D1PW"))
glist.wo7 <- prep_gsea(subset(cellmarkers, cluster == "D7PW"))
glist.wo30 <- prep_gsea(subset(cellmarkers, cluster == "D30PW"))

# Run GSEA
gsea.uwo <- gseGO(geneList = glist.uwo, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID", pvalueCutoff = 0.2)
gsea.wo1 <- gseGO(geneList = glist.wo1, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID", pvalueCutoff = 0.2)
gsea.wo7 <- gseGO(geneList = glist.wo7, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID", pvalueCutoff = 0.2)
gsea.wo30 <- gseGO(geneList = glist.wo30, OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENTREZID", pvalueCutoff = 0.2)

# Optional: Convert to readable for display
gsea.uwo <- setReadable(gsea.uwo, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
gsea.wo1 <- setReadable(gsea.wo1, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
gsea.wo7 <- setReadable(gsea.wo7, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
gsea.wo30 <- setReadable(gsea.wo30, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Plotting: Example dotplots or ridgeplots for each
p1 <- dotplot(gsea.uwo, showCategory = 10, font.size = 14, title = "Unwounded")
p2 <- dotplot(gsea.wo1, showCategory = 10, font.size = 14, title = "D1PW")
p3 <- dotplot(gsea.wo7, showCategory = 10, font.size = 14, title = "D7PW")
p5 <- dotplot(gsea.wo30, showCategory = 10, font.size = 14, title = "D30PW")

ggsave("GSEA_Unwounded.png", p1, width = 10, height = 6)
ggsave("GSEA_D1PW.png", p2, width = 10, height = 6)
ggsave("GSEA_D7PW.png", p3, width = 10, height = 6)
ggsave("GSEA_D30PW.png", p5, width = 10, height = 6)

# Save results to CSV
write.csv(as.data.frame(gsea.uwo), "GSEA_Unwounded.csv", row.names = FALSE)
write.csv(as.data.frame(gsea.wo1), "GSEA_D1PW.csv", row.names = FALSE)
write.csv(as.data.frame(gsea.wo7), "GSEA_D7PW.csv", row.names = FALSE)
write.csv(as.data.frame(gsea.wo30), "GSEA_D30PW.csv", row.names = FALSE)


library(dplyr)
library(ggplot2)
library(forcats)


# Convert results to dataframes
df.uwo <- as.data.frame(gsea.uwo) %>% slice_max(order_by = abs(NES), n = 5) %>% mutate(Timepoint = "UWO")
df.wo1 <- as.data.frame(gsea.wo1) %>% slice_max(order_by = abs(NES), n = 5) %>% mutate(Timepoint = "D1PW")
df.wo7 <- as.data.frame(gsea.wo7) %>% slice_max(order_by = abs(NES), n = 5) %>% mutate(Timepoint = "D7PW")
df.wo30 <- as.data.frame(gsea.wo30) %>% slice_max(order_by = abs(NES), n = 5) %>% mutate(Timepoint = "D30PW")

# Combine into one dataframe
gsea_combined <- bind_rows(df.uwo, df.wo1, df.wo7, df.wo30)

# Set timepoint order
gsea_combined$Timepoint <- factor(gsea_combined$Timepoint, levels = c("UWO", "D1PW", "D7PW", "D30PW"))

# Make unique pathway labels to avoid overlap
gsea_combined$Description <- make.unique(gsea_combined$Description)

top_terms <- gsea_combined %>%
  count(Description, sort = TRUE) %>%
  pull(Description)

gsea_combined$Description <- factor(gsea_combined$Description, levels = rev(top_terms))

# Plot with facetting by Timepoint
p_facet <- ggplot(gsea_combined, aes(x = NES, y = Description, size = setSize, color = p.adjust)) +
  geom_point() +
  facet_wrap(~ Timepoint, nrow = 1) +
  scale_color_gradient(low = "red", high = "blue", name = "adj. p-value") +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  labs(x = "Normalized Enrichment Score (NES)", y = "GO Term", size = "Gene Set Size",
       title = "Top 5 GSEA Terms per Timepoint")

# Save the faceted plot
ggsave("Top5_GSEA_Faceted.png", p_facet, width = 16, height = 8)

library(dplyr)
library(ggplot2)
library(forcats)

# Make sure Timepoint factor has the right order
gsea_combined$Timepoint <- factor(gsea_combined$timepoint, 
                                  levels = c("UWO", "D1PW", "D7PW", "D30PW"))

# Find the timepoint of max NES for each Description
term_max_timepoint <- gsea_combined %>%
  group_by(Description) %>%
  slice_max(order_by = NES, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Timepoint, desc(NES))

# Create a factor for Description ordered first by Timepoint, then NES
gsea_combined$Description <- factor(gsea_combined$Description,
                                    levels = term_max_timepoint$Description)

# Now plot heatmap with reordered Description factor
p_heatmap_ordered <- ggplot(gsea_combined, aes(x = Timepoint, y = Description, fill = NES)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(x = "Timepoint", y = "GO Term",
       title = "GSEA NES Heatmap Ordered by Max NES Timepoint")

# Save the plot
ggsave("GSEA_NES_Heatmap_Ordered.png", p_heatmap_ordered, width = 10, height = 8)

# Dot plot: x = Timepoint, y = Description, size = setSize, color = NES
p_dotplot <- ggplot(gsea_combined, aes(x = Timepoint, y = Description, color = NES, size = setSize)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
  scale_size(range = c(2, 8), name = "Gene Set Size") +
  theme_minimal(base_size = 16, base_family = "Arial") +   # increased base_size
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 16),  # bigger x axis text
    axis.text.y = element_text(face = "bold", size = 14),                        # bigger y axis text
    axis.title = element_text(face = "bold", size = 18),                         # bigger axis titles
    legend.title = element_text(face = "bold", size = 16),                       # bigger legend title
    legend.text = element_text(size = 14)                                        # bigger legend text
  ) +
  labs(x = "Timepoint", y = "GO Term",
       title = "GSEA Dotplot Ordered by Max NES Timepoint")
# Save plot
ggsave("GSEA_Dotplot_Ordered.png", p_dotplot, width = 15, height = 8)

```
```{r}
cross_seu = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 4/CrossSpeciesFineCellTypeLabels_2025-08-05.rds")


species_fine_table <- table(cross_seu$species, cross_seu$fine_names)

write.csv(
  species_fine_table,
  file = "Species_by_FineCellTypeRawCounts.csv",
  row.names = TRUE
)





Idents(cross_seu) <- cross_seu$fine_names
VlnPlot(cross_seu,"Sox6",split.by = "species")
ggsave("Sox6VlnPlotSpecies.png", width = 20)

cross_seu <- SplitObject(cross_seu, split.by = "species")


DotPlot_scCustom(cross_seu$Mouse, features = c("Sox6","Cdh13","Sema3c","Krt16")) +
  theme_ggprism_mod() +
  customtheme +
  scale_y_discrete(limits = rev) +  # Reverse y-axis
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
ggsave("BasalIVMarkersDotPlotInMouse.png", width = 5, height = 20)


FeaturePlot_scCustom(cross_seu$Human, "Sema3c", reduction = "umap.harmony") +
  theme_ggprism_mod() +
  customtheme + NoAxes()
ggsave("FeaturePlotSox6InHuman.png", width = 4, height = 4)


FeaturePlot_scCustom(cross_seu, "Sema3c", reduction = "umap.harmony", split.by = "species", na_cutoff = 2) +
  theme_ggprism_mod() &
  customtheme + NoAxes()
ggsave("FeaturePlotSema3cInHuman.png", width = 6, height = 3)


# Subset object for Basal I–IV
basal_clusters <- c("Basal I", "Basal II", "Basal III", "Basal IV")
alldata_basal <- subset(alldata.1$Human, idents = basal_clusters)

# Plot
p <- DotPlot_scCustom(
  alldata_basal,
  features = c("Sox6", "Cdh13", "Sema3c", "Krt16")
) +
  theme_ggprism_mod() +
  customtheme +
  scale_y_discrete(limits = rev) +  # Reverse y-axis order
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) & NoLegend()  # Rotate x-axis labels

# Save
ggsave("BasalIVMarkersDotPlotInHuman.png", plot = p, width = 2, height = 2)
ggsave("BasalIVMarkersDotPlotInHuman.svg", plot = p, width = 2, height = 2)



basal_clusters <- c("Basal I", "Basal II", "Basal III", "Basal IV")
alldata_basal <- subset(alldata.1$Mouse, idents = basal_clusters)

# Plot
p <- DotPlot_scCustom(
  alldata_basal,
  features = c("Sox6", "Cdh13", "Sema3c", "Krt16")
) +
  theme_ggprism_mod() +
  customtheme +
  scale_y_discrete(limits = rev) +  # Reverse y-axis order
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) & NoLegend()  # Rotate x-axis labels

# Save
ggsave("BasalIVMarkersDotPlotInMouse.png", plot = p, width = 2, height = 2)
ggsave("BasalIVMarkersDotPlotInMouse.svg", plot = p, width = 2, height = 2)


p <- DotPlot_scCustom(
  alldata_basal,
  features = c("Sox6", "Cdh13", "Sema3c", "Krt16")
) +
  theme_ggprism_mod() +
  customtheme +
  scale_y_discrete(limits = rev) +  # Reverse y-axis order
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave("BasalIVMarkersDotPlotInMouseLegend.png", plot = p, width = 5, height = 5)
ggsave("BasalIVMarkersDotPlotInMouseLegend.png", plot = p, width = 5, height = 5)
```


```{r Heatmap of Mouse Marker Genes}

alldata.1 = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/scRNA SEQ of wounds/Jonathan/Single-Cell RNA seq/Liu 2024/LiuHumanData_Labeled2025-01-24.rds")


alldata.1$metaclusters <- Idents(alldata.1)
alldata.1$timepoint <- as.character(alldata.1$timepoint)  # convert to character if factor
alldata.1$timepoint[alldata.1$timepoint == "Unwounded"] <- "UW"
alldata.1$timepoint <- factor(alldata.1$timepoint, levels = c("UW", "D1PW", "D7PW", "D30PW"))


library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(grid)
library(ggprism)  # assuming you have customtheme from this or replace

# Your human gene list
human_genes <- c(
  "KRT15",       # Krt15
  "CCL27",       # Ccl27a :contentReference[oaicite:0]{index=0}
  "ANXA8",       # Anxa8 :contentReference[oaicite:1]{index=1}
  "KRT5",        # Krt5
  "FXYD3",       # Fxyd3
  "DMKN",        # Dmkn
  "LY6D",        # Ly6d
  "CILP2",       # Fam25c — may map to CILP2 (fictional; needs validation)
  "KRT77",       # Krt77
  "URAH",        # Urah
  "PERP",        # Perp
  "KRT10",       # Krt10
  "DGCR6",       # Dgcr6
  "TACSTD2",     # Tacstd2
  "COL17A1",     # Col17a1
  "SFN",         # Sfn
  "KRT14",       # Krt14
  "APOC1",       # Apoc1
  "CXCL14",      # Cxcl14
  "SBSN",        # Sbsn
  "KRT1",        # Krt1
  "LGALS7",      # Lgals7
  "DAPL1",       # Dapl1
  "S100A14",     # S100a14
  # "Gm10076",    # no human ortholog
  "KRTDAP",      # Krtdap
  "KLF5",        # Klf5
  "CALM4",       # Calm4
  "AHNAK2",      # Ahnak2
  "TRIM29",      # Trim29
  # "Chil3",      # no human ortholog
  "TREM2",       # Treml2
  "RETNLB",      # Retnlg
  # "SLFN4",      # no human ortholog
  "MMP8",        # Mmp8
  "LCN2",        # Lcn2
  "CXCR2",       # Cxcr2
  "FCGR1A",      # Fcgr1
  "MS4A4A",      # Ms4a4c
  "IL36G",       # Il1f9
  # "Wfdc21",     # no human ortholog
  "IFI47",       # Ifi47
  "OAS3",        # Oas3
  # "Lilr4b",     # no human ortholog
  "RSAD2",       # Rsad2
  "S100A8",      # S100a8
  # "Phf11b",     # ambiguous
  # "Fcgr4",      # no human ortholog
  # "Trim30b",    # no human ortholog
  # "B430306N03Rik", # no human ortholog
  # "Pira2",      # no human ortholog
  "SLC9A3R1",    # Slc9a3r1
  "MMP9",        # Mmp9
  "DDX58",       # Ddx58
  "HELZ2",       # Helz2
  "CXCL3",       # Cxcl3
  "CCL6",        # Ccl6
  "IFIT1",       # Ifit1
  "S100A9",      # S100a9
  # "Ifi209",     # no human ortholog
  # "Ankrd37",    # no human ortholog
  "ARG1",        # Arg1
  "ERO1L",       # Ero1l
  "CEACAM1",     # Ceacam1
  "SPP1",        # Spp1
  # "AA467197",   # no human ortholog
  "CLEC5A",      # Clec5a
  "BST1",        # Bst1
  "TMEM189",     # Tmem189
  "BNIP3",       # Bnip3
  # "Stfa2l1",    # no human ortholog
  "ITGB2",       # Itgb2
  "OLR1",        # Olr1
  # "Cstdc4",     # no human ortholog
  "CCR1",        # Ccr1
  # "Ndufb1-ps",  # no human ortholog
  # "Tlr13",      # no human ortholog
  "SLC16A3",     # Slc16a3
  "MPEG1",       # Mpeg1
  # "Mirt2",      # no human ortholog
  "C3AR1",       # C3ar1
  "PGK1",        # Pgk1
  "ITGAM",       # Itgam
  "PGAM1",       # Pgam1
  "BASP1",       # Basp1
  # "Cyrib",      # no human ortholog
  # "Lrmda",      # no human ortholog
  "GADD45B",     # Gadd45b
  # "ENSMUSG00000121513", # no human ortholog
  "QKI",         # Qki
  "SEPTIN9",     # Septin9
  # "Morrbid",    # no human ortholog
  "IGFBP7",      # Igfbp7
  "CCND3",       # Ccnd3
  # "Nav1",       # Nav1 partial
  "ZEB1",        # Zeb1
  "FNDC3B",      # Fndc3b
  "TNFAIP6",     # Tnfaip6
  "COL6A2",      # Col6a2
  "CCL7",        # Ccl7
  "GADD45G",     # Gadd45g
  # "Gm26917",    # no human ortholog
  "MFAP2",       # Mfap2
  "GPX3",        # Gpx3
  "COL6A1",      # Col6a1
  "MFAP4",       # Mfap4
  "LGALS1",      # Lgals1
  "BGN",         # Bgn
  "THBS2",       # Thbs2
  "POSTN",       # Postn
  "SERPINF1",    # Serpinf1
  "RCN3",        # Rcn3
  "LUM",         # Lum
  "CXCL1",       # Cxcl1
  "CCL2",        # Ccl2
  "PPIC",        # Ppic
  "CTSK",        # Ctsk
  "PCOLCE",      # Pcolce
  "MMP2",        # Mmp2
  # "Selenom",    # no human ortholog
  "TMSB10",      # Tmsb10
  "HSPA1B",      # Hspa1b
  "HSPA1A",      # Hspa1a
  "RPL35",       # Rpl35
  "ASPN",        # Aspn
  # "Gm15564",    # no human ortholog
  "TCHH",        # Tchh
  "SLC39A10",    # Slc39a10
  "FHOD3",       # Fhod3
  "NHS",         # Nhs
  # "2610307P16Rik", # no human ortholog
  "TENM2",       # Tenm2
  "MICAL2",      # Mical2
  "CCSER1",      # Ccser1
  "DACH1",       # Dach1
  "SOX5",        # Sox5
  "NCKAP5",      # Nckap5
  "HEXB",        # Hexb
  "NAV2",        # Nav2
  "BNC2",        # Bnc2
  "NEO1",        # Neo1
  "MEGF9",       # Megf9
  "CUX1",        # Cux1
  "PDZD2",       # Pdzd2
  "PLEKHG1",     # Plekhg1
  "H2AFJ",       # H2afj
  "SKP1",        # Skp1a
  "ID1",         # Id1
  "MARCH3"       # March3
)


# Intersect with genes present in the Seurat object to avoid errors
human_genes <- intersect(human_genes, rownames(alldata.1))

# Scale data for these genes
alldata.1 <- ScaleData(alldata.1, features = human_genes)

# Extract scaled data matrix
mat <- GetAssayData(alldata.1, slot = "scale.data")[human_genes, ]

# Extract metadata and order cells by timepoint and metacluster (or replace metacluster with a relevant column in alldata.1)
meta <- alldata.1@meta.data[colnames(mat), ]

# If you don't have "metaclusters", replace it with a suitable metadata column name
# e.g., meta$metaclusters <- meta$celltype or meta$cluster if needed

meta$order <- paste(meta$timepoint, meta$metaclusters, sep = "_")

ordered_cells <- rownames(meta[order(meta$timepoint, meta$metaclusters), ])

mat <- mat[, ordered_cells, drop = FALSE]

# Define colors for annotations (replace or customize as needed)
metacluster_colors <- c(
  "IFE Keratinocytes" = "#364B9A",
  "HF Keratinocytes"  = "#56B4E9",
  "Sebocytes"         = "#8DD3C7",
  "Fibroblasts"       = "#F95B3C",
  "Immune Cells"      = "#198F1B",   # check casing in metadata
  "Mast Cells"        = "chartreuse",
  "Endothelial cells" = "#F6C800",
  "Endothelial"       = "#F6C800",   # added exact name from error
  "Pericytes"         = "#D12E82",
  "Muscle cells"      = "#EC6F2D",
  "Schwann cells"     = "#7A4FCF",
  "Melanocytes"       = "#FA5B77",
  "Lymphatic Endothelial" = "#EC6F2D",  # fixed spelling
  "Red Blood Cells"   = "#A7222B",
  "Keratinocytes"     = "#364B9A"         # added missing key
)

nightfall_pal <- color("nightfall")(4)  # get 7 colors

# Assign 4 colors for your 4 timepoints (example: take 1st, 3rd, 5th, 7th colors)
timepoint_colors <- c(
  "UW"   = nightfall_pal[1],
  "D1PW" = nightfall_pal[2],
  "D7PW" = nightfall_pal[3],
  "D30PW"= nightfall_pal[4]
)

# Create HeatmapAnnotation
col_annot <- HeatmapAnnotation(
  df = meta[ordered_cells, c("metaclusters", "timepoint")],
  col = list(
    metaclusters = metacluster_colors,
    timepoint = timepoint_colors
  )
)

# Color function for heatmap values
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Draw heatmap and save as PNG
png(filename = "Human_HeatmapByTimepoint-Metacluster_TopGenes.png", width = 12, height = 14, units = "in", res = 300, bg = "transparent")

hm <- Heatmap(
  mat,
  name = "Scaled Expression (Z-score)",
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
```{r Making Human Heatmap}
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)  # ⬅️ Mouse annotation database
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(khroma)

# Set identity and assay
Idents(alldata.1) <- "timepoint"

alldata.1$metaclusters <- Idents(alldata.1)
alldata.1$timepoint <- as.character(alldata.1$timepoint)  # convert to character if factor
alldata.1$timepoint[alldata.1$timepoint == "Unwounded"] <- "UW"
alldata.1$timepoint <- factor(alldata.1$timepoint, levels = c("UW", "D1PW", "D7PW", "D30PW"))



alldata.1@meta.data$timepoint <- factor(alldata.1@meta.data$timepoint,
                                        levels = c("UW", "D1PW", "D7PW","D30PW"))
Idents(alldata.1) <- "timepoint"
# Find markers
cellmarkers_mouse <- FindAllMarkers(alldata.1, min.pct = 0.15, recorrect_umi = FALSE, only.pos = TRUE)

# Subset DE genes per timepoint
timepoints <- levels(alldata.1@meta.data$timepoint)
de_genes_list <- lapply(timepoints, function(tp) {
  subset(cellmarkers_mouse, cluster == tp & avg_log2FC > 1)
})
names(de_genes_list) <- timepoints

# Convert to ENTREZ (mouse)
gene_lists_entrez <- lapply(de_genes_list, function(df) {
  bitr(df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
})

# GO enrichment
ck <- compareCluster(geneCluster = gene_lists_entrez, fun = enrichGO, OrgDb = org.Hs.eg.db,
                     pvalueCutoff = 0.05, ont = "BP")
ck <- setReadable(ck, OrgDb =org.Hs.eg.db, keyType = "ENTREZID")
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
ggsave("Human_GeneOntology_GlobalPathways_AllTimepoints.png", width = 10, height = 7)

# Save GO results
write.csv(as.data.frame(ck), file = "Human_GlobalGOcompareCluster_results.csv", row.names = FALSE)

# Heatmap section
# Top 30 genes per timepoint
top_genes_list <- lapply(names(de_genes_list), function(tp) {
  top <- head(de_genes_list[[tp]][order(-de_genes_list[[tp]]$avg_log2FC), ], 30)
  transform(top, condition = tp)
})
top.genes <- do.call(rbind, top_genes_list)
mouse_genes <- unique(top.genes$gene)
mouse_genes <- setdiff(mouse_genes, "GM42418")

# Scale and extract matrix
alldata.1 <- ScaleData(alldata.1, features = mouse_genes)
mouse_genes <- intersect(mouse_genes, rownames(alldata.1))
matDE1 <- GetAssayData(alldata.1, slot = "scale.data")[mouse_genes, ]

# Order meta
meta <- alldata.1@meta.data[colnames(matDE1), ]
meta$order <- paste(meta$timepoint, meta$metaclusters, sep = "_")
ordered_cells <- rownames(meta[order(meta$timepoint, meta$metaclusters), ])
matDE1 <- matDE1[, ordered_cells, drop = FALSE]

# Annotations
metacluster_colors <- c(
  "IFE Keratinocytes" = "#364B9A",
  "HF Keratinocytes"  = "#56B4E9",
  "Sebocytes"         = "#8DD3C7",
  "Fibroblasts"       = "#F95B3C",
  "Immune Cells"      = "#198F1B",   # check casing in metadata
  "Mast Cells"        = "chartreuse",
  "Endothelial cells" = "#F6C800",
  "Endothelial"       = "#F6C800",   # added exact name from error
  "Pericytes"         = "#D12E82",
  "Muscle cells"      = "#EC6F2D",
  "Schwann cells"     = "#7A4FCF",
  "Melanocytes"       = "#FA5B77",
  "Lymphatic Endothelial" = "#EC6F2D",  # fixed spelling
  "Red Blood Cells"   = "#A7222B",
  "Keratinocytes"     = "#364B9A",
  "Immune cells"      = "#198F1B"# added missing key
)
timepoint_colors <- color("nightfall")(4)
names(timepoint_colors) <- c("UW", "D1PW", "D7PW","D30PW")

col_annot <- HeatmapAnnotation(
  df = meta[ordered_cells, c("metaclusters", "timepoint")],
  col = list(
    metaclusters = metacluster_colors,
    timepoint = timepoint_colors
  )
)

col_fun <- colorRamp2(c(-1, 0, 2), c("blue", "white", "red"))

# Save heatmap as PNG
png(
  file = "Human_HeatmapByTimepoint-Metacluster_TopGenes.png",
  width = 12, height = 15, units = 'in', bg = "transparent", res = 300
)
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

# Save heatmap as SVG
svg(
  filename = "Human_HeatmapByTimepoint-Metacluster_TopGenes.svg",
  width = 12,
  height = 15,
  bg = "transparent"
)
draw(hm, merge_legend = TRUE)
dev.off()

```

```{r}
alldata  = readRDS("CrossSpeciesFineCellTypeLabels_2025-08-05.rds")

color_pals = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 1/color_palette.rds")

DimPlot(alldata, cols = color_pals, label = T, reduction = "umap.harmony", repel = T) + NoLegend()

ggsave("LabeledUAMPCrossSpecies.png", height = 20, width = 20)
```

```{r}
library(dplyr)
library(ggplot2)
library(ggalluvial)

metacluster_colors <- c(
  "IFE Keratinocytes" = "#364B9A",
  "HF Keratinocytes"  = "#56B4E9",
  "Sebocytes"         = "#8DD3C7",
  "Fibroblasts"       = "#F95B3C",
  "Immune Cells"      = "#198F1B",   # check casing in metadata
  "Mast Cells"        = "chartreuse",
  "Endothelial cells" = "#F6C800",
  "Endothelial"       = "#F6C800",   # added exact name from error
  "Pericytes"         = "#D12E82",
  "Muscle cells"      = "#EC6F2D",
  "Schwann cells"     = "#7A4FCF",
  "Melanocytes"       = "#FA5B77",
  "Lymphatic Endothelial" = "#EC6F2D",  # fixed spelling
  "Red Blood Cells"   = "#A7222B",
  "Keratinocytes"     = "#364B9A",   
  "Immune cells"      = "#198F1B",
  "Adipocytes"        = "deeppink"
  # added missing key
)


# Optional: remove plyr if loaded
if("package:plyr" %in% loadedNamespaces()) detach("package:plyr", unload = TRUE)

# ----------------------------
# 1. (Optional) Subsample metadata
# ----------------------------
meta_sub <- cross_seu@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  as_tibble() %>%
  group_by(species) %>%
  sample_n(size = min(80000, n()), replace = FALSE) %>%
  ungroup()

# ----------------------------
# 2. Separate mouse / human metadata
# ----------------------------
mouse_meta <- meta_sub %>%
  filter(species == "Mouse") %>%
  mutate(
    cluster = as.character(harmony_clusters),
    mouse_label = as.character(metaclusters)
  ) %>%
  select(cluster, mouse_label)

human_meta <- meta_sub %>%
  filter(species == "Human") %>%
  mutate(
    cluster = as.character(harmony_clusters),
    human_label = as.character(metaclusters)
  ) %>%
  select(cluster, human_label)

# ----------------------------
# 3. Count cells per cluster/label
# ----------------------------
mouse_counts <- dplyr::count(mouse_meta, cluster, mouse_label, name = "mouse_n")
human_counts <- dplyr::count(human_meta, cluster, human_label, name = "human_n")

# ----------------------------
# 4. Join counts for Sankey plot
# ----------------------------
sankey_df <- full_join(mouse_counts, human_counts, by = "cluster") %>%
  filter(!is.na(mouse_label), !is.na(human_label)) %>%
  mutate(n = pmin(mouse_n, human_n))  # balanced flow width

# ----------------------------
# 5. Plot Sankey
# ----------------------------
ggplot(sankey_df,
       aes(axis1 = mouse_label,
           axis2 = paste0("Cl_", cluster),
           axis3 = human_label,
           y = n)) +
  geom_alluvium(aes(fill = mouse_label), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/12, fill = "grey90", color = "black", na.rm = TRUE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.5) +
  scale_x_discrete(limits = c("Mouse metacluster", "Integrated cluster", "Human metacluster"),
                   expand = c(0.1, 0.05)) +
  theme_minimal() +
  labs(title = "Mouse → Integrated Cluster → Human Cell Flow",
       y = "Cell Count",
       x = "") +
  theme(axis.text.x = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.title = element_text(size = 14, face = "bold"))

# ----------------------------
# 6. Save plot
# ----------------------------
ggsave("ggaluvialThreeAxisSankeyPlots.svg", width = 10, height = 10)



library(ggalluvial)
library(ggplot2)

# Make sure the data is alluvial-compatible
sankey_df <- sankey_df %>%
  mutate(
    axis1 = mouse_label,
    axis2 = paste0("Cl_", cluster),
    axis3 = human_label
  )

ggplot(sankey_df,
       aes(axis1 = axis1, axis2 = axis2, axis3 = axis3, y = n)) +
  geom_alluvium(aes(fill = mouse_label), width = 1/12, alpha = 0.8) +    # Color by mouse metacluster
  geom_stratum(width = 1/12, fill = "grey90", color = "black") +          # Blocks
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) + # Stratum labels
  scale_fill_manual(values = metacluster_colors) +                        # Apply custom colors
  scale_x_discrete(limits = c("Mouse metacluster", "Integrated cluster", "Human metacluster"),
                   expand = c(0.1, 0.05)) +
  theme_minimal() +
  labs(title = "Mouse → Integrated Cluster → Human Cell Flow",
       y = "Cell Count", x = "") +
  theme(axis.text.x = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.title = element_text(size = 14, face = "bold"))
ggsave("ggaluvialThreeAxisSankeyPlots.svg", width = 20, height = 20)

```


```{r}
# ==========================
# Libraries
# ==========================
library(dplyr)
library(ggplot2)
library(ggalluvial)

# ==========================
# Colors
# ==========================
metacluster_colors <- c(
  "IFE Keratinocytes" = "#364B9A",
  "HF Keratinocytes"  = "#56B4E9",
  "Sebocytes"         = "#8DD3C7",
  "Fibroblasts"       = "#F95B3C",
  "Immune Cells"      = "#198F1B",
  "Mast Cells"        = "chartreuse",
  "Endothelial cells" = "#F6C800",
  "Endothelial"       = "#F6C800",
  "Pericytes"         = "#D12E82",
  "Muscle cells"      = "#EC6F2D",
  "Schwann cells"     = "#7A4FCF",
  "Melanocytes"       = "#FA5B77",
  "Lymphatic Endothelial" = "#EC6F2D",
  "Red Blood Cells"   = "#A7222B",
  "Keratinocytes"     = "#364B9A",
  "Immune cells"      = "#198F1B",
  "Adipocytes"        = "deeppink"
)

# ==========================
# Optional cleanup
# ==========================
if("package:plyr" %in% loadedNamespaces()) detach("package:plyr", unload = TRUE)

# ==========================
# 1. Subsample metadata
# ==========================
meta_sub <- cross_seu@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  as_tibble() %>%
  group_by(species) %>%
  sample_n(size = min(80000, n()), replace = FALSE) %>%
  ungroup()

# ==========================
# 2. Separate mouse / human metadata
# ==========================
mouse_meta <- meta_sub %>%
  filter(species == "Mouse") %>%
  mutate(
    cluster = as.character(harmony_clusters),
    mouse_label = as.character(metaclusters)
  ) %>%
  select(cluster, mouse_label)

human_meta <- meta_sub %>%
  filter(species == "Human") %>%
  mutate(
    cluster = as.character(harmony_clusters),
    human_label = as.character(metaclusters)
  ) %>%
  select(cluster, human_label)

# ==========================
# 3. Count cells per cluster/label
# ==========================
mouse_counts <- dplyr::count(mouse_meta, cluster, mouse_label, name = "mouse_n")
human_counts <- dplyr::count(human_meta, cluster, human_label, name = "human_n")

# ==========================
# 4. Join counts for Sankey plot
# ==========================
sankey_df <- full_join(mouse_counts, human_counts, by = "cluster") %>%
  filter(!is.na(mouse_label), !is.na(human_label)) %>%
  mutate(
    n = pmin(mouse_n, human_n), # balanced flow width
    cluster = factor(cluster, levels = as.character(0:30))  # ensure numeric order
  )

# ==========================
# 5. Prepare for ggalluvial
# ==========================
sankey_df <- sankey_df %>%
  mutate(
    axis1 = mouse_label,
    axis2 = paste0("Cl_", cluster),
    axis3 = human_label
  )

# ==========================
# 6. Plot Sankey
# ==========================
p <- ggplot(sankey_df,
            aes(axis1 = axis1,
                axis2 = axis2,
                axis3 = axis3,
                y = n)) +
  geom_alluvium(aes(fill = mouse_label), width = 1/12, alpha = 0.85) +
  geom_stratum(width = 1/12, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4.5, fontface = "bold") +
  scale_fill_manual(values = metacluster_colors) +
  scale_x_discrete(
    limits = c("Mouse metacluster", "Integrated cluster", "Human metacluster"),
    expand = c(0.1, 0.05)
  ) +
  theme_minimal(base_size = 15) +
  labs(
    title = "Mouse → Integrated Cluster → Human Cell Flow",
    y = "Cell Count",
    x = ""
  ) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  )

# Preview plot (in RStudio or Jupyter)
p

```
```{r}
# ==========================
# Libraries
# ==========================
library(dplyr)
library(ggplot2)
library(ggalluvial)

# ==========================
# Colors
# ==========================
metacluster_colors <- c(
  "IFE Keratinocytes" = "#364B9A",
  "HF Keratinocytes"  = "#56B4E9",
  "Sebocytes"         = "#8DD3C7",
  "Fibroblasts"       = "#F95B3C",
  "Immune Cells"      = "#198F1B",
  "Mast Cells"        = "chartreuse",
  "Endothelial cells" = "#F6C800",
  "Endothelial"       = "#F6C800",
  "Pericytes"         = "#D12E82",
  "Muscle cells"      = "#EC6F2D",
  "Schwann cells"     = "#7A4FCF",
  "Melanocytes"       = "#FA5B77",
  "Lymphatic Endothelial" = "#EC6F2D",
  "Red Blood Cells"   = "#A7222B",
  "Keratinocytes"     = "#364B9A",
  "Immune cells"      = "#198F1B",
  "Adipocytes"        = "deeppink",
  "Small clusters"    = "grey80"  # color for grouped clusters
)

# ==========================
# Optional cleanup
# ==========================
if ("package:plyr" %in% loadedNamespaces()) detach("package:plyr", unload = TRUE)

# ==========================
# 1. Subsample metadata
# ==========================
meta_sub <- cross_seu@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  as_tibble() %>%
  group_by(species) %>%
  sample_n(size = min(80000, n()), replace = FALSE) %>%
  ungroup()

# ==========================
# 2. Separate mouse / human metadata
# ==========================
mouse_meta <- meta_sub %>%
  filter(species == "Mouse") %>%
  mutate(
    cluster = as.character(harmony_clusters),
    mouse_label = as.character(metaclusters)
  ) %>%
  select(cluster, mouse_label)

human_meta <- meta_sub %>%
  filter(species == "Human") %>%
  mutate(
    cluster = as.character(harmony_clusters),
    human_label = as.character(metaclusters)
  ) %>%
  select(cluster, human_label)

# ==========================
# 3. Count cells per cluster/label
# ==========================
mouse_counts <- count(mouse_meta, cluster, mouse_label, name = "mouse_n")
human_counts <- count(human_meta, cluster, human_label, name = "human_n")

# ==========================
# 4. Join counts for Sankey plot
# ==========================
sankey_df <- full_join(mouse_counts, human_counts, by = "cluster") %>%
  filter(!is.na(mouse_label), !is.na(human_label)) %>%
  mutate(
    n = pmin(mouse_n, human_n),
    cluster = as.character(cluster)
  )

# ==========================
# 5. Filter and group small clusters/flows
# ==========================

# Identify small clusters
cluster_totals <- sankey_df %>%
  group_by(cluster) %>%
  summarise(total_n = sum(n, na.rm = TRUE))

small_clusters <- cluster_totals %>%
  filter(total_n < 400) %>%
  pull(cluster)

# Apply grouping and filter out small flows
sankey_df <- sankey_df %>%
  mutate(
    cluster = ifelse(cluster %in% small_clusters, "Small clusters", cluster)
  ) %>%
  group_by(mouse_label, cluster, human_label) %>%
  summarise(n = sum(n, na.rm = TRUE), .groups = "drop") %>%
  filter(n >= 20) %>%  # remove small edges
  mutate(
    cluster = factor(cluster,
                     levels = c(as.character(0:30), "Small clusters"))
  )

# ==========================
# 6. Prepare for ggalluvial
# ==========================
sankey_df <- sankey_df %>%
  mutate(
    axis1 = mouse_label,
    axis2 = cluster,  # numeric clusters or "Small clusters"
    axis3 = human_label
  )

# ==========================
# 7. Plot Sankey
# ==========================
p <- ggplot(
  sankey_df,
  aes(axis1 = axis1, axis2 = axis2, axis3 = axis3, y = n)
) +
  geom_alluvium(aes(fill = mouse_label), width = 1/10, alpha = 0.9) +
  geom_stratum(width = 1/10, fill = "grey95", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 5, fontface = "bold") +
  scale_fill_manual(values = metacluster_colors, na.value = "grey80") +
  scale_x_discrete(
    limits = c("Mouse metacluster", "Integrated cluster", "Human metacluster"),
    expand = c(0.1, 0.05)
  ) +
  theme_minimal(base_size = 18) +
  labs(
    title = "Mouse → Integrated Cluster → Human Cell Flow",
    y = "Cell Count",
    x = ""
  ) +
  theme(
    axis.text.x = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 14),
    legend.title = element_blank()
  )

# ==========================
# 8. Save plot
# ==========================
ggsave("ggalluvialThreeAxisSankey_filtered.svg", plot = p, width = 20, height = 20)

# Preview
p


```

```{r}
# ========================== #
# 1. Libraries
# ========================== #
library(dplyr)
library(ggplot2)
library(ggalluvial)

# ========================== #
# 2. Define color palette
# ========================== #
metacluster_colors <- c(
  "IFE Keratinocytes" = "#364B9A",
  "HF Keratinocytes"  = "#56B4E9",
  "Sebocytes"         = "#8DD3C7",
  "Fibroblasts"       = "#F95B3C",
  "Immune Cells"      = "#198F1B",
  "Immune cells"      = "#198F1B",
  "Mast Cells"        = "chartreuse",
  "Endothelial cells" = "#F6C800",
  "Endothelial"       = "#F6C800",
  "Pericytes"         = "#D12E82",
  "Muscle cells"      = "#EC6F2D",
  "Schwann cells"     = "#7A4FCF",
  "Melanocytes"       = "#FA5B77",
  "Red Blood Cells"   = "#A7222B",
  "Adipocytes"        = "deeppink"
)

if("package:plyr" %in% loadedNamespaces()) detach("package:plyr", unload = TRUE)

# ========================== #
# 3. Prepare metadata (subsample)
# ========================== #
meta_sub <- cross_seu@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  as_tibble() %>%
  group_by(species) %>%
  sample_n(size = min(80000, n()), replace = FALSE) %>%
  ungroup()

mouse_meta <- meta_sub %>%
  filter(species == "Mouse") %>%
  mutate(cluster = as.character(harmony_clusters),
         mouse_label = as.character(metaclusters)) %>%
  select(cluster, mouse_label)

human_meta <- meta_sub %>%
  filter(species == "Human") %>%
  mutate(cluster = as.character(harmony_clusters),
         human_label = as.character(metaclusters)) %>%
  select(cluster, human_label)

# ========================== #
# 4. Count cells per cluster/label
# ========================== #
mouse_counts <- mouse_meta %>%
  count(cluster, mouse_label, name = "mouse_n")

human_counts <- human_meta %>%
  count(cluster, human_label, name = "human_n")

# ========================== #
# 5. Join counts and clean
# ========================== #
sankey_df <- full_join(mouse_counts, human_counts, by = "cluster") %>%
  filter(!is.na(mouse_label), !is.na(human_label)) %>%
  mutate(
    mouse_n = ifelse(is.na(mouse_n), 0, mouse_n),
    human_n = ifelse(is.na(human_n), 0, human_n),
    n = pmin(mouse_n, human_n, na.rm = TRUE)
  ) %>%
  filter(n >= 10) %>%  # remove small flows
  mutate(cluster = as.character(cluster))

# ========================== #
# 6. Group small clusters together
# ========================== #
sankey_df <- sankey_df %>%
  group_by(cluster) %>%
  mutate(total_cluster_n = sum(n)) %>%
  ungroup() %>%
  mutate(cluster_grouped = ifelse(total_cluster_n < 200, "Small clusters", cluster))

# ========================== #
# 7. Define consistent cell order for both sides
# ========================== #
cell_order <- c(
  "IFE Keratinocytes",
  "HF Keratinocytes",
  "Keratinocytes",
  "Fibroblasts",
  "Immune Cells",
  "Immune cells",
  "Mast Cells",
  "Lymphatic Endothelial",
  "Endothelial cells",
  "Endothelial",
  "Pericytes",
  "Sebocytes",
  "Schwann cells",
  "Muscle cells",
  "Melanocytes",
  "Adipocytes",
  "Red Blood Cells"
)

sankey_df <- sankey_df %>%
  mutate(
    mouse_label = factor(mouse_label, levels = cell_order),
    human_label = factor(human_label, levels = cell_order),
    cluster_grouped = factor(cluster_grouped,
                             levels = c(as.character(0:30), "Small clusters"))
  )

# ========================== #
# 8. Prepare for ggalluvial
# ========================== #
sankey_df <- sankey_df %>%
  mutate(
    axis1 = mouse_label,
    axis2 = cluster_grouped,
    axis3 = human_label
  )

# ========================== #
# 9. Plot Sankey
# ========================== #
p <- ggplot(
  sankey_df,
  aes(axis1 = axis1, axis2 = axis2, axis3 = axis3, y = n)
) +
  geom_alluvium(aes(fill = mouse_label), width = 1/12, alpha = 0.85) +
  geom_stratum(width = 1/12, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 4.5, fontface = "bold") +
  scale_fill_manual(values = c(metacluster_colors, "Small clusters" = "grey80")) +
  scale_x_discrete(
    limits = c("Mouse metacluster", "Integrated cluster", "Human metacluster"),
    expand = c(0.1, 0.05)
  ) +
  theme_minimal(base_size = 15) +
  labs(
    title = "Mouse → Integrated Cluster → Human Cell Flow",
    y = "Cell Count",
    x = ""
  ) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  )

# ========================== #
# 10. Preview or save
# ========================== #
print(p)
ggsave("ggalluvial_three_axis_sankey.svg", p, width = 14, height = 10)


```

