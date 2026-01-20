
  This script will load, process, and analyze data output from a Visium HD experiment. 
To begin, we load the necessary packages for analysis and set some parameters:
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
# Sets root directory for chunks to the directory where the file is located
knitr::opts_knit$set(root.dir = "~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Shared Data/Spatial Data")

set.seed(123)
future.seed = TRUE

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(dplyr)
library(tidyr)
library(SeuratWrappers)
library(Banksy)
library(hdf5r)
library(spacexr)
library(sfarrow)
library(Rfast2)
library(rlang)
#install.packages("ggsignif")  # only once
library(ggsignif)
library(khroma)
library(ggbreak)
library(extrafont)

#font_import() # Needed to import fonts from computer 
#Load needed functions
source("~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/scRNA SEQ of wounds/Skin Wound Healing Atlas/Useful Functions/GGsave_multidevice.R")

# needs to be set for large dataset analysis
options(future.globals.maxSize = 20*1024^3)

#Save a custom theme for ggplots, including setting Arial as default font and
#making background transparent
custom_theme <- theme(
  text = element_text(family = "Arial", face = "bold"),
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  legend.background = element_rect(fill = "transparent")
)
```
If needed, install necessary packages for analysis.  
```{r install packages, include=FALSE}
## Install packages required for Visium HD analysis

# Set timeout for package downloads
options(timeout = 600) 

# List of CRAN packages to check and install
cran_packages <- c("hdf5r", "arrow", "sfarrow","Rfast2","R.utils","khroma")

# Install missing CRAN packages
missing_cran <- cran_packages[!sapply(cran_packages, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

bioc_packages <- c("MAST", "SpatialExperiment")

for(package in bioc_packages){
  if (!requireNamespace(package, quietly = TRUE)){
    BiocManager::install(package)
  }
}

# Install GitHub packages
github_packages <- list(
  Banksy = "prabhakarlab/Banksy@devel",
  spacexr = "dmcable/spacexr",
  SeuratWrappers = "satijalab/seurat-wrappers"
)
for (pkg in names(github_packages)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    remotes::install_github(github_packages[[pkg]], build_vignettes = FALSE)
  }
}

# Clear space
rm("bioc_packages","cran_packages","github_packages","missing_cran","pkg")
```
## Load and QC data
Load Visium HD data, grouped by 8 um bin. Here we can also specify the dataset we are analyzing
```{r Load Visium HD data}
#Set which dataset we are going to analyze
dataset <- "YWV04_D4PW"
# Set directory for data - this should be the SpaceRanger output folder including
# the expression matrix and associated images
localdir <- paste0("~/SRSP Laboratory Dropbox/SRSP Lab/Team/Immunodynamics/Simon Van Deursen/DATA/Spatial Transcriptomics/", dataset, "/", paste0(dataset,"_Realignment"), "/binned_outputs/square_008um")

#Load in the spatial object
object <- Load10X_Spatial(data.dir = localdir,
                          filename = "filtered_feature_bc_matrix.h5")

# Save the dataset as orig.ident
object$orig.ident <- dataset

#ESSENTIAL:Change the scaling factors to only use hires - without this, the image will be very low resolution and blurry
#Discussed here: https://github.com/satijalab/seurat/issues/5614
object@images$slice1@scale.factors$lowres = object@images$slice1@scale.factors$hires
```
Plot some general QC plots
```{r QC plotting}
vln.plot <- VlnPlot(object, features = c("nCount_Spatial","nFeature_Spatial"), pt.size = 0)
count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial") + theme(legend.position = "right")
feature.plot <- SpatialFeaturePlot(object, features = "nFeature_Spatial") + theme(legend.position = "right")

# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions

vln.plot & custom_theme
# ggsave_multidevice(filename_base = paste0(dataset,"_QCVlnplots"), 
#                    width = 10, height = 5)

count.plot + feature.plot & custom_theme & theme(panel.grid = element_blank())
# ggsave_multidevice(filename_base = paste0(dataset,"_QCFeatplots"), 
#                    width = 10, height = 5)
```
```{r Plot QC Scatter plot}
FeatureScatter(object, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial",
               group.by = "orig.ident", cols = "grey25") & custom_theme
ggsave_multidevice(filename_base = paste0(dataset,"nCount_nFeature_Scatter"), height = 7, width = 8)
```
To view the distributions of counts, we can plot a histogram
```{r nCount Histogram}
ggplot(object@meta.data, aes(x = nCount_Spatial)) +
  geom_histogram(binwidth = 50, fill = "grey50", color = "grey40") +
  geom_vline(xintercept = c(10), color = "red") +
  theme_minimal() +            
  # scale_x_break(c(2000,4000), scales = 0.25, ticklabels = c(4500, 5500)) + # Add axis breaks if necessary
  theme(panel.border = element_rect(colour = "black", fill = NA)) +  # Box the plot
  labs(x = "nCount_Spatial", y = "Count") +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(face = "bold", family = "Arial"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        plot.margin = margin(0, 20, 0, 0)) +
  custom_theme

ggsave_multidevice(filename_base = paste0(dataset,"nCountHistogram"), height = 7, width = 7)
```
Show spot counts before subset
```{r Show spot count presubset}
ncol(object)
```

As we can see, there are many spots where there are fewer than 10 counts detected. We can remove these from the object before downstream analysis
```{r Subset spots with few counts}
object <- subset(object, subset = nCount_Spatial > 10)
```
Show spot counts after subset
```{r Show spot count post-subset}
ncol(object)
```
## Normalization and Clusterings
Normalize data using default log normalization - could also run SCTransform here
```{r Normalize Spatial Assay}
DefaultAssay(object) <- "Spatial"
object <- NormalizeData(object) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(assay = "Spatial", reduction.name = "pca")
```
Plot spatial gene expression plot for marker genes of interest
```{r Marker Featplot check}
DefaultAssay(object) <- "Spatial"
marker1 <- "Krt14"
marker2 <- "Lyz2"

SpatialFeaturePlot(object, features = marker1) + ggtitle(paste(marker1, "expression", sep = " "))

SpatialFeaturePlot(object, features = marker2) + ggtitle(paste(marker2, "expression", sep = " "))
```
### Segment Regions
We can then use the 10X Loupe browser to select regions of interest - in our case,
we will separate the two halves of the wound and can treat them as replicates
```{r Add Region Metadata}
# Save determined coordinates bounding region of interest, or if using Loupe browser save the csv of cell barcodes and their associated regions
region.coordinates <- read.csv(paste0("~/SRSP Laboratory Dropbox/SRSP Lab/Team/Immunodynamics/Simon Van Deursen/DATA/Spatial Transcriptomics/",dataset,"/",dataset,"_Realignment/Regions.csv"))

# Save as a named vector
region.assignments <- region.coordinates$Regions
names(region.assignments) <- region.coordinates$Barcode

# View the updated vector to ensure its properly formatted
print(region.assignments)

# Add the region info to Seurat object's metadata
object <- AddMetaData(object, region.assignments, col.name = "Region")

# Check that it was added correctly - all spots in one tissue should be the same
SpatialDimPlot(object, group.by = "Region")
```
We can then separate the images by Region using GetTissueCoordinates
```{r Segment Region Images}
# Extract coordinates from the Seurat object
coords <- GetTissueCoordinates(object)

# Add Region assignment to coordinates
coords$Region <- object$Region[rownames(coords)]

# Split the coordinates by region
region_coords_list <- split(coords, coords$Region)

# Loop through each region and create the segmentation
for (region_name in names(region_coords_list)) {
  # Extract the coordinates for the current region
  region_coords <- region_coords_list[[region_name]]
  
  # Compute convex hull (this gives us the boundary points)
  hull_indices <- chull(region_coords$x, region_coords$y)
  
  # Extract the hull coordinates
  hull_coords <- region_coords[hull_indices, c("x", "y")]
  
  # Add the 'cell' column for region label - essential to tell CreateSegmentation that the coordinates are all for the same object
  hull_coords$cell <- region_name
  
  # View the hull coordinates to ensure it's correct
  print(hull_coords)
  
  # Now you can pass formatted_coords to CreateSegmentation
  region_segmentation <- CreateSegmentation(hull_coords)
  
  # Optionally, overlay onto the tissue image
  object[[region_name]] <- Overlay(object[["slice1"]], region_segmentation)
}

# Plot the segmented regions
SpatialDimPlot(object, images = c("Region1","Region2"))
```
### Add tissue annotations
Add in anatomical annotations
```{r Add anatomical annotations}
# Save determined coordinates bounding region of interest, or if using Loupe browser save the csv of cell barcodes and their associated regions
anatomical.coordinates <- read.csv(paste0("~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Sequencing Repository/Data/Spatial Data/",dataset,"/",dataset,"_Realignment/AnatomicalSites.csv"))

# Save as a named vector
anatomical.assignments <- anatomical.coordinates$AnatomicalSites
names(anatomical.assignments) <- anatomical.coordinates$Barcode

# View the updated vector to ensure its properly formatted
print(anatomical.assignments)

# Add the region info to Seurat object's metadata
object_d30 <- AddMetaData(object_d30, anatomical.assignments, col.name = "AnatomicalSite")

# Check that it was added correctly - all spots in one tissue should be the same
SpatialDimPlot(object_d30, group.by = "AnatomicalSite")
```
Add woundsites where applicable
```{r Add woundsite annotation}
# Save determined coordinates bounding region of interest, or if using Loupe browser save the csv of cell barcodes and their associated regions
wound.coordinates <- read.csv(paste0("~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Sequencing Repository/Data/Spatial Data/",dataset,"/",dataset,"_Realignment/WoundSite.csv"))

# Save as a named vector
wound.assignments <- wound.coordinates$WoundSite
names(wound.assignments) <- wound.coordinates$Barcode

# View the updated vector to ensure its properly formatted
print(wound.assignments)

# Add the region info to Seurat object's metadata
object_d30 <- AddMetaData(object_d30, wound.assignments, col.name = "WoundSite")

# Check that it was added correctly - all spots in one tissue should be the same
SpatialDimPlot(object_d30, group.by = "WoundSite")
```
### Add 10X clustering results
We can try importing the clustering results from the 10X alignment to be used on our object
```{r Add 10X graph clusters}
clusters <- read.csv(paste0("~/SRSP Laboratory Dropbox/SRSP Lab/Team/Immunodynamics/Simon Van Deursen/DATA/Spatial Transcriptomics/", dataset, "/", paste0(dataset,"_Realignment"), "/binned_outputs/square_008um/analysis/clustering/gene_expression_graphclust/clusters.csv"))

cluster <- clusters$Cluster
names(cluster) <- clusters$Barcode
object <- AddMetaData(object, metadata = cluster, col.name = "alignment_clusters")
```

```{r Check alignment clusters}
table(object$alignment_clusters)
```
With the alignment clusters added, we can now save a dedicated plot with a defined coloring
```{r Plot alignment clusters}
Idents(object) <- "alignment_clusters"
levels(object) <- as.character(sort(as.numeric(levels(object))))
object$alignment_clusters <- Idents(object)

pal <- viridis(length(levels(object)), option = "H")
names(pal) <- levels(object)

SpatialDimPlot(object, label = F, repel = T, label.size = 4, images = "slice1") & 
  scale_fill_manual(values = pal) &
  custom_theme + theme(panel.grid = element_blank())
ggsave_multidevice(filename_base = paste0(dataset,"_10xClustering_Legend"), width = 12, height = 10)

SpatialDimPlot(object, pt.size.factor = 1.5, images = "Region1") & 
  scale_fill_manual(values = pal) & 
  custom_theme + theme(panel.grid = element_blank()) + NoLegend()
ggsave_multidevice(filename_base = paste0(dataset,"_10xClustering_Region1"), width = 8, height = 8)

SpatialDimPlot(object, pt.size.factor = 1.5, images = "Region2") & 
  scale_fill_manual(values = pal) &
  custom_theme + theme(panel.grid = element_blank()) + NoLegend()
ggsave_multidevice(filename_base = paste0(dataset,"_10xClustering_Region2"), width = 8, height = 8)
```
This will generate individual UMAP plots for each cluster, highlighted in yellow
```{r Alignment individual cluster plot}
cells <- CellsByIdentities(object)
p <- SpatialDimPlot(object, cells.highlight = cells[setdiff(names(cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
p
```
With clusters added, find markers for each compared to the rest
```{r Find Markers for alignment clusters}
Idents(object) <- "alignment_clusters"

marks <- FindAllMarkers(object, min.pct = 0.1, assay = "Spatial")
#Pull the top 50 positive markers for each cluster
top_marks <- marks %>%
  filter(avg_log2FC > 0) %>%
  group_by(cluster) %>%
  slice_head(n = 50)

DoHeatmap(object, unique(top_marks$gene), assay = "Spatial")
```
For D7PW, there are some interesting clustering results in the migrating tongue
```{r D7PW migrating tongue clusters}
# Save clusters associated with migrating tongue populations
migrating_tongue <- c("7","10","11","13","17","18","21")

# Determine cell ids for these groups
cells <- CellsByIdentities(object, idents = migrating_tongue)

# Plot
p <- SpatialDimPlot(object, images = "slice1", 
                    cells.highlight = cells[setdiff(names(cells), "NA")], 
                    cols.highlight = c(pal, "grey75"), 
                    facet.highlight = F, combine = T) 
p & custom_theme

migrating_tongue_marks <- top_marks %>%
  filter(cluster %in% migrating_tongue)

DoHeatmap(object, cells = WhichCells(object, idents = migrating_tongue), unique(migrating_tongue_marks$gene), assay = "Spatial")
```

```{r}
# Fast vectorized function to add annotations
add_spatial_annotations <- function(merged_object, datasets, annotation_type = "AnatomicalSites") {
  
  # Define column name based on annotation type
  col_name <- ifelse(annotation_type == "AnatomicalSites", "AnatomicalSite", "WoundSite")
  csv_file <- ifelse(annotation_type == "AnatomicalSites", "AnatomicalSites.csv", "WoundSite.csv")
  csv_column <- ifelse(annotation_type == "AnatomicalSites", "AnatomicalSites", "WoundSite")
  
  # Get all cell barcodes from merged object
  merged_barcodes <- colnames(merged_object)
  
  # Extract barcode prefix (everything before the last underscore and number)
  # For "s_008um_00738_00232-1_1", this extracts "s_008um_00738_00232-1"
  barcode_prefixes <- sub("_[0-9]+$", "", merged_barcodes)
  
  # Create a lookup data frame
  barcode_lookup <- data.frame(
    full_barcode = merged_barcodes,
    prefix = barcode_prefixes,
    orig.ident = merged_object$orig.ident,
    stringsAsFactors = FALSE
  )
  
  # Initialize annotation vector
  all_annotations <- rep(NA, length(merged_barcodes))
  names(all_annotations) <- merged_barcodes
  
  # Loop through each dataset
  for(dataset in datasets) {
    
    # Construct file path
    file_path <- paste0("~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Sequencing Repository/Data/Spatial Data/",
                        dataset, "/", dataset, "_Realignment/", csv_file)
    
    # Check if file exists
    if(!file.exists(file_path)) {
      warning(paste("File not found for dataset:", dataset, "- Skipping"))
      next
    }
    
    # Read coordinates
    coordinates <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Filter lookup table for this dataset
    dataset_lookup <- barcode_lookup[barcode_lookup$orig.ident == dataset, ]
    
    # Merge annotations with lookup table
    dataset_lookup <- merge(dataset_lookup, coordinates, 
                            by.x = "prefix", by.y = "Barcode", 
                            all.x = TRUE, sort = FALSE)
    
    # Assign annotations
    all_annotations[dataset_lookup$full_barcode] <- dataset_lookup[[csv_column]]
  }
  
  # Remove NA names (spots without annotations)
  all_annotations <- all_annotations[!is.na(all_annotations)]
  
  # Add annotations to merged object
  merged_object <- AddMetaData(merged_object, all_annotations, col.name = col_name)
  
  return(merged_object)
}

# Add anatomical annotations
merged_st <- add_spatial_annotations(
  merged_object = merged_st,
  datasets = unique(merged_st$orig.ident),
  annotation_type = "AnatomicalSites"
)

# Verify anatomical annotations
SpatialDimPlot(merged_st, group.by = "AnatomicalSite")

# Add wound site annotations
merged_st <- add_spatial_annotations(
  merged_object = merged_st,
  datasets = unique(merged_st$orig.ident),
  annotation_type = "WoundSite"
)

# Verify wound site annotations
SpatialDimPlot(merged_st, group.by = "WoundSite")
```

### Seurat Clustering
Seurat recommends sketching the data to downsize it and make analysis easier
```{r Sketch-based analysis for seurat clustering}
# note that data is already normalized
DefaultAssay(object) <- "Spatial"

# we select 50,0000 cells and create a new 'sketch' assay
object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch",
  features = VariableFeatures(object)
)

# switch analysis to sketched cells
DefaultAssay(object) <- "sketch"

# perform clustering workflow
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:30)
object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = 1)
object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:30)
```

```{r Project sketched data to all cells}
object <- ProjectData(
  object = object,
  assay = "Spatial",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:30,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)
```

```{r Plot Seurat clustering on UMAP}
DefaultAssay(object) <- "sketch"
Idents(object) <- "seurat_cluster.sketched"
p1 <- DimPlot(object, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(object) <- "Spatial"
Idents(object) <- "seurat_cluster.projected"
p2 <- DimPlot(object, reduction = "full.umap.sketch", label = F, raster = F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p1 | p2
```

```{r Save plots of seurat clustering}
Idents(object) <- "seurat_cluster.projected"
levels(object) <- as.character(sort(as.numeric(levels(object))))
object$seurat_cluster.projected <- Idents(object)

SpatialDimPlot(object, label = F, repel = T, label.size = 4, images = "slice1") & 
  scale_fill_manual(values = viridis(length(levels(object)), option = "H")) &
  custom_theme
ggsave(paste0(dataset,"_SeuratClustering_res1.png"), width = 8, height = 7)
```
### BANKSY clustering
Lastly, we will try the BANKSY implementation in Seurat: Before running BANKSY, there are two important model parameters that users should consider:
  
  k_geom : Local neighborhood size. Larger values will yield larger domains. A k_geom of 24 corresponds to the second order spots surrounding the index spot in Visium HD data, as oligos are arranged in 2 um squares and binned into 8 μm bins.

lambda : Influence of the neighborhood. Larger values yield more spatially coherent domains. Must be within [0,1], and a value of 0 will result in clusters entirely by cell type with no spatial influence. Default values implemented in the associated paper were 0.2 for cell type clustering and 0.8 for spatial clustering.
The RunBanksy function creates a new BANKSY assay, which can be used for dimensional reduction and clustering:
  ```{r Run BANKSY 0.2 Lambda}
DefaultAssay(object) <- "Spatial"
object <- RunBanksy(object,
                    lambda = 0.2, verbose = TRUE,
                    assay = "Spatial", slot = "data", features = "variable",
                    k_geom = 24
)
```
```{r 0.2 Lambda Clustering}
DefaultAssay(object) <- "BANKSY"
object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy.0.2", features = rownames(object), npcs = 30)
object <- FindNeighbors(object, reduction = "pca.banksy.0.2", dims = 1:30)
object <- FindClusters(object, cluster.name = "banksy_cluster_0.2", resolution = 0.5)
```
In order to test BANKSY clustering, we will run it twice with two different lambda values: 0.2 to generate cell type oriented clusters, and 0.8 to generate spatial clusters. 

```{r Run BANKSY with Lambda 0.8}
DefaultAssay(object) <- "Spatial"
object <- RunBanksy(object,
                    lambda = 0.8, verbose = TRUE,
                    assay = "Spatial", slot = "data", features = "variable",
                    k_geom = 24
)

DefaultAssay(object) <- "BANKSY"
object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy.0.8", features = rownames(object), npcs = 30)
object <- FindNeighbors(object, reduction = "pca.banksy.0.8", dims = 1:30)
object <- FindClusters(object, cluster.name = "banksy_cluster_0.8", resolution = 0.8)
```

```{r Save plots of BANKSY clusterings with different lambda values}
#Idents(object) <- "banksy_cluster"
levels(object$banksy_cluster_0.2) <- as.character(sort(as.numeric(levels(object$banksy_cluster_0.2))))
p2 <- SpatialDimPlot(object, group.by = "banksy_cluster_0.2", label = F,
                     repel = T, label.size = 4, images = "slice1") + 
  labs(title = "lambda = 0.2 and k_geom = 24") & 
  scale_fill_manual(values = viridis(length(unique(object$banksy_cluster_0.2)), 
                                     option = "H")) &
  custom_theme
p2
ggsave_multidevice(filename_base = paste0(dataset,"_BANKSYClustering_0.2Lambda"),
                   width = 8, height = 7)

p3 <- SpatialDimPlot(object, group.by = "banksy_cluster_0.8", label = F, 
                     repel = T, label.size = 4, images = "slice1") + 
  labs(title = "lambda = 0.8 and k_geom = 24") & 
  scale_fill_manual(values = viridis(length(unique(object$banksy_cluster_0.8)), 
                                     option = "H")) &
  custom_theme
p3

ggsave_multidevice(filename_base = paste0(dataset,"_BANKSYClustering_0.8Lambda"),
                   width = 8, height = 7)
```
We can test additional lambda values
```{r Run BANKSY with lambda 0.5}
lambda <- 0.5

DefaultAssay(object) <- "Spatial"
object <- RunBanksy(object,
                    lambda = lambda, verbose = TRUE,
                    assay = "Spatial", slot = "data", features = "variable",
                    k_geom = 24
)

DefaultAssay(object) <- "BANKSY"
object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy.0.5", features = rownames(object), npcs = 30)
object <- FindNeighbors(object, reduction = "pca.banksy.0.5", dims = 1:30)
object <- FindClusters(object, cluster.name = "banksy_cluster_0.5", resolution = 1)

p2 <- SpatialDimPlot(object, group.by = "banksy_cluster_0.5", label = F, 
                     repel = T, label.size = 4, images = "slice1") + 
  labs(title = "lambda = 0.5 and k_geom = 24") & 
  scale_fill_manual(values = viridis(length(unique(object$banksy_cluster_0.5)), 
                                     option = "H")) & 
  custom_theme
p2

ggsave(paste0(dataset,"_BANKSYClustering_0.5Lambda.png"), width = 8, height = 7)
```
We can also use the BANKSY clustering to generate a UMAP plot
```{r Use BANKSY to run UMAP}
# switch to full dataset
Idents(object) <- "banksy_cluster_0.5"

object <- RunUMAP(object, reduction = "pca.banksy.0.5", reduction.name = "umap.banksy", dims = 1:30)

DimPlot(object, reduction = "umap.banksy", label = F, raster = F) + ggtitle("BANKSY Clustering") + theme(legend.position = "bottom")
```
Lastly, with the object divided, we can generate new BANKSY clusters specifying each region - this will treat each tissue piece as its own batch
```{r Run BANKSY on each Region}
# Run BANKSY
DefaultAssay(object) <- "Spatial"
object <- RunBanksy(object,
                    lambda = 0.8, k_geom = 24, verbose = TRUE,
                    assay = "Spatial", slot = "data", features = "variable",
                    group = "Region"
)

DefaultAssay(object) <- "BANKSY"
object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(object), npcs = 30)
object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
object <- FindClusters(object, cluster.name = "banksy_cluster_split", resolution = 0.5)
```

```{r Plot region-specific BANKSY clustering}
SpatialDimPlot(object, images = c("Region.1","Region.2"), group.by = "banksy_cluster_split", label = F, repel = T) & scale_fill_manual(values = viridis(length(levels(object$banksy_cluster_split)), begin = 0.05, option = "H")) & NoLegend()
ggsave(paste0(dataset,"_SeperatedRegions_SplitBanksyClustering.png"), width = 16, height = 8)
```
```{r Plot BANKSY split clustering on slice1}
SpatialDimPlot(object, images = c("slice1"), group.by = "banksy_cluster_split", label = T, repel = T) & scale_fill_manual(values = viridis(length(levels(object$banksy_cluster_split)), begin = 0.05, option = "H")) 

ggsave(paste0(dataset,"_LabeledSplitBanksyClustering.png"), width = 12, height = 8)
```
Find markers for the split banksy clusters
```{r Find markers for banksy clustering}
Idents(object) <- "banksy_cluster_split"

#Find markers for all banksy clusters
spatial_markers <- FindAllMarkers(object, assay = "Spatial", min.pct = 0.1)

top_marks <- spatial_markers %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 1) %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  ungroup()

object_subset <- subset(object, cells = Cells(object[["Spatial"]]), downsample = 5000)
object_subset <- ScaleData(object_subset, assay = "Spatial", features = top_marks$gene)

p <- DoHeatmap(object_subset, assay = "Spatial", features = top_marks$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p
```

This is a good stopping point before proceeding to annotation, so save the current object
```{r Save RDS pre-annotation}
saveRDS(object, paste0(dataset,"_Unannotated_051225.rds"))

object <- readRDS(paste0(dataset,"_Unannotated_031725.rds"))
```
## RCTD Annotation
Instead of trying to cluster the spots on their own, we can try integrating them with the existing atlas data to annotate spots. Seurat has a pipeline for this that uses this package: https://www.nature.com/articles/s41587-021-00830-w
```{r Sketch analysis}
# sketch the cortical subset of the Visium HD dataset
DefaultAssay(object) <- "Spatial"
object <- FindVariableFeatures(object)
object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch",
  features = VariableFeatures(object)
)

DefaultAssay(object) <- "sketch"
object <- ScaleData(object)
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch", verbose = T)
object <- FindNeighbors(object, reduction = "pca.sketch", dims = 1:50)
object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50, verbose = T)
```

Now we load in the reference data to be used - since we only are analyzing unwounded spatial data, we can subset that out from the larger dataset
```{r Load RCTD reference}
# load in the reference scRNA-seq dataset
ref <- readRDS("~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/scRNA SEQ of wounds/Skin Wound Healing Atlas/IntegratedMultimodal_FineCellTypes_061225.rds")

Idents(ref) <- "timepoint"
#Subset unwounded cells
#ref <- subset(ref, idents = "Unwounded")

ref[["RNA"]] <- JoinLayers(ref[["RNA"]])
```
Run RCTD on full object, first by converting the reference dataset loaded to the format needed for RCTD, and extracting count and object data from our Visium object
```{r Create RCTD reference and query objects}
# Set identities to fine cell type annotation in the atlas reference
Idents(ref) <- "celltype_fine"
# Save RNA counts
counts <- ref[["RNA"]]$counts
# Factor cell type
cluster <- as.factor(ref$celltype_fine)
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)

counts_hd <- object[["Spatial"]]$counts
object_cells_hd <- colnames(object[["Spatial"]])
coords <- GetTissueCoordinates(object)[object_cells_hd, 1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

rm(ref)
```
With the objects formatted, we can run RCTD. This step takes a long time (~2-4 hours) for each object - when annotating with fine cell type labels, will take significantly longer - as of 3/28/25, has taken 36 hours and not completed
```{r Run RCTD}
# run RCTD - the min UMI needs to be set very low, as the Visium HD is fairly sparse in some regions
RCTD <- create.RCTD(query, reference, max_cores = 12, UMI_min = 10,CELL_MIN_INSTANCE = 20)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# add results back to Seurat object
object <- AddMetaData(object, metadata = RCTD@results$results_df)
```
Save the RCTD output
```{r Save RCTD results}
saveRDS(RCTD, paste0(dataset,"_RCTDOutput.rds"))
```
RCTD will fail to annotate certain spots depending on the number of features present, among other things. In these cases, we can overwrite the annotation so NA values are annotated as unknown
```{r Make NA first_type Unknown}
# project RCTD labels from sketched cortical cells to all cortical cells
object$first_type <- as.character(object$first_type)
object$first_type[is.na(object$first_type)] <- "Unknown"
```
Examine how RCTD called spots, with singlets, doublets, and others
```{r Examine RCTD doublet calling}
SpatialDimPlot(object, group.by = "spot_class", images = "slice1") & 
  custom_theme + theme(panel.grid = element_blank())
ggsave(paste0(dataset,"_spacexrRCTD_doublet.png"), width = 9, height = 7)
```
To check the cell type labeling, we can plot markers of different populations to see if they are correlated
```{r Check RCTD with marker gene feature plots}
DefaultAssay(object) <- "Spatial"

SpatialFeaturePlot(object, features = "Krt72", images = "slice1") + ggtitle("IRS Keratinocytes")

SpatialFeaturePlot(object, features = "Lgr5", images = "slice1") + ggtitle("ORS Keratinocytes")

SpatialFeaturePlot(object, features = "Krt79", images = "slice1") + ggtitle("uHF Keratinocytes")

SpatialFeaturePlot(object, features = "Adipoq", images = "slice1") + ggtitle("Adipocytes")

SpatialFeaturePlot(object, features = "Krt14", images = "slice1") + ggtitle("Keratinocytes")

SpatialFeaturePlot(object, features = "Mylpf", images = "slice1") + ggtitle("Muscle")

SpatialFeaturePlot(object, features = "Mmp9", images = "slice1") + ggtitle("OC-like")
```

```{r Set levels of RCTD annotation}
# We can use the levels from the RCTD reference object to create the cell type annotation
object$first_type <- factor(object$first_type, levels = c(levels(reference@cell_types),"Unknown"))
Idents(object) <- "first_type"
```
Now that levels are set appropriately, we can generate a color palette based on metaclusters
```{r Save metacluster and fine cell type color palette}
Idents(object) <- "first_type"

# Set color palette for main metacluster colors
metacluster_colors <- c(
  "IFE Keratinocytes" = "#364B9A",   # deep blue
  "HF Keratinocytes"  = "#56B4E9",   # vivid turquoise
  "Sebocytes"         = "#8DD3C7",   # strong teal
  "Fibroblasts"       = "#F95B3C",   # bright red-orange
  "Immune cells"      = "#198F1B",   # saturated green
  "Endothelial cells" = "#F6C800",   # bright yellow
  "Pericytes"         = "#D12E82",   # strong magenta
  "Muscle cells"      = "#EC6F2D",   # vibrant orange
  "Schwann cells"     = "#7A4FCF",   # bold violet
  "Melanocytes"       = "#FA5B77",   # punchy coral-pink
  "Adipocytes"        = "#19A1D1",   # bright cyan-blue
  "Red Blood Cells"   = "#A7222B",   # vivid blood red
  "Unknown" = "grey70"
)
# Load fine cell type color palette
cell_colors <- readRDS("~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 1/color_palette.rds")
# Add Unknown to cell palettte
cell_colors <- c(cell_colors, "Unknown" = "grey70")
```

```{r Plot RCTD annotation with metacluster coloring}
# Generate SpatialDimPlot with the custom color palette
p <- SpatialDimPlot(object, images = "slice1") +
  scale_fill_manual(values = cell_colors)

# Display the plot
print(p)
```
Save plots of the annotations overlayed on the spatial dimplot with new color palette
```{r Plot RCTD annotations on spatial DimPlot}
Idents(object) <- "first_type"
p <- SpatialDimPlot(object, images = "slice1", pt.size.factor = 1) +
  scale_fill_manual(values = cell_colors)

p & custom_theme + theme(panel.grid = element_blank())
ggsave_multidevice(filename_base = paste0(dataset,"_FineCellType_Legend"), width = 18, height = 10)

p + NoLegend() & custom_theme + theme(panel.grid = element_blank())
ggsave_multidevice(filename_base = paste0(dataset,"_FineCellType"), width = 12, height = 10)
```
```{r Test point size factors}
Idents(object) <- "first_type"

sizes <- c(0.5,0.75,1,1.5,2)
for(size in sizes){
  SpatialDimPlot(object, images = "Region.1", pt.size.factor = size) +
    scale_fill_manual(values = cell_colors) + NoLegend() & custom_theme + 
    theme(panel.grid = element_blank())
  
  ggsave(paste0(dataset,"_FineCellType_Section1_ptSize",size,".png"),
         width = 10, height = 10)
}
```

Save plots for each region - for YWV01 and 2, image is listed as "Region.X",
whereas for YWV03 and 4 it is "RegionX"
```{r Plot cell types by region}
Idents(object) <- "first_type"

SpatialDimPlot(object, images = "Region1", pt.size.factor = 1.25) + # 1.5 for YWV01, 1 for YWV02
  scale_fill_manual(values = cell_colors) + NoLegend() & custom_theme + theme(panel.grid = element_blank())

ggsave_multidevice(filename_base = paste0(dataset,"_FineCellType_Section1"),
                   width = 10, height = 10)

SpatialDimPlot(object, images = "Region2", pt.size.factor = 1.25) +
  scale_fill_manual(values = cell_colors) + NoLegend() & custom_theme + theme(panel.grid = element_blank())

ggsave_multidevice(filename_base = paste0(dataset,"_FineCellType_Section2"),
                   width = 10, height = 10)
```
Show cell types on UMAP
```{r Plot cell types on UMAP}
Idents(object) <- "first_type"

p1 <- DimPlot(object, reduction = "full.umap.sketch", label = F, 
              repel = T, raster = F, 
              cols = cell_colors) 

p1 + NoLegend() + NoAxes() & custom_theme
ggsave(paste0(dataset,"_CellTypeUMAP_NoLegend.png"), width = 8, height = 8)

p1 & custom_theme
ggsave(paste0(dataset,"_CellTypeUMAP_Legend.png"), width = 16, height = 8)
```
We can also group the fine cell types by metacluster
```{r Add metacluster annotations}
# Create a named vector for renaming cell types to metaclusters
metacluster_levels <- c(
  # IFE Keratinocytes
  "Basal I" = "IFE Keratinocytes", "Basal II" = "IFE Keratinocytes", 
  "Basal III" = "IFE Keratinocytes", "Basal IV" = "IFE Keratinocytes",
  "uHF I" = "IFE Keratinocytes","uHF II" = "IFE Keratinocytes",
  "Spinous I" = "IFE Keratinocytes", "Spinous II" = "IFE Keratinocytes", 
  "Spinous III" = "IFE Keratinocytes", "Spinous IV" = "IFE Keratinocytes", 
  "Cornified" = "IFE Keratinocytes",
  # HF Keratinocytes
  "Ker.Cycling I" = "HF Keratinocytes", "Ker.Cycling II" = "HF Keratinocytes",
  "Ker.Cycling III" = "HF Keratinocytes","IRS I" = "HF Keratinocytes",
  "IRS II" = "HF Keratinocytes","Medulla" = "HF Keratinocytes", 
  "Cortex" = "HF Keratinocytes", "Outer Bulge I" = "HF Keratinocytes", 
  "Outer Bulge II" = "HF Keratinocytes", "Outer Bulge III" = "HF Keratinocytes", 
  "HF I" = "HF Keratinocytes", "HF II" = "HF Keratinocytes", 
  # Sebocytes
  "Transitional Basal Cells" = "Sebocytes", "Sebocyte I" = "Sebocytes", 
  "Sebocyte II" = "Sebocytes", "Sebocyte III" = "Sebocytes",  
  # Fibroblasts
  "Papillary I" = "Fibroblasts", "Papillary II" = "Fibroblasts", 
  "Papillary III" = "Fibroblasts", "Papillary IV" = "Fibroblasts", 
  "Reticular I" = "Fibroblasts", "Reticular II" = "Fibroblasts", 
  "Reticular III" = "Fibroblasts", "Reticular IV" = "Fibroblasts", "Reticular V" = "Fibroblasts", 
  "Fascia" = "Fibroblasts", "Myofibroblasts I" = "Fibroblasts", 
  "Myofibroblasts II" = "Fibroblasts", "Myofibroblasts III" = "Fibroblasts", 
  "Lef1+ Fibroblasts" = "Fibroblasts", "Tnc+ Fibroblasts" = "Fibroblasts", 
  "Myoc+ Fibroblasts" = "Fibroblasts", "Lox+ Fibroblasts" = "Fibroblasts", 
  "Nfkb+ Fibroblasts" = "Fibroblasts",
  # Immune cells
  "Monocytes I" = "Immune cells", "Monocytes II" = "Immune cells","Monocytes III" = "Immune cells",
  "Macrophages I" = "Immune cells", "Macrophages II" = "Immune cells", "Macrophages III" = "Immune cells",
  "Macrophages IV" = "Immune cells", "Cycling Macrophages" = "Immune cells",
  "OC-like cells" = "Immune cells", "moDCs" = "Immune cells", "cDC1" = "Immune cells", 
  "cDC2" = "Immune cells", "Migratory DCs" = "Immune cells", "Langerhans cells" = "Immune cells",
  "pDCs" = "Immune cells", "Neutrophils I" = "Immune cells", "Neutrophils II" = "Immune cells", 
  "Neutrophils III" = "Immune cells","Basophils" = "Immune cells", "Mast cells" = "Immune cells",
  "NK cells" = "Immune cells", "γδT cells" = "Immune cells", "T cells" = "Immune cells", 
  "Tregs" = "Immune cells", "ILCs" = "Immune cells", "B cells" = "Immune cells",  
  # Endothelial cells
  "Tacr1+ ECs" = "Endothelial cells","Capillary ECs" = "Endothelial cells", 
  "Mature ECs" = "Endothelial cells","LEC I" = "Endothelial cells", 
  "LEC II" = "Endothelial cells", "Artery ECs" = "Endothelial cells", 
  "Proliferating ECs" = "Endothelial cells", "Progenitor ECs" = "Endothelial cells", 
  "EC I" = "Endothelial cells", "EC II" = "Endothelial cells", 
  "Vein ECs" = "Endothelial cells","PVMs" = "Endothelial cells",
  "Proliferating LECs" = "Endothelial cells","LEC III" = "Endothelial cells", 
  "EC III" = "Endothelial cells", "VWF-hi ECs" = "Endothelial cells",
  "Pericytes" = "Endothelial cells","VSMCs" = "Endothelial cells",
  # Pericytes
  "Pericyte I" = "Pericytes", "Pericyte II" = "Pericytes", "Pericyte III" = "Pericytes", 
  "Pericyte IV" = "Pericytes", 
  # Muscle
  "Muscle Progenitors" = "Muscle cells", "Skeletal Muscle" = "Muscle cells", 
  "Myonuclei" = "Muscle cells", 
  # Schwann cells
  "Myelinating SC I" = "Schwann cells", "Myelinating SC II" = "Schwann cells", 
  "Repair SC I" = "Schwann cells", "Repair SC II" = "Schwann cells",  
  "Repair SC III" = "Schwann cells",
  # Melanocytes
  "Oca2- Melanocytes" = "Melanocytes", "Oca2-lo Melanocytes" = "Melanocytes",
  "Oca2+ Melanocytes" = "Melanocytes",
  # Adipocytes
  "Adipocyte I" = "Adipocytes", "Adipocyte II" = "Adipocytes",  
  # RBCs
  "Red Blood Cells" = "Red Blood Cells",
  # Unknown
  "Unknown" = "Unknown"  # Unknown
)

# Assign the new identifiers
Idents(object) <- "first_type"

# Rename the identifiers based on the new levels
object <- RenameIdents(object, metacluster_levels)

# Optionally, add the new metacluster information to the object
object$metaclusters <- Idents(object)

SpatialDimPlot(object, images = "slice1", pt.size.factor = 1) +
  scale_fill_manual(values = metacluster_colors) & custom_theme + theme(panel.grid = element_blank())

ggsave_multidevice(filename_base = paste0(dataset,"_Metaclusters"), width = 12, height = 10)
```
Save metacluster plots for each section
```{r Plot metaclusters for each tissue section}
SpatialDimPlot(object, images = "Region1", pt.size.factor = 1.5) + # 1.5 for YWV01, 1 for YWV02
  scale_fill_manual(values = metacluster_colors) + NoLegend() & custom_theme + theme(panel.grid = element_blank())

ggsave_multidevice(filename_base = paste0(dataset,"_Metaclusters_Section1"), width = 10, height = 10)

SpatialDimPlot(object, images = "Region2", pt.size.factor = 1.5) +
  scale_fill_manual(values = metacluster_colors) + NoLegend() & custom_theme + theme(panel.grid = element_blank())

ggsave_multidevice(filename_base = paste0(dataset,"_Metaclusters_Section2"), width = 10, height = 10)
```
Plot metacluster annotations on the UMAP
```{r Plot metaclusters on UMAP}
Idents(object) <- "metaclusters"

p <- DimPlot(object, reduction = "full.umap.sketch", label = F, 
             repel = T, raster = F, 
             cols = metacluster_colors) & custom_theme

p + NoLegend() + NoAxes() 
ggsave(paste0(dataset,"_MetaclusterUMAP_NoLegend.png"), width = 8, height = 8)

p
ggsave(paste0(dataset,"_MetaclusterUMAP_Legend.png"), width = 10, height = 8)
```
Now make a plot in which each metacluster is highlighted separately, with metacluster colors
```{r Plot each metacluster highlighted in one facet plot}
Idents(object) <- "metaclusters"
# Get cells grouped by metacluster
metacluster_cells <- CellsByIdentities(object)

# Build a full palette including “Unselected” - necessary to change background color
full_pal <- c(
  metacluster_colors,
  Unselected = "grey65"
)

# Make SpatialDimPlot
p <- SpatialDimPlot(
  object,
  images          = "slice1",
  cells.highlight = metacluster_cells,
  cols.highlight  = c(metacluster_colors, "grey65"),  # placeholder
  facet.highlight = TRUE,
  combine         = TRUE
) + NoLegend() & custom_theme

# Override the fill scale
p <- p & scale_fill_manual(
  values = full_pal,
  breaks = names(full_pal)  # ensures both the metacluster and “Unselected” levels are known
)

p & theme(title = element_text(face = "bold"),
          panel.grid = element_blank())

ggsave_multidevice(filename_base = paste0(dataset,"_MetaclusterSplitDimPlot"), width = 10, height = 15)
```
To highlight specific groups on the DimPlot, starting with metaclusters
```{r Highlight immune metacluster}
cells <- CellsByIdentities(object, idents = c("Immune cells"))
p <- SpatialDimPlot(object, cells.highlight = cells[setdiff(names(cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
p & custom_theme

ggsave(paste0(dataset,"_ImmuneHighlight.png"), width = 12, height = 10)
```
We can also highlight specific cell types, starting with all cell types within a metacluster
```{r Highlight all immune cell types}
# Set idents back to fine annotation
Idents(object) <- "first_type"

# Set list of cells to highlight
cells_highlight <- levels(droplevels(factor(object$first_type[object$metaclusters == "Immune cells"], 
                                            levels = levels(immune$celltype))))

# Determine cell ids for these groups
cells <- CellsByIdentities(object, idents = cells_highlight)
cells <- cells[cells_highlight]  # Reorder the list by your desired factor levels

pal <- color("smoothrainbow")(length(cells_highlight), range = c(0.1,0.9))
names(pal) <- cells_highlight
pal <- c(pal, "Unselected" = "grey75")

# Plot
p <- SpatialDimPlot(object, images = "slice1", 
                    cells.highlight = cells[setdiff(names(cells), "NA")], 
                    cols.highlight = c(pal, "grey75"),
                    facet.highlight = F, combine = T) 
p & custom_theme & theme(panel.grid = element_blank())

ggsave_multidevice(filename_base = paste0(dataset,"_ImmuneCellTypeHighlight"), 
                   width = 10, height = 8)

# Plot
p1 <- SpatialDimPlot(object, images = "Region1", 
                     cells.highlight = cells[setdiff(names(cells), "NA")], 
                     cols.highlight = c(pal, "grey75"), 
                     facet.highlight = F, combine = T) 
p1 & custom_theme & theme(panel.grid = element_blank()) + NoLegend()

ggsave_multidevice(filename_base = paste0(dataset,"_ImmuneCellTypeHighlight_Region1"), 
                   width = 8, height = 8)

# Plot
p2 <- SpatialDimPlot(object, images = "Region2", 
                     cells.highlight = cells[setdiff(names(cells), "NA")], 
                     cols.highlight = c(pal, "grey75"), 
                     facet.highlight = F, combine = T) 
p2 & custom_theme & theme(panel.grid = element_blank()) + NoLegend()

ggsave_multidevice(filename_base = paste0(dataset,"_ImmuneCellTypeHighlight_Region2"), 
                   width = 8, height = 8)
```

```{r Highlight all immune cell types}
# Set idents back to fine annotation
Idents(object) <- "first_type"

# Set list of cells to highlight
cells_highlight <- levels(droplevels(factor(object$first_type[object$metaclusters == "IFE Keratinocytes"], 
                                            levels = levels(object$first_type))))

# Determine cell ids for these groups
cells <- CellsByIdentities(object, idents = cells_highlight)
cells <- cells[cells_highlight]  # Reorder the list by your desired factor levels

pal <- c(cell_colors, "Unselected" = "grey75")

# Plot
p <- SpatialDimPlot(object, images = "slice1", 
                    cells.highlight = cells[setdiff(names(cells), "NA")], 
                    cols.highlight = c(pal, "grey75"),
                    facet.highlight = T, combine = T) 
p & custom_theme & theme(panel.grid = element_blank())

# Plot
p1 <- SpatialDimPlot(object, images = "Region1", 
                     cells.highlight = cells[setdiff(names(cells), "NA")], 
                     cols.highlight = c(pal, "grey75"), 
                     facet.highlight = F, combine = T) 
p1 & custom_theme & theme(panel.grid = element_blank()) + NoLegend()

# Plot
p2 <- SpatialDimPlot(object, images = "Region2", 
                     cells.highlight = cells[setdiff(names(cells), "NA")], 
                     cols.highlight = c(pal, "grey75"), 
                     facet.highlight = F, combine = T) 
p2 & custom_theme & theme(panel.grid = element_blank()) + NoLegend()

```

```{r Highlight specific cell types}
# Set idents back to fine annotation
Idents(object) <- "first_type"

# Set list of cells to highlight
cells_highlight <- c("Migratory DCs","LEC I","LEC II","LEC III")

# Determine cell ids for these groups
cells <- CellsByIdentities(object, idents = cells_highlight)

# Plot
p <- SpatialDimPlot(object, images = "slice1", 
                    cells.highlight = cells[setdiff(names(cells), "NA")], 
                    cols.highlight = c(pal,"grey75"), 
                    facet.highlight = F, combine = T) 
p & custom_theme
```

Lastly we can determine broader cell type labels that are still more specific than the metacluster labels
```{r Add broad cell type labels}
# Create a named vector for renaming cell types to metaclusters
levels <- c(
  "Spinous I" = "Spinous", "Spinous II" = "Spinous", "Spinous III" = "Spinous", "Spinous IV" = "Spinous", "Basal I" = "Basal", 
  "Basal II" = "Basal", "Basal III" = "Basal", "Basal IV" = "Basal", "Cornified" = "Spinous",  # Keratinocytes
  
  "HF I" = "HF", "HF II" = "HF", "HF III" = "HF", "HF IV" = "HF", 
  "Outer Bulge I" = "Outer Bulge", "Outer Bulge II" = "Outer Bulge", "Outer Bulge III" = "Outer Bulge", 
  "Cycling I" = "HF Cycling", "Cycling II" = "HF Cycling", "Medulla I" = "Medulla", 
  "Medulla II" = "Medulla", "Cortex" = "Cortex", "Inner Rooth Sheath I" = "IRS",  # HF Keratinocytes
  
  "Transitional Basal Cells" = "Sebocytes", "Sebocyte I" = "Sebocytes", 
  "Sebocyte II" = "Sebocytes", "Sebocyte III" = "Sebocytes",  # Sebocytes
  
  "Papillary I" = "Papillary Fib", "Papillary II" = "Papillary Fib", 
  "Papillary III" = "Papillary Fib", "Papillary IV" = "Papillary Fib", 
  "Reticular I" = "Reticular Fib", "Reticular II" = "Reticular Fib", 
  "Fascia" = "Fascia Fib", "Myofibroblasts I" = "Myofibroblasts", 
  "Myofibroblasts II" = "Myofibroblasts", "Myofibroblasts III" = "Myofibroblasts", 
  "Myofibroblasts IV" = "Myofibroblasts", "Cycling" = "Cycling Fibroblasts", 
  "Lox+ Fibroblasts" = "Other Fibroblasts", "Lgr5+ Fibroblasts" = "Other Fibroblasts", 
  "Tnc+ Fibroblasts" = "Other Fibroblasts", "Myoc+ Fibroblasts" = "Other Fibroblasts",  # Fibroblasts
  
  "Mono-I" = "Monocytes", "Mono-II" = "Monocytes", "Mac-I" = "Macrophages", 
  "Mac-II" = "Macrophages", "Mac-III" = "Macrophages", "Mac-IV" = "Macrophages", 
  "Mac-V" = "Macrophages", "OC-like" = "Macrophages", "Unclear Mac-I" = "Macrophages", 
  "Unclear Mac-II" = "Macrophages", "Mono-III" = "Monocytes", "moDCs" = "DCs", 
  "cDC1" = "DCs", "cDC2" = "DCs", "migDCs" = "DCs", 
  "LCs" = "LCs", "Neu-I" = "Neutrophils", "Neu-II" = "Neutrophils", 
  "Baso" = "Granulocytes", "Mast" = "Granulocytes", "NK" = "ILCs", 
  "γδT" = "γδT", "T cells" = "T cells", "Tregs" = "T cells", 
  "ILC2" = "ILCs", "B cells" = "B cells", "AP1-hi" = "AP1-hi Immune",  # Immune cells
  
  "Proliferating ECs" = "Endothelial cells", "Intermediate ECs" = "Endothelial cells", 
  "Tacr1+ ECs" = "Endothelial cells", "Vwf-hi ECs" = "Endothelial cells", 
  "Artery ECs" = "Endothelial cells", "Capillary ECs" = "Endothelial cells", 
  "Vein ECs" = "Endothelial cells", "Pdgfra+ ECs" = "Endothelial cells", 
  "Mural cells" = "Endothelial cells", "Vascular muscle" = "Endothelial cells",  # Endothelial cells
  
  "Foxp2-hi Lymphatic ECs" = "Lymphatic Endothelial cells", 
  "Lyve1+ Lymphatic ECs" = "Lymphatic Endothelial cells", 
  "Cldn11+ Lymphatic ECs" = "Lymphatic Endothelial cells",  # Lymphatic Endothelial cells
  
  "Pericyte I" = "Pericytes", "Pericyte II" = "Pericytes", "Pericyte III" = "Pericytes", 
  "Sox6+ Pericytes" = "Pericytes", "Mt1+ Pericytes" = "Pericytes", "Actg2+ Pericytes" = "Pericytes", 
  "Pericyte Cycling" = "Pericytes",  # Pericytes
  
  "Skeletal Muscle" = "Muscle", "Myonuclei" = "Muscle",  # Muscle
  
  "Oca2-" = "Melanocytes", "Oca2-lo" = "Melanocytes", "Oca2+" = "Melanocytes",  # Melanocytes
  
  "Myl SC-I" = "Schwann cells", "Myl SC-II" = "Schwann cells", 
  "Repair SC-I" = "Schwann cells", "Repair SC-II" = "Schwann cells",  # Schwann cells
  
  "Adipocyte I" = "Adipocytes", "Adipocyte II" = "Adipocytes",  # Adipocytes
  
  "Red Blood Cells" = "Red Blood Cells",  # RBCs
  
  "Unknown" = "Unknown"  # Unknown
)

# Assign the new identifiers
Idents(object) <- "first_type"

# Rename the identifiers based on the new levels
object <- RenameIdents(object, levels)

# Optionally, add the new metacluster information to the object
object$broad_celltype <- Idents(object)

# Define the metacluster names and corresponding cell counts for broad cell types
metaclusters <- c("Keratinocytes" = 8,
                  "Sebocytes" = 1,
                  "Fibroblasts" = 6,
                  "Immune cells" = 11,
                  "Endothelial cells" = 2,
                  "Pericytes" = 1,
                  "Muscle cells" = 1,
                  "Melanocytes" = 1,
                  "Schwann cells" = 1,
                  "Adipocytes" = 1,
                  "RBCs" = 1,
                  "Unknown" = 1)

broad_celltypes <- levels(object$broad_celltype)

broadcell_colors <- unlist(mapply(function(base_color, n, cluster_name) {
  if (n == 1) return(setNames(base_color, cluster_name))
  pal <- colorway::build_palette(central_color = base_color, n_colors = n, val_range = 0.2, sat_range = -0.3)$palette
  setNames(pal, rep(cluster_name, n))
}, base_color = metacluster_colors, n = metaclusters, cluster_name = names(metaclusters)))

# Ensure cell colors are correctly named for plotting
names(broadcell_colors) <- broad_celltypes

SpatialDimPlot(object, images = "slice1") +
  scale_fill_manual(values = broadcell_colors) & custom_theme
```
Save plots for each region
```{r Plot broad cell types by region}
Idents(object) <- "broad_celltype"

SpatialDimPlot(object, images = "slice1", pt.size.factor = 1) + 
  scale_fill_manual(values = broadcell_colors) & custom_theme + theme(panel.grid = element_blank())

ggsave_multidevice(filename_base = paste0(dataset,"_BroadCellType"),
                   width = 14, height = 10)

SpatialDimPlot(object, images = "Region.1", pt.size.factor = 1) + # 1.5 for YWV01, 1 for YWV02
  scale_fill_manual(values = broadcell_colors) + NoLegend() & custom_theme + theme(panel.grid = element_blank())

ggsave_multidevice(filename_base = paste0(dataset,"_BroadCellType_Section1"),
                   width = 10, height = 10)

SpatialDimPlot(object, images = "Region.2", pt.size.factor = 1) +
  scale_fill_manual(values = broadcell_colors) + NoLegend() & custom_theme + theme(panel.grid = element_blank())

ggsave_multidevice(filename_base = paste0(dataset,"_BroadCellType_Section2"),
                   width = 10, height = 10)
```
Now that we have cell types identified, we can find marker genes in the spatial data that distinguish these groups
```{r Find fine cell type markers}
Idents(object) <- "first_type"
celltype_markers <- FindAllMarkers(object)

top_celltype_markers <- celltype_markers %>%
  group_by(cluster) %>%
  slice_head(n = 20)

feats <- top_celltype_markers %>%
  group_by(cluster) %>%
  slice_head(n = 2)

#We can reset the levels
DotPlot(object, features = unique(feats$gene), cols = "RdBu") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave(paste0(dataset,"_FineCellTypeMarker_DotPlot.png"), width = 20, height = 16)
```
Plot the two images side by side
```{r Plot seperate regions with clustering and annotation}
SpatialDimPlot(object, images = c("Region.1","Region.2"), group.by = "banksy_cluster_0.8", label = F, repel = T) & scale_fill_manual(values = viridis(length(levels(object$banksy_cluster_0.8)), begin = 0.05, option = "H")) & NoLegend()

ggsave(paste0(dataset,"_SeperatedRegions_Banksy.png"), width = 16, height = 8)

SpatialDimPlot(object, images = c("Region.1","Region.2"), group.by = "first_type", label = F, repel = T) & scale_fill_manual(values = viridis(length(levels(object$first_type)), begin = 0.05, option = "H")) & NoLegend()

ggsave(paste0(dataset,"_SeperatedRegions_Celltype.png"), width = 16, height = 8)
```

```{r Save annotated RDS}
saveRDS(object, paste0(dataset,"_FineAnnotation_",Sys.Date(),".rds"))
```

```{r}
object_uw <- plot_info$YWV01_UW$object
object_d4 <- plot_info$YWV04_D4PW$object
object_d7 <- plot_info$YWV03_D7PW$object
object_d30 <- plot_info$YWV02_D30PW$object

setwd("~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Sequencing Repository/Users/Simon/Spatial Data/RDS Files")
saveRDS(object_uw, "YWV01_UW_FineAnnotation_092525.rds")
saveRDS(object_d4, "YWV04_D4PW_FineAnnotation_092525.rds")
saveRDS(object_d7, "YWV03_D7PW_FineAnnotation_092525.rds")
saveRDS(object_d30, "YWV02_D30PW_FineAnnotation_092525.rds")
```

Load the object if desired
```{r Load RDS}
dataset <- "YWV02_D30PW"
object <- readRDS(paste0("~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Shared Data/Spatial Data/RDS Files/",dataset,"_FineAnnotation_2025-05-19.rds"))

object <- readRDS("~/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Shared Data/Spatial Data/RDS Files/YWV04_D4PW_FineAnnotation_2025-07-14.rds")
```

## Compare clusterings and annotations
View the overlaps between RCTD annotation and the various clusterings we have done
```{r Assess similarity between clustering and annotation}
Idents(object) <- "first_type"

p <- SpatialDimPlot(object) + scale_fill_manual(values = viridis(length(levels(object)), begin = 0.05, option = "H"))
p

Idents(object) <- "seurat_cluster.projected"

p1 <- SpatialDimPlot(object) + scale_fill_manual(values = viridis(length(levels(object)), begin = 0.05, option = "H"))
p1

Idents(object) <- "banksy_cluster_0.5"

p2 <- SpatialDimPlot(object) + scale_fill_manual(values = viridis(length(levels(object)), begin = 0.05, option = "H"))
p2

Idents(object) <- "banksy_cluster_split"

p2 <- SpatialDimPlot(object) + scale_fill_manual(values = viridis(length(levels(object)), begin = 0.05, option = "H"))
p2

Idents(object) <- "alignment_clusters"

p3 <- SpatialDimPlot(object) + scale_fill_manual(values = viridis(length(levels(object)), begin = 0.05, option = "H"))
p3
```
We can also use a plotting function to show how the cell types are grouped into different clusterings
```{r}
plot_cell_types <- function(data, label) {
  p <- ggplot(data, aes(x = get(label), y = n, fill = first_type)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = ifelse(n >= min_count_to_show_label, first_type, "")), position = position_stack(vjust = 0.5), size = 2) +
    xlab(label) +
    ylab("# of Spots") +
    ggtitle(paste0("Distribution of Cell Types across ", label)) +
    theme_minimal() + 
    scale_fill_manual(values = viridis(length(unique(data$first_type)), begin = 0.05, option = "H"))
}

min_count_to_show_label <- 100

cell_type_banksy0.2_counts <- object[[]] %>%
  dplyr::count(first_type, banksy_cluster_0.2)

p <- plot_cell_types(cell_type_banksy0.2_counts, "banksy_cluster_0.2")
p
ggsave(paste0(dataset,"_Banksy0.2ClusterCelltypeBarplot.png"), width = 10, height = 6)

cell_type_banksy0.5_counts <- object[[]] %>%
  dplyr::count(first_type, banksy_cluster_0.5)

p <- plot_cell_types(cell_type_banksy0.5_counts, "banksy_cluster_0.5")
p
ggsave(paste0(dataset,"_Banksy0.5ClusterCelltypeBarplot.png"), width = 10, height = 6)

cell_type_banksy0.8_counts <- object[[]] %>%
  dplyr::count(first_type, banksy_cluster_0.8)

p <- plot_cell_types(cell_type_banksy0.8_counts, "banksy_cluster_0.8")
p
ggsave(paste0(dataset,"_Banksy0.8ClusterCelltypeBarplot.png"), width = 10, height = 6)

cell_type_seurat_counts <- object[[]] %>%
  dplyr::count(first_type, seurat_cluster.projected)

p <- plot_cell_types(cell_type_seurat_counts, "seurat_cluster.projected")
p
ggsave(paste0(dataset,"_SeuratClusterCelltypeBarplot.png"), width = 10, height = 6)

cell_type_10x_counts <- object[[]] %>%
  dplyr::count(first_type, alignment_clusters)

p <- plot_cell_types(cell_type_10x_counts, "alignment_clusters")
p
ggsave(paste0(dataset,"_10XClusterCelltypeBarplot.png"), width = 10, height = 6)
```
Plot a heatmap of marker genes by predicted cell type 
```{r}
# Crete downsampled object to make visualization easier
DefaultAssay(object) <- "Spatial"
Idents(object) <- "first_type"
object_subset <- subset(object, cells = Cells(object[["Spatial"]]), downsample = 1000)

# # Order clusters by similarity
# DefaultAssay(object_subset) <- "Spatial"
# Idents(object_subset) <- "first_type"
# object_subset <- BuildClusterTree(object_subset, assay = "Spatial", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "Spatial", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p
```
Using the same marker genes, we can show how well the different clustering methods group together the cell type markers
```{r}
DefaultAssay(object) <- "Spatial"
Idents(object) <- "first_type"
celltype_markers <- FindAllMarkers(object)

top_celltype_markers <- celltype_markers %>%
  group_by(cluster) %>%
  slice_head(n = 20)

feats <- top_celltype_markers %>%
  group_by(cluster) %>%
  slice_head(n = 4)

# List of clustering methods
clustering_methods <- c("alignment_clusters", "banksy_cluster", "seurat_cluster.projected")

# Iterate over each clustering method
for (cluster_method in clustering_methods) {
  
  # Set the active identity to the current clustering method
  Idents(object) <- cluster_method
  
  # Subset the object for ease of plotting
  object_subset <- subset(object, cells = Cells(object[["Spatial"]]), downsample = 10000)
  object_subset <- ScaleData(object_subset, assay = "Spatial", features = feats$gene)
  
  # Print the clustering method being used
  print(paste("Plotting heatmap for:", cluster_method))
  
  # Generate the heatmap
  p <- DoHeatmap(object_subset, features = feats$gene, size = 4, assay = "Spatial") +
    ggtitle(paste("Marker Gene Heatmap -", cluster_method)) +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
  
  # Define the file name dynamically
  file_name <- paste0(dataset,"_", cluster_method, "_CellTypeMarkerHeatmap.png")
  
  # Save the plot using ggsave
  ggsave(filename = file_name, plot = p, width = 15, height = 12, dpi = 300)
  
  # Print confirmation
  print(paste("Saved heatmap to:", file_name))
  # Print the plot
  print(p)
}
```
To quantify the clusterings, we can use normalized mutual information (NMI) to determine how much agreement there is between the two groupings
```{r}
# For alignment_clusters
alignment_table <- table(object$alignment_clusters, object$first_type)

# For banksy clusterings
banksy0.2_table <- table(object$banksy_cluster_0.2, object$first_type)

banksy0.5_table <- table(object$banksy_cluster_0.5, object$first_type)

banksy0.8_table <- table(object$banksy_cluster_0.8, object$first_type)

# For seurat_cluster.projected
seurat_projected_table <- table(object$seurat_cluster.projected, object$first_type)

alignment_prop <- prop.table(alignment_table, margin = 1)  # Row-wise proportions
```

```{r}
library(ggplot2)

# Melt the data for visualization
library(reshape2)

# Alignment Clusters
alignment_melt <- as.data.frame(as.table(alignment_table))
colnames(alignment_melt) <- c("Cluster", "CellType", "Count")

ggplot(alignment_melt, aes(x = Cluster, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Cell Type Distribution in Alignment Clusters",
       x = "Alignment Clusters", y = "Proportion") +
  theme_minimal() + 
  scale_fill_manual(values = viridis(length(unique(alignment_melt$CellType)), begin = 0.05, option = "H"))

ggsave(paste0(dataset,"_10XClusterCelltypeProportion.png"), width = 10, height = 8)

# Repeat for Banksy Clusters
banksy0.2_melt <- as.data.frame(as.table(banksy0.2_table))
colnames(banksy0.2_melt) <- c("Cluster", "CellType", "Count")

ggplot(banksy0.2_melt, aes(x = Cluster, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Cell Type Distribution in Banksy Clusters - λ = 0.2",
       x = "Banksy Clusters - λ = 0.2", y = "Proportion") +
  theme_minimal() + 
  scale_fill_manual(values = viridis(length(unique(banksy0.2_melt$CellType)), begin = 0.05, option = "H"))

ggsave(paste0(dataset,"_Banksy0.2ClusterCelltypeProportion.png"), width = 10, height = 8)

banksy0.5_melt <- as.data.frame(as.table(banksy0.5_table))
colnames(banksy0.5_melt) <- c("Cluster", "CellType", "Count")

ggplot(banksy0.5_melt, aes(x = Cluster, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Cell Type Distribution in Banksy Clusters - λ = 0.5",
       x = "Banksy Clusters - λ = 0.5", y = "Proportion") +
  theme_minimal() + 
  scale_fill_manual(values = viridis(length(unique(banksy0.5_melt$CellType)), begin = 0.05, option = "H"))

ggsave(paste0(dataset,"_Banksy0.5ClusterCelltypeProportion.png"), width = 10, height = 8)

banksy0.8_melt <- as.data.frame(as.table(banksy0.8_table))
colnames(banksy0.8_melt) <- c("Cluster", "CellType", "Count")

ggplot(banksy0.8_melt, aes(x = Cluster, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Cell Type Distribution in Banksy Clusters - λ = 0.8",
       x = "Banksy Clusters - λ = 0.8", y = "Proportion") +
  theme_minimal() + 
  scale_fill_manual(values = viridis(length(unique(banksy0.8_melt$CellType)), begin = 0.05, option = "H"))

ggsave(paste0(dataset,"_Banksy0.8ClusterCelltypeProportion.png"), width = 10, height = 8)

# Repeat for Seurat Projected Clusters
seurat_projected_melt <- as.data.frame(as.table(seurat_projected_table))
colnames(seurat_projected_melt) <- c("Cluster", "CellType", "Count")

ggplot(seurat_projected_melt, aes(x = Cluster, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Cell Type Distribution in Seurat Projected Clusters",
       x = "Seurat Projected Clusters", y = "Proportion") +
  theme_minimal() + 
  scale_fill_manual(values = viridis(length(unique(seurat_projected_melt$CellType)), begin = 0.05, option = "H"))

ggsave(paste0(dataset,"_SeuratClusterCelltypeProportion.png"), width = 10, height = 8)

```
To calculate agreement between clustering and annotation, we will use Normalized Mutual Information (NMI) analysis.
```{r}
if(!requireNamespace("aricode")){
  install.packages("aricode")
}

library(aricode)

# Extract cluster labels and cell types as vectors
alignment_vector <- as.character(object$alignment_clusters)
banksy0.2_vector <- as.character(object$banksy_cluster_0.2)
banksy0.5_vector <- as.character(object$banksy_cluster_0.5)
banksy0.8_vector <- as.character(object$banksy_cluster_0.8)
seurat_vector <- as.character(object$seurat_cluster.projected)
cell_type_vector <- as.character(object$first_type)

# Calculate NMI for each clustering method
alignment_NMI <- NMI(alignment_vector, cell_type_vector)
banksy0.2_NMI <- NMI(banksy0.2_vector, cell_type_vector)
banksy0.5_NMI <- NMI(banksy0.5_vector, cell_type_vector)
banksy0.8_NMI <- NMI(banksy0.8_vector, cell_type_vector)
seurat_NMI <- NMI(seurat_vector, cell_type_vector)

# Print results
print(paste("NMI for Alignment Clusters vs Cell Types:", round(alignment_NMI, 3)))
print(paste("NMI for Banksy Clusters (0.2 λ) vs Cell Types:", round(banksy0.2_NMI, 3)))
print(paste("NMI for Banksy Clusters (0.5 λ) vs Cell Types:", round(banksy0.5_NMI, 3)))
print(paste("NMI for Banksy Clusters (0.8 λ) vs Cell Types:", round(banksy0.8_NMI, 3)))
print(paste("NMI for Seurat Projected Clusters vs Cell Types:", round(seurat_NMI, 3)))

nmi_results <- data.frame(
  Method = c("Alignment Clusters", "Banksy Clusters - λ = 0.2", "Banksy Clusters - λ = 0.5","Banksy Clusters - λ = 0.8", "Seurat Projected Clusters"),
  NMI = c(alignment_NMI, banksy0.2_NMI, banksy0.5_NMI, banksy0.8_NMI, seurat_NMI)
)

# Plot the NMI results
ggplot(nmi_results, aes(x = Method, y = NMI, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = round(NMI, 3)), vjust = -0.5, size = 4) +
  labs(title = "Normalized Mutual Information (NMI) for Clusterings vs Cell Types",
       x = "Clustering Method", y = "NMI") +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +  # Optional color palette
  theme(legend.position = "none")

#ggsave(paste0(dataset, "_ClusteringAnnotationNMI.png"), width = 10, height = 6)
```
We can also examine the cell type annotations overlayed on the UMAP generated using BANKSY
```{r}
DimPlot(data.ref, reduction = "umap.banksy", group.by = "first_type", label = F, raster = F) + ggtitle("BANKSY Clustering") + theme(legend.position = "bottom")
```
