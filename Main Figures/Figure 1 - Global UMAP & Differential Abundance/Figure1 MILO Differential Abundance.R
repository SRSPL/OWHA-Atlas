#This code runs differntial abundance for the OWHA atlas.
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r setup, include=FALSE}
#Set seed and directories
set.seed(123)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir ="/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 2")
#Load Libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Rfast2)
library(viridis)
library(scCustomize)
library(harmony)
library(miloR)
library(SingleCellExperiment)
library(scater)                    
```
Load dataset to analyze
```{r Load data}
#alldata = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 1-2/PSUUnit.rds")

alldata = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/scRNA SEQ of wounds/Skin Wound Healing Atlas/2. Integration of Data/RDS files/IntegratedMultimodal_FineCelltypes_102124.rds")
```
Rename existing dataset if desired
```{r, include = FALSE}
library(dplyr)
#Add new modality group 
alldata$modality_group <- ifelse(alldata$modality %in% c("scRNAseq", "CITEseq"), "whole cell", "nuclei")

# Add a new metadata column called 'condition' based on 'timepoint'
alldata@meta.data <- alldata@meta.data %>%
  mutate(condition = case_when(
    timepoint == "Unwounded" ~ "Unwounded",
    timepoint %in% c("Wounded_D1PW", "Wounded_D2PW", "Wounded_D4PW", "Wounded_D7PW", "Wounded_D15PW", "Wounded_D30PW") ~ "Wounded",
    TRUE ~ as.character(timepoint)  # Keeps any other potential values unchanged
  ))
```
Format data as singleCellExperiment object and create the milo object
```{r}
DefaultAssay(alldata) <- "SCT"
#Subset only necessary data
alldata.diet <- DietSeurat(alldata, assay = "SCT", dimreducs = c("pca","umap.rpca"), graph = "pca")

alldata.diet$nCount_ADT <- NULL
alldata.diet$nFeature_ADT <- NULL
alldata <- as.SingleCellExperiment(alldata.diet)

alldata.milo <- Milo(alldata)
reducedDim(alldata.milo, "UMAP") <- reducedDim(alldata, "UMAP.RPCA")
alldata.milo
```
```{r}
plotUMAP(alldata.milo,  colour_by="metaclusters")
```

Once we have defined neighbourhoods, it’s good to take a look at how big the neighbourhoods are (i.e. how many cells form each neighbourhood). This affects the power of DA testing. We can check this out using the plotNhoodSizeHist function. Empirically, we found it’s best to have a distribution peaking above 20. Otherwise you might consider rerunning makeNhoods increasing k and/or prop.
```{r}
alldata.milo <- buildGraph(alldata.milo, k = 30, d = 30)
alldata.milo <- makeNhoods(alldata.milo, prop = 0.10, k = 30, d=30, refined = TRUE, refinement_scheme = "graph")
plotNhoodSizeHist(alldata.milo)
```

```{r}
alldata.milo <- countCells(alldata.milo, meta.data = data.frame(colData(alldata.milo)), sample = "orig.ident")
head(nhoodCounts(alldata.milo))
```
Defining experimental design
Now we are all set to test for differential abundance in neighbourhoods. We implement this hypothesis testing in a generalized linear model (GLM) framework, specifically using the Negative Binomial GLM implementation in edgeR.

We first need to think about our experimental design. The design matrix should match each sample to the experimental condition of interest for DA testing. In this case, we want to detect DA between embryonic stages, stored in the stage column of the dataset colData. We also include the sequencing.batch column in the design matrix. This represents a known technical covariate that we want to account for in DA testing.

```{r}
#Extract metadata
alldata.design <- data.frame(colData(alldata))[,c("orig.ident", "timepoint", "condition")]
alldata.design <- distinct(alldata.design)
## Convert batch info from integer to factor
rownames(alldata.design) <- alldata.design$orig.ident
#Reorder rownames to match columns of nhoodCounts
alldata.design <- alldata.design[colnames(nhoodCounts(alldata.milo)), , drop=FALSE]
table(alldata.design$timepoint)
table(alldata.design$condition)
```
Computing neighbourhood connectivity
Milo uses an adaptation of the Spatial FDR correction introduced by cydar, which accounts for the overlap between neighbourhoods. Specifically, each hypothesis test P-value is weighted by the reciprocal of the kth nearest neighbour distance. To use this statistic we first need to store the distances between nearest neighbors in the Milo object. This is done by the calcNhoodDistance function (N.B. this step is the most time consuming of the analysis workflow and might take a couple of minutes for large datasets).

This might not need to be used, as the calculation of neighborhoods using graph correction may make this step unnecessary: https://github.com/MarioniLab/miloR/issues/314
```{r, results='hide'}
alldata.milo <- calcNhoodDistance(alldata.milo, d=30)
```
Testing
Now we can do the DA test, explicitly defining our experimental design. In this case, we want to dest for differences between experimental stages, while accounting for the variability between technical batches (You can find more info on how to use formulas to define a testing design in R here)

First, we can compare conditions specified in the contrast.1 object
```{r}
#To specify comparisons for design - must be written as <ColumnName><Condition1> - <ColumnName><Condition2> and will plot as what is up in condition 1
contrast.1 <- c("conditionWounded - conditionUnwounded")
da_results <- testNhoods(alldata.milo, design = ~ 0 + condition, design.df = alldata.design, model.contrasts = contrast.1, fdr.weighting = "graph-overlap")
#This calculates a Fold-change and corrected P-value for each neighbourhood, which indicates wheather there is significant differential abundance between conditions

da_results %>%
  arrange(SpatialFDR) %>%
  head()
```
We can also compare individual timepoints
```{r}
#To specify comparisons for design - must be written as <ColumnName><Condition1> - <ColumnName><Condition2> and will plot as what is up in condition 1
contrast.1 <- c("timepointWounded_D7PW - timepointWounded_D4PW")
da_results <- testNhoods(alldata.milo, design = ~ 0 + timepoint, design.df = alldata.design, model.contrasts = contrast.1, fdr.weighting = "graph-overlap")
#This calculates a Fold-change and corrected P-value for each neighbourhood, which indicates wheather there is significant differential abundance between conditions

da_results %>%
  arrange(SpatialFDR) %>%
  head()

#To specify comparisons for design - must be written as <ColumnName><Condition1> - <ColumnName><Condition2> and will plot as what is up in condition 1
contrast.2 <- c("timepointWounded_D15PW - timepointWounded_D7PW")
da_results.2 <- testNhoods(alldata.milo, design = ~ 0 + timepoint, design.df = alldata.design, model.contrasts = contrast.2, fdr.weighting = "graph-overlap")
#This calculates a Fold-change and corrected P-value for each neighbourhood, which indicates wheather there is significant differential abundance between developmental stages.




da_results.2 %>%
  arrange(SpatialFDR) %>%
  head()
```
We can also specify multiple contrasts to be made simultaneously
```{r}
#Specify contrasts using same format, now just in a list
contrast.all <- c("timepointUnwounded - timepointWounded_D4PW", 
                  "timepointWounded_D4PW - timepointWounded_D7PW",
                  "timepointWounded_D7PW - timepointWounded_D15PW",
                  "timepointWounded_D15PW - timepointWounded_D30PW",
                  "timepointWounded_D30PW - timepointUnwounded")

contrast.res <- testNhoods(alldata.milo, design= ~0 + timepoint, design.df=alldata.design, fdr.weighting="graph-overlap", model.contrasts = contrast.all)

head(contrast.res)
```

Inspecting DA testing results
We can start inspecting the results of our DA analysis from a couple of standard diagnostic plots. We first inspect the distribution of uncorrected P values, to verify that the test was balanced.
```{r}
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
```
To visualize results, build an abstracted graph of neighbourhoods which can be superimposed on single cell embedding
```{r}
alldata.milo <- buildNhoodGraph(alldata.milo)
```

Plot the results on a UMAP
```{r}
plotUMAP(alldata.milo) + plotNhoodGraphDA(alldata.milo, da_results, alpha=0.05 ) +
  plot_layout(guides="collect")

ggsave("UMAPGlobalMiloDifferentialExpressionRPCA.png", width = 9)
```

```{r}
da_results <- annotateNhoods(alldata.milo, da_results, coldata_col = "metaclusters")
head(da_results)
```

```{r}

second_level_terms <- unlist(metacluster_groups)
second_level_terms <- unlist(second_level_terms)

#Set the factor levels for plotting
da_results$metaclusters <- factor(da_results$metaclusters,
                                  levels = rev(alldata$metaclusters))

da_results$metaclusters <- factor(da_results$metaclusters)
plotDAbeeswarm(da_results, group.by = "metaclusters", alpha = 0.05) + scale_color_gradient2(low = "blue", high = "red")


ggsave("GlobalDifferentialAbundanceUnwoundedvsWounded.png", height = 7)


#plotDAbeeswarm(da_results, group.by = "metaclusters", alpha = 0.2) +  scale_color_grey(
na.value = "grey50")   # color for NA values
```