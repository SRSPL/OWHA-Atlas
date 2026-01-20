
  
  ```{r}
set.seed(123)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 3")
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

Load Data Sets

```{r}
alldata = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/scRNA SEQ of wounds/Skin Wound Healing Atlas/3. Subclustering/Fibroblast Subclustering/FibroblastsAnnotated_2025-06-09.rds")


alldata$subclusters <- Idents(alldata)
table(alldata$timepoint)
Idents(alldata) <- "timepoint"
alldata <- RenameIdents(alldata,
                        "Unwounded" = "UW",
                        "Wounded_D1PW" = "D1PW",
                        "Wounded_D2PW" = "D2PW",
                        "Wounded_D4PW" = "D4PW",
                        "Wounded_D7PW" = "D7PW",
                        "Wounded_D15PW" = "D15PW",
                        "Wounded_D30PW" = "D30PW")
alldata$timepoint <- Idents(alldata)
Idents(alldata) <- alldata$subclusters
```

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 3/Cell Chat Data")
#Install required packages as needed
# Check and install CRAN package
if (!requireNamespace("NMF", quietly = TRUE)) {
  install.packages("NMF")
}
# Check and install GitHub packages
if (!requireNamespace("CellChat", quietly = TRUE)) {
  devtools::install_github("jinworks/CellChat")
}
if (!requireNamespace("circlize", quietly = TRUE)) {
  devtools::install_github("jokergoo/circlize")
}
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  devtools::install_github("jokergoo/ComplexHeatmap")
}

library(CellChat)
library(Seurat)
library(patchwork)
library(viridis)


options(future.globals.maxSize = 10000*1024^2)
```

Prep the cellchat database

```{r Prep Data and load cellchat DB}
CellChatDB <- CellChatDB.mouse

showDatabaseCategory(CellChatDB)
#This removes non-protein interactions that are contained in the CellChat DB such as synaptic interactions and metabolic signals
CellChatDB.use <- subsetDB(CellChatDB)
```

Subset the seurat object by timepoint for comparisons

```{r pressure, echo=FALSE}

seuratobj <- readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/scRNA SEQ of wounds/Skin Wound Healing Atlas/IntegratedMultimodal_FineCellTypes_061225.rds")

colors_pals = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 1/color_palette.rds")


Idents(seuratobj) <- seuratobj$metaclusters

seuratobj$subclusters <- Idents(seuratobj)


seuratobj@meta.data$timepoint <- recode(seuratobj@meta.data$timepoint,
                                        "Unwounded" = "UW",
                                        "Wounded_D1PW" = "D1PW",
                                        "Wounded_D2PW" = "D2PW",
                                        "Wounded_D4PW" = "D4PW",
                                        "Wounded_D7PW" = "D7PW",
                                        "Wounded_D15PW" = "D15PW",
                                        "Wounded_D30PW" = "D30PW"
)


Idents(seuratobj) <- seuratobj$timepoint
seuratobj <- subset(seuratobj, idents = c("D1PW","D2PW"), invert = T)
Idents(seuratobj) <- seuratobj$metaclusters)

seuratobj$subclusters <- Idents(seuratobj)

cell_type <- "Global"  # Change this variable as needed
Idents(seuratobj) <- "subclusters"
#Cellchat requires that the batches are stored in a column called samples - for ease of cellchat object generation we can just create a new column with this info
seuratobj <- AddMetaData(seuratobj, seuratobj$orig.ident, col.name = "samples")
seuratobj.uw <- subset(seuratobj, timepoint == 'Unwounded')
seuratobj.d4 <- subset(seuratobj, timepoint == 'Wounded_D4PW')
seuratobj.d7 <- subset(seuratobj, timepoint == 'Wounded_D7PW')
seuratobj.d15 <- subset(seuratobj, timepoint == 'Wounded_D15PW')
seuratobj.d30 <- subset(seuratobj, timepoint == 'Wounded_D30PW')

seuratobj.uw <- subset(seuratobj, timepoint == 'UW')
seuratobj.d1 <- subset(seuratobj, timepoint == 'D1PW')
seuratobj.d2 <- subset(seuratobj, timepoint == 'D2PW')
seuratobj.d4 <- subset(seuratobj, timepoint == 'D4PW')
seuratobj.d7 <- subset(seuratobj, timepoint == 'D7PW')
seuratobj.d15 <- subset(seuratobj, timepoint == 'D15PW')
seuratobj.d30 <- subset(seuratobj, timepoint == 'D30PW')
```

Run cellchat on each timepoint object

```{r}
grouping <- "subclusters"
#Set up parallelization - there is currently an issue with the futures interaction with CellChat that will produce RNG errors
#future::plan("multisession", workers = 4)
#future.seed=TRUE
#Set the length for how many objects to make
times <- 1:length(unique(seuratobj$timepoint))
#times <-1:5
#Provide timepoint names which will be used to pull data from seurat objects
names(times) <- c('uw','d4','d1','d2','d7',"d15","d30")
for (j in times) {
  cellchat <- createCellChat(get(paste0('seuratobj.', names(times)[times == j])), group.by = grouping, assay = "RNA")
  levels(cellchat@idents)
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  #projects data onto a protein protein interaction network to smooth gene expression based on neighbors
  cellchat <- smoothData(cellchat, adj = PPI.mouse) 
  #Alter factor levels for cell type labels based on presence in particular dataset 
  cellchat@meta[[grouping]] <- droplevels(cellchat@meta[[grouping]], exclude = setdiff(levels(cellchat@meta[[grouping]]), unique(cellchat@meta[[grouping]])))
  cellchat <- setIdent(cellchat, ident.use = grouping)
  #runs the cellchat communication algorithm
  cellchat <- computeCommunProb(cellchat, population.size = TRUE, type = "triMean")
  #filters out communications between groups of less than 10 cells
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  #summarizes ligand/receptor interactions to compute communication probability
  cellchat <- computeCommunProbPathway(cellchat)
  #creates aggregated communication netwrok by counting links resulting from above
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  assign(paste0('cellchat',names(times)[times == j]), cellchat)
}
#If these do not all match, we need to lift the objects that are missing celltypes
levels(cellchatuw@idents)
levels(cellchatd4@idents)
levels(cellchatd7@idents)
levels(cellchatd15@idents)
levels(cellchatd30@idents)
```

```{r}
#create a list of the cellchat objects - this is needed for downstream applications
object.list <- list(UW = cellchatuw, D4PW = cellchatd4, D7PW =cellchatd7, D15PW = cellchatd15, D30PW = cellchatd30)

object.list <- list(UW = cellchatuw, D1PW = cellchatd1, D2PW = cellchatd2, D4PW = cellchatd4, D7PW =cellchatd7, D15PW = cellchatd15, D30PW = cellchatd30)

#merge them together into one object
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#rm(cellchatuw,cellchatd4,cellchatd7,cellchatd15,cellchatd30, seuratobj, seuratobj.uw, seuratobj.d4, seuratobj.d7, seuratobj.d15, seuratobj.d30)

# Get unique timepoint labels


```

```{r}
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4,5))
gg2 <- compareInteractions(cellchat, show.legend = F, measure = "weight", group = c(1,2,3,4,5))
gg1 + gg2
ggsave("IntegratedMultimodal_Metacluster_Interaction_strengths_barplot.tiff", width = 12, height = 7)
```

```{r}
idents_used <- levels(object.list[[1]]@idents)  # Ensure color vector matches only the identities used
color.use <- colors_pals[idents_used]
cell_type <- "BasalIVGlobal"


color.use <- metacluster_colors[idents_used]


weight.max <- getMaxWeight(object.list, attribute = c("idents", "weight"))

for (i in 1:length(object.list)) {
  # Define file name including the cell type
  filename <- paste0("CellChat_CirclePlot_", cell_type, "_", names(object.list)[i], ".svg")
  
  # Start PNG device
  png(filename, width = 1600, height = 1600, res = 300)
  
  # Set plot parameters
  par(mfrow = c(1, 1), xpd = TRUE)
  
  # Get matching colors for this object
  idents_i <- levels(object.list[[i]]@idents)
  color_i <- metacluster_colors[idents_i]
  
  # Draw the plot
  netVisual_circle(
    object.list[[i]]@net$weight,
    color.use = color_i,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight.max[2],
    edge.width.max = 12,
    arrow.width = 0.75,
    arrow.size = 0.05,
    title.name = paste0("Weight of interactions - ", names(object.list)[i])
  )
  
  # Close PNG device
  dev.off()
}

```


```{r}
# Define your metacluster colors
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

# Get the max weights and max interaction counts across all objects
weight.max <- getMaxWeight(object.list, attribute = c("idents", "weight"))
count.max  <- getMaxWeight(object.list, attribute = c("idents", "count"))

# Loop over objects to generate plots
for (i in seq_along(object.list)) {
  
  # Define file-safe label
  name_i <- gsub(" ", "_", names(object.list)[i])
  celltype_label <- gsub(" ", "_", cell_type)
  
  # Colors for the identities present in this object
  idents_i <- levels(object.list[[i]]@idents)
  color_i <- metacluster_colors[idents_i]
  
  ## ---- Plot 1: Weight of interactions ---- ##
  
  # PNG
  png(paste0("CellChat_CirclePlot_Weight_", celltype_label, "_", name_i, ".png"),
      width = 1600, height = 1600, res = 300)
  par(mfrow = c(1, 1), xpd = TRUE)
  netVisual_circle(
    object.list[[i]]@net$weight,
    color.use = color_i,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight.max[2] * 1.1,   # slightly higher max
    edge.width.max = 15,                     # thicker strong edges
    arrow.width = 0.5,                       # thinner weak edges
    arrow.size = 0.03,                       # smaller arrows
    title.name = paste0("Weight of interactions - ", names(object.list)[i])
  )
  dev.off()
  
  # SVG
  svg(paste0("CellChat_CirclePlot_Weight_", celltype_label, "_", name_i, ".svg"),
      width = 8, height = 8)
  par(mfrow = c(1, 1), xpd = TRUE)
  netVisual_circle(
    object.list[[i]]@net$weight,
    color.use = color_i,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight.max[2] * 1.1,
    edge.width.max = 15,
    arrow.width = 0.5,
    arrow.size = 0.03,
    title.name = paste0("Weight of interactions - ", names(object.list)[i])
  )
  dev.off()
  
  ## ---- Plot 2: Number of interactions ---- ##
  
  # PNG
  png(paste0("CellChat_CirclePlot_Count_", celltype_label, "_", name_i, ".png"),
      width = 1600, height = 1600, res = 300)
  par(mfrow = c(1, 1), xpd = TRUE)
  netVisual_circle(
    object.list[[i]]@net$count,
    color.use = color_i,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = count.max[2] * 1.1,
    edge.width.max = 15,
    arrow.width = 0.5,
    arrow.size = 0.03,
    title.name = paste0("Number of interactions - ", names(object.list)[i])
  )
  dev.off()
  
  # SVG
  svg(paste0("CellChat_CirclePlot_Count_", celltype_label, "_", name_i, ".svg"),
      width = 8, height = 8)
  par(mfrow = c(1, 1), xpd = TRUE)
  netVisual_circle(
    object.list[[i]]@net$count,
    color.use = color_i,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = count.max[2] * 1.1,
    edge.width.max = 15,
    arrow.width = 0.5,
    arrow.size = 0.03,
    title.name = paste0("Number of interactions - ", names(object.list)[i])
  )
  dev.off()
}
```


```{r}
library(patchwork)
library(ggplot2)

# Initialize plot lists
gg_with_labels <- list()
gg_no_labels <- list()

# Loop over each CellChat object
for (i in seq_along(object.list)) {
  # Plot with labels
  p1 <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    title = names(object.list)[i],
    label.size = 5,
    font.size = 12
  ) + scale_color_manual(values = colors_pals)
  
  # Plot without labels
  p2 <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    title = names(object.list)[i],
    label.size = 0,
    font.size = 12
  ) + scale_color_manual(values = colors_pals)
  
  # Optionally remove text geoms if still present (safety)
  p2 <- p2 + theme(legend.text = element_text(size = 10)) +
    guides(label = "none")
  
  # Save plots
  gg_with_labels[[i]] <- p1
  gg_no_labels[[i]] <- p2
}

# Combine with patchwork
g_with_labels <- wrap_plots(plots = gg_with_labels, ncol = 5)
g_no_labels   <- wrap_plots(plots = gg_no_labels,   ncol = 5)

# Save both versions
ggsave(
  paste0("IntegratedMultimodal_CellChat2DSignaling_", cell_type, "_with_labels.png"),
  g_with_labels,
  width = 30, height = 5, dpi = 300, limitsize = F
)

ggsave(
  paste0("IntegratedMultimodal_CellChat2DSignaling_", cell_type, "_no_labels.png"),
  g_no_labels,
  width = 30, height = 5, dpi = 300
)



```

```{r}
library(ComplexHeatmap)
library(CellChat)

# Define cell type you're analyzing (this gets inserted into filenames)
cell_type <- ""

# Define sample names
sample_names <- c("UW", "D4PW", "D7PW", "D15PW", "D30PW")


# Loop over each sample in the object list
for (sample_name in names(object.list)) {
  message("Processing: ", sample_name)
  
  # Get the CellChat object for the current sample
  obj <- object.list[[sample_name]]
  
  # Create the outgoing and incoming signaling role heatmaps
  h1 <- netAnalysis_signalingRole_heatmap(obj, pattern = "outgoing", height = 20, width = 10)
  h2 <- netAnalysis_signalingRole_heatmap(obj, pattern = "incoming", height = 20, width = 10)
  
  # Define the output filename, inserting sample name and cell_type
  outfile <- paste0(sample_name, "_", cell_type, "_SignalingRoleHeatmaps.png")
  
  # Save the combined heatmap to a PNG file
  png(outfile, width = 4000, height = 3000, res = 300)
  draw(h1 + h2, annotation_legend_side = "right")
  dev.off()
}
# Save the full object list with the cell type in the name
saveRDS(object.list, file = paste0("cellchat_", cell_type, ".rds"))

```

```{r}
# Save ligand-receptor pairs for each sample
lr_pairs <- cellchatuw@LR
write.csv(lr_pairs, paste0(cell_type, "_UW_ligand_receptor_pairs.csv"), row.names = FALSE)

lr_pairs <- cellchatd4@LR
write.csv(lr_pairs, paste0(cell_type, "_D4PW_ligand_receptor_pairs.csv"), row.names = FALSE)

lr_pairs <- cellchatd7@LR
write.csv(lr_pairs, paste0(cell_type, "_D7PW_ligand_receptor_pairs.csv"), row.names = FALSE)

lr_pairs <- cellchatd15@LR
write.csv(lr_pairs, paste0(cell_type, "_D15PW_ligand_receptor_pairs.csv"), row.names = FALSE)

lr_pairs <- cellchatd30@LR
write.csv(lr_pairs, paste0(cell_type, "_D30PW_ligand_receptor_pairs.csv"), row.names = FALSE)

```


```{r}
library(dplyr)
library(scatterplot3d)

# Load CellChat objects
cellchat = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 3/Cell Chat Data/cellchat_IFEEndoSubclusters.rds")
cellchat = readRDS("/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 3/Cell Chat Data/RDS Files/To Use/UPDATE/cellchat_GlobalD1PWD2PW.rds")


object.list <- list(
  UW = cellchat$UW, 
  D4PW = cellchat$D4PW, 
  D7PW = cellchat$D7PW, 
  D15PW = cellchat$D15PW, 
  D30PW = cellchat$D30PW
)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# --- Extract centrality function (same as before)
extract_signaling_roles_df <- function(object,
                                       signaling = NULL,
                                       slot.name = "netP",
                                       x.measure = "outdeg",
                                       y.measure = "indeg") {
  centr <- slot(object, slot.name)$centr
  if (length(centr) == 0) {
    stop("Run `netAnalysis_computeCentrality` first.")
  }
  if (!all(c(x.measure, y.measure) %in% names(centr[[1]]))) {
    stop(paste0("x.measure and y.measure must be one of: ",
                paste(names(centr[[1]]), collapse = ", ")))
  }
  
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  
  for (i in seq_along(centr)) {
    outgoing[, i] <- centr[[i]][[x.measure]]
    incoming[, i] <- centr[[i]][[y.measure]]
  }
  
  if (!is.null(signaling)) {
    signaling <- signaling[signaling %in% names(centr)]
    outgoing <- outgoing[, signaling, drop = FALSE]
    incoming <- incoming[, signaling, drop = FALSE]
  }
  
  outgoing.cells <- rowSums(outgoing)
  incoming.cells <- rowSums(incoming)
  
  num.link <- aggregateNet(object, signaling = signaling, return.object = FALSE, remove.isolate = FALSE)$count
  num.link <- rowSums(num.link) + colSums(num.link) - diag(num.link)
  
  df <- data.frame(
    x = outgoing.cells,
    y = incoming.cells,
    labels = names(incoming.cells),
    Count = num.link
  )
  df$labels <- factor(df$labels, levels = names(incoming.cells))
  return(df)
}

all_roles_df <- lapply(seq_along(object.list), function(i) {
  df <- extract_signaling_roles_df(object.list[[i]])
  df$sample <- names(object.list)[i]
  return(df)
}) %>% bind_rows()


write.csv(all_roles_df, "Global3DScatterData.csv")



# Map numeric timepoints
time_map <- c("UW" = 0, "D1PW" = 1, "D2PW" = 2, "D4PW" = 4,
              "D7PW" = 7, "D15PW" = 15, "D30PW" = 30)
all_roles_df$time_numeric <- time_map[all_roles_df$sample]

# Remove unwanted labels if needed
all_roles_df <- all_roles_df[!(all_roles_df$labels == "Fascia"), ]

# ============================
# Assign custom colors
# ============================
my_colors <- c(
  "#FF0000", "#00CED1", "#0000FF", "#FFFF00",
  "#8B008B", "#008B8B", "#FFA500", "#8000FF",
  "#006400", "#C71585", "#7FFF00", "#1E90FF",
  "#FF69B4", "#FFD700", "lightpink", "#A9A9A9"
)

# Map labels to colors (cycling through if more labels than colors)
label_levels <- unique(all_roles_df$labels)
color_map <- setNames(rep(my_colors, length.out = length(label_levels)), label_levels)

all_roles_df$col <- color_map[all_roles_df$labels]
write.csv(all_roles_df, "Global3DScatterData.csv")
# ============================
# Plot function
# ============================
plot_cellchat3d <- function(add_legend = FALSE) {
  s3d <- scatterplot3d(
    x = all_roles_df$time_numeric,
    y = all_roles_df$x,
    z = all_roles_df$y,
    color = all_roles_df$col,
    pch = 19,
    cex.symbols = 2,
    angle = 55,
    scale.y = 1,
    xlab = "Time",
    ylab = "Sender Score",
    zlab = "Receiver Score",
    type = "h",
    lty.hplot = 2,
    box = FALSE,
    cex.lab = 1.5,
    cex.axis = 1.2,
    cex.main = 1.5
  )
  
  if (add_legend) {
    legend("topleft", legend = names(color_map),
           col = color_map, pch = 19, pt.cex = 1.5,
           cex = 1, bty = "n")
  }
}

# ============================
# Save plots
# ============================
png("FibroblastsCellChat3DScatterPlot_nolegend.png", width = 1200, height = 1000, res = 150)
plot_cellchat3d(add_legend = FALSE)
dev.off()

png("FibroblastsCellChat3DScatterPlot_withlegend.png", width = 1200, height = 1000, res = 150)
plot_cellchat3d(add_legend = TRUE)
dev.off()

svg("FibroblastsCellChat3DScatterPlot_nolegend.svg", width = 7, height = 7)
plot_cellchat3d(add_legend = FALSE)
dev.off()

svg("FibroblastsCellChat3DScatterPlot_withlegend.svg", width = 7, height = 7)
plot_cellchat3d(add_legend = TRUE)
dev.off()
```


```{r}
library(dplyr)
library(ggplot2)
library(readr)

# Load your data (already saved from your previous step)
all_roles_df <- read_csv("Global3DScatterData.csv")

# Add numeric time
time_map <- c("UW" = 0, "D1PW" = 1, "D2PW" = 2, "D4PW" = 4,
              "D7PW" = 7, "D15PW" = 15, "D30PW" = 30)
all_roles_df$time_numeric <- time_map[all_roles_df$sample]

# Compute total signaling (sender + receiver)
all_roles_df <- all_roles_df %>%
  mutate(total_signaling = x + y)

# Summarize by timepoint and cell group
signaling_summary <- all_roles_df %>%
  group_by(labels, sample, time_numeric) %>%
  summarise(total_signaling = mean(total_signaling, na.rm = TRUE), .groups = "drop")

# (Optional) Add metacluster info if available
# For now, if 'labels' are already descriptive of metaclusters, skip this
# Otherwise, you can merge() or mutate() a mapping table here
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


p_facet <- ggplot(signaling_summary, aes(x = time_numeric, y = total_signaling)) +
  geom_area(aes(fill = labels), alpha = 0.6, color = "black", linewidth = 0.6) +
  geom_smooth(
    aes(color = labels, y = total_signaling),
    method = "loess",
    se = FALSE,
    linewidth = 1
  ) +
  facet_wrap(~ labels, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = metacluster_colors, na.value = "grey80") +
  scale_color_manual(values = metacluster_colors, na.value = "grey80") +
  scale_x_continuous(breaks = c(0, 1, 2, 4, 7, 15, 30)) +
  theme_ggprism_mod(base_size = 14) +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0),
    panel.spacing = unit(0.2, "lines"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 7),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    legend.position = "none"  # hides redundant legend since faceted
  ) +
  labs(
    x = "Time (Days Post Wounding)",
    y = "Total Signaling Strength",
    fill = "Cell Type",
    color = "Cell Type"
  )

# ===============================================
# Save versions (PNG + SVG)
# ===============================================
ggsave("facet_signaling_metacluster.png", p_facet, width = 4, height = 9, dpi = 300, bg = "white")
ggsave("facet_signaling_metacluster.svg", p_facet, width = 4, height = 9, bg = "white")




p_overlay <- ggplot(signaling_summary, aes(x = time_numeric, y = total_signaling,
                                           fill = labels, color = labels)) +
  geom_area(alpha = 0.25, position = "identity") +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  scale_color_manual(values = metacluster_colors, na.value = "grey80") +
  scale_fill_manual(values = metacluster_colors, na.value = "grey80") +
  scale_x_continuous(breaks = c(0, 1, 2, 4, 7, 15, 30)) +
  theme_ggprism_mod(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    axis.title = element_text(face = "bold")
  ) +
  labs(
    x = "Time (Days Post Wounding)",
    y = "Total Signaling Strength",
    color = "Cell Type",
    fill = "Cell Type"
  )

p_overlay

ggsave("facetoverlay_signaling_metacluster.png", p_overlay, width = 10, height = 5, dpi = 300, bg = "white")
ggsave("facetoverlay_signaling_metacluster.svg", p_overlay, width = 10, height = 5, bg = "white")
```

