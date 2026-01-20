
#This code takes made CellChat objects in a rds folder and loops through them to produce plots

library(CellChat)
library(ComplexHeatmap)
library(grid)
library(dplyr)

# =========================
# Step 1: List CellChat RDS files
# =========================
rds_folder <- "/Users/jchincheong/SRSP Laboratory Dropbox/SRSP Lab/Resources/Town Square - Data Share/Manuscripts/Wound Healing Atlas Paper/Hair Follicle Paper/Figure 3/Cell Chat Data/RDS Files/To Use/UPDATE"
rds_files <- list.files(rds_folder, pattern = "cellchat_.*\\.rds$", full.names = TRUE)
celltypes <- gsub("^cellchat_|\\.rds$", "", basename(rds_files))
timepoints <- c("UW","D1PW", "D2PW","D4PW","D7PW","D15PW","D30PW")

# =========================
# Step 2: Define metacluster colors
# =========================
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

# =========================
# Step 3: Output directories
# =========================
outdir_chord   <- "CellChat_chord_plots"
outdir_heat    <- "CellChat_heatmaps"
outdir_global_heat <- "CellChat_global_heatmaps"
outdir_circle  <- "CellChat_circle_plots"
outdir_csv     <- "CellChat_topN_csv"

dir.create(outdir_chord, showWarnings = FALSE)
dir.create(outdir_heat, showWarnings = FALSE)
dir.create(outdir_global_heat, showWarnings = FALSE)
dir.create(outdir_circle, showWarnings = FALSE)
dir.create(outdir_csv, showWarnings = FALSE)

# =========================
# Step 4: Loop through each RDS file
# =========================
for(i in seq_along(rds_files)){
  
  # Read CellChat object
  cellchat <- readRDS(rds_files[i])
  ct <- celltypes[i]  # cell type name
  cat("Processing:", ct, "\n")
  
  # Extract timepoint objects
  object.list <- list(
    UW    = cellchat$UW,
    D1PW = cellchat$D1PW,
    D2PW = cellchat$D2PW,
    D4PW  = cellchat$D4PW,
    D7PW  = cellchat$D7PW,
    D15PW = cellchat$D15PW,
    D30PW = cellchat$D30PW
  )
  
  # Merge into one object for differential plots
  cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))
  
  # =========================
  # 1. Global signaling role heatmaps (per timepoint)
  # =========================
  for(sample_name in names(object.list)){
    obj <- object.list[[sample_name]]
    
    h_out <- netAnalysis_signalingRole_heatmap(obj, pattern = "outgoing", height = 20, width = 10)
    h_in  <- netAnalysis_signalingRole_heatmap(obj, pattern = "incoming", height = 20, width = 10)
    
    outfile_png <- file.path(outdir_global_heat, paste0(sample_name, "_", ct, "_SignalingRoleHeatmaps.png"))
    
    png(outfile_png, width = 4000, height = 3000, res = 300)
    draw(h_out + h_in, annotation_legend_side = "right")
    dev.off()
  }
  
  # =========================
  # 2. Chord plots per timepoint with CSV
  # =========================
  topN <- 70
  for(j in seq_along(object.list)){
    obj <- object.list[[j]]
    tp <- timepoints[j]
    
    # Get top interactions
    net.full <- subsetCommunication(obj)
    net.top <- net.full[order(net.full$prob, decreasing = FALSE), ][1:min(topN, nrow(net.full)), ]
    if(nrow(net.top) == 0) next
    
    # Save CSV
    csv_file <- file.path(outdir_csv, paste0(ct, "_", tp, "_Top", topN, "_interactions.csv"))
    write.csv(net.top, csv_file, row.names = FALSE)
    
    # Match colors to identities present in this object
    idents_i <- levels(obj@idents)
    color_i <- metacluster_colors[idents_i]
    
    # Filenames
    outfile_png <- file.path(outdir_chord, paste0("Chord_", ct, "_Top", topN, "_", tp, ".png"))
    outfile_svg <- file.path(outdir_chord, paste0("Chord_", ct, "_Top", topN, "_", tp, ".svg"))
    
    # PNG
    png(outfile_png, width = 3000, height = 3000, res = 300)
    netVisual_chord_gene(obj, sources.use = NULL, targets.use = NULL,
                         slot.name = "net", net = net.top,
                         color.use = color_i,
                         lab.cex = 1.1, small.gap = 1.7,
                         annotationTrackHeight = c(0.03,0.03),
                         title.name = paste0("Top ", topN, " interactions - ", tp))
    dev.off()
    
    # SVG
    svg(outfile_svg, width = 10, height = 10)
    netVisual_chord_gene(obj, sources.use = NULL, targets.use = NULL,
                         slot.name = "net", net = net.top,
                         color.use = color_i,
                         lab.cex = 1.1, small.gap = 1.7,
                         annotationTrackHeight = c(0.03,0.03),
                         title.name = paste0("Top ", topN, " interactions - ", tp))
    dev.off()
  }
  
  # =========================
  # 3. Circle plots (global per object)
  # =========================
  # PNG
  
  comparisons <- list(
    c(1,4),  # UW vs D4PW
    c(4,5),  # D4PW vs D7PW
    c(5,6),  # D7PW vs D15PW
    c(6,7)   # D15PW vs D30PW
  )
  
  for (comp in comparisons) {
    comp_name <- paste(group.new[comp], collapse = "_vs_")
    
    # PNG
    outfile_png <- file.path(outdir_circle, paste0("CirclePlot_", ct, "_", comp_name, ".png"))
    png(outfile_png, width = 4000, height = 2000, res = 300)
    par(mfrow = c(1,2), xpd = TRUE)
    netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, comparison = comp, color.use = metacluster_colors)
    netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight", comparison = comp, color.use = metacluster_colors)
    dev.off()
    
    # SVG
    outfile_svg <- file.path(outdir_circle, paste0("CirclePlot_", ct, "_", comp_name, ".svg"))
    svg(outfile_svg, width = 14, height = 7)
    par(mfrow = c(1,2), xpd = TRUE)
    netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, comparison = comp, color.use = metacluster_colors)
    netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight", comparison = comp, color.use = metacluster_colors)
    dev.off()
  }
  
  # =========================
  # 4. Differential heatmaps across consecutive comparisons
  # =========================
  comparisons <- list(c(1,2), c(2,3), c(3,4), c(4,5))
  
  for(comp in comparisons){
    ht1 <- netVisual_heatmap(cellchat_merged, comparison = comp)
    ht2 <- netVisual_heatmap(cellchat_merged, measure = "weight", comparison = comp)
    
    fname_base1 <- file.path(outdir_heat, paste0("heatmap_default_", ct, "_", comp[1], "_vs_", comp[2]))
    fname_base2 <- file.path(outdir_heat, paste0("heatmap_weight_", ct, "_", comp[1], "_vs_", comp[2]))
    
    # PNG
    png(paste0(fname_base1, ".png"), width = 2800, height = 2000, res = 300)
    draw(ht1)
    dev.off()
    
    png(paste0(fname_base2, ".png"), width = 2800, height = 2000, res = 300)
    draw(ht2)
    dev.off()
    
    # SVG
    svg(paste0(fname_base1, ".svg"), width = 14, height = 11)
    draw(ht1)
    dev.off()
    
    svg(paste0(fname_base2, ".svg"), width = 14, height = 11)
    draw(ht2)
    dev.off()
  }
  
  cat("Finished processing", ct, "\n\n")
}


```
