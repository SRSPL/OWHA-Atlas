RunDoubletFinder <- function(seur.list) {
  # Determine the estimated number of doublets based on the number of cells in the individual objects
  doublet_rates <- estimateDoubletRate.DWM(seur.list)
  
  # Iterate over each Seurat object in the list
  seur.list <- lapply(names(seur.list), function(seurat_name) {
    seurat_obj <- seur.list[[seurat_name]]
    
    sweep.res <- paramSweep(seurat_obj, PCs = 1:20, sct = TRUE)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    # Select the pK that corresponds to max bcmvn to optimize doublet detection
    pK <- bcmvn %>%
      filter(BCmetric == max(BCmetric)) %>%
      select(pK) %>%
      as.numeric(as.character(.))
    
    annotations <- table(seurat_obj[["seurat_clusters"]])
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi <- round(doublet_rates[seurat_name] * nrow(seurat_obj@meta.data) / 100)  # Use doublet rate for the current Seurat object
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    
    seurat_obj <- doubletFinder(seurat_obj, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    
    # Find the pANN column name dynamically
    pANN_colname <- grep("^pANN", names(seurat_obj@meta.data), value = TRUE)
    
    seurat_obj <- doubletFinder(seurat_obj, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = pANN_colname, sct = TRUE)
    
    return(seurat_obj)
  })
  
  return(seur.list)
}
