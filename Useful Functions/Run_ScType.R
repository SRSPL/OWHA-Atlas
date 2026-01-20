Run_ScType <- function(seurat_object, path_to_db_file, tissue_type, assay = "SCT", use_negative_markers = TRUE, cluster_column = NULL, reduction = NULL, annotation_column = NULL) {
  #Load necessary packages 
  library(openxlsx)
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  # Define the gene_sets_prepare function
  # A modified gene_sets_prepare function to work with mouse gene names (custom DB)
  gene_sets_prepare <- function(path_to_db_file, tissue_type){
    cell_markers = openxlsx::read.xlsx(path_to_db_file)
    cell_markers = cell_markers[cell_markers$tissueType == tissue_type,] 
    cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
    # correct gene symbols from the given DB (up-genes)
    cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
      markers_all = gsub(" ","",unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
      #markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""]) # Only for human, removed for mouse analysis
      markers_all = sort(markers_all)
      
      if(length(markers_all) > 0){
        # suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all, species = "mouse")$Suggested.Symbol))}) # This is used if the symbols are not in the format of HGNC
        paste0(markers_all, collapse = ",")
      } else {
        ""
      }
    })
    
    # correct gene symbols from the given DB (down-genes)
    cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
      markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
      # markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""]) # Only for human
      markers_all = sort(markers_all)
      
      if (length(markers_all) > 0){
        #suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all, species = "mouse")$Suggested.Symbol))}) # This is used if the symbols are not in the format of HGNC
        paste0(markers_all, collapse = ",")
      } else {
        ""
      }
    })
    
    cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1)
    cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
    cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2)
    cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
    
    gs = lapply(1:nrow(cell_markers), function(j) {
      gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),","))
      )
    })
    names(gs) = cell_markers$cellName
    
    gs2 = lapply(1:nrow(cell_markers), function(j) {
      gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),","))
      )
    })
    names(gs2) = cell_markers$cellName
    
    list(gs_positive = gs, gs_negative = gs2)
  }
  # Prepare gene sets
  gs_list <- gene_sets_prepare(path_to_db_file, tissue_type)
  
  if(use_negative_markers) {
    # Get cell-type by cell matrix
    es.max <- sctype_score(scRNAseqData = seurat_object[[assay]]$scale.data, scaled = TRUE, 
                           gs = gs_list$gs_positive, gs2 = gs_list$gs_negative, gene_names_to_uppercase = FALSE)
  } else {
    es.max <- sctype_score(scRNAseqData = seurat_object[[assay]]$scale.data, scaled = TRUE, 
                           gs = gs_list$gs_positive, gs2 = NULL, gene_names_to_uppercase = FALSE)
  }
  
  # Merge by cluster
  cL_results <- do.call("rbind", lapply(unique(seurat_object@meta.data[[cluster_column]]), function(cl) {
    es.max.cl <- sort(rowSums(es.max[, rownames(seurat_object@meta.data[seurat_object@meta.data[[cluster_column]] == cl, ])]), decreasing = TRUE)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_object@meta.data[[cluster_column]] == cl)), 15)
  }))
  
  sctype_scores <- cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 4] <- "Unknown"
  
  # Annotate Seurat object with predicted cell types
  seurat_object@meta.data[annotation_column] <- ""
  for (j in unique(sctype_scores$cluster)) {
    cl_type <- sctype_scores[sctype_scores$cluster == j, ]
    seurat_object@meta.data[[annotation_column]][seurat_object@meta.data[[cluster_column]] == j] <- as.character(cl_type$type[1])
  }
  
  # Overlay the predicted cell types on the UMAP
  #DimPlot(seurat_object, group.by = "customclassif", reduction = reduction, label = T, raster = F, repel = T)
  
  return(seurat_object)
} 
