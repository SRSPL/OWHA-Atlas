#For comparisons, the input should be formatted as "Condition 1 vs. Condition 2". 
#For example, if we want to compare unwounded and D1PW, comparisons should be 
#formatted as "Unwounded vs. Wounded_DPW1"
runGO <- function(srat, group, comparisons, sample, enrich_function = c("enrichGO", "gseGO"), make_plot = TRUE, plot_type = c("dotplot", "ridgeplot"), merge_results = FALSE) {
 #Load necessary packages
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  
  # Set up needed variables
  Idents(srat) <- group
  celltypes <- levels(srat)
  
  # Initialize a list to store results for each cluster
  result_list <- list()
  
  # Iterate through the individual cell types or groups to be compared
  for (celltype in celltypes) {
    # Initialize a list to store results for each comparison
    comparison_results <- list()
    
    # Iterate through comparisons to be made, running differential expression for each comparison
    for (comparison in comparisons) {
      # Extract conditions for comparison
      conditions <- strsplit(comparison, " vs. ")[[1]]
      
      # Perform DE analysis for the current comparison
      df <- FindMarkers(srat, 
                        ident.1 = conditions[2], 
                        ident.2 = conditions[1], 
                        group.by = sample,
                        subset.ident = celltype,
                        method = 'deseq')
      
      # Convert results to a dataframe
      genelist <- data.frame(gene = rownames(df),
                             FC = df$avg_log2FC,
                             group = rep(comparison, nrow(df)))
      
      # Set DE status
      genelist$DE <- ifelse(genelist$FC > 0, "upregulated", "downregulated")
      
      # Append results to the comparison_results list
      comparison_results[[comparison]] <- genelist
    }
    
    # Combine results for all comparisons into a single dataframe if merge_results is TRUE
    if (merge_results) {
      genelist_combined <- do.call(rbind, comparison_results)
      
      # Perform enrichment analysis
      if (enrich_function == "enrichGO") {
        gse <- compareCluster(gene ~ group + DE, 
                              universe = genelist_combined$gene, 
                              data = genelist_combined, 
                              fun = enrichGO, 
                              OrgDb = org.Mm.eg.db,
                              keyType = "SYMBOL", 
                              pvalueCutoff = 0.05, 
                              pAdjustMethod = "BH", 
                              ont = "BP")
      } else if (enrich_function == "gseGO") {
        gse <- compareCluster(gene | FC ~ group + DE, 
                              data = genelist_combined, 
                              fun = gseGO, 
                              OrgDb = org.Mm.eg.db,
                              keyType = "SYMBOL", 
                              pvalueCutoff = 0.05, 
                              pAdjustMethod = "BH", 
                              ont = "BP", 
                              nPermSimple = 10000)
      } else {
        stop('No enrichment function provided')
      }
      
      # Set group levels for plotting
      gse@compareClusterResult$group <- factor(gse@compareClusterResult$group, levels = comparisons)
      
      # Append results to the result_list
      result_list[[celltype]] <- gse
      
      # Plot if required
      if (make_plot) {
        if (plot_type == "dotplot") {
          tiff(paste0(enrich_function, '_dotplot_', celltype, '.tiff'), width = 4000, height = 3000, res = 300)
          print(dotplot(gse, x = "DE") + facet_grid(~group)) + ggtitle(paste0(celltype, ' GO'))
          dev.off()       
        } else {
          stop("Plot type not specified/included")
        } 
      } else {
        print("No plots made")
      }
    } else {
      # Perform enrichment analysis separately for each comparison
      for (comparison in comparisons) {
        genelist <- comparison_results[[comparison]]
        genelist <- genelist %>%
          arrange(desc(FC))
        
        # Perform enrichment analysis
        if (enrich_function == "enrichGO") {
          gse <- compareCluster(gene ~ DE, 
                                universe = genelist$gene, 
                                data = genelist, 
                                fun = enrichGO, 
                                OrgDb = org.Mm.eg.db,
                                keyType = "SYMBOL", 
                                pvalueCutoff = 0.05, 
                                pAdjustMethod = "BH", 
                                ont = "BP")
        } else if (enrich_function == "gseGO") {
          genes <- genelist$FC
          names(genes) <- genelist$gene
          gse <- gseGO(geneList = genes,
                       OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL", 
                       pvalueCutoff = 0.05, 
                       pAdjustMethod = "BH", 
                       ont = "BP", 
                       nPermSimple = 10000)
        } else {
          stop('No enrichment function provided')
        }
        
        # Append results to the result_list
        result_list[[paste0(celltype, "_", comparison)]] <- gse
        
        # Plot if required
        if (make_plot) {
          if (plot_type == "dotplot") {
            tiff(paste0(enrich_function, '_dotplot_', celltype, '_', comparison, '.tiff'), width = 4000, height = 3000, res = 300)
            print(dotplot(gse, x = "NES") + ggtitle(paste0(celltype, ' ', comparison, ' GO')))
            dev.off()       
          } else if (plot_type == "ridgeplot") {
            tiff(paste0(enrich_function, '_ridgeplot_', celltype, '_', comparison, '.tiff'), width = 4000, height = 3000, res = 300)
            print(ridgeplot(gse) + ggtitle(paste0(celltype, ' ', comparison, ' GSEA')))
            dev.off()    
          } else {
            stop("Plot type not specified/included")
          }
        } else {
          print("No plots made")
        }
      }
    }
  }
  
  return(result_list)
}
