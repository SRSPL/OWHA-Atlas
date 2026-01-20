npcs <- function(
    seu,
    var.total=0.95,
    reduction="pca"
){
  if(is.null(seu@reductions[[reduction]])){
    cat("Reduction", reduction, "not found!")
    return(NULL)
  }
  
  tmp.var <- (seu@reductions[[reduction]]@stdev)^2
  var.cut <- var.total*sum(tmp.var)
  n.pcs=0
  var.sum = 0
  while(var.sum < var.cut){
    n.pcs = n.pcs + 1
    var.sum <- var.sum + tmp.var[n.pcs]
  }
  
  return(n.pcs)
}

