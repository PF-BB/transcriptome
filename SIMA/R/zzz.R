#' Internal utils functions.
#' @title Internal utils functions


eSet2lumiBatch <- function(eSet) {
  lumiBatch = new('LumiBatch', 
                  exprs = eSet@assayData$exprs, 
                  phenoData = new("AnnotatedDataFrame", eSet@phenoData@data), 
                  featureData = new("AnnotatedDataFrame", eSet@featureData@data),
                  annotation = eSet@annotation
  )
  return(lumiBatch)
}
eSet2affyBatch <- function(eSet){
  affyBatch = new('AffyBatch', 
                  exprs = eSet@assayData$exprs, 
                  phenoData = new("AnnotatedDataFrame", eSet@phenoData@data), 
                  featureData = new("AnnotatedDataFrame", eSet@featureData@data),
                  annotation = eSet@annotation
}

batch2eSet <- function(batch){
  eSet = new("ExpressionSet", 
             exprs = batch@assayData$exprs, 
             phenoData = new("AnnotatedDataFrame", batch@phenoData@data), 
             featureData = new("AnnotatedDataFrame", batch@featureData@data),
             annotation = batch@annotation
  )
  return(eSet)
}
