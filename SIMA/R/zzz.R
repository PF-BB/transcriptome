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
  )
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

thresholdIntensity <- function(eSet, seuil=NULL){
  
  if(is.null(seuil))
    stop("L'argulent threshold est manquant avec aucune valeur par défaut",call.=FALSE)
  
  if ( !is.numeric(seuil))
    stop("L'argument threshold doit être une valeur numérique",call.=FALSE)
  
  data = exprs(eSet)
  if(max(data) < 20 & seuil >= 10)
    warning("Les valeurs d'intensités sont en log2. Vérifiez votre seuil ! ")
  
  data[data < seuil] = seuil
  
  eSetnew = new('ExpressionSet',
                exprs = data, 
                phenoData = new("AnnotatedDataFrame", eSet@phenoData@data),
                featureData = new("AnnotatedDataFrame", eSet@featureData@data),
                annotation = eSet@annotation
  )
  
  return(eSetnew)
  
}
