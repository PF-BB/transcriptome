normalisation <-
function(noms_fichiers){
  GCRMA <- justGCRMA(filenames=noms_fichiers, optimize.by="memory")
  logGCRMA <- exprs(GCRMA)
  foo <- function(chaine) substr(chaine,1,nchar(chaine)-4)
  colnames(logGCRMA) <- sapply(colnames(logGCRMA),foo)
  logGCRMA
}
