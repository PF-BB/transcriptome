#' Normalization of Affymetrix datasets.
#' @param inputFolder A character string: where to find the CEL files and the phenotype data.
#' @return \item{eset}{An object of classe Expression-set.}
#' @references FARMS.
#' @title Normalization of Affymetrix
#' @export affy_normalization

affy_normalization <- function(inputFolder){
  require(farms)
  affybatch <- ReadAffy(celfile.path = inputFolder, phenoData=new("AnnotatedDataFrame")) ; gc()
  eset <- expFarms(affybatch, bgcorrect.method = "none", pmcorrect.method = "pmonly", 
                   normalize.method = "constant", weight=0.5) ; gc()
  return(eset)
}




