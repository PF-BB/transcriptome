#' Normalization of Affymetrix or Illumina datasets.
#' @param input A character string: either the name of the folder where to find the CEL files and the phenotype data OR the name of the raw Illumina file obtained with BeadStudio.
#' @param phenoFile A character string: where to find the phenotypes of interest. It is required to be a tab-separated file.
#' @param type A character string: must be either "affy" or "illumina".
#' @return \item{eset}{An object of classe Expression-set.}
#' @references FARMS.
#' @title Normalization of Affymetrix
#' @export affy_normalization

#' @param dataFile A character string for BeadStudio file.
#' @return \item{eset}{An object of classe Expression-set.}
#' @references lumi.
#' @title Normalization of Illumina
#' @export lumi_normalization

#' @param inputFolder A character string: where to find the CEL files and the phenotype data.
#' @param bg.method A character string: the method used to deals with microarray background noise. If the specified method is not supported, an error is generated
#' @param norm.method A character string: the method used perform an inter-array normalization. If the specified method is not supported, an error is generated
#' @return \item{eset}{An object of class Expression-set.}
#' @references FARMS.
#' @title Normalization of Affymetrix
#' @export affy_normalization

normalization <- function(input, phenoFile, type, bg.method="none", norm.method="quantile"){
  
  # 1. tests
  
  if (missing(input)) stop("\n\tL'argument input est manquant !\n")
  if (missing(phenoFile)) stop("\n\tL'argument phenoFile est manquant !\n")
  if (missing(type)) stop("\n\tL'argument type est manquant !\n")
  
  if (!( type %in% c("affy","lumi"))) 
    stop("Le type de puces ", type, " n'est pas reconnu.")
  
  # 2. Call the proper normalization method
  if (type == "affy") {
    eset <- affy_normalization(inputFolder=input,  phenoFile=phenoFile,  bg.method=bg.method, norm.method=norm.method)
  } else if (type == "lumi") {
    eset <- lumi_normalization(dataFile=input, phenoFile=phenoFile, bg.method=bg.method, norm.method=norm.method){
      
  } else {
    stop("\n\t Probleme avec l'argument 'type'\n")
  }
  
  return(eset)
}




