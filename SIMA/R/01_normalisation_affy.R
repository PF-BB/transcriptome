#' Normalization of Affymetrix datasets.
#' @param inputFolder A character string: where to find the CEL files and the phenotype data.
#' @param phenoFile A character string: where to find the phenotypes of interest. It is required to be a tab-separated file.
#' @param bg.method A character string: the method used to deals with microarray background noise. If the specified method is not supported, an error is generated
#' @param norm.method A character string: the method used perform an inter-array normalization. If the specified method is not supported, an error is generated
#' @return \item{eset}{An object of class Expression-set.}
#' @references FARMS.
#' @title Normalization of Affymetrix
#' @export affy_normalization

affy_normalization <- function(inputFolder, phenoFile, bg.method="none", norm.method="quantile"){
  
  # 0. Load packages (move somewhere else?)
  require(farms)
  require(affy)
  
  # 1. Test on arguments 
  if (missing (inputFolder)) 
    stop("\n\tL'argument inputFolder est manquant.\n",call.=FALSE)
  if (missing (phenoFile)) 
    stop("\n\tL'argument phenoFile est manquant.\n",call.=FALSE)
  
  if (! file.exists(inputFolder)) 
    stop( "\n\Le dossier ",inputFolder," n'existe pas !\n" )  
  if (! file.exists(phenoFile)) 
    stop( "\n\tLe fichier ",phenoFile," n'existe pas !\n" )  
  
  if(!( bg.method %in% bgcorrect.methods() ))
    stop("\n\tLa methode de correction du bruit de fond ",bg.method," n'est pas supportee !\n",call.=FALSE)
  
  # 2. Read phenotype data file
  pheno <- read.table(phenoFile,sep="\t",header=T)
  
  # 3. Create affyBatch
  affybatch <- ReadAffy(celfile.path = inputFolder, phenoData = pheno) ; gc()
  
  if(all(norm.method != normalize.methods(affybatch) ))
    stop("\n\tLa methode de normalisation ",norm.method," n'est pas supportee !\n", call.=FALSE)
  
  # 4. Apply normalization with wrapper function "farms::expFarms"
  eset <- expFarms(affybatch, bgcorrect.method = bg.method, pmcorrect.method = "pmonly", 
                   normalize.method = norm.method, weight=0.5) ; gc()
  
  return(eset)
}