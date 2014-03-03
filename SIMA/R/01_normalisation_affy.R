#' Normalization of Affymetrix datasets.
#' @param inputFolder A character string: where to find the CEL files and the phenotype data.
#' @param pheno An oject of class AnnotatedDataFrame.
#' @param bg.method A character string: the method used to deals with microarray background noise. If the specified method is not supported, an error is generated
#' @param norm.method A character string: the method used perform an inter-array normalization. If the specified method is not supported, an error is generated
#' @return \item{eset}{An object of class Expression-set.}
#' @references FARMS.
#' @title Normalization of Affymetrix
#' @export affy_normalization

affy_normalization <- function(inputFolder, pheno, bg.method="none", norm.method="quantile"){
  
  # 0. Load packages (move somewhere else?)
  require(farms)
  require(affy)
  require(annotate)
  
  # 1. Test on arguments 
  if (missing (inputFolder)) 
    stop("\n\tL'argument inputFolder est manquant.\n",call.=FALSE)

  if (! file.exists(inputFolder)) 
    stop( "\n\tLe dossier ",inputFolder," n'existe pas !\n" )  

  if(!( bg.method %in% bgcorrect.methods() ))
    stop("\n\tLa methode de correction du bruit de fond ",bg.method," n'est pas supportee !\n",call.=FALSE)
  
  # 2. Create affyBatch
  affybatch <- ReadAffy(celfile.path = inputFolder, phenoData = pheno)
  
  if(all(norm.method != normalize.methods(affybatch) ))
    stop("\n\tLa methode de normalisation ",norm.method," n'est pas supportee !\n", call.=FALSE)
  
  # 3. Apply normalization with wrapper function "farms::expFarms"
  eset <- expFarms(affybatch, bgcorrect.method = bg.method, pmcorrect.method = "pmonly", 
                   normalize.method = norm.method, weight=0.5)
  annot_lib <- paste0(annotation(eset),".db")
  library(annot_lib, character.only = TRUE)
  
  eset@featureData = new("AnnotatedDataFrame", data.frame(ID = getSYMBOL(rownames(eset), annotation(eset)) ) )
  
  gc()
  return(eset)
}
