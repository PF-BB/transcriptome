#' Normalization of Affymetrix or Illumina datasets.
#' @param input A character string: either the name of the folder where to find the CEL files and the phenotype data OR the name of the raw Illumina file obtained with BeadStudio.
#' @param phenoFile A character string: where to find the phenotypes of interest. It is required to be a tab-separated file.
#' @param type A character string: must be either "affy" or "illumina".
#' @param bg.method A character string: the method used to deals with microarray background noise. If the specified method is not supported, an error is generated
#' @param norm.method A character string: the method used perform an inter-array normalization. If the specified method is not supported, an error is generated
#' @param bwrite A boolean to indicate if the data should be written
#' @param bplot A boolean to indicate if plots whould be saved. In the affirmative, the plots are saved as BMP images in the QC sub-folder.
#' @param output A character string, the path of the working directory.
#' @return \item{eset}{An object of classe Expression-set.}
#' @references FARMS and lumi.
#' @title Normalization of Affymetrix or Illumina datasets.
#' @export normalization

normalization <- function(input, phenoFile, type, bg.method="none", norm.method="quantile", bplot=TRUE, bwrite=TRUE, output="./"){
  
  # 1. Tests
  
  if (missing(input)) stop("\n\tL'argument input est manquant !")
  if (missing(phenoFile)) stop("\n\tL'argument phenoFile est manquant !")
  if (missing(type)) stop("\n\tL'argument type est manquant !")
  
  if (!( type %in% c("affy","lumi"))) 
    stop("Le type de puces ", type, " n'est pas reconnu.")
  
  # 2. Read phenotype
  
  if (missing (phenoFile))
    stop("\n\tL'argument phenoFile est manquant.",call.=FALSE)
  
  if(!file.exists(phenoFile))
    stop("\n\tLe fichier ",phenoFile," n'existe pas.",call.=FALSE)
  
  pheno <- read.AnnotatedDataFrame(phenoFile)
  
  # 3. Call the proper normalization method
  if (type == "affy") {
    eset <- affy_normalization(inputFolder=input,  pheno=pheno,  bg.method=bg.method, norm.method=norm.method)
  } else if (type == "lumi") {
    eset <- lumi_normalization(dataFile=input, pheno=phenoFile, bg.method=bg.method, norm.method=norm.method)
  } else {
    stop("\n\tProbleme avec l'argument 'type'")
  }
  
  # 4. Write outputs if required
  if (bplot)  {
    if(!file.exists("QC/")){
      dir.create("QC")
      cat("Creation du repertoire QC/\n")
    }
    bmp(file.path(output,"QC/norm_plot.bmp"))
    plot(eset, what="boxplot", main="Apres normalisation")
    dev.off()
  }

  if (bwrite) writeExprs(eset, output=paste0(output,"/Data"), type=type, bg.method=bg.method, norm.method=norm.method)
  
  return(eset)
}




