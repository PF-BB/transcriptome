#' Normalization of Illumina datasets.
#' @param dataFile A character string for BeadStudio file.
#' @return \item{eset}{An object of classe Expression-set.}
#' @references lumi.
#' @title Normalization of Illumina
#' @export lumi_normalization

lumi_normalization <- function(dataFile, phenoFile, bg.method="none", norm.method="quantile"){
  
  # 0. Load packages (move somewhere else?)
  require(lumi)
  
  # 1. Test on arguments 
  if (missing (dataFile)){
    stop("\n\t L'argument dataFile est manquant.",call.=FALSE)
  }
  
  if(!file.exists(dataFile)){
    stop("\n\t Le fichier ",dataFile," n'existe pas.",call.=FALSE)
  }
  
  if (missing (phenoFile)){
    stop("\n\t L'argument phenoFile est manquant.",call.=FALSE)
  }
  
  if(!file.exists(phenoFile)){
    stop("\n\t Le fichier ",phenoFile," n'existe pas.",call.=FALSE)
  }
  
  if(all(bg.method != c('none', 'bgAdjust', 'forcePositive', 'bgAdjust.affy'))){
    stop("\n\t La méthode de correction du bruit de fond ",bg.method," n'est pas supportée !",call.=FALSE)
  }
  
  if(all(norm.method != c("quantile", "rsn", "ssn", "loess", "vsn", "rankinvariant","average"))){
    stop("\n\t La méthode de normalisation ",norm.method," n'est pas supportée !",call.=FALSE)
  }
  
  #2. Read phenotype data file
  pheno = read.table(phenoFile, sep="\t", header=T)
  
  # 3. Create lumiBatch
  data = tryCatch(
    lumiR.batch(fileList=dataFile,sampleInfoFile=pheno),
    warning=function(w){
      stop("\n\t Les IDs du fichier ",phenoFile," ne correspondent pas à ceux de ",dataFile,".\nVérifier votre fichier de phenotype ! ",call.=FALSE)
    },
    error=function(e) {
      stop("\n\t Le fichier ",dataFile," n'est pas un fichier BeadStudio.\n Vérifiez votre fichier !",call.=FALSE)
    }
  )  
  
  # 4. Apply normalization with wrapper function "lumi::lumiExpresso"
  data.norm = lumiExpresso(data,bgcorrect.param=list(method=bg.method),normalize.param=list(method=norm.method),varianceStabilize.param=list(method="log2"))
  
  #Normalisation "average" à implémenter
  #if(norm.method=="average"){
  #  
  #}
  
  eSet = lumiBatch2eSet(data.norm)
  
  return(eSet)
}
