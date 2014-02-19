#' Normalization of Illumina datasets.
#' @param dataFile A character string for BeadStudio file.
#' @return \item{eset}{An object of classe Expression-set.}
#' @references lumi.
#' @title Normalization of Illumina
#' @export lumi_normalization

lumi_normalization <- function(dataFile,phenoFile,bg.method="none",norm.method="quantile"){
  
  require(lumi)
  
  #Tests
  if (missing (dataFile)){
    stop("\n\tL'argument dataFile est manquant.\n",call.=FALSE)
  }
  
  if(!file.exists(dataFile)){
    stop(paste0("\n\tLe fichier ",dataFile," n'existe pas.\n"),call.=FALSE)
  }
  
  if (missing (phenoFile)){
    stop("\n\tL'argument phenoFile est manquant.\n",call.=FALSE)
  }
  
  if(!file.exists(phenoFile)){
    stop(paste0("\n\tLe fichier ",phenoFile," n'existe pas.\n"),call.=FALSE)
  }
  
  if(all(bg.method != c('none', 'bgAdjust', 'forcePositive', 'bgAdjust.affy'))){
    stop(paste0("\n\tLa méthode de correction du bruit de fond ",bg.method," n'est pas supportée !\n"),call.=FALSE)
  }
  
  if(all(norm.method != c("quantile", "rsn", "ssn", "loess", "vsn", "rankinvariant","average"))){
    stop(paste0("\n\tLa méthode de normalisation ",norm.method," n'est pas supportée !\n"),call.=FALSE)
  }
  
  #Lecture du fichier phénotype
  pheno = read.table(phenoFile,sep="\t",header=T)
  
  #lecture du fichier avec la fonction lumiR.batch
  data = tryCatch(
    lumiR.batch(fileList=dataFile,sampleInfoFile=pheno),
    warning=function(w){
      stop(paste0("\n\tLes IDs du fichier ",phenoFile," ne correspondent pas à ceux de ",dataFile,".\nVérifier votre fichier de phenotype ! \n"),call.=FALSE)
    },
    error=function(e) {
      stop(paste0("\n\tLe fichier ",dataFile," n'est pas un fichier BeadStudio.\n Vérifiez votre fichier !\n"),call.=FALSE)
    }
  )  
  
  #Normalisation avec la fonction lumiExpresso
  data.norm = lumiExpresso(data,bgcorrect.param=list(method=bg.method),normalize.param=list(method=norm.method),varianceStabilize.param=list(method="log2"))
  
  #Normalisation "average" à implémenter
  #if(norm.method=="average"){
  #  
  #}
  
  eSet = lumiBatch2eSet(data.norm)
  
  return(eSet)
  
}
