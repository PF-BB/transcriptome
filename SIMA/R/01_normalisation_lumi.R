#' Normalization of Illumina datasets.
#' @param dataFile A character string for BeadStudio file.
#' @param pheno An oject of class AnnotatedDataFrame.
#' @param bg.method A character string: the method used to deals with microarray background noise. If the specified method is not supported, an error is generated
#' @param norm.method A character string: the method used perform an inter-array normalization. If the specified method is not supported, an error is generated
#' @return \item{eset}{An object of classe Expression-set.}
#' @references lumi.
#' @title Normalization of Illumina
#' @export lumi_normalization

lumi_normalization <- function(dataFile, pheno, bg.method="none", norm.method="quantile"){
  
  # 0. Load packages (move somewhere else?)
  require(lumi)
  
  # 1. Test on arguments 
  if (missing (dataFile))
    stop("\n\tL'argument dataFile est manquant.",call.=FALSE)
  
  if(!file.exists(dataFile))
    stop(paste0("\n\tLe fichier ",dataFile," n'existe pas."),call.=FALSE)
  
  
  if(all(bg.method != c('none', 'bgAdjust', 'forcePositive', 'bgAdjust.affy')))
    stop("\n\tLa méthode de correction du bruit de fond ",bg.method," n'est pas supportée !",call.=FALSE)
  
  if(all(norm.method != c("quantile", "rsn", "ssn", "loess", "vsn", "rankinvariant","none")))
    stop("\n\tLa méthode de normalisation ",norm.method," n'est pas supportée !",call.=FALSE)
      
  # 2. Create lumiBatch
  data=lumiR.batch(fileList=dataFile,sampleInfoFile=pheno,annotationColumn=c("PROBE_ID","SYMBOL","ACCESSION"),lib.mapping="lumiHumanIDMapping")
  cat("\n")
  
  # 3. Normalization
  ## If already performed with BeadStudio
  if(bg.method=="none" & norm.method=="none"){
    cat("Perform thresholding ...\n")
    cat("Perform log2 transformation ...\ndone.\n")
    data = thresholdIntensity(data, seuil=10)
    data.log2 = log2(exprs(data))
    eSet = new('ExpressionSet',
                  exprs = data.log2, 
                  phenoData = new("AnnotatedDataFrame", data@phenoData@data),
                  featureData = new("AnnotatedDataFrame", data@featureData@data),
                  annotation = data@annotation
    )
  }else{
    ## Else: apply normalization with wrapper function "lumi::lumiExpresso"
    data.norm = lumiExpresso(data,bgcorrect.param=list(method=bg.method),normalize.param=list(method=norm.method),varianceStabilize.param=list(method="log2"))
    eSet = batch2eSet(data.norm)
  }
  
  if(length(grep("SYMBOL",colnames(eSet@featureData@data),ignore.case=TRUE))==0){
    eSet@featureData@data=cbind(eSet@featureData@data,nuID2RefSeqID(rownames(eSet), lib.mapping='lumiHumanIDMapping',returnAllInfo=TRUE))
  }
  
  return(eSet)
}

