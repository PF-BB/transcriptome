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
  if (missing (dataFile))
    stop("\n\tL'argument dataFile est manquant.",call.=FALSE)
  
  if(!file.exists(dataFile))
    stop(paste0("\n\tLe fichier ",dataFile," n'existe pas."),call.=FALSE)
  
  if (missing (phenoFile))
    stop("\n\tL'argument phenoFile est manquant.",call.=FALSE)
    
  if(!file.exists(phenoFile))
    stop("\n\tLe fichier ",phenoFile," n'existe pas.",call.=FALSE)
  
  if(all(bg.method != c('none', 'bgAdjust', 'forcePositive', 'bgAdjust.affy')))
    stop("\n\tLa méthode de correction du bruit de fond ",bg.method," n'est pas supportée !",call.=FALSE)
  
  if(all(norm.method != c("quantile", "rsn", "ssn", "loess", "vsn", "rankinvariant","average")))
    stop("\n\tLa méthode de normalisation ",norm.method," n'est pas supportée !",call.=FALSE)
    
  #2. Read phenotype data file
  pheno = read.table(phenoFile, sep="\t", header=T)
  
  # 3. Create lumiBatch
  data = tryCatch(
    lumiR.batch(fileList=dataFile,sampleInfoFile=pheno),
    warning=function(w){
      stop("\n\tLes IDs du fichier ",phenoFile," ne correspondent pas à ceux de ",dataFile,".\nVérifier votre fichier de phenotype !",call.=FALSE)
    },
    error=function(e) {
      stop("\n\tLe fichier ",dataFile," n'est pas un fichier BeadStudio.\n Vérifiez votre fichier !",call.=FALSE)
    }
  )  
  cat("\n")
  # 4. Normalisation
  if(norm.method=="average"){
    
        eSet = lumi_normalization_average(data)
        
  }else{
    # 4. Apply normalization with wrapper function "lumi::lumiExpresso"
    data.norm = lumiExpresso(data,bgcorrect.param=list(method=bg.method),normalize.param=list(method=norm.method),varianceStabilize.param=list(method="log2"))
    eSet = batch2eSet(data.norm)
    
  }
  
  return(eSet)
}


lumi_normalization_average <- function(rawLumiBatch){
  cat("Background Correction: none\n")
  cat("Variance Stabilizing Transform method: log2\n")
  cat("Normalisation method: average\n\n")
  cat("Choose reference array(s) in the table\n")
  #Modifier la ou les arrays de reference (1=reference, 0=autre array)
  M = cbind(rawLumiBatch@phenoData@data,refArray=c(1,rep(0, ncol(rawLumiBatch)-1)))
  M=fix(M)
  cat(paste0("Number of reference arrays = ",length(which(M$refArray=="1")),"\n\n"))
  cat("Perform average normalization ...\n")
  
  dataraw.mat = exprs(rawLumiBatch)
  dataraw.mat[dataraw.mat < 1] = 1.01
  dataraw.mat = log2(dataraw.mat)
  datanorm.mat= dataraw.mat
    
  
  ref.mean.array = mean(dataraw.mat[,M$refArray==1])
  scaling.factor=list()
  for(i in 1:ncol(dataraw.mat)){
    scaling.factor[[i]] = ref.mean.array/mean(dataraw.mat[,i])
    datanorm.mat[,i]=dataraw.mat[,i]*scaling.factor[[i]]
  }
  
  eSet = new('ExpressionSet',
             exprs = datanorm.mat, 
             phenoData = new("AnnotatedDataFrame", rawLumiBatch@phenoData@data),
             featureData = new("AnnotatedDataFrame", rawLumiBatch@featureData@data),
             annotation = rawLumiBatch@annotation
  )
  cat("done.\n")
  return(eSet)
  
}