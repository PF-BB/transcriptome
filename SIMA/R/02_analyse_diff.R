#' Declare the contrasts of interests for the differential analysis.
#' @param X An Expression Set.
#' @param y A vector contaning the class whose each sample blongs to.
#' @param contr A vector of character strongs describing the contrasts. 
#' @return \item{top}{An object topTable obtained after using limma's fit function.}
#' @references limma.
#' @title Perform differential analysis.
#' @export analyse_diff

analyse_diff <- function(X, y, contr){
	design <- model.matrix(~ 0 + y)
	colnames(design) <- levels(y)
	cont <- makeContrasts(contrasts=contr, levels=design)
	fit  <- lmFit(X, design)
	fit  <- contrasts.fit(fit, cont)
	fit  <- eBayes(fit)
	return( fit )
}

#' Modify the original topTable function to allow non-log fold changes and average expressions.
#' @param fit .
#' @param FC .
#' @param PV . 
#' @param p.val . 
#' @return \item{top}{An object topTable with a twist obtained after using limma's fit function.}
#' @references limma.
#' @title Top table reducer.
#' @export topTableReducer

topTableReducer <- function(fit, coef=NULL,  FC = 2, PV=0.05, p.val = c("adj","raw"), tlog2=FALSE, eset,bWrite=T){
  p.val <- match.arg(p.val)
  sort.by.pval <- ifelse(ncol(fit) > 1 & is.null(coef), "F", "P")
  top  <- topTable(fit, coef=coef, adjust.method="BH", sort.by=sort.by.pval, number=nrow(fit), 
                   p.value=ifelse(p.val=="adj",PV,1), lfc=log2(FC))
  if (p.val=="raw") {
    top <- subset(top, P.Value < PV)
  }
  foo <- function(A) sign(A) * 2^( abs(A) )
  if (nrow(top) > 0){
    if (tlog2 == FALSE) {
      N1 <- ifelse( is.null(ncol(fit$genes)) ,0,ncol(fit$genes))
      N2 <- ifelse(is.null(coef), ncol(fit), 1)
      top[, N1+(1:N2)] <- foo(top[, N1+(1:N2)])
      names(top)[names(top)=="logFC"] <- "FC"
    }
    top$AveExpr <- 2^(top$AveExpr)
    
    if(bWrite){
      idxPdata = grep(colnames(fit$design)[1],pData(eset))
      class1 = which(pData(eset)[idxPdata]==rownames(fit$contrasts)[which(fit$contrasts==1)])
      class2 = which(pData(eset)[idxPdata]==rownames(fit$contrasts)[which(fit$contrasts==-1)])
      if(file.exists("AnalyseDifferentielle/"))
        output="AnalyseDifferentielle/"
      else
        output="./"
      write.table(cbind(top,2^(exprs(eset[top[,1],c(class1,class2)]))),
                  paste0(output,"/topTable_",colnames(fit$contrasts),"_PV",p.val,PV,"_FC",FC,".xls"),
                  sep="\t", row.names=F,quote=F)
    }
    
    return(top)
  }
  else
    stop("Aucun gene differentiellement exprime selon ces critères ! ")
}
