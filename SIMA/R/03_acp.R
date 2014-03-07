#' Principal Component Analysis.
#' @param fichier .
#' @param eset .
#' @title PCA.
#' @export acp

acp <-function(fichier,eset) {
  logValues <- exprs(eset)
  classes <- pData(eset)
  couleurs <- as.numeric(as.factor(classes))
  PRC       <- prcomp(t(logValues),scale=TRUE)
  PC_prc    <- PRC$x
  noms      <- colnames(logValues)
  cc  <- palette()[as.numeric(levels(as.factor(couleurs)))]
  cc2 <- palette()[couleurs]
  pourcent <- xyplot( 100*(summary(PRC)$importance[3,])~1:length(PRC$sdev),main="Cumulative proportion of Variance",xlab="PC number",ylab="Percentage", panel=function(x,y){llines(x,y);lpoints(x,y,pch=5)})
  prcplot <-xyplot(PC_prc[,2]~PC_prc[,1],groups=couleurs,col=cc,xlab="PC1",ylab="PC2",main = "2D-PCA", pch=16, key=list(space="right",points=list(tab=classes, pch=16),text=list(tab=classes, pch=16),col=cc))
  prctextplot<-xyplot(PC_prc[,2]~PC_prc[,1],groups = noms,panel = function(x, y, subscripts, groups)ltext(x = x, y = y, label = groups[subscripts], cex=1,fontfamily = "HersheySans",col=cc2))
  pdf(fichier,sep="",collapse="")
  plot(prcplot)
  biplot(PRC)
  plot(update(prctextplot,aspect = "xy",xlab="PC1",ylab="PC2",main = "2D-PCA"))
  plot(pourcent)
  dev.off()
}
