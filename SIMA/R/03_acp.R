#' Principal Component Analysis.
#' @param fichier .
#' @param eset .
#' @param bWrite Boolean, whether to save figure in a file. Default to FALSE.
#' @title PCA.S
#' @export acp

acp <-function(fichier, eset, bWrite=FALSE) {
  require(ggplot2)
  require(reshape2)
  logValues <- exprs(eset)
  Classes <- unlist(pData(eset))
  couleurs <- as.numeric(as.factor(Classes))
  PRC       <- prcomp(t(logValues),scale=TRUE)
  PC_prc    <- PRC$x
  noms      <- colnames(logValues)
  
  pca <- prcomp(t(logValues), scale=T)
  scores <- data.frame(Classes, pca$x[,1:3])
  pc1.2 <- qplot(x=PC1, y=PC2, data=scores, colour=Classes, size=2) + scale_size(guide = 'none') + 
    labs(x="1st principal component", y="2nd principal component") + 
    geom_text(aes(label=colnames(eset)),vjust=-1,show_guide  = F) +
    xlim(range(scores$PC1)*1.8) # used for unusually long sample names --> pass as parameter?
  
  if (bWrite) {
    pdf(fichier)
    pc1.2
    dev.off()
  }
#   logValues <- exprs(eset)
#   classes <- unlist(pData(eset))
#   couleurs <- as.numeric(as.factor(classes))
#   PRC       <- prcomp(t(logValues),scale=TRUE)
#   PC_prc    <- PRC$x
#   noms      <- colnames(logValues)
#   cc  <- palette()[as.numeric(levels(as.factor(couleurs)))]
#   cc2 <- palette()[couleurs]
#   pourcent <- xyplot( 100*(summary(PRC)$importance[3,])~1:length(PRC$sdev),main="Cumulative proportion of Variance",xlab="PC number",ylab="Percentage", panel=function(x,y){llines(x,y);lpoints(x,y,pch=5)})
#   prcplot <-xyplot(PC_prc[,2]~PC_prc[,1],groups=couleurs,col=cc,xlab="PC1",ylab="PC2",main = "2D-PCA", pch=16, key=list(space="right",points=list(tab=classes, pch=16),text=list(tab=classes, pch=16),col=cc))
#   prctextplot<-xyplot(PC_prc[,2]~PC_prc[,1],groups = noms,panel = function(x, y, subscripts, groups)ltext(x = x, y = y, label = groups[subscripts], cex=1,fontfamily = "HersheySans",col=cc2))
#   pdf(fichier)
#   plot(prcplot)
#   biplot(PRC)
#   plot(update(prctextplot,aspect = "xy",xlab="PC1",ylab="PC2",main = "2D-PCA"))
#   plot(pourcent)
#   dev.off()
  pc1.2
}
