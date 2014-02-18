acplight <-
function(titre="",logValues,classes,couleurs,formes=NULL) {
  if (is.null(formes)) {
    formes <- 16 ; ff <- 16
  }else{
    ff  <- as.numeric(levels(as.factor(formes)))
  }
  PRC       <- prcomp(t(logValues),scale=TRUE)
  PC_prc    <- PRC$x
  cc  <- palette()[as.numeric(levels(as.factor(couleurs)))]
  prcplot <-xyplot(PC_prc[,2]~PC_prc[,1], groups=couleurs,
                    col=cc, xlab="PC1",ylab="PC2",main = titre,
                    pch=ff, key=list(space="right",points=list(tab=classes, pch=ff),
                    text=list(tab=classes, pch=ff),col=cc))
  plot(prcplot)
}
