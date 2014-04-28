#' Volcano plot performed after a differential analysis. It represents the p-value obtained after a two sample (limma) test as a function of the Fold Change.
#' @param fit an an ‘MArrayLM’ fitted linear model object.
#' @param coef, an Integer value corresponding to the desired contrast. Default to NULL.
#' @param FC, numeric, Fold Change threshold.
#' @param PV, numeric, p-value threshold.
#' @param p.val, character string, compute raw p-value ("raw") or adjusted p-value (adj"). Default to "adj". The p-values are adjusted with the procedure defined by Benjamini and Hochberg.
#' @param name, character string, title for the graphic. Default is an empty title ("").
#' @title Volcano plot
#' @export VolcanoPlot

VolcanoPlot <- function(fit, coef=NULL, FC=2, PV=0.05, p.val = c("adj","raw"), name = ""){  # , plotInFile=TRUE
	
	titre = NULL
	pval  = NULL
	ylab = ""
	p.val <- match.arg(p.val)
  
	if (p.val == "adj"){
	  topT = topTable(fit, coef=coef, number=Inf)
		if (is.element("adj.P.Val", colnames(topT))) {
			pval  = -log10(topT$adj.P.Val)
			if (name == ""){
				titre = paste("Volcanoplot","(FDR)", sep=" - ")
				
			}
			else {
				titre = paste("Volcanoplot",name,"(FDR)", sep=" - ")
			}
			ylab  = "-log10(Adj.P.Val)"
		}
		else{
			warning("The topTable object does not contain column 'adj.P.Val'")
			return()
		}
	}
	else if (p.val == "raw"){
	  topT = topTable(fit, coef=coef, adjust.method="none",number=Inf)
		if (is.element("P.Value", colnames(topT))) {
			pval = -log10(topT$P.Value)
			if (name == ""){
				titre = paste("Volcanoplot","(noFDR)", sep=" - ")
			}
			else {
				titre = paste("Volcanoplot",name,"(noFDR)", sep=" - ")
			}
			ylab  = "-log10(P.Value)"		
		}
		else{
			warning("The topTable object does not contain column 'P.Value'")
			return()
		}
	}
	else{
		warning("The option 'p.val' should be either 'adj' or 'raw'")
		return()
	}	
	
	
	gauche = intersect( which(pval >= -log10(PV))  , which(topT$logFC <= -log2(FC)))
	droite = intersect( which(pval >= -log10(PV))  , which(topT$logFC >= log2(FC)))
	reste  = union( which(pval< -log10(PV)), which(abs(topT$logFC)< log2(FC) ))
	
	# --- definition de l'echelle des abscisses
	min.log = min(topT$logFC,na.rm=T)
	max.log = max(topT$logFC,na.rm=T)
	max.xax = floor(max(abs(min.log), abs(max.log)))+1
	xaxis = c(-max.xax, max.xax)
	
	# --- definition de l'echelle des ordonnÈe
	max.pv = max(pval, na.rm=TRUE)
	max.yax = 10*(max.pv%/%10) + 10
	yaxis = c(0,max.yax)
	
	# --- Plot
	layout( matrix(c(1,1,2,3), ncol=2,nrow=2, byrow=T), heights=c(0.9,0.1) )
	#layout.show(3)
	
	par(mar = c(5,5,5,5))
	xlab = "log2(FC)"
	plot  (topT$logFC[reste] , pval[reste] , col="grey" , pch=16, cex=0.7,main=titre, xlab=xlab,ylab=ylab,xlim=xaxis, ylim=yaxis)
	points(topT$logFC[gauche], pval[gauche], col="chartreuse4", pch=16, cex=1)
	points(topT$logFC[droite], pval[droite], col="red", pch=16, cex=1)
	
	abline(v=c(-log2(FC), 0, log2(FC)), col=c("blue","black","blue"), lty=2)
	abline(h= -log10(PV), col="blue", lty=2)
	
	
	labs1 = c("Threshold:",paste(" - FC:", FC, sep=" "), paste(" - PV:", PV))
	par(mar = c(0.1,0.1,0.1, 1.1))
	plot( c(0,1), c(0,1), type="n", ylab="", xlab = "", axes = F)
	text( x = c(0,0,0), y=c(0.8,0.5,0.1), labels = labs1 , adj = c(0,0),font=2)
	
	labs2 = c(paste("Underexpressed:", length(gauche), sep=" "),
			paste("Overexpressed :" , length(droite), sep=" "))
	par(mar = c(0.1,0.1,0.1, 1.1))
	plot( c(0,1), c(0,1), type="n", ylab="", xlab = "", axes = F)
	text( x = c(0,0), y=c(0.8,0.5), labels = labs2 , adj = c(0,0),col=c("chartreuse4","red"),font=2)
}



