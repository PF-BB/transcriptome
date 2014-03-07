#' Heatmap with annotations.
#' @param X numeric matrix of the values to be plotted
#' @param dist.method function used to compute the distance (dissimilarity) between both rows and columns
#' @param hclust.method function used to compute the hierarchical clustering
#' @param cols colors of the heatmap
#' @param annot a vector or matrix to annotate the samples
#' @param labRow character vectors with row/column labels to use 
#' @param labCol character vectors with row/column labels to use 
#' @param title character string for the title of the cluster
#' @param scale character indicating if the values should be centered and scaled ("row" or "none")
#' @param Rowv booleen. TRUE plots samples/genes dendrogram. FALSE otherwise
#' @param Colv booleen. TRUE plots samples/genes dendrogram. FALSE otherwise
#' @param hcr object of class hclust if dendrogram already calculated
#' @param hcc object of class hclust if dendrogram already calculated
#' @param saturate.cols boolen.
#' @return \item{list}{An list with 2 objects of class hclust.}
#' @title Draw a Heat Map with Annotations
#' @export sima.heatmap
sima.heatmap = function(	X, 
		dist.method=c("pearson","euclidean","manhattan","L1"), 
		hclust.method=c("ward","complete"),
		cols=c("RedGreen","OrangeBlue","RedBlue"),
		annot,
		labRow=rownames(X),
		labCol=colnames(X),
		title="Heatmap",
		scale="row",
		Rowv=TRUE,
		Colv=TRUE,
		hcr=NULL,
		hcc=NULL,
		saturate.cols=F){
  library(RColorBrewer)
	nbSample = ncol(X)
	
	# --- Calcul de la matrice de distance
	if(dist.method=="pearson"){
		dd = distPearson
	}else if(dist.method=="spearman"){
		dd= distSpearman
	}else{
		dd= function(x){dist(x)}
	}
	#else if(dist.method="L1"){
	#	dd = distL1
	#}
	
	#Colorisation de la heatmap
	if(cols=="RedGreen")
		Zcol = colz(low = "green", high = "red", mid="black", n=100)
	if(cols=="OrangeBlue")
		Zcol = colz(low = "blue", high = "orange", mid="white", n=100)
	if(cols=="RedBlue")
		Zcol = colz(low = "blue", high = "red", mid="white", n=100)
	
	#Annotations
	nbLabel = 0
	if(!missing(annot)){
		annot=as.matrix(annot)
		if(nrow(annot) != nbSample)
			stop("Number of arrays in annotation matrix is different from expression matrix")
		else{
			nbLabel = ncol(annot)
			if(nbLabel > 4) stop("Can't display more than 4 annotations")
		}
	}
	
	
	
	#ID Genes
	if(is.null(labRow))
		labRow=rep("",nrow(X))
	#ID Arrays
	if(is.null(labCol))
		labCol=rep("",ncol(X))
	
	#Calcul dendrogram genes
	if(is.null(hcr))
		hcr = hclust(dd(X), method=hclust.method)
	
	#calcul dendrogram arrays
	if(is.null(hcc))
		hcc = hclust(dd(t(X)), method=hclust.method)
	
	#ordonner matrice
	if(!Rowv )
		Xorder = X[,hcc$order]
	else if (!Colv)
		Xorder = X[hcr$order,]
	else if (!Colv & !Rowv)
		Xorder = X
	else
		Xorder = X[hcr$order, hcc$order]
	
	#scaling
	if (scale == "row") {
		Xorder <- sweep(Xorder, 1L, rowMeans(Xorder, na.rm = T), check.margin = FALSE)
		sXorder <- apply(Xorder, 1L, sd, na.rm = T)
		Xorder <- sweep(Xorder, 1L, sXorder, "/", check.margin = FALSE)
		hm.legend = "Z-score"
	}else if (scale == "column") {
		Xorder <- sweep(Xorder, 2L, colMeans(Xorder, na.rm = T), check.margin = FALSE)
		sXorder <- apply(Xorder, 2L, sd, na.rm = T)
		Xorder <- sweep(Xorder, 2L, sXorder, "/", check.margin = FALSE)
		hm.legend = "Z-score"
	} else if (scale == "none"){
		hm.legend = "Log2"
	}
	
	#modifier échelle de couleur
	if (saturate.cols ) {
		centile = quantile(Xorder, probs=seq(0,1,by=0.01))
		centile2 = centile[3]
		centile98 = centile[99]
		Xorder[which(Xorder > centile98)] = centile98
		Xorder[which(Xorder < centile2)] = centile2
	}
	
	display.layout(nbLabel)
	
	## -1- Plot titre
	par(mar = c(1,0.1,0.05, 1.1))
	plot( c(0,1), c(0,1), type="n", ylab="", xlab = "", axes = F)
	text(0.5,0.5,title,font=2,cex=1.5)
	
	## -2- Plot dendrogram
	# arrays
	a=marge.array.gene(labCol,labRow)
	par(mar=c(0.1, 0.5, 0.05, a$right))
	if(Colv)
		plot(as.dendrogram(hcc), axes=FALSE, horiz=F, yaxs="i", xaxs="i",leaflab = "none")
	else
		plot(0, bty="n", xaxt="n", yaxt="n", type="n", xlab="", ylab="")
	
	# genes
	par(mar=c(a$bottom, 1, 0.1, 0))
	if(Rowv)
		plot(as.dendrogram(hcr), axes=FALSE, horiz=T, yaxs="i", xaxs="i",leaflab = "none")
	else
		plot(0, bty="n", xaxt="n", yaxt="n", type="n", xlab="", ylab="")
	
	## -3- Plot heatmap
	par(mar=c(a$bottom, 0.5, 0.1, a$right))
	image(	x=1:ncol(Xorder),
			y=1:nrow(Xorder),
			z=t(Xorder), 
			col=Zcol, 
			axes=F, xaxs="i", yaxs="i",
			xlab="")
	
	#nom genes 
	label.row = labRow[hcr$order]
	if(!Rowv )
		label.row = rownames(X)
	
	label.col = labCol[hcc$order]
	if (!Colv)
		label.col = colnames(X)
	
	axis(	4, 1:nrow(Xorder), labels = label.row, 
			las = 2, line = -0.5, tick = 0, cex.axis = cex.array.gene(Xorder)["cexRow"])
	
	#nom arrays
	axis(	1, 1:ncol(Xorder), labels = label.col, 
			las = 2, line = -0.5, tick = 0, cex.axis = cex.array.gene(Xorder)["cexCol"])
	
	
	## -5- Plot legende bar de couleur
	min.raw <- min(Xorder, na.rm = TRUE)
	max.raw <- max(Xorder, na.rm = TRUE)
	z <- seq(min.raw, max.raw, length = length(Zcol))
	
	par(mar=c(0.5, 1, 2, 1))
	image(z = matrix(z, ncol = 1), col = Zcol, xaxt = "n", yaxt = "n", xlim=c(0, 1)) 
	mtext(hm.legend, side=3, cex=0.7)
	mtext(c("down", "up"), side=1,at=c(0.1, 0.9), cex=0.7)
	mtext(c(round(min.raw,1),round(max.raw,1)),side=3,line=-2,col="white",at=c(0.1,0.9),cex=0.7,font=2)
	
	## -6- Plot Annotations
	myColors = heatmap.colors()
	if(nbLabel != 0){
		for(i in 1:nbLabel){
			# Annotations
			par(mar = c(0,0.5,0, a$righ))
			myCol = myColors[i,getColor(annot[,i])]
			myCol = myCol[hcc$order]
			
			image(cbind(1:ncol(Xorder)),col=myCol,axes=F)
			
			#Legende Annotations
			nbLevel = length(levels(factor(annot[,i])))
			
			
			par(mar = c(0.2,5,0, 0.5))
			
			legendColor = myColors[i,getLegendColor(annot[,i])]
			if (length(which(is.na(annot[,i])))>0){
				nbLevel = nbLevel + 1
				legendColor = c("white",legendColor)
			}
			
			
			image(	x=c(1:nbLevel),y=c(0,0.5),z=matrix(seq(1,nbLevel,length=nbLevel),ncol=1),
					col= legendColor, ylim=c(0,1),
					axes=F, xaxs="i", yaxs="i",xlab="",ylab="")
			
			axis(side=2, at=0.4 ,labels=colnames(annot)[i],las=2, , tick=0, cex.axis=1)
			legendLab = getLegendLabel(annot[,i])
			
			axis(side=3, at=1:nbLevel ,labels=legendLab, , tick=0, cex.axis=0.8,line=-1.5)		
		}
	}
	
	return(list(hcc=hcc,hcr=hcr))
}



### Fonction avec possiblité de specifier deux type de distance pour le biclustering (echantillons / genes)

ugf.heatmap.biclust = function(	X, 
		gene.dist.method=c("pearson","euclidean","manhattan","L1"), 
		gene.hclust.method=c("ward","complete"),
		sample.dist.method=c("pearson","euclidean","manhattan","L1"), 
		sample.hclust.method=c("ward","complete"),
		
		cols=c("RedGreen","OrangeBlue","RedBlue"),
		annot,
		labRow=rownames(X),
		labCol=colnames(X),
		title="Heatmap",
		scale="row",
		Rowv=TRUE,
		Colv=TRUE,
		hcr=NULL,
		hcc=NULL,
		saturate.cols=F){
	
	nbSample = ncol(X)
	
	# --- Calcul de la matrice de distance
	if(gene.dist.method=="pearson"){
		g.dd = distPearson
	}else if(gene.dist.method=="spearman"){
		g.dd= distSpearman
	}else{
		g.dd= function(x){dist(x)}
	}
	#else if(dist.method="L1"){
	#	g.dd = distL1
	#}
	if(sample.dist.method=="pearson"){
		s.dd = distPearson
	}else if(sample.dist.method=="spearman"){
		s.dd= distSpearman
	}else{
		s.dd= function(x){dist(x)}
	}
	#else if(dist.method="L1"){
	#	g.dd = distL1
	#}
	
	
	#Colorisation de la heatmap
	if(cols=="RedGreen")
		Zcol = colz(low = "green", high = "red", mid="black", n=100)
	if(cols=="OrangeBlue")
		Zcol = colz(low = "blue", high = "orange", mid="white", n=100)
	if(cols=="RedBlue")
		Zcol = colz(low = "blue", high = "red", mid="white", n=100)
	
	#Annotations
	nbLabel = 0
	if(!missing(annot)){
		annot=as.matrix(annot)
		if(nrow(annot) != nbSample)
			stop("Number of arrays in annotation matrix is different from expression matrix")
		else{
			nbLabel = ncol(annot)
			if(nbLabel > 4) stop("Can't display more than 4 annotations")
		}
	}
	
	
	
	#ID Genes
	if(is.null(labRow))
		labRow=rep("",nrow(X))
	#ID Arrays
	if(is.null(labCol))
		labCol=rep("",ncol(X))
	
	#Calcul dendrogram genes
	if(is.null(hcr))
		hcr = hclust(g.dd(X), method=gene.hclust.method)
	
	#calcul dendrogram arrays
	if(is.null(hcc))
		hcc = hclust(s.dd(t(X)), method=sample.hclust.method)
	
	#ordonner matrice
	if(!Rowv )
		Xorder = X[,hcc$order]
	else if (!Colv)
		Xorder = X[hcr$order,]
	else if (!Colv & !Rowv)
		Xorder = X
	else
		Xorder = X[hcr$order, hcc$order]
	
	#scaling
	if (scale == "row") {
		Xorder <- sweep(Xorder, 1L, rowMeans(Xorder, na.rm = T), check.margin = FALSE)
		sXorder <- apply(Xorder, 1L, sd, na.rm = T)
		Xorder <- sweep(Xorder, 1L, sXorder, "/", check.margin = FALSE)
		hm.legend = "Z-score"
	}else if (scale == "column") {
		Xorder <- sweep(Xorder, 2L, colMeans(Xorder, na.rm = T), check.margin = FALSE)
		sXorder <- apply(Xorder, 2L, sd, na.rm = T)
		Xorder <- sweep(Xorder, 2L, sXorder, "/", check.margin = FALSE)
		hm.legend = "Z-score"
	}else if (scale == "none"){
		hm.legend = "Log2"
	}
	
	
	#modifier échelle de couleur
	if (saturate.cols ) {
		centile = quantile(Xorder, probs=seq(0,1,by=0.01))
		centile2 = centile[3]
		centile98 = centile[99]
		Xorder[which(Xorder > centile98)] = centile98
		Xorder[which(Xorder < centile2)] = centile2
	}
	
	display.layout(nbLabel)
	
	## -1- Plot titre
	par(mar = c(1,0.1,0.05, 1.1))
	plot( c(0,1), c(0,1), type="n", ylab="", xlab = "", axes = F)
	text(0.5,0.5,title,font=2,cex=1.5)
	
	## -2- Plot dendrogram
	# arrays
	a=marge.array.gene(labCol,labRow)
	par(mar=c(0.1, 0.5, 0.05, a$right))
	if(Colv)
		plot(as.dendrogram(hcc), axes=FALSE, horiz=F, yaxs="i", xaxs="i",leaflab = "none")
	else
		plot(0, bty="n", xaxt="n", yaxt="n", type="n", xlab="", ylab="")
	
	# genes
	par(mar=c(a$bottom, 1, 0.1, 0))
	if(Rowv)
		plot(as.dendrogram(hcr), axes=FALSE, horiz=T, yaxs="i", xaxs="i",leaflab = "none")
	else
		plot(0, bty="n", xaxt="n", yaxt="n", type="n", xlab="", ylab="")
	
	## -3- Plot heatmap
	par(mar=c(a$bottom, 0.5, 0.1, a$right))
	image(	x=1:ncol(Xorder),
			y=1:nrow(Xorder),
			z=t(Xorder), 
			col=Zcol, 
			axes=F, xaxs="i", yaxs="i",
			xlab="")
	
	#nom genes 
	label.row = labRow[hcr$order]
	if(!Rowv )
		label.row = rownames(X)
	
	label.col = labCol[hcc$order]
	if (!Colv)
		label.col = colnames(X)
	
	axis(	4, 1:nrow(Xorder), labels = label.row, 
			las = 2, line = -0.5, tick = 0, cex.axis = cex.array.gene(Xorder)["cexRow"])
	
	#nom arrays
	axis(	1, 1:ncol(Xorder), labels = label.col, 
			las = 2, line = -0.5, tick = 0, cex.axis = cex.array.gene(Xorder)["cexCol"])
	
	
	## -5- Plot legende bar de couleur
	min.raw <- min(Xorder, na.rm = TRUE)
	max.raw <- max(Xorder, na.rm = TRUE)
	z <- seq(min.raw, max.raw, length = length(Zcol))
	
	par(mar=c(0.5, 1, 2, 1))
	image(z = matrix(z, ncol = 1), col = Zcol, xaxt = "n", yaxt = "n", xlim=c(0, 1)) 
	mtext(hm.legend, side=3, cex=0.7)
	mtext(c("down", "up"), side=1,at=c(0.1, 0.9), cex=0.7)
	mtext(c(round(min.raw,1),round(max.raw,1)),side=3,line=-2,col="white",at=c(0.1,0.9),cex=0.7,font=2)
	
	## -6- Plot Annotations
	myColors = heatmap.colors()
	if(nbLabel != 0){
		for(i in 1:nbLabel){
			# Annotations
			par(mar = c(0,0.5,0, a$righ))
			myCol = myColors[i,getColor(annot[,i])]
			myCol = myCol[hcc$order]
			
			image(cbind(1:ncol(Xorder)),col=myCol,axes=F)
			
			#Legende Annotations
			nbLevel = length(levels(factor(annot[,i])))
			
			
			par(mar = c(0.2,5,0, 0.5))
			
			legendColor = myColors[i,getLegendColor(annot[,i])]
			if (length(which(is.na(annot[,i])))>0){
				nbLevel = nbLevel + 1
				legendColor = c("white",legendColor)
			}
			
			
			image(	x=c(1:nbLevel),y=c(0,0.5),z=matrix(seq(1,nbLevel,length=nbLevel),ncol=1),
					col= legendColor, ylim=c(0,1),
					axes=F, xaxs="i", yaxs="i",xlab="",ylab="")
			
			axis(side=2, at=0.4 ,labels=colnames(annot)[i],las=2, , tick=0, cex.axis=1)
			legendLab = getLegendLabel(annot[,i])
			
			axis(side=3, at=1:nbLevel ,labels=legendLab, , tick=0, cex.axis=0.8,line=-1.5)		
		}
	}
	
	return(list(hcc=hcc,hcr=hcr))
}



### Fonction avec plot pvalue
###Arguments : 

### cols = type de colorisation de la heatmap
### annot = matrice d'annotations
### labRow = vecteur d'annotations des genes. NULL si on ne veut pas plotter de noms
### labCol = vecteur d'annotations des samples. NULL si on ne veut pas plotter de noms
### scale = "row" ou "none" si pas de scaling
### Rowv = booleen. TRUE affiche le dendrogramme des gènes. FALSE sinon
### Colwv = booleen. TRUE affiche le dendrogramme des samples. FALSE sinon
### hcr = objet de classe hclust (si le dendrogramme des gènes est déjà calculé)
### hcc = objet de classe hclust (si le dendrogramme des samples est déjà calculé)
### pvalue = vecteur de pvalue
### pvalue.seuil = seuil de significativité de la pvalue

ugf.heatmap.pvalue = function(	X, 
		dist.method=c("pearson","euclidean","manhattan","L1"), 
		hclust.method=c("ward","complete"),
		cols=c("RedGreen","OrangeBlue"),
		annot,
		labRow=NULL,
		labCol=colnames(X),
		title="Heatmap",
		scale="row",
		Rowv=TRUE,
		Colv=TRUE,
		hcr=NULL,
		hcc=NULL,
		pvalue,
		pvalue.seuil=0.0001){
	
	
	nbSample = ncol(X)
	
	# --- Calcul de la matrice de distance
	if(dist.method=="pearson"){
		dd = distPearson
	}else if(dist.method=="spearman"){
		dd= distSpearman
	}else{
		dd= function(x){dist(x)}
	}
	#else if(dist.method="L1"){
	#	dd = distL1
	#}
	
	#Calcul de la pvalue en -log10
	pvalue = -log10(pvalue)
	
	#Colorisation de la heatmap
	if(cols=="RedGreen")
		Zcol = colz(low = "green", high = "red", mid="black", n=100)
	if(cols=="OrangeBlue")
		Zcol = colz(low = "blue", high = "orange", mid="white", n=100)
	
	#Annotations
	nbLabel = 0
	if(!missing(annot)){
		annot=as.matrix(annot)
		if(nrow(annot) != nbSample)
			stop("Number of arrays in annotation matrix is different from expression matrix")
		else{
			nbLabel = ncol(annot)
			if(nbLabel > 4) stop("Can't display more than 4 annotations")
		}
	}
	
	#ID Arrays
	if(is.null(labCol))
		labCol=rep("",ncol(X))
	
	#Calcul dendrogram genes
	if(is.null(hcr))
		hcr = hclust(dd(X), method=hclust.method)
	
	#calcul dendrogram arrays
	if(is.null(hcc))
		hcc = hclust(dd(t(X)), method=hclust.method)
	
	#ordonner matrice
	if(!Rowv ){
		Xorder = X[,hcc$order]
	}else if (!Colv){
		Xorder = X[hcr$order,]
	}else if (!Colv & !Rowv){
		Xorder = X
	}else{
		Xorder = X[hcr$order, hcc$order]
	}
	pvalueOrder =pvalue[hcr$order] 
	
	#scaling
	if (scale == "row") {
		Xorder <- sweep(Xorder, 1L, rowMeans(Xorder, na.rm = T), check.margin = FALSE)
		sXorder <- apply(Xorder, 1L, sd, na.rm = T)
		Xorder <- sweep(Xorder, 1L, sXorder, "/", check.margin = FALSE)
		hm.legend = "Z-score"
	}else if (scale == "column") {
		Xorder <- sweep(Xorder, 2L, colMeans(Xorder, na.rm = T), check.margin = FALSE)
		sXorder <- apply(Xorder, 2L, sd, na.rm = T)
		Xorder <- sweep(Xorder, 2L, sXorder, "/", check.margin = FALSE)
		hm.legend = "Z-score"
	} else if (scale == "none"){
		hm.legend = "Log2"
	}
	
	display.layout.pvalue(nbLabel)
	
	## -1- Plot titre
	par(mar = c(1,0.1,0.05, 1.1))
	plot( c(0,1), c(0,1), type="n", ylab="", xlab = "", axes = F)
	mtext(title,side=1,at=0.5,font=2,cex=0.9,line=-1)
	
	## -2- Plot dendrogram
	# arrays
	a=marge.array.gene(labCol,labRow)
	par(mar=c(0.1, 0.5, 0.05, 0.3))
	if(Colv){
		plot(as.dendrogram(hcc), axes=FALSE, horiz=F, yaxs="i", xaxs="i",leaflab = "none")
	}else{
		plot(0, bty="n", xaxt="n", yaxt="n", type="n", xlab="", ylab="")
	}
	# genes
	par(mar=c(a$bottom, 1, 0.1, 0))
	if(Rowv){
		plot(as.dendrogram(hcr), axes=FALSE, horiz=T, yaxs="i", xaxs="i",leaflab = "none")
	}else{
		plot(0, bty="n", xaxt="n", yaxt="n", type="n", xlab="", ylab="")
	}
	## -3- Plot heatmap
	par(mar=c(a$bottom, 0.5, 0.1, 0.3))
	image(	x=1:ncol(Xorder),
			y=1:nrow(Xorder),
			z=t(Xorder), 
			col=Zcol, 
			axes=F, xaxs="i", yaxs="i",
			xlab="")
	
	#nom arrays
	axis(	1, 1:ncol(Xorder), labels = labCol[hcc$order], 
			las = 2, line = -0.5, tick = 0, cex.axis = cex.array.gene(Xorder)["cexCol"])
	
	
	## -5- Plot legende bar de couleur
	min.raw <- min(Xorder, na.rm = TRUE)
	max.raw <- max(Xorder, na.rm = TRUE)
	z <- seq(min.raw, max.raw, length = length(Zcol))
	
	par(mar=c(4, 1, 2, 1))
	image(z = matrix(z, ncol = 1), col = Zcol, xaxt = "n", yaxt = "n", xlim=c(0, 1)) 
	mtext(hm.legend, side=3, cex=0.7)
	mtext(c("down", "up"), side=1,at=c(0.1, 0.9), cex=0.7)
	mtext(c(round(min.raw,1),round(max.raw,1)),side=3,line=-2,col="white",at=c(0.1,0.9),cex=0.7,font=2)
	
	## -6 plot pvalue
	par(mar=c(a$bottom, 0.5, 0.1, 2))
	plot(0,0,xlim=c(0,5),ylim=c(0,nrow(X)),axes=F, xaxs="i", yaxs="i", xaxt="n", yaxt="n", type="n",xlab="")
	rect(0, 0:nrow(X)-1, pvalueOrder, 1:nrow(X), col="gray",border=NA)
	axis(3)
	box()
	mtext("-log10(pvalue)",side=1,cex=0.5)
	abline(v=-log10(pvalue.seuil),col="red",lty=2)
	
	## -7- Plot Annotations
	myColors = heatmap.colors()
	if(nbLabel != 0){
		for(i in 1:nbLabel){
			# Annotations
			par(mar = c(0,0.5,0, 0.3))
			myCol = myColors[i,getColor(annot[,i])]
			myCol = myCol[hcc$order]
			
			image(cbind(1:ncol(Xorder)),col=myCol,axes=F)
			
			#Legende Annotations
			nbLevel = length(levels(factor(annot[,i])))
			par(mar = c(0.2,5,0, 0.5))
			
			
			legendColor = myColors[i,getLegendColor(annot[,i])]
			if (length(which(is.na(annot[,i])))>0){
				nbLevel = nbLevel + 1
				legendColor = c("white",legendColor)
			}
			
			image(	x=c(1:nbLevel),y=c(0,0.5),z=matrix(seq(1,nbLevel,length=nbLevel),ncol=1),
					col=legendColor , ylim=c(0,1),
					axes=F, xaxs="i", yaxs="i",xlab="",ylab="")
			axis(side=2, at=0.4 ,labels=colnames(annot)[i],las=2, , tick=0, cex.axis=1)
			legendLab = getLegendLabel(annot[,i])
			
			axis(side=3, at=1:nbLevel ,labels=legendLab, , tick=0, cex.axis=0.8,line=-1.5)		
		}
	}
	
	
	
	return(list(hcc=hcc,hcr=hcr))
}


#####################################################
#						FONCTIONS					#			
#####################################################


# Coloration des annotations
heatmap.colors = function(){
	bleu = c("blue","cyan","cornflowerblue","midnightblue","blueviolet","aliceblue","magenta","darkolivegreen1")
	rouge = c("red","lightpink1","orange2","yellow1","deeppink1","orange4","pink","wheat")
	vert = c("darkgreen","green","darkolivegreen1","darkslategrey","seagreen","turquoise2","yellow2","wheat")
	pink = c("mediumvioletred","purple","palevioletred","pink4","rosybrown1","magenta","lightslateblue","wheat")
	myColors= as.matrix(rbind(
					bleu = bleu,
					rouge = rouge,
					vert = vert,
					pink = pink))
}


### Generation d'un vecteur de couleurs pour la heatmap
colz <- function(low = lowcol, high = highcol, mid=midcol, n=100) {
	if (is.character(low)) low <- col2rgb(low)/255
	if (is.character(high)) high <- col2rgb(high)/255
	if (is.character(mid)) mid <- col2rgb(mid)/255
	col <- rgb(c(seq(low[1], mid[1], len=n/2), seq(mid[1], high[1], len = n/2)), c(seq(low[2], mid[2], len=n/2), seq(mid[2], high[2], len = n/2)), c(seq(low[3], mid[3], len=n/2), seq(mid[3], high[3], len = n/2)))
	return(col)
}

#distance de pearson
distPearson=function(x){as.dist((1-cor(t(x),use="complete.obs"))/2)}

#distance de Spearman
distSpearman=function(x){as.dist((1-cor(X, method="spearman"))/2)}

### distance L1 pour les profils
### formule = mean(abs(X2 - X2))
distL1 = function(X, do.logtransformation = TRUE){
	
	### calcul de l'array de reference
	refArray = apply(X$E, 1, median)
	
	res = matrix(nrow=ncol(X),ncol=ncol(X))
	
	colnames(res) = rownames(res) = colnames(X)
	for (i in 1:ncol(X)){
		for (k in i:ncol(X)){
			
			if (do.logtransformation){
				v1 = log2(X$E[,i]/refArray)
				v2 = log2(X$E[,k]/refArray)
			}
			else{
				v1 = X$E[,i]/refArray
				v2 = X$E[,k]/refArray
			}
			val = mean(abs(v1-v2))
			res[i,k] = val
			res[k,i] = val	
		}
	}
	return(res)
}



# Fonction qui retourne 
getColor = function(x){
	
	res = x
	f = as.factor(res)
	
	res[which(is.na(x))] = 0
	
	for (i in 1 : length(levels(f))){
		l = levels(f)[i]
		res[which(x==l)] = i
	}
	as.numeric(res)
	
}


#retourne le numéro de colonne de la matrice myColors
getLegendColor = function(x){
	
	f = as.factor(x)
	res = c()
	if (length(which(is.na(x))) > 0 ){
		res =  c(0)
	}
	
	for (i in 1 : length(levels(f))){
		colF = i
		res = c(res,colF)	  	
	}
	as.numeric(res)
}

getLegendLabel = function(x){
	
	f = as.factor(x)
	res = c()
	if (length(which(is.na(x))) > 0 ){
		res =  c("NA")
	}
	
	for (i in 1 : length(levels(f))){
		l = levels(f)[i]
		res = c(res,l)	  	
	}
	res   
}


#x = nb d'annotation
display.layout = function(x){
	
	if(is.null(x)) x=0
	nr = x + 3       # nombre de ligne du layout
	nc = 2          # nombre de colonne du layout
	
	myStepL = 0.02 #hauteur de chaque annot
	myStepB = if(x==0) 0.8 else 0.75 #heatmap
	myStepH = (1-myStepB)-x*myStepL #hauteur dendrogram array
	
	zone=NULL
	zone = matrix(,nrow=nr , ncol=nc ,byrow=T)
	zone[1,1]=5
	zone[1,2]=1
	zone[2,1]=0
	zone[2,2]=2
	zone[nr,1]=3
	zone[nr,2]=4
	
	if(x >0){
		for (i in 1:x){
			zone[i+2,2] = max(zone,na.rm=T)+1
			zone[i+2,1] = max(zone,na.rm=T)+1
		}
	}

	nrz = nrow(zone)
	zone = cbind(zone, rep(0, nrz))
	
	layW = c(0.20,0.75, 0.05)
	layH = c(0.1,myStepH, rep(myStepL,x), myStepB)
	
	
	layout(zone, widths=layW, heights = layH)
#	layout.show(max(zone))
}

#x = nb d'annotation
display.layout.pvalue = function(x){
	
	if(is.null(x)) x=0
	nr = x + 3       # nombre de ligne du layout
	nc = 3          # nombre de colonne du layout
	
	myStepL = 0.02 #hauteur de chaque annot
	myStepB = if(x==0) 0.8 else 0.75 #heatmap
	myStepH = (1-myStepB)-x*myStepL #hauteur dendrogram array
	
	zone=NULL
	zone = matrix(,nrow=nr , ncol=nc ,byrow=T)
	zone[1,1]=0
	zone[1,2]=1
	zone[1,3]=0
	zone[2,1]=5
	zone[2,2]=2
	zone[2,3]=0
	zone[nr,1]=3
	zone[nr,2]=4
	zone[nr,3]=6
	
	
	if(x >0){
		for (i in 1:x){
			zone[i+2,2] = max(zone,na.rm=T)+1
			zone[i+2,1] = max(zone,na.rm=T)+1
			zone[i+2,3] = 0 
		}
	}
	
	layW = c(0.25,0.65,0.1)
	layH = c(0.05,myStepH, rep(myStepL,x), myStepB)
	
	
	layout(zone, widths=layW, heights = layH)
#	layout.show(max(zone))
}


# Calcul la taille des caracteres nom de genes et nom arrays
# en fonction de la taille de la matrice
cex.array.gene = function(X){
	return(c(cexRow=0.2 + 1/log10(nrow(X)), cexCol = 0.2 + 1/log10(ncol(X))))
}

# Ajustement des marges de la heatmap en fonction de
# la longueur des noms de colonnes et noms de gènes
marge.array.gene = function(labCol, labRow){
	col=as.numeric(max(nchar(labCol)))
	row=as.numeric(max(nchar(labRow)))
	
	if(all(col <=10 & row <=10))
		return(list(bottom=5,right=5))
	if(all(col <=10 & row >10 & row <=20))
		return(list(bottom=5,right=7.5))
	if(all(col <=10 & row > 20)){
		cat("Le nom des genes est trop long pour apparaitre en entier\n")
		return(list(bottom=5,right=7.5))
	}
	
	if(all(col > 10 & col <=20 & row <=10))
		return(list(bottom=7.5,right=5))	
	if(all(col > 10 & col <=20 & row > 10 & row <=20))
		return(list(bottom=7.5,right=7.5))
	if(all((col > 10 & col <=20) & row > 20)){
		cat("Le nom des genes est trop long pour apparaitre en entier\n")
		return(list(bottom=7.5,right=7.5))
	}
	
	
	if(all(col > 20  & row <=10)){
		cat("Le nom des arrays est trop long pour apparaitre en entier\n")
		return(list(bottom=7.5,right=5))	
	}	
	if(all((col > 20  & row > 10) & row <=20)){
		cat("Le nom des arrays est trop long pour apparaitre en entier\n")
		return(list(bottom=7.5,right=7.5))
	}			
	if(all(col > 20  & row > 20)){
		cat("Le nom des arrays est trop long pour apparaitre en entier\n")
		cat("Le nom des genes est trop long pour apparaitre en entier\n")
		return(list(bottom=7.5,right=7.5))
	}
	
}



