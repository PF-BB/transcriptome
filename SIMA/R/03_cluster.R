#' Heatmap with annotations.
#' @param X numeric matrix of the values to be plotted
#' @param pheno a vector or matrix to annotate the samples
#' @param dist.method function used to compute the distance (dissimilarity) between both rows and columns
#' @param hclust.method function used to compute the hierarchical clustering
#' @param title character string for the title of the cluster
#' @param sample.name vector for names of the sample
#' @param rect.hclust integer for numbers of groups : draw k rectangles
#' @return \item{hclust object}{An object of class hclust.}
#' @references 
#' @title Draw a Cluster with Annotations
#' @export clusterWithParameter


library(RColorBrewer)


### ------------------------------------------------------------------------ ###
clusterWithParameter = function(X, pheno, title="",sample.name=colnames(X),dist.method=c("pearson","euclidean","manhattan","spearman"), hclust.method = "ward", rect.hclust.k=0 ){
	
	bleu   = c("blue","cyan","cornflowerblue","midnightblue","blueviolet","aliceblue","magenta","darkolivegreen1")
	rouge  = c("red","lightpink1","orange2","yellow1","deeppink1","orange4","pink","wheat")
	vert   = c("darkgreen","green","darkolivegreen1","darkslategrey","seagreen","turquoise2","yellow2","wheat")
	pink   = c("mediumvioletred","purple","palevioletred","pink4","rosybrown1","magenta","lightslateblue","wheat")
	marron = c("#35201D","#8B7355","#DEB887","#8B4500","#CD661D","#8B3E2F","#380700","#5B0000")
	violet = rev(c("#F7FCFD","#BFD3E6","#9EBCDA","#8C96C6","#8C6BB1","#88419D","#810F7C","#4D004B"))
	set4 = rev(brewer.pal(8,"Paired"))
	set5 = rev(brewer.pal(8,"Set1"))
	set6 = rev(brewer.pal(8,"Set2"))
	set7 = rev(brewer.pal(8,"Set3"))
	set8 = rev(brewer.pal(8,"Pastel1"))
	cluster.colors = as.matrix(rbind(
					bleu = bleu,
					rouge = rouge,
					vert = vert,
					pink = pink,
					marron = marron,
					violet = violet, 
					set4 = set4[c(1,4,7,2,5,8,3,6)],
					set5 = set5[c(1,4,7,2,5,8,3,6)],
					set6 = set6[c(1,4,7,2,5,8,3,6)],
					set7 = set7[c(1,4,7,2,5,8,3,6)],
					set8 = set8[c(1,4,7,2,5,8,3,6)]))
	cluster.colors = cbind(rep("white",nrow(cluster.colors)),cluster.colors)
	
	myXcol  = c(0,10,20,30,0,10,20,30)
	myYcol1 = c(-0.3,-0.3,-0.3,-0.3)
	myYcol2 = c(0.2,0.2,0.2,0.2,-0.8,-0.8,-0.8,-0.8,-0.8)
	myYlab1 = 0
	myYlab2 = c(0.5,0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5)
	
	#check arguments
	if(is.factor(pheno))
		pheno= as.vector(pheno)
	
	if (length(sample.name) != ncol(X)){
		if (is.null(colnames(X))){
			warning("spurious length of sample.name : setting sample.name to empty labels")
			sample.name = rep("",ncol(X))
		}
		else{
			warning("spurious length of sample.name : using colnames(X) to set sample.name")
			sample.name = colnames(X)
		}
	}
	
	# --- Calcul de la matrice de distance
	if(dist.method=="pearson"){
		dd= as.dist((1-cor(X))/2)
	}
	else if (dist.method=="spearman"){
		dd= as.dist((1-cor(X, method="spearman"))/2)
	}
	else {
		dd = dist(t(X), method=dist.method)
	}
	
	# --- Calcul du cluster
	hc = hclust(dd,method=hclust.method)
	coph.d = cophenetic(hc)
	hcLabel = hc$labels
	hc$labels = ""
	coph.cor = cor(dd,coph.d)
	
	# --- Nombre de sample
	nbSample = ncol(X)
	print (paste("Nb of sample: ",nbSample,sep=""))
	
	
	# Verification du label / definition des labels 'utils'
	myLabelName = cluster.getLabel(pheno)
	
	nbLabel     = length(myLabelName)
	print (paste("Nb of labels: ",nbLabel,sep=""))
	print (paste("Names       : ",myLabelName,sep=""))
	
	# Definition du layout
	cluster.displayLayout(nbLabel)
	
	# Plot Dendrogram
	par(mar = c(0.1,0.1,2.1, 1.1))
	dendroName = title
	plot(as.dendrogram(hc),xlim=c(1,nbSample), main="")
	
	#plot rectangle
	print(rect.hclust.k)
	if (rect.hclust.k>1){
		hc$rect = rect.hclust(hc, k=rect.hclust.k)
	}
	
	# Plot MetaData
	myMaxLen = max(strwidth(myLabelName,cex=1))
	
	for (i in 1:nbLabel){
		myLab = myLabelName[i]
		myLabelData = c()
		if (nbLabel == 1){
			myLabelData = pheno
		}
		else{
			myLabelData = pheno[,myLab]
		}
		
		myLabF = as.numeric(as.factor(myLabelData))
		
		# Plot n : modalit�s du facteur
		par(mar = c(0.1,0.1,0.1, 1.1))
		plot( c(1,nbSample), c(-1,1), type="n", ylab="", xlab = "", axes = F)
		myCol = cluster.colors[i,cluster.getColor(myLabF)]
		myCol = myCol[hc$order]
		rect(
				xleft   = (1:nbSample)-0.5,
				ybottom = -0.5,
				xright  = (1:nbSample)+0.5,
				ytop    = 0.5,
				border  = "black",
				col     =myCol
		)
		#abline(h=c(-1,0,1))
		
		
		# Plot n+1 : La legende
		par(mar = c(0.1,0.1,0.1,0.1),adj = 1)
		plot(x=0,y=0,type = "n", axes=F, ylab="",xlab="", xlim=c(0,40))
		
		legendCol = cluster.colors[i,cluster.getLegendColor(myLabF)]
		legendLab = cluster.getLegendLabel(myLabelData)
		
		xCol = myXcol[1:length(legendCol)]
		yCol =c()
		xLab = xCol + 4
		yLab = c()
		if (length(legendCol) > 4){
			yCol = myYcol2[1:length(legendCol)]
			yLab = myYlab2[1:length(legendCol)]
		}
		else if(length(legendCol) > 8){
			print ("Sorry, can't display legend with more than 8 levels")
		}
		else{
			yCol = myYcol1[1:length(legendCol)]
			yLab = myYlab1
		}
		
		rect( xleft   = xCol ,
				ybottom = yCol,
				xright  = xCol+3 ,
				ytop    = yCol + 0.6,
				col=legendCol,
				border="black")
		
		text( x   = xLab ,
				y   = yLab,
				labels  = legendLab,
				cex  = 0.9,
				adj  = c(0,0.5))
		
		#abline(v=c(0,10,20,30,40),h=c(-1,0,1))
		
		# Plot n+2 : Le nom du facteur      
		par(mar = c(0.1,1.1,0.1,0.1))
		plot(x=0,y=0,type = "n", axes=F, ylab="",xlab="",xlim=c(0,14))
		text(x = 0, y = 0, labels =  myLab, cex = 1.1, adj=0 )
		#abline(v=c(0,14),h=c(-1,0,1))
	}
	
	# Plotting cluster sample name
	par(mar = c(0.1,0.1,0.1, 1.1))
	plot( c(1,nbSample), c(0,1), type="n", ylab="", xlab = "", axes = F)
	text(  x   = c(1:nbSample) ,
			y   = 1,
			labels  = sample.name[hc$order],
			cex  = 1.3,
			srt = 90,
			adj  = c(1,0.5))
	#abline(v=c(1:nbSample))
	
	# plotting info up-left
	par(mar = c(1.1,0.1,2.1, 1.1))
	plot( c(0,1), c(0,1), type="n", ylab="", xlab = "", axes = F)
	#abline(h=c(0,1), v=c(0,1), col="green")
	text(x=0,y=1,    cex=2.5, labels=dendroName, adj=0, font=2)
	text(x=0,y=0.90, cex=2, labels=paste("Distance:",dist.method), adj=0)
	text(x=0,y=0.80, cex=2, labels=paste("Cluster  :",hc$method)  , adj=0)
	text(x=0,y=0.70, cex=2, labels=paste("Cor.cophenetic:",format(coph.cor,digits=3))  , adj=0)
	
	hc$labels = hcLabel  
	return(hc)
}

#retourne le numéro de colonne de la matrice cluster.colors
cluster.getLegendColor = function(x){
	
	f = as.factor(x)
	res = c()
	if (length(which(is.na(x))) > 0 ){
		res =  c(1)
	}
	
	for (i in 1 : length(levels(f))){
		colF = i+1
		res = c(res,colF)	  	
	}
	as.numeric(res)
}


cluster.getLegendLabel = function(x){
	
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

#retourne le num�ro de colonne de la matrice cluster.colors
cluster.getColor = function(x){
	res = x
	f = as.factor(res)
	
	res[which(is.na(x))] = 1
	
	for (i in 1 : length(levels(f))){
		l = levels(f)[i]
		res[which(x==l)] = i+1
	}
	as.numeric(res)
}

cluster.displayLayout = function(x){
	nr = x +2       # nombre de ligne du layout
	nc = 3          # nombre de colonne du layout
	
	myStepL = 0.05
	myStepB = 0.2
	myStepH = (1-myStepB)-x*myStepL
	
	
	zone = matrix(nrow=nr , ncol=nc ,byrow=T)
	zone[1,1] = 3*(nr-2)+3
	zone[1,2] = 3*(nr-2)+3
	zone[1,3] = 1
	
	for (i in 1:x){
		zone[i+1,] = rep((3*(i-1)+4):(3*(i-1)+2),1)
	}
	zone[x+2,1] = 3*(nr-2)+4
	zone[x+2,2] = 3*(nr-2)+4
	zone[x+2,3] = 3*(nr-2)+2
	
	layW = c(0.15,0.2,0.65)
	layH = c(myStepH, rep(myStepL,x), myStepB)
	layout(zone, widths=layW, heights = layH)
	#layout.show(3*(nr-2)+4)
}

cluster.getLabel = function(pheno){
	res = c()
	if (is.null(ncol(pheno))){
	  res = c("")
	}
	else{
		res = colnames(pheno)
	}
	res
}



