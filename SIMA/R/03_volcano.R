volcano <-
function(top,seuilFC,seuilPval,correction = correction, titre = "Volcano Plot") {
	# Calcul de l'abscisse et de l'ordonnee
	if (correction) {B <- -log10(top$adj.P.Val); lab <- "-log10(adjusted p-value)"; titre <- paste(c(titre," (FDR)"),collapse="",sep="")}
	else {B <- -log10(top$P.Value);lab <- "-log10(p-value)"}
	
	FC <- top$logFC
	# Trois groupes de genes : UR, DR et le reste.
	gauche <- intersect( which(B>-log10(seuilPval)) , which(FC < -seuilFC))
	droite <- intersect( which(B>-log10(seuilPval)) , which(FC > seuilFC))

	textseuil1 <- paste("Seuil FC = ",seuilFC,collapse="",sep="")
	textseuil2 <- paste("p-value = ",seuilPval,collapse="",sep="")
	UR <- paste("UR = ",length(droite),collapse="",sep="")
	DR <- paste("DR = ",length(gauche),collapse="",sep="")

	g <- rep(1,length(FC)); g[gauche] <- 2; g[droite] <- 3
	donnees <- data.frame(FC=FC,B=B,g=g)
	print(donnees[sample(1:length(FC),10),])

	graphes <- xyplot(B ~ FC,groups=g, panel = function(x,y,subscripts,groups) { panel.abline(h=-log10(seuilPval),v=c(-seuilFC,seuilFC),lty=2); panel.superpose(x,y,groups=g,subscripts=subscripts,col=c("black","green","red"))} ,xlab="log2(FC)",ylab=lab, main = titre ,key=list(space="top",col=c("black","black","red","green"), text=list(lab=c(textseuil1,textseuil2,UR,DR)),columns = 2))	
	plot(graphes)
}
