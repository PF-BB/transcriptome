ecrire <-
function(logValues,top,titre="Tableau.tsv") {
	Ncol  <- dim(logValues)[2]
	g     <- dim(logValues)[1]
	Ncomp <- length(top)
	print(c(Ncol,g,Ncomp))
	tableau <- as.data.frame(matrix(0,g,1+Ncol+4*Ncomp))
	ordreID <- rownames(logValues)
	noms_tab <- vector(length=1+Ncol+4*Ncomp)
	tableau[,1] <- ordreID
	noms_tab[1] <- "ID"
	tableau[,1+1:Ncol] <- logValues
	noms_tab[1+1:Ncol] <- colnames(logGCRMA)
	
	for (ind in 1:Ncomp){
		print(dim(tableau))
		temp <- top[[ind]]
		print(names(temp))
		temp <- temp[order(temp$ID),]
		DF <- sign(temp$logFC)*2**abs(temp$logFC)
		tableau[,1+Ncol+4*(ind-1)+1:4] <- cbind(temp$logFC,DF,temp$P.Value,temp$adj.P.Val)
		noms_tab[1+Ncol+4*(ind-1)+1:4] <- paste(names(top[ind]), c("logFC","DF","p-value","adjusted p-value"))
	}
	colnames(tableau)<-noms_tab
	print(tableau[1:10,])
	write.table(tableau,titre,dec=",",sep="\t",row.names=FALSE)
}
