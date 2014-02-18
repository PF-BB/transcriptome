ecrireLight <-
function(logValues,top,fichier="tableau.tsv") {
  smbls <-as.list(hgu133plus2SYMBOL)
	Ncol  <- dim(logValues)[2]
	g     <- dim(logValues)[1]
	ordreID <- rownames(logValues)
  temp <- top
  print(names(temp))
	temp <- temp[order(temp$ID),]
	DF <- sign(temp$logFC)*2**abs(temp$logFC)
  tableau <- dataframe(ID=ordreID,logValues,"logFC"=temp$logFC,"DF"=DF,"p-value"=temp$P.Value,"adjusted p-value"=temp$adj.P.Val)
	print(tableau[1:10,])
	write.table(tableau,titre,dec=",",sep="\t",row.names=FALSE)
}
