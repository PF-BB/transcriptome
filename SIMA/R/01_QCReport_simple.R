QCReport_simple <-
function(projet,nom_Etude){
	# Le tableau des ratio 3'/5'
	QCproj <- qc(projet) 
	ind1 <- 1; ind2 <- 3; ind3 <- 4; ind4 <- 6
	prob1 <- substr(colnames(QCproj@qc.probes)[ind1], 1 , nchar(colnames(QCproj@qc.probes )[ind1])-5)
	prob2 <- substr(colnames(QCproj@qc.probes )[ind3], 1, nchar(colnames(QCproj@qc.probes )[ind3])-5)
	ratios <- cbind(QCproj@qc.probes[,ind1] - QCproj@qc.probes[,ind2], QCproj@qc.probes[,ind3] - QCproj@qc.probes[,ind4]) 
	rownames(ratios) <- rownames(QCproj@qc.probes); colnames(ratios) <- c(prob1,prob2)
	expratios <- 2**ratios
	
	pdf(paste(c(nom_Etude,"_QC.pdf"),sep="",collapse=""))
	afficheTable3col(2**ratios)
	plotQC2(QCproj,cex=0.5)
	plot.new();titlePage(projet)
	signalDist2(projet)
	dev.off()
}
