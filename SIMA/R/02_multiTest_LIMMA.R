multiTest_LIMMA <-
function(X,y){
	g <- nrow(X)
	design <- cbind(y,1-y)
	colnames(design) <- c("un","deux")
	cont <- makeContrasts(un-deux,levels=design)
	fit  <- lmFit(X,design)
	fit <- contrasts.fit(fit,cont)
	ordinary.t <- fit$coef / fit$stdev.unscaled / fit$sigma
	fit  <- eBayes(fit)
	top  <- topTable(fit,number=g,adjust.method="BH",sort.by="P")
	list("top"=top,"ordinary.t"=ordinary.t)
}

#Fonction qui reduit un topTable selon seuil FC et pvalue et qui enleve les log2
topTableReducer = function(top, FC = 2, PV=0.05, p.val = c("adj","raw")){
  FC = log2(FC)
  
  if (p.val == "adj")
    idxPV = which(top$adj.P.Val  <= PV)
  
  else if (p.val == "raw")
    idxPV = which(top$P.Value  <= PV)
  
  #quelles colonnes correspondent à des ratios ?
  # si t-test il y a une seule colonne, si ANOVA on ne sait pas combien de contrast il y a dans l'objet top
  idxContrast = c((which(colnames(top)== "ID")+1):(which(colnames(top)== "AveExpr")-1))
  
  if(idxContrast > 1)
    idxFC = which(apply(abs(top[,idxContrast]),1,max) >= FC)
  else 
    idxFC = which(abs(top[,idxContrast]) >= FC)
  
  red = intersect(idxFC,idxPV)
  
  if (length(red) > 0){
    topred = top[red,]
    ID = as.vector(topred[,1])
    PValue = topred[,"P.Value"]
    adj.PValue = topred[,"adj.P.Val"]
    ratio = 2^topred[,idxContrast]
    if(idxContrast > 1)
      fc = apply(ratio,c(1,2), function(r){ if(r < 1){ r = -1/r}else{r}})
    else
      fc = sapply(ratio, function(r){ if(r < 1){ r = -1/r}else{r}})
    AveExpr = 2^(topred$AveExpr)
    
    topred = as.data.frame(cbind(ID, round(fc,3), round(AveExpr,3), sprintf("%0.3e",PValue), sprintf("%0.3e",adj.PValue)))
    colnames(topred) = c("ID",colnames(top)[idxContrast],"AveExpr","P.Value","adj.P.Val")
    return(topred)
  }
  else
    stop("Aucun gene differentiel selon ces critères ! ")
  
}
