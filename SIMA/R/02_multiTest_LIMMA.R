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
