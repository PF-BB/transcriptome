# script pour shiny
# l'utilisateur choisit une table 


# l'utilisateur choisit un Rdata contenant un eset
# puis il choisi un gene

# 1) transformer en data frame
require(ggplot2)
geneData = data.frame(int=2^(data@assayData$exprs[10,]),ID = colnames(data),y = data@phenoData@data$temps)

# 2) plot
p<-ggplot(geneData,aes(x=ID,y=int,fill=ID))
p <- p + geom_bar(position=dodge,stat='identity')
p<-p + scale_fill_hue(l=50,h.start=200) + scale_color_hue(l=50,h.start=200)
p

# idem mais les couleurs correspondent à un facteur
p<-ggplot(geneData,aes(x=ID,y=int,fill=y))
p <- p + geom_bar(position=dodge,stat='identity')
p<-p + scale_fill_hue(l=50,h.start=200) + scale_color_hue(l=50,h.start=200)
p <- p + scale_x_discrete("",breaks=tracking_ids,labels=gene_labels) + 
  theme(axis.text.x=element_text(hjust=0,angle=-90))
p

#si on veut plotter en fonction d'un facteur sans voir tous les arrays
showErrorbars=TRUE
conf_low=conf_high=NULL
geneData.mean = geneData
for(i in levels(geneData$y)){
  t=t.test(geneData$int[geneData$y==i])
  geneData.mean$int[geneData$y==i] = t$estimate
  conf_low = c(conf_low, t$conf.int[1])
  conf_high = c(conf_high, t$conf.int[2])
}
geneData.mean = unique(geneData.mean[,c(1,3)])
rownames(geneData.mean) = geneData.mean$y
geneData.mean = geneData.mean[levels(geneData$y),]
geneData.mean$conf_l = conf_low
geneData.mean$conf_h = conf_high
p<-ggplot(geneData.mean,aes(x=y,y=int,fill=y))
p <- p + geom_bar(position=dodge,stat='identity')
if (showErrorbars)
  p <- p + geom_errorbar(aes(ymin=conf_l,ymax=conf_h),position=dodge,width=0.5)
p<-p + scale_fill_hue(l=50,h.start=200) + scale_color_hue(l=50,h.start=200)
p
