library(shiny)
require(ggplot2)

ggplot2_barplot <- function(X,phenoValue){
  
  if(missing(phenoValue) | is.null(phenoValue) | phenoValue=="")
    phenoValue="ID"
  
  dodge <- position_dodge(width=0.9) 
  p<-eval(parse(text=paste0("ggplot(X,aes(x=",phenoValue,",y=int,fill=ID))")))
  p <- p + geom_bar(position=dodge,stat='identity')
  p <-p + scale_fill_hue(l=50,h.start=200) + scale_color_hue(l=50,h.start=200)
  print(p)
  
}

 ggplot2_barplot_rep <- function(X,phenoValue){

  if(missing(phenoValue) | is.null(phenoValue) | phenoValue=="")
    phenoValue="ID"
  
  dodge <- position_dodge(width=0.9) 
  p <-eval(parse(text=paste0("ggplot(X,aes(x=",phenoValue,",y=int,fill=ID))")))
  p <- p + geom_bar(position=dodge,stat='identity')
  p <-p + scale_fill_hue(l=50,h.start=200) + scale_color_hue(l=50,h.start=200)
  p <- p + facet_grid(probe~.)
  print(p)

}


shinyServer(function(input, output) {
  
  ##### GET DATA #######
  
  datasetInput <- reactive({
    datasetInput = read.table(paste0("../",input$dataset,"_data.txt"),sep="\t", header=T)
  })

  phenoInput <- reactive({
    phenoInput = read.table(paste0("../",input$dataset,"_pheno.txt"),sep="\t", header=T)
  })
    
  output$pheno <- renderUI({
    radioButtons("phenotype",
                  label="Choose phenotype",
                  c(NoPhenotype="",colnames(phenoInput())[-1])
    )
  })
  
  getGene <- reactive({
    if (!nzchar(input$gene))
      return()
    gene=grep(paste0("^",input$gene,"$"),datasetInput()[,1],ignore.case=TRUE)
    if(length(gene> 0))
      return(gene)
    else
      stop("Gene non connu")
  })
  
  
  #### OUTPUT ######
  myheight <- function(){
    myGene = getGene()
    210*length(myGene)
  }

  output$view <- renderPlot({
    myGene = getGene()
    if(is.null(myGene) ) {
      #fudge to stop shiny trying to plot NULL, which results in an ugly black square:
      return(plot(1,type="n",bty="n",yaxt="n",xaxt="n",ylab="",xlab=""))
    }
    
    else if(length(myGene) == 1){
      mat = data.frame(t(datasetInput()[myGene,-1]), colnames(datasetInput())[-1], phenoInput()[,-1])
      colnames(mat) = c("int","ID",colnames(phenoInput())[-1])
      ggplot2_barplot(mat,input$phenotype)  
    }
    else if (length(myGene) > 1){
      data = data.frame(c(t(datasetInput()[myGene,-1])),
                        rep(colnames(datasetInput())[-1],times=length(myGene)),
                        rep(c(1:length(myGene)),each=ncol(datasetInput()[,-1])))
      tmp = NULL
      for(i in 1:(length(myGene))){
        tmp = rbind(tmp,as.matrix(phenoInput()[,-1]))
       }
      mat = cbind(data, tmp)
      colnames(mat) = c("int","ID","probe",colnames(phenoInput())[-1])
      ggplot2_barplot_rep(mat,input$phenotype)
    }
    
  }, height=myheight)

})
