library(shiny)


f.loadData <- function(rdata) {
  myData = new.env()  # create environment for the data
  load(rdata, envir = myData)
  myData  # return the data
} 


shinyServer(function(input, output) {
  
  datasetInput <- reactive({
    paste0("../",input$dataset,".RData")

  })
  

  output$view <- renderPlot({
    load(datasetInput())
    mat.env = f.loadData(datasetInput())
    #mat.env$mat #ne pas utliser [[]] pour les variables environnement
    #     barplot(subset(mat.env$mat,Gene==input$gene))
        barplot(mat.env$mat[input$gene,])
    
    })
})
