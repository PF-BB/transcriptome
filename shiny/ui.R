library(shiny)


shinyUI(pageWithSidebar(
  headerPanel("Gene Profile"),
  sidebarPanel(
    selectInput(inputId="dataset",
                label="Choose a dataset:", 
                choices =gsub("_data\\.txt","",dir("/mnt/ihu-3-analyses-users/PROJET/LOBSIGER_GEO_ARRAYEXPRESS/data4shiny",pattern="*_data.txt"))
                ),
    
    uiOutput("pheno"),
    textInput("gene","Gene"),
    submitButton(text = "Submit")
    
  ),
  
  mainPanel(
    conditionalPanel(condition="$('html').hasClass('shiny-busy')",img(src="http://i1215.photobucket.com/albums/cc508/nawel12/loading.gif")),
    plotOutput("view")
  )
))



#