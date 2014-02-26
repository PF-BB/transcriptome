library(shiny)


shinyUI(pageWithSidebar(
  headerPanel("Loading data"),
  sidebarPanel(
    selectInput(inputId="dataset",
                label="Choose a dataset:", 
                choices =gsub("_data\\.txt","",dir("/mnt/ihu-3-analyses-users/PROJET/LOBSIGER_GEO_ARRAYEXPRESS/data4shiny",pattern="*_data.txt"))
                ),
    
    textInput("gene","Gene"),
          
    # Display this only if the phenotype is shown
    conditionalPanel(condition="input.gene.length > 0",
                     uiOutput("pheno")
                    ),
    
    submitButton(text = "Submit")
  ),
  
  mainPanel(
    conditionalPanel(condition="$('html').hasClass('shiny-busy')",img(src="http://i1215.photobucket.com/albums/cc508/nawel12/loading.gif")),
    plotOutput("view")
  )
))



#