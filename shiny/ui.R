library(shiny)


shinyUI(pageWithSidebar(
  headerPanel("Loading data"),
  sidebarPanel(
    selectInput("dataset", "Choose a dataset:", 
                choices =gsub("\\.RData","",dir("/home/guegan/Documents/ICM/Dev/shiny/test_vincent/",pattern="*.RData"))),
    textInput("gene","Gene"),    
    submitButton("Submit")

  ),
  mainPanel(
    plotOutput("view")
  )
))



#