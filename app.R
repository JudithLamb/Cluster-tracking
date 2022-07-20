#############
# Libraries #
#############

if(!requireNamespace("easyalluvial", quietly = TRUE)) install.packages("easyalluvial")
if(!requireNamespace("DiagrammeR", quietly = TRUE)) install.packages("DiagrammeR")
if(!requireNamespace("parcats", quietly = TRUE)) install.packages("parcats")
if(!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(easyalluvial)
library(DiagrammeR)
library(parcats)
library(dplyr)

##################
# Functions file #
##################

source("Code/functions.R")

#####################################
# Arrow thickness of the flowchart #
#####################################

pixlink <- function(x){
  #x ["integer vector"]: number of common patients associated to all arrows in the flowchart
  
  a <- min(x) 
  b <- max(x)
  
  #Conversion of the number of common patients in pixels
  pix <- lapply(x, function(x){
    round(( ((x-a) * (5-0.4)) / (b-a) ) + 0.4, 2)
  })
  return(unlist(pix))
}

###########################
# User Interface of Shiny #
###########################

ui <- fluidPage(
  
  # App title
  titlePanel("Cluster-tracking approach from 60 to 70 years old"),
  
  # Sidebar layout
  sidebarPanel(
    #Method selection
    selectInput("method", "Clustering strategy", c("Network-based", "Raw-data-based"), "Network-based"),
    
    #cluster size selection
    sliderInput("size", "Cluster limit size", 0, 100 , 30, 10),
    
    hr(),
    
    #Display button
    actionButton("display", "Display"),
    
    #Clear button
    actionButton("reset", "Clear"),
    
    width = 2
  ),
  
  # Main panel
  tabsetPanel(type = "pills",
              tabPanel("Alluvial Plot", parcatsOutput("Alluvial", height = "900px")),
              tabPanel("Cluster Trajectories", DiagrammeROutput(outputId = "top2", width = "2500px", height = "2500px"))
  )
)



############################
# Server function of Shiny #
############################

server <- function(input, output, session) {
  
  value <- reactiveValues()
  
  observeEvent(input$display, {
    allu_tab <- alluvial(60:70, input$method, as.numeric(input$size))
    value$allu <- parcats(alluvial_wide(allu_tab), marginal_histograms = FALSE, hoverinfo = "count")
    value$traj <- flowchart(60:70, input$method, as.numeric(input$size))
  })
  
  observeEvent(input$reset, {
    session$reload()
  })
  
  output$Alluvial <- render_parcats(value$allu)
  
  output$top2 <- renderDiagrammeR(mermaid(value$traj))
  
}

shinyApp(ui, server)