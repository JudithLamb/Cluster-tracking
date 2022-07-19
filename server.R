#############
# Libraries #
#############

if(!requireNamespace("easyalluvial", quietly = TRUE)) install.packages("easyalluvial")
if(!requireNamespace("DiagrammeR", quietly = TRUE)) install.packages("DiagrammeR")
if(!requireNamespace("parcats", quietly = TRUE)) install.packages("parcats")

library(easyalluvial)
library(DiagrammeR)
library(parcats)



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
    sliderInput("size", "Limit cluster size", 0, 100 , 30, 10),
    
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