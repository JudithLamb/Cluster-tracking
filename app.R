library(shiny)
library(easyalluvial)
library(DiagrammeR)
library(parcats)

#User Interface of Shiny
ui <- fluidPage(
  
  # App title
  titlePanel("Cluster Tracking"),
  
  # Sidebar layout
  sidebarPanel(
    #Selection of ange range
    h5(HTML("<b>Age range</b>")),
    
    splitLayout(cellWidths = c("50%","50%"),
                selectInput("age1", NULL, c("from",seq(40,60,10)), width = "200px"),
                selectInput("age2", NULL, c("to",seq(50,70,10)), width = "200px")
                ),
    
    #selection of cluster size
    selectInput("size", "Cluster size", c(100), 100),
    
    hr(),
    
    actionButton("display", "Display"), #display button
    
    actionButton("reset", "Clear"), #clear button
    
    width = 2
  ),
  
  # Main panel
  tabsetPanel(type = "pills",
              tabPanel("Alluvial Plot", parcatsOutput("Alluvial", height = "900px")), #alluvial tab
              tabPanel("Cluster Trajectories", DiagrammeROutput(outputId = "Trajectory", width = "2500px", height = "2500px")) #trajectory tab
  )
)

#Server function of Shiny
server <- function(input, output, session) {
  
  value <- reactiveValues()
  
  #Update of age range upper limit
  observeEvent(input$age1, {if(input$age1!="from"){
    updateSelectInput(session, "age1", choices = c(seq(40,60,10)) , selected = input$age1)
    updateSelectInput(session, "age2", choices = c(seq(50,70,10)) , selected = as.numeric(input$age1)+10)
  }
  })
  
  #Update of age range lower limit
  observeEvent(input$age2, {if(input$age2!="to"){
    updateSelectInput(session, "age2", choices = c(seq(50,70,10)) , selected = input$age2)
    updateSelectInput(session, "age1", choices = c(seq(40,60,10)) , selected = as.numeric(input$age2)-10)
  }
  })
  
  #Action of display button
  observeEvent(input$display, {if(input$age1=="from" | input$age2=="to"){
    showNotification("Incorrect age range selected", duration = NULL, id = "ntf", type = "error")
  }else{
    removeNotification("ntf")
    all_data <- read.csv(paste0("Data/alluvial_", input$age1, ".", input$age2, "_t" , input$size, ".csv"), sep = ";", header = TRUE, check.names = FALSE)
    value$all <- parcats(alluvial_wide(all_data), marginal_histograms = FALSE, hoverinfo = "count")
    value$traj <- read.csv(paste0("Data/traj_clust_", input$age1, ".", input$age2, "_t" , input$size, ".txt"), header = FALSE)
  }
  })
  
  #Action of clear button
  observeEvent(input$reset, {
    session$reload()
  })
  
  #Output in alluvial tab
  output$Alluvial <- render_parcats(value$all)
  
  #Output in trajectory tab
  output$Trajectory <- renderDiagrammeR(mermaid(value$traj[1,]))
}

shinyApp(ui, server)