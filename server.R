#############
# Libraries #
#############

if(!requireNamespace("easyalluvial", quietly = TRUE)) install.packages("easyalluvial")
if(!requireNamespace("DiagrammeR", quietly = TRUE)) install.packages("DiagrammeR")
if(!requireNamespace("parcats", quietly = TRUE)) install.packages("parcats")

library(easyalluvial)
library(DiagrammeR)
library(parcats)

#############
# Functions #
#############


# Interactive alluvial plot
###########################

alluvial <- function(age, method, size){
  
  #age ["integer vector"]: follow-up period 
  #method ["character"]: clustering strategy used
  #size ["integer"]: limit size of clusters that will be displayed
  
  #selected method
  if(method == "Network-based"){
    tab_method <- "net"
  }else{
    tab_method <- "raw"
  }
  
  #Total of patients belonging to clusters >= 'size' from 60 to 70 years old 
  pat_tot <- lapply(age[1]:tail(age,1), function(x){
    tab <- read.csv(paste0("Data/clusters_", tab_method, "_", x, ".csv"), sep = ";")
    tab[tab[,2] %in% names(which(table(tab[,2])>=size)),1]})
  pat_tot <- as.character(unique(unlist(pat_tot)))
  
  #Patient characteristics
  pat_ch <- read.csv("Data/patient_characteristics.csv", sep = ";", header = TRUE, colClasses = c("character", NA, NA, NA))[, c(1,3,4)] #patient Ids + year of  birth + year of death
  rownames(pat_ch) <- pat_ch[,1]
  pat_ch <- pat_ch[pat_tot,] #Characteristics of selected patients
  pat_ch[,3] <- ifelse(is.na(pat_ch[,3]), Inf, pat_ch[,3]) #replacing NA in year of death by Inf
  
  for(i in age){
    #Clusters identified at age i with the selected 'method'
    clust_tab <- read.csv(paste0("Data/clusters_", tab_method, "_", x, ".csv"), sep = ";", colClasses = c("character", NA))
    
    #Keep only clusters >= 'size'
    isolate <- names(which(table(clust_tab[,2])<size))
    if(length(isolate) != 0){
      clust_tab <- clust_tab[-which(clust_tab[,2] %in% isolate),]
    }
    
    #Ranking clusters from the largest to the smallest
    id_ord <- names(sort(table(clust_tab[,2]), decreasing = TRUE))
    clust_tab[,2] <- as.character(factor(clust_tab[,2], id_ord, 1:length(id_ord)))
    
    
    #Identifying patients with truncated data
    tr_pat <- setdiff(pat_tot, clust_tab[,1]) #patients with truncated data
    tr_pat_ch <- pat_ch[tr_pat,] #characteristics of patients with truncated data
    tr_pat_ch$Y <- tr_pat_ch[,2]+i #Year they turned i
    
    #Identifying the types of truncated data
    D_pat <- tr_pat_ch[tr_pat_ch[,4]>tr_pat_ch[,3],1] #deceased patients
    S_pat <- tr_pat_ch[tr_pat_ch[,4]<2008,1] #patients aged 60 after 2008
    E_pat <- setdiff(tr_pat_ch[tr_pat_ch[,4]>2018,1], D_pat) #patients aged 60 before 2008
    
    #Adding patients with truncated data in the table of clusters ('clust_tab') 
    if(length(D_pat)!=0 | length(S_pat)!=0 | length(E_pat)!=0){
      
      #Deceased patients
      if(length(D_pat)!=0){
        clust_tab <- rbind(clust_tab, data.frame(Patient = D_pat, cluster = "TT3"))
        tr_pat <- setdiff(tr_pat, D_pat)
      }
      
      #Patients aged 60 after 2008
      if(length(S_pat)!=0){
        clust_tab <- rbind(clust_tab, data.frame(Patient = S_pat, cluster = "TT1"))
        tr_pat <- setdiff(tr_pat, S_pat)
      }
      
      #Patients aged 60 before 2008
      if(length(E_pat)!=0){
        clust_tab <- rbind(clust_tab, data.frame(Patient = E_pat, cluster = "TT2"))
        tr_pat <- setdiff(tr_pat, E_pat)
      }
      
      #Patients without prescriptions at age i
      clust_tab <- rbind(clust_tab, data.frame(Patient = tr_pat, cluster = "0")) 
      
      
    }else{
      
      #Patients without prescriptions at age i
      clust_tab <- rbind(clust_tab, data.frame(Patient = tr_pat, cluster = "0"))
      
    }
    
    colnames(clust_tab)[2] = i
    
    if(i == age[1]){
      allu_tab <- clust_tab
    }else{
      allu_tab <- merge(allu_tab, clust_tab)
    }
  }
  rownames(allu_tab) <- allu_tab[,1]
  return(allu_tab[,-1])
}



# Cluster-trajectories flowchart 
#################################

flowchart <- function(age, method, size){
  
  #age ["integer vector"]: follow-up period
  #method ["character"]: clustering strategy used
  #size ["integer"]: limit size of clusters that will be displayed
  
  #Table obtained from the alluvial function
  allu_tab <- allu_inter_simu(age, size, method)
  
  #Patient characteristics
  pat_ch <- read.csv("Data/patient_characteristics.csv", sep = ";", header = TRUE, colClasses = c("character", NA, NA, NA))[, c(1,2)] #patient ids + sex
  rownames(pat_ch) <- pat_ch[,1]
  
  for(i in head(age, -1)){
    clust_label <- unique(allu_tab[,paste(i)]) #label of clusters at age i
    clust_label <- clust_label[-which(clust_label%in%c("TT1", "TT2", "TT3"))] #removing type of truncated data
    pat_clust1 <- lapply(clust_label, function(x,tab)rownames(tab)[which(tab[,paste(i)]%in%x)], allu_tab) #getting all patient ids per cluster at age i
    clust_label <- paste0(i, ".", clust_label) #adding age before the cluster label ("age.label")
    names(pat_clust1) <- clust_label
    
    clust_label <- unique(allu_tab[,paste(i+1)]) #label of clusters at age i+1
    clust_label <- clust_label[-which(clust_label%in%c("TT1", "TT2", "TT3"))] #removing type of truncated data
    pat_clust2 <- lapply(clust_label, function(x,tab)rownames(tab)[which(tab[,paste(i+1)]%in%x)], allu_tab) #getting all patient ids per cluster at age i+1
    clust_label <- paste0(i+1, ".", clust_label) #adding age before the cluster label ("age.label")
    names(pat_clust2) <- clust_label
    
    #Tables of prescriptions
    pres_tab1 <- read.csv(paste0("Data/pres_", i, ".csv"), sep = ";", check.names = FALSE) #at age i
    pres_tab2 <- read.csv(paste0("Data/pres_", i+1, ".csv"), sep = ";", check.names = FALSE) #at age i+1
    
    #Combinations of clusters at age i and i+1
    pat_combi <- expand.grid(pat_clust1,pat_clust2)
    
    #Name of clusters associated to combinations of clusters
    name_combi <- expand.grid(names(pat_clust1),names(pat_clust2))
    
    #Computing the number of common patients between all combinations of clusters at age i and i+1
    name_combi$nbpat <- apply(pat_combi, 1, function(x){
      length(intersect(x[[1]], x[[2]]))
    })
    name_combi <- name_combi[name_combi[,3]>0,] #number of common patients > 0
    
    #For each cluster at age i, we keep the cluster at age i+1 with which it has the greatest number of common patients
    name_combi <- name_combi %>% group_by(Var1) %>%
      filter(nbpat >= max(nbpat))
    name_combi <- name_combi[order(name_combi[,3], decreasing = TRUE),] #ordering table from the greatest to the smallest number of common patients 
    name_combi <- as.data.frame(name_combi)
    name_combi[,1] <- as.character(name_combi[,1])
    name_combi[,2] <- as.character(name_combi[,2])
    
    #Sex ratio in each cluster of age i
    SR1 <- unlist(lapply(pat_clust1, function(x){
      return(round(mean(pat_ch[x,2]==1),2))
    }))
    
    #Sex ratio in each cluster of age i+1
    SR2 <- unlist(lapply(pat_clust2, function(x){
      return(round(mean(pat_ch[x,2]==1),2))
    }))
    
    #Characterizing each clusters at age i and i+1 with the two most prescribed drugs, the sex ratio and  the total number of patients 
    cluster_cha <- lapply(unique(c(name_combi[,1], name_combi[,2])), function(x){
      if(substr(x, 1, 2)==i){
        
        #characterizing clusters at age i
        if(gsub("^.{3}", "", x)=="0"){
          #cluster containing patients without prescriptions
          res <- paste0(x, "[", x, " <br> SR=", SR1[x], " <br> n=", length(pat_clust1[[x]]), "]")
        }else{
          pres_clust <- pres_tab1[pat_clust1[[x]],] #prescriptions of patients belonging to the cluster x
          top_pres <- round(sort(colMeans(pres_clust!=0),decreasing = TRUE)[1:2],2) #the two most prescribed drugs in cluster x
          res <- paste0(x, "[", x, " <br> ", names(top_pres[1]), " : ", top_pres[1], " <br> ", names(top_pres[2]), " : ", top_pres[2], " <br> SR=", SR1[x], " <br> n=", length(pat_clust1[[x]]), "]")
        }
        
      }else{
        
        #characterizing clusters at age i+1
        if(gsub("^.{3}", "", x)=="0"){
          #cluster containing patients without prescriptions
          res <- paste0(x, "[", x, " <br> SR=", SR2[x], " <br> n=", length(pat_clust2[[x]]), "]")
        }else{
          pres_clust <- pres_tab2[pat_clust2[[x]],] #prescriptions of patients belonging to the cluster x
          top_pres <- round(sort(colMeans(pres_clust!=0),decreasing = TRUE)[1:2],2) #the two most prescribed drugs in cluster x
          res <- paste0(x, "[", x, " <br> ", names(top_pres[1]), " : ", top_pres[1], " <br> ", names(top_pres[2]), " : ", top_pres[2], " <br> SR=", SR2[x], " <br> n=", length(pat_clust2[[x]]), "]")
        }
        
      }
      return(res)
    })
    names(cluster_cha) <- unique(c(name_combi[,1], name_combi[,2]))
    
    #cluster-trajectories
    traj <- apply(name_combi, 1, function(x){
      return(paste0(cluster_cha[[ x[1] ]], " --> |", x[3], "|", cluster_cha[[ x[2] ]]))
    })
    
    if(i == age[1]){
      clust_traj <- traj
      link <- name_combi[,3] #defining arrows in flowchart by the number of common patients
    }else{
      clust_traj <- append(clust_traj, traj)
      link <- append(link, name_combi[,3])
    }
  }
  
  #Arrow thickness
  link <- pixlink(link)
  link <- paste0("linkStyle ", 0:(length(link)-1), " stroke-width:", link)
  
  clust_traj <- paste(c(clust_traj,link), collapse = ";")
  
  return(paste("graph LR; linkStyle default interpolate basis", clust_traj, sep = ";"))
}

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