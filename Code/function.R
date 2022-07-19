#############################
# Interactive alluvial plot #
#############################

alluvial = function(age, method, size){
  
  #age: follow-up period 
  #method: clustering strategy used
  #size: limit size of clusters that will be displayed
  
  if(method == "Network-based"){
    tab_method = "net"
  }else{
    tab_method = "raw"
  }
  
  #Total of patients belonging to clusters >= 'size' from 60 to 70 years old 
  pat_tot = lapply(age[1]:tail(age,1), function(x){
    tab=read.csv(paste0("Data/clusters_", tab_method, "_", x, ".csv"), sep = ";")
    tab[tab[,2] %in% names(which(table(tab[,2])>=size)),1]})
  pat_tot = as.character(unique(unlist(pat_tot)))
  
  #Patient characteristics
  pat_ch = read.csv("patient_simu.csv", sep = ";", header = TRUE, colClasses = c("character", NA, NA, NA))[, c(1,3,4)] #patient Ids + year of  birth + year of death
  rownames(pat_ch) = pat_ch[,1]
  pat_ch = pat_ch[pat_tot,] #Characteristics of selected patients
  pat_ch[,3] = ifelse(is.na(pat_ch[,3]), Inf, pat_ch[,3]) #replacing NA in year of death by Inf
  
  for(i in age){
    #Clusters identified at age i with the selected 'method'
    clust_tab = read.csv(paste0("Data/clusters_", tab_method, "_", x, ".csv"), sep = ";", colClasses = c("character", NA))
    
    #Keep only clusters >= 'size'
    isolate = names(which(table(clust_tab[,2])<size))
    if(length(isolate) != 0){
      clust_tab = clust_tab[-which(clust_tab[,2] %in% isolate),]
    }
    
    #Ranking clusters from the largest to the smallest
    id_ord = names(sort(table(clust_tab[,2]), decreasing = TRUE))
    clust_tab[,2] = as.character(factor(clust_tab[,2], id_ord, 1:length(id_ord)))
    
    
    #Identifying patients with truncated data
    tr_pat = setdiff(pat_tot, clust_tab[,1]) #patients with truncated data
    tr_pat_ch = pat_ch[tr_pat,] #characteristics of patients with truncated data
    tr_pat_ch$Y = tr_pat_ch[,2]+i #Year they turned i
    
    #Identifying the types of truncated data
    D_pat = tr_pat_ch[tr_pat_ch[,4]>tr_pat_ch[,3],1] #deceased patients
    S_pat = tr_pat_ch[tr_pat_ch[,4]<2008,1] #patients aged 60 after 2008
    E_pat = setdiff(tr_pat_ch[tr_pat_ch[,4]>2018,1], D_pat) #patients aged 60 before 2008
    
    #Adding patients with truncated data in the table of clusters ('clust_tab') 
    if(length(D_pat)!=0 | length(S_pat)!=0 | length(E_pat)!=0){
      
      #Deceased patients
      if(length(D_pat)!=0){
        clust_tab = rbind(clust_tab, data.frame(Patient = D_pat, cluster = "TT3"))
        tr_pat = setdiff(tr_pat, D_pat)
      }
      
      #Patients aged 60 after 2008
      if(length(S_pat)!=0){
        clust_tab = rbind(clust_tab, data.frame(Patient = S_pat, cluster = "TT1"))
        tr_pat = setdiff(tr_pat, S_pat)
      }
      
      #Patients aged 60 before 2008
      if(length(E_pat)!=0){
        clust_tab = rbind(clust_tab, data.frame(Patient = E_pat, cluster = "TT2"))
        tr_pat = setdiff(tr_pat, E_pat)
      }
      
      #Patients with no prescription at age i
      clust_tab = rbind(clust_tab, data.frame(Patient = tr_pat, cluster = "0")) 
      
      
    }else{
      
      #Patients with no prescription at age i
      clust_tab = rbind(clust_tab, data.frame(Patient = tr_pat, cluster = "0"))
      
    }
    
    colnames(clust_tab)[2] = i
    
    if(i == age[1]){
      allu_tab = clust_tab
    }else{
      allu_tab = merge(allu_tab, clust_tab)
    }
  }
  rownames(allu_tab) = allu_tab[,1]
  return(allu_tab[,-1])
}




##################################
# Cluster-trajectories flowchart #
##################################

traj_simu = function(age, traj_clust, size, method){
  
  #age: follow-up period
  #method: clustering strategy used
  #size ["integer"]: limit size of clusters that will be displayed
  
  if(method == "Network-based"){
    tab_method = "B01_per_age/clust_simu_B01_"
  }else{
    tab_method = "B01_per_age/clustKM_simu_B01_"
  }
  
  #patient characteristics
  pat = read.csv("patient_simu.csv", sep = ";", header = TRUE, colClasses = c("character", NA, NA, NA))[, c(1,2)] #we extract Id patient and sex
  rownames(pat) = pat[,1]
  #ald = read.csv("ALD2.csv", header = TRUE, sep = ";") #ALD
  
  clust_tab1 = read.csv(paste0(tab_method, age[1], ".csv"), sep = ";", colClasses = c("character", NA))
  isolate = names(which(table(clust_tab1[,2])<size))
  if(length(isolate) != 0){
    clust_tab1 = clust_tab1[-which(clust_tab1[,2] %in% isolate),]
  }
  
  #We order name cluster from 1, the bigger cluster to the last, the smaller
  id_ord = names(sort(table(clust_tab1[,2]), decreasing = TRUE)) #cluster names ordered by size
  clust_tab1[,2] = as.character(factor(clust_tab1[,2], id_ord, 1:length(id_ord))) #renaming
  
  num_clust1 = paste0(age[1], ".", unique(clust_tab1[,2]))
  ben1 = lapply(unique(clust_tab1[,2]), function(x,tab)tab[tab[,2]==x,1], clust_tab1) #we get all patient id per cluster -> return a list 
  names(ben1) = num_clust1
  
  #Table of prescription
  pres_tab1 = read.csv(paste0("B01_per_age/pres_simu_B01_age_", age[1], ".csv"), sep = ";", check.names = FALSE)
  
  #Sex ratio
  SR1 = unlist(lapply(ben1, function(x){
    return(round(mean(pat[x,2]==1),2))
  }))
  names(SR1) = num_clust1
  
  for(i in age[-1]){
    clust_tab2 = read.csv(paste0(tab_method, i, ".csv"), sep = ";", colClasses = c("character", NA))
    isolate = names(which(table(clust_tab2[,2])<size))
    if(length(isolate) != 0){
      clust_tab2 = clust_tab2[-which(clust_tab2[,2] %in% isolate),]
    }
    
    #We order name cluster from 1, the bigger cluster to the last, the smaller
    id_ord = names(sort(table(clust_tab2[,2]), decreasing = TRUE)) #cluster names ordered by size
    clust_tab2[,2] = as.character(factor(clust_tab2[,2], id_ord, 1:length(id_ord))) #renaming
    
    #We add patient in cluster 0 identified from allu_inter function
    clust_tab2 = rbind(clust_tab2, data.frame(Patient = rownames(traj_clust)[traj_clust[,paste(i)]=="0"], cluster = 0))
    
    num_clust2 = paste0(i, ".", unique(clust_tab2[,2]))
    ben2 = lapply(unique(clust_tab2[,2]), function(x,tab)tab[tab[,2]==x,1], clust_tab2)
    names(ben2) = num_clust2
    
    #Table of prescription
    pres_tab2 = read.csv(paste0("B01_per_age/pres_simu_B01_age_", i, ".csv"), sep = ";", check.names = FALSE)
    
    #Sex ratio
    SR2 = unlist(lapply(ben2, function(x){
      return(round(mean(pat[x,2]==1),2))
    }))
    names(SR2) = num_clust2
    
    #combination of clusters t and t+1 to compute similarity between them
    ben_pair = expand.grid(ben1,ben2)
    
    #combination of name of cluster associated (name cluster = age.numberOfCluster)
    name_pair = expand.grid(num_clust1,num_clust2)
    
    #Jaccard index calculation between consecutive cluster
    name_pair$nbpat = apply(ben_pair, 1, function(x){
      length(intersect(x[[1]], x[[2]]))
    })
    
    #Keep only cluster couples with Jaccard index >0
    name_pair = name_pair[name_pair[,3]>0,]
    name_pair = name_pair %>% group_by(Var1) %>%
      filter(nbpat >= max(nbpat))
    name_pair = name_pair[order(name_pair[,3], decreasing = TRUE),]
    name_pair = as.data.frame(name_pair)
    name_pair[,1] = as.character(name_pair[,1])
    name_pair[,2] = as.character(name_pair[,2])
    
    top_pres = lapply(unique(c(name_pair[,1], name_pair[,2])), function(x){
      if(substr(x, 1, 2)==i){
        if(gsub("^.{3}", "", x)=="0"){
          res = paste0(x, "[", x, "<br> SR=", SR2[x], "<br> n=", nrow(clust_tab2[clust_tab2[,2]=="0",]), "]")
        }else{
          pat_clust = clust_tab2[clust_tab2[,2]==gsub("^.{3}", "", x),1]
          pres_clust = pres_tab2[pat_clust,]
          top = round(sort(colMeans(pres_clust!=0),decreasing = TRUE)[1:2],2)
          res = paste0(x, "[", x, " <br> ", names(top[1]), " : ", top[1], " <br> ", names(top[2]), " : ", top[2], "<br> n=", length(pat_clust), "]")
        }
      }else{
        if(gsub("^.{3}", "", x)=="0"){
          res = paste0(x, "[", x, "<br> SR=", SR1[x], "<br> n=", nrow(clust_tab1[clust_tab1[,2]=="0",]), "]")
        }else{
          pat_clust = clust_tab1[clust_tab1[,2]==gsub("^.{3}", "", x),1]
          pres_clust = pres_tab1[pat_clust,]
          top = round(sort(colMeans(pres_clust!=0),decreasing = TRUE)[1:2],2)
          res = paste0(x, "[", x, " <br> ", names(top[1]), " : ", top[1], " <br> ", names(top[2]), " : ", top[2], "<br> SR=", SR1[x], "<br> n=", length(pat_clust), "]")
        }
      }
      return(res)
    })
    names(top_pres) = unique(c(name_pair[,1], name_pair[,2]))
    
    trj = apply(name_pair, 1, function(x){
      return(paste0(top_pres[[ x[1] ]], " --> |", x[3], "|", top_pres[[ x[2] ]]))
    })
    
    
    ben1 = ben2
    num_clust1 = num_clust2
    clust_tab1 = clust_tab2
    pres_tab1 = pres_tab2
    SR1 = SR2
    
    if(i == age[2]){
      res = trj
      link = name_pair[,3] #edge weights
    }else{
      res = append(res, trj)
      link = append(link, name_pair[,3])
    }
  }
  link = pixedge(link)
  link = paste0("linkStyle ", 0:(length(link)-1), " stroke-width:", link)
  res = paste(c(res,link), collapse = ";")
  return(paste("graph LR; linkStyle default interpolate basis", res, sep = ";"))
}