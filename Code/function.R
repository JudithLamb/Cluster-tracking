#############################
# Interactive alluvial plot #
#############################

alluvial = function(method, size){
  
  #method: clustering strategy used
  #size: limit size of clusters that will be displayed
  
  if(method == "Network-based"){
    tab_method = "B01_per_age/clust_simu_B01_"
  }else{
    tab_method = "B01_per_age/clustKM_simu_B01_"
  }
  
  #Total of patients belonging to clusters greater than or equal to 'size' 
  pat_tot = lapply(60:70, function(x){
    tab=read.csv(paste0(tab_method, x, ".csv"), sep = ";")
    tab[tab[,2] %in% names(which(table(tab[,2])>=size)),1]})
  pat_tot = as.character(unique(unlist(pat_tot)))
  
  #patient characteristics
  pat = read.csv("patient_simu.csv", sep = ";", header = TRUE, colClasses = c("character", NA, NA, NA))[, c(1,3,4)] #Patient Ids + year of  birth + year of death
  rownames(pat) = pat[,1]
  pat = pat[pat_tot,] #Characterics of selected patients
  
  #we replace NA in year of death by Inf
  pat[,3] = ifelse(is.na(pat[,3]), Inf, pat[,3])
  
  for(i in age){
    clust_tab = read.csv(paste0(tab_method, i, ".csv"), sep = ";", colClasses = c("character", NA))
    isolate = names(which(table(clust_tab[,2])<size))
    if(length(isolate) != 0){
      clust_tab = clust_tab[-which(clust_tab[,2] %in% isolate),]
    }
    
    #We order name cluster from 1, the bigger cluster to the last, the smaller
    id_ord = names(sort(table(clust_tab[,2]), decreasing = TRUE)) #cluster names ordered by size
    clust_tab[,2] = as.character(factor(clust_tab[,2], id_ord, 1:length(id_ord))) #renaming
    
    #We add patients clustered during the period but not during this given age = they belongs to the "0" cluster if they haven't prescription and "D" cluster if they are died
    ex_ben = setdiff(ben_tot, clust_tab[,1])
    
    pat_clust = pat[ex_ben,]
    pat_clust$Y = pat_clust[,2]+i #Year of a i-th birthday
    D_pat = pat_clust[pat_clust[,4]>pat_clust[,3],1] #patients died
    S_pat = pat_clust[pat_clust[,4]<2008,1] #patients before the start of the follow-up
    E_pat = setdiff(pat_clust[pat_clust[,4]>2018,1], D_pat) #patients at the end of the follow-up and not died
    
    if(length(D_pat)!=0 | length(S_pat)!=0 | length(E_pat)!=0){
      
      
      if(length(D_pat)!=0){
        clust_tab = rbind(clust_tab, data.frame(Patient = D_pat, cluster = "TT3"))
        ex_ben = setdiff(ex_ben, D_pat)
      }
      
      if(length(S_pat)!=0){
        clust_tab = rbind(clust_tab, data.frame(Patient = S_pat, cluster = "TT1"))
        ex_ben = setdiff(ex_ben, S_pat)
      }
      
      if(length(E_pat)!=0){
        clust_tab = rbind(clust_tab, data.frame(Patient = E_pat, cluster = "TT2"))
        ex_ben = setdiff(ex_ben, E_pat)
      }
      
      clust_tab = rbind(clust_tab, data.frame(Patient = ex_ben, cluster = "0"))
      
      
    }else{
      
      clust_tab = rbind(clust_tab, data.frame(Patient = ex_ben, cluster = "0"))
      
    }
    
    colnames(clust_tab)[2] = i #New name of cluster column
    
    if(i == age[1]){
      res = clust_tab
    }else{
      res = merge(res, clust_tab)
    }
  }
  rownames(res) = res[,1]
  return(res[,-1])
}




##################################
# Cluster-trajectories flowchart #
##################################

traj_simu = function(age, traj_clust, size, method){
  
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