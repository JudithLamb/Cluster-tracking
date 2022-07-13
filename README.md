# Cluster-tracking approach
We propose novel approaches based on cluster-tracking for clustering patients from longitudinal data extracted from medico-administrative databases. These approaches start by identifying clusters of patients at each considered age. To this goal, we used two different clustering strategies: Markov Cluster algorithm (MCL) applied on patient networks built from raw data and Kmeans applied directly on raw data. Clusters are then tracked over ages to define cluster-trajectories. We applied our approaches to the analysis of antithrombotic drug prescriptions extracted from the Echantillon Généraliste des Bénéficiaires (EGB, a French cohort) between 2008 and 2018 in patients aged from 60 to 70 years old. From this data, we simulated 5594 patients with their drug prescriptions. This simulated sample is used in the following to apply our two cluster-tracking approaches.

## Identifying clusters of patients from patient networks
We started by constructing a patient network for each age considered. We then applied the MCL clustering algorithm on each network.

### Constructing patient networks
A patient network is a graph $G = (V,E)$ with $V$ patient nodes and $E$ edges representing interactions between patient nodes. We built a network for each patient age. Each network is constructed using a similarity matrix. In this similarity matrix, we computed the similarity between patients of the same age using the Cosine similarity.

```python
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity

# Construction of similarity matrices from 60 to 70 years old
for i in np.arange(60,71):
    #Table of prescriptions at age i
    pres_tab = pd.read_csv(("Data/pres_%d.csv" %(i)), sep = ";")
    
    #Computation of the Cosine similarity between patients 
    cos_DF = pd.DataFrame(cosine_similarity(pres_tab), columns = pres_tab.index.astype("str"),
                          index = pres_tab.index.astype("str")) 
    
    #Saving
    cos_DF.to_parquet(("Data/cosine_%d.gzip" %(i)), compression="gzip") 
```

We then filtered the similaritry matrices according to a threshold. This threshold is chosen over all the matrices in order to reduce the number of edges in the networks while obtaining a minimum number of isolated patients. From the matrices constructed with the simulated sample, we choose a threshold of 0.7 because this is where we observe the fastest decrease in the number of edges and there is only a small number of isolated patients (see figure below). 

![example visualization](Figure/cosine_threshold.png)

From each filtered matrix, we obtain a patient network in which patients are connected only if they have a Cosine similarity $\ge 0.7$. The patient network at 60 years of age is represented in the figure below.

![example visualization](Figure/network_60.png)

### Identifying clusters of patients from raw data

