# Explanation of data files

## patient_characteristics.csv
The file contains the patients IDs, the patient sex, the years of birth and the years of death.

## pres_60.csv
The file contains the total number of prescriptions per drug that patients have had at age 60. This information is available from age 60 to 70 in separate files.

## cosine_60.gzip
The file contains the similarity matrix at age 60 constructed by computing the Cosine similarity between all the 60-year-old patients in the network-based cluster-tracking approach. The similarity matrices are available from age 60 to 70 in separate files.

## clusters_net_60.csv
The file contains the clusters identified with MCL applied on the 60-year-old patient network. Each patient is associated to their cluster label. Clusters identified at each age (from 60 to 70) in the network-based cluster-tracking approach are available in separate files.

## clusters_raw_60.csv
The file contains the clusters identified with Kmeans applied on the drug prescriptions of 60-year-old patients. Each patient is associated to their cluster label. Clusters identified at each age (from 60 to 70) in the raw-data-based cluster-tracking approach are available in separate files.
