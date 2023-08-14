#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script Name: select_clus_rep.py

Description:
The code reads in two input files, vp_clusters_cluster.tsv(clusters' reps and clusters' members in 
 tabular format) and cluster_sizes.csv,and performs several operations on them. The first file
 contains information about protein clusters produced by MMSeq, including a cluster 
 representative and cluster members. The second file contains information about 
 cluster sizes (the number of members in each cluster). 
 The code then generates statistics about the clusters and plots 
 a histogram of the data. It also replaces some of the cluster representatives 
 produced by MMSeq with new ones from PDB or AFDB, and outputs the top 
 clusters into a tsv file.

Procedure:

1- Import the pandas library.
2- Read the vp_clusters_cluster.tsv file into a dataframe.
3- Read the cluster_sizes.csv file into a dataframe df.
4- Generate statistics about the clusters using the describe() method and 
   plot a histogram of the data using the plot() method.
5- Sort the dataframe by the number of members in each cluster.
6- Get the top 10,000 cluster representatives from the dataframe.
7- Calculate the percentage of proteins covered in the top clusters.
8- Create a new dataframe newClusterReps for the new cluster representatives.
9- Iterate through the top cluster representatives and check if they have a
   refseq representative.
10- If a refseq representative is found, get all the members of the cluster and 
   check if they have any PDB or AFDB members.
11- Once a PDB or AFDB member is found, make it the new cluster representative 
   and add a row to the newClusterReps dataframe.
12- Augment the topClustersReps dataframe with a new column showing the final 
   cluster representatives that will be used for folding.
13- Iterate through the newClusterReps dataframe and assign the new cluster 
   representative to the corresponding top cluster.
14- Write the topClustersReps dataframe into a tsv file.

Inputs:

1- vp_clusters_cluster.tsv: A tab-separated file containing information about 
   protein clusters produced by MMSeq, including a cluster representative and cluster members.
cluster_sizes.csv: A comma-separated file containing information about cluster 
   sizes (#of members).
          
Outputs:

1- A histogram of the clusters statistics.
2- top_viral_protein_cluster_reps.tsv: A tab-separated file containing the 
   top 10,000 cluster representatives and their associated size.

Description Author: ChatGPT, reviewed and corrected by Tamim AlMurad.
Code author: Tamim AlMurad

Version 1.0
last review Date: July 11th 2023
"""

import pandas as pd
import sys

if len(sys.argv) != 3:
    print("Usage: python process_clusters.py <clusters_file> <cluster_size_file>")
    exit()
    
    
clustersFile = sys.argv[1]
clusterSizeFile = sys.argv[2]

#Read clusters file into a dataframe. File is produced by MMSeq with cluster rep and cluster members.
clustersDF = pd.read_table(clustersFile)
clustersDF.columns=['Clus_Rep','Clus_Mem']

#Read clusters stats file with cluster rep and the number of members in the cluster.
df = pd.read_table(clusterSizeFile)
df.columns=['ACCN', '#of members']

#Get statistics of the clusters and plot a histogram of the data.
df.plot(kind='hist', xlabel='Cluster Size',logx=True, logy=True, bins=5000,title='Histogram of Cluster Sizes (log)').figure.savefig('cluster_sizes_Hist.png')
stat=df.describe()

#Sort the clusters by the size.
df.sort_values(by=['#of members'],inplace=True,ascending=False)



#Get the top 10,000 clusters reps.
topClustersReps = df.iloc[0:10000]

#%%
#Percentage of proteins covered in the top clusters. 
print('10k clusters coverage of the whole viral proteins is %f' 
      %(topClustersReps['#of members'].sum()/len(clustersDF)*100))
#%%
#This section is to replace cluster reps produced by MMSeq to be reps from PDB ot AFDB.
#Dataframe for the new cluster reps.
newClusterReps = pd.DataFrame(columns=['old_clus_rep', 'new_clus_rep'])



for rep in topClustersReps['ACCN']:
    #Check if the cluster has a refseq rep.
    if rep.find('.') > 0:
        
        #Get all the members of the cluster.
        tempDF = clustersDF[clustersDF['Clus_Rep'] == rep]
      
        for mem in tempDF['Clus_Mem']:
            #Check if the cluster has any AFDB or PDB members.
            if mem.find('.') < 0:
                #Once a PDF or AFDB member is found, make it the cluster rep and break.
                row = {'old_clus_rep': rep,'new_clus_rep': mem}
                rowDF=pd.DataFrame.from_dict(row,orient='index').T
                print("cluster rep is changed \n " + rowDF.to_string())
                #Add a row showing old cluster rep and new cluster rep.
                newClusterReps = pd.concat([newClusterReps, rowDF],ignore_index=True)       
                break
                        


#%%
#In this section, the df with top cluster reps is augmented with a new column showing the 
# the final cluster reps that will be used.
topClustersReps['final_clus_rep']=topClustersReps['ACCN']
topClustersReps.reset_index()

i =0
for rep in newClusterReps['old_clus_rep']:
    #Get the index of the top clusters with new cluster rep.
    index = topClustersReps[topClustersReps['ACCN']==rep].index.values
    #Assign the new cluster rep to a top cluster.
    topClustersReps.loc[index,'final_clus_rep'] = newClusterReps.at[i,'new_clus_rep']
    i=i+1

print('###Number of clusters rep to be folded is: '+str(len(topClustersReps[topClustersReps['final_clus_rep'].str.contains('\.')])))
#Write the top clusters dataframe into a tsv file.
topClustersReps.to_csv('top_viral_protein_cluster_reps_.tsv',sep='\t')
