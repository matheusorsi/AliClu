# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:56:45 2018

@author: kisha_000
"""

# Validation of Hierarchical Clustering

import traceback as tb
from scipy.cluster.hierarchy import cut_tree
from fastcluster import linkage
from clustering_scores import sequence_indices_in_clusters, cluster_external_index
from statistics import mean, stdev
import pandas as pd
import numpy as np
import itertools

#Function that performs validation of hierarchical clustering as on the paper
#'On validation of Hierarchical Clustering'
#INPUT:
#       bootstrapSamples (M): number of bootstrap samples
#       df_encoded: dataframe containing the temporal sequences
#       main_results: dataframe with all pairwise alignments
#       main_partition (Z): result of hierarchical clustering on results
#       method: distance metric used on hierarchical clustering of results, this will be used
#               again for hierarchical clustering of the bootstrap samples
#       min_K: minimum number of clusters that we want to analyze
#       max_K: maximum number of clusters that we want to analyze
def validation(bootstrapSamples, df_encoded, main_results, main_partition, method, min_K, max_K):
    try:
        ##############################################################################
        # HOW MANY CLUSTERS?
        ###############################################################################
        # bootstrap method - sampling without replacement

        #dictionary to store all computed indexes for each number of clusters K=min_K,...max_K
        dicio_statistics = {k:{} for k in range(min_K, max_K)}
        for k in range(min_K, max_K):
            dicio_statistics[k]['rand'] = []
            dicio_statistics[k]['adjusted'] = []
            dicio_statistics[k]['FM'] = []
            dicio_statistics[k]['jaccard'] = []
            dicio_statistics[k]['adjusted_wallace'] = []


        #for each bootstrap sample
        for i in range(bootstrapSamples):
            # sampling rows of the original data
            idx = np.random.choice(len(df_encoded), int((3/4)*len(df_encoded)), replace = False)
            idx = np.sort(idx)
            #get all the possible combinations between the sampled patients
            patient_comb_bootstrap = list(itertools.combinations(df_encoded.loc[idx,'id_patient'], 2))
            patient_comb_bootstrap = pd.DataFrame(patient_comb_bootstrap,columns = ['patient1','patient2'])
            #extract the scores regarding the previous sampled combinations to be used in hierarchical clustering
            results_bootstrap = pd.merge(main_results, patient_comb_bootstrap, how='inner', on=['patient1','patient2'])
            # Hierarchical Clustering of the bootstrap sample
            partitionBootstrap = linkage(results_bootstrap['score'], method)

            #for each number of clusters k=min_K,...,max_K
            for k in range(min_K, max_K):
                cluster_assignments_original = cut_tree(main_partition, k)
                cluster_assignments_bootstrap = cut_tree(partitionBootstrap, k)
                #list of clusters for the clustering result with the original data
                cluster_indices_original = sequence_indices_in_clusters(cluster_assignments_original,df_encoded.index.tolist())
                #list of clusters for the clustering result with the bootstrap sample
                cluster_indices_bootstrap = sequence_indices_in_clusters(cluster_assignments_bootstrap,idx)

                #compute 4 different cluster external indexes between the partitions
                computed_indexes = cluster_external_index(cluster_indices_original, cluster_indices_bootstrap)
                #print(computed_indexes)
                dicio_statistics[k]['rand'].append(computed_indexes[0])
                dicio_statistics[k]['adjusted'].append(computed_indexes[1])
                dicio_statistics[k]['FM'].append(computed_indexes[2])
                dicio_statistics[k]['jaccard'].append(computed_indexes[3])
                dicio_statistics[k]['adjusted_wallace'].append(computed_indexes[4])


        ###########################################################################
        #  DECISION ON THE NUMBER OF CLUSTERS
        # The correct number of clusters is the k that yield most maximum average values of
        # clustering indices.
        # Also the k found before needs to have a low value of standard deviation - it has to
        # be the minimum between all k's or a value that is somehow still low compared to others
        ###########################################################################

        #dataframe that stores the clustering indices averages for each k
        df_avgs = pd.DataFrame(index = range(min_K, max_K), columns = ['k','Rand','Adjusted Rand','Fowlkes and Mallows','Jaccard','Adjusted Wallace','k_score_avg'], dtype='float')
        #dataframe that stores the AR and AW indices standard deviations for each k
        df_stds = pd.DataFrame(index = range(min_K, max_K), columns = ['k','Rand','Adjusted Rand','Fowlkes and Mallows','Jaccard','Adjusted Wallace'], dtype = 'float')

        #computing the means and standard deviations
        for k in range(min_K, max_K):
            df_avgs.loc[k]['k'] = k
            df_avgs.loc[k]['Rand'] = mean(dicio_statistics[k]['rand'])
            df_avgs.loc[k]['Adjusted Rand'] = mean(dicio_statistics[k]['adjusted'])
            df_avgs.loc[k]['Fowlkes and Mallows']= mean(dicio_statistics[k]['FM'])
            df_avgs.loc[k]['Jaccard']= mean(dicio_statistics[k]['jaccard'])
            df_avgs.loc[k]['Adjusted Wallace'] = mean(dicio_statistics[k]['adjusted_wallace'])
            df_avgs.loc[k]['k_score_avg'] = 0

            df_stds.loc[k]['k'] = k
            df_stds.loc[k]['Rand'] = stdev(dicio_statistics[k]['rand'])
            df_stds.loc[k]['Adjusted Rand'] = stdev(dicio_statistics[k]['adjusted'])
            df_stds.loc[k]['Fowlkes and Mallows']  =stdev(dicio_statistics[k]['FM'])
            df_stds.loc[k]['Jaccard'] = stdev(dicio_statistics[k]['jaccard'])
            df_stds.loc[k]['Adjusted Wallace'] = stdev(dicio_statistics[k]['adjusted_wallace'])
            #df_stds.loc[k]['k_score_std'] = 0
            #df_stds.loc[k]['k_score_std_2'] = 0

        #weights given to each clustering indice, Rand Index does not value as much as the other indices
        weights = {'Adjusted Rand': 1/4, 'Fowlkes and Mallows': 1/4,
                    'Jaccard':1/4, 'Adjusted Wallace':1/4}
        #found the maximum value for each clustering index and locate in which k it happens
        # compute the scores for each k as being the sum of weights whenever that k has maximums of clustering indices
        columns = df_avgs.columns
        analyzed_columns = columns[2:-1]
        for column in analyzed_columns:
            idx_max = df_avgs[column].idxmax()
            df_avgs.loc[idx_max]['k_score_avg'] = df_avgs.loc[idx_max]['k_score_avg'] + weights[column]

        #final number of clusters chosen by analysing df_avgs
        final_k = df_avgs['k_score_avg'].idxmax()

        return [df_avgs,df_stds,final_k]
    except:
            print(tb.format_exc(), flush=True)
            #k=input("press close to exit") 

#Function to retrieve the final number of cluster
#This function is used after having retrieved all the information when the method was used
#for different variations of the gap penalty
def final_decision(df_final_decision):

    final_k = df_final_decision['k'].value_counts().idxmax()
    df_final_decision = df_final_decision[df_final_decision['k']==final_k]

    df_aux = pd.DataFrame(0,index = range(0,len(df_final_decision)),columns = ['k_score'],dtype = 'float')
    df_final_decision = df_final_decision.reset_index(drop = 'True')
    #weights given to each clustering indice, Rand Index does not value as much as the other indices
    weights = {'Adjusted Rand': 1/4, 'Fowlkes and Mallows': 1/4,
                   'Jaccard':1/4, 'Adjusted Wallace':1/4}
    #found the maximum value for each clustering index and locate in which k it happens
    # compute the scores for each k as being the sum of weights whenever that k has maximums of clustering indices
    for column in df_final_decision.drop(columns = ['k','Rand','k_score_avg','gap']).columns:
        idx_max = df_final_decision[column].idxmax()
        df_aux.loc[idx_max]['k_score'] = df_aux.loc[idx_max]['k_score'] + weights[column]

    #final number of clusters and best results
    final_k_results = df_final_decision.loc[df_aux['k_score'].idxmax()]

    return final_k_results
