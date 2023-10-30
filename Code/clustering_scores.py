# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 11:26:07 2018

@author: kisha_000
"""


import numpy as np
from statistics import mean
import math

def sequence_indices_in_clusters(cluster_assignments,idx):
    n = cluster_assignments.max()
    clusters = []
    for cluster_number in range(0, n + 1):
        aux = np.where(cluster_assignments == cluster_number)[0].tolist()
        cluster = list(idx[i] for i in aux )
        clusters.append(cluster)
    return clusters

def cluster_external_index(partition_a, partition_b):
    #size of contigency table
    R = len(partition_a)
    C = len(partition_b)
    #contigency table
    ct = np.zeros((R+1,C+1))
    #fill the contigency table
    for i in range(0,R+1):
        for j in range(0,C):
            if(i in range(0,R)):  
                n_common_elements = len(set(partition_a[i]).intersection(partition_b[j]))
                ct[i][j] = n_common_elements
            else:
                ct[i][j] = ct[:,j].sum()
                      
        ct[i][j+1] = ct[i].sum()  
    
    N = ct[R][C]
    #condensed information of ct into a mismatch matrix (pairwise agreement)
    sum_all_squared = np.sum(ct[0:R][:,range(0,C)]**2)   
    sum_R_squared = np.sum(ct[0:R,C]**2)
    sum_R = np.sum(ct[0:R,C])
    sum_C_squared = np.sum(ct[R,0:C]**2)
    sum_C = np.sum(ct[R,0:C])
    #computing the number of pairs that are in the same cluster both in partition A and partition B
    a = 0
    for i in range(0,R):
        for j in range(0,C):
            a = a + ct[i][j]*(ct[i][j]-1)
    a = a/2
    #computing the number of pair in the same cluster in partition A but in different cluster in partition B
    b = (sum_R_squared- sum_all_squared)/2
    #computing the number of pair in different cluster in partition A but in the same cluster in partition B
    c = (sum_C_squared - sum_all_squared)/2
    #computing the number of pairs in different cluster both in partition A and partition B
    d = (N**2 + sum_all_squared - (sum_R_squared + sum_C_squared))/2
    
    #Rand Index
    rand_index = (a+d)/(a+b+c+d)

    #Adjusted Rand Index
    nc = ((sum_R_squared - sum_R)*(sum_C_squared -sum_C))/(2*N*(N-1))
    nd = (sum_R_squared - sum_R + sum_C_squared - sum_C)/4
    if(nd==nc):
        adjusted_rand_index = 0
    else:      
        adjusted_rand_index = (a-nc)/(nd - nc)
   
    #Fowlks and Mallows
    if((a+b)==0 or (a+c)==0):
        FM = 0
    else:     
        FM = a/math.sqrt((a+b)*(a+c))
    
    #Jaccard
    if(a+b+c == 0):
        jaccard = 1
    else:
        jaccard = a/(a+b+c)
        
    #Adjusted Wallace
    if((a+b)==0):
        wallace = 0
    else:
        wallace = a/(a+b)
    SID_B = 1-((sum_C_squared-sum_C)/(N*(N-1)))
    if((SID_B)==0):
        adjusted_wallace = 0
    else:
        adjusted_wallace = (wallace-(1-SID_B))/(1-(1-SID_B))

    return [rand_index, adjusted_rand_index, FM, jaccard, adjusted_wallace]


def cluster_internal_indices(partition, distance_matrix, patients_id):
    #print(partition)
    #print(patients_id)

    sum_within_cluster = 0
    sum_between_cluster = 0
    total_distances_within = 0
    term = []
    max_cluster_diameter = 0
    for cluster in partition:
        #print(cluster)
        dist_within_cluster = []
        dist_between_cluster = []
        if len(cluster) != 1:
            for pat1, pat2 in zip(distance_matrix['patient1'], distance_matrix['patient2']):
                idxPat1 = patients_id.index[patients_id == pat1].tolist()[0]
                idxPat2 = patients_id.index[patients_id == pat2].tolist()[0]
                
                if idxPat1 not in cluster and idxPat2 not in cluster:
                    continue

                pair_distance = distance_matrix.loc[(distance_matrix['patient1'] == pat1) &
                                                                (distance_matrix['patient2'] == pat2)]['score']
                if pair_distance is None or len(pair_distance) == 0:
                    pair_distance = distance_matrix.loc[(distance_matrix['patient1'] == pat2) &
                                                                (distance_matrix['patient2'] == pat1)]['score']
                pair_distance = float(pair_distance)
                # print(pair_distance)
                # print('pats', pat1, pat2)
                # print('idxpats', idxPat1, idxPat2)
                
                if all(dupla in cluster for dupla in [idxPat1, idxPat2]):
                    #print("both")
                    dist_within_cluster.append(pair_distance)
                if idxPat1 in cluster and idxPat2 not in cluster:
                    #print("soh primeiro")
                    dist_between_cluster.append(pair_distance)
                if idxPat1 not in cluster and idxPat2 in cluster:
                    #print("soh segundo")
                    dist_between_cluster.append(pair_distance)
            if(len(dist_within_cluster) == 0):
                print(cluster)
                print(distance_matrix)
            number_distances_within = int(len(cluster) * (len(cluster) - 1)/2)
            total_distances_within += number_distances_within
            number_distances_between = len(cluster) * (len(patients_id) - len(cluster))
            #print('vamos verificar', number_distances_between, len(dist_between_cluster), number_distances_within, len(dist_within_cluster))
            sum_between_cluster += sum(dist_between_cluster)
            sum_within_cluster += sum(dist_within_cluster)
            if max(dist_within_cluster) > max_cluster_diameter:
                max_cluster_diameter = max(dist_within_cluster)
            #print("c",cluster, number_distances_within, dist_between_cluster, sum_between_cluster, flush=True)
            term.append((number_distances_between * sum(dist_within_cluster))/(number_distances_within * sum(dist_between_cluster)))

    #McClain
    total_distances = len(distance_matrix)
    number_distances_between = total_distances - total_distances_within
    if(total_distances_within == 0 or sum_between_cluster == 0):
        mcclainN = 10
    else:
        mcclainN = (number_distances_between * sum_within_cluster)/(total_distances_within * sum_between_cluster)
    #print("mcclain novo", mcclainN)

    mcclainV = float(sum(term)/len(partition))
    #print("mcclain velho", mcclainV)

    #C-index
    smin = sum(distance_matrix['score'].sort_values().head(total_distances_within))
    smax = sum(distance_matrix['score'].sort_values(ascending = False).head(total_distances_within))
    if smax == smin:
        cindex = 1
    else:
        cindex = round((sum_within_cluster - smin)/(smax - smin),10)
    if cindex < 0:
        print("cindex", cindex, smin, smax, sum_within_cluster, distance_matrix['score'].sort_values())
    #print("cindex", cindex)

    s_clusters = []
    s_general = []
    dist_i_j = [[10000]*len(partition) for _ in range(len(partition))]

    for i, cluster_i in enumerate(partition):
        s_cluster_term = []
        for pat1 in cluster_i:
            object_to_cluster_mean_list = []
            for j, cluster_j in enumerate(partition):
                #silhouette
                if cluster_i == cluster_j:
                    if len(cluster_i) == 1:
                        continue
                    a = mean(object_to_cluster_distances(pat1, cluster_i, distance_matrix, patients_id))
                    continue
                object_to_cluster_dist = object_to_cluster_distances(pat1, cluster_j, distance_matrix, patients_id)
                object_to_cluster_mean_list.append(mean(object_to_cluster_dist))
                #dunn 
                if min(object_to_cluster_dist) < dist_i_j[i][j]:
                    dist_i_j[i][j] = min(object_to_cluster_dist)
            #silhouette
            if len(cluster_i) == 1:
                s_obj = 0
            else:   
                b = min(object_to_cluster_mean_list)
                s_obj = (b - a)/max(a, b)
            s_cluster_term.append(s_obj)
            s_general.append(s_obj)
        #silhouette
        s_cluster = mean(s_cluster_term)
        s_clusters.append(s_cluster)

    #print('silhouette_C', s_clusters)
    silhouette = mean(s_general)
    #print('silhouette', silhouette)

    #dunn
    min_dist_between = np.matrix(dist_i_j).min()
    if max_cluster_diameter == 0:
        dunn = -10000
    else:
        dunn = min_dist_between/max_cluster_diameter
    #print('dunn', dunn)

    return[mcclainN, mcclainV, cindex, silhouette, dunn]

            
def object_to_cluster_distances(patient_i, cluster, distance_matrix, patients_id):
    dist = []
    idPatient_i = patients_id[patient_i]

    for patient in cluster:
        patient = patients_id[patient]
        if patient == idPatient_i:
            continue
        pair_distance = distance_matrix.loc[(distance_matrix['patient1'] == patient) &
                                                                (distance_matrix['patient2'] == idPatient_i)]['score']
        if pair_distance is None or len(pair_distance) == 0:
            pair_distance = distance_matrix.loc[(distance_matrix['patient1'] == idPatient_i) &
                                                        (distance_matrix['patient2'] == patient)]['score']
        pair_distance = float(pair_distance)
        dist.append(pair_distance)

    return dist


def cluster_validation_indexes(cluster_a,cluster_b):
    #jaccard index
    num_jaccard = len(set(cluster_a).intersection(cluster_b))
    den_jaccard = len(set(cluster_a).union(cluster_b))
    jaccard = num_jaccard/den_jaccard
    
    #the asymmetric measure gama - rate of recovery
    gama = num_jaccard/len(cluster_a)
    
    #the symmetric Dice coefficient
    dice = num_jaccard/(len(cluster_a)+len(cluster_b))
    
    return [jaccard, gama, dice]


