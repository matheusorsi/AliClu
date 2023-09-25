import numpy as np
import pandas as pd
from synthetic_data import ctmc_sequences

def generate_dataset(n_sequences, n_dataset, filename="time_experiment.csv"):
#Generates differents datasets, configurated by parameter n_dataset. Each configuration defines the number and types of
#sequence clusters, each one containing n_sequences
#Dataset 2 generates a bigger variety of clusters
    
    #initialize list that will contain the auxiliary dataframes to be concataneted
    concat = [] 

    if n_dataset > 3:
        return n_events_datasets(n_sequences, n_dataset, filename)

    #maximum number of different events
    n_events = 4

    ###############################################################
    #CLUSTER 1 and 2 - A->B
    ##############################################################
    
    clusters = 2
    #rates of state changes for the clusters
    rates = [1000, 10]
    
    #generate sequences
    for i in range(0, clusters):
    
        rate = rates[i] #rate of the transition

        # initial distribution for the states
        iDistribution = [0]*n_events
        iDistribution[0] = 1 

        genMatrix = np.zeros((n_events, n_events)) # generator matrix
        genMatrix[0][0:2] = [-rate, rate]

        df_aux = ctmc_sequences(n_events, iDistribution, genMatrix, n_sequences) # temporal sequences
        concat.append(df_aux)

    ################################################################
    #CLUSTER 3 e 4 - B->C
    ################################################################
    if n_dataset == 2 or n_dataset == 3:
        #generate sequences
        for i in range(0,clusters):
            
            rate = rates[i] #rate of the transition

            # initial distribution for the states
            iDistribution = [0]*n_events
            iDistribution[1] = 1 

            genMatrix = np.zeros((n_events, n_events)) # generator matrix
            genMatrix[1][1:3] = [-rate, rate]
            
            df_aux = ctmc_sequences(n_events, iDistribution, genMatrix, n_sequences) # temporal sequences
            concat.append(df_aux)

    ################################################################
    #CLUSTER 5, 6, 7 e 8 - C -> D, A -> C
    ################################################################
    if n_dataset == 3:
        #generate sequences C -> D
        for i in range(0,clusters):
            
            rate = rates[i] #rate of the transition

            # initial distribution for the states
            iDistribution = [0]*n_events
            iDistribution[2] = 1 

            genMatrix = np.zeros((n_events, n_events)) # generator matrix
            genMatrix[2][2:4] = [-rate, rate]
            
            df_aux = ctmc_sequences(n_events, iDistribution, genMatrix, n_sequences) # temporal sequences
            concat.append(df_aux)

        #generate sequences A -> C
        for i in range(0,clusters):
            
            rate = rates[i] #rate of the transition

            # initial distribution for the states
            iDistribution = [0]*n_events
            iDistribution[0] = 1 

            genMatrix = np.zeros((n_events, n_events)) # generator matrix
            genMatrix[0][0] = -rate
            genMatrix[0][2] = rate
            
            df_aux = ctmc_sequences(n_events, iDistribution, genMatrix, n_sequences) # temporal sequences
            concat.append(df_aux)

    df_encoded = pd.concat(concat, ignore_index = True)
    #numerate patients from 0 to N-1, where N is the number patients
    df_encoded['id_patient'] = df_encoded.index.tolist()
    df_encoded.to_csv(filename, index=False)
    return df_encoded

def n_events_datasets(n_sequences, n_dataset, filename):
    #initialize list that will contain the auxiliary dataframes to be concataneted
    concat = [] 
    #maximum number of different events
    match n_dataset:
        case 4:
            n_events = 2
        case 5:
            n_events = 4
        case 6:
            n_events = 8

    clusters = 2
    #rates of state changes for the clusters
    rates = [1000, 10]
    
    #generate sequences
    for i in range(0, clusters):
    
        rate = rates[i] #rate of the transition

        # initial distribution for the states
        iDistribution = [0]*8
        iDistribution[0] = 1 

        genMatrix = generator_matrix(rate, 8)

        df_aux = ctmc_sequences(n_events, iDistribution, genMatrix, n_sequences) # temporal sequences
        concat.append(df_aux)

    df_encoded = pd.concat(concat, ignore_index = True)
    #numerate patients from 0 to N-1, where N is the number patients
    df_encoded['id_patient'] = df_encoded.index.tolist()
    df_encoded.to_csv(filename, index=False)
    return df_encoded

def generator_matrix(rate, size):

    genMatrix = np.zeros((size, size)) # generator matrix

    for i in range(0, size):
        genMatrix[i][:] = rate
        genMatrix[i][i] = -genMatrix.sum(axis=1)[i] + rate
    #print (genMatrix)
    return genMatrix