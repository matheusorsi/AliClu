import pandas as pd
import numpy as np
from collections import ChainMap
from pathos.multiprocessing import ProcessingPool as Pool



def ctmc_sequences(n_events, iDistribution, genMatrix, n_sequences):
#Obtain a sample path with a maximum of n_events for a continuous-time Markov chain with
#initial distribution iDistribution (pi) and generator matrix genMatrix (Q). The procedure is repetead n_sequences 
#times to obtain sequences 

    #the jump matrix P gives the probability P[m][n] of getting to the n_th state given the current state m
    jumpMatrix = compute_jump_matrix(genMatrix)

    #dictionary used to replace state by alphabetical order letters
    encodeState = {0: "A", 1: "B", 2:"C", 3:"D", 4:"E", 5:"F", 6:"G", 7:"H", 8:"I", 9:"J"}

    n_states = len(genMatrix) #number of possible states   

    def create_sequence(seqIter):

        sequence = {}

        #initialize the process at accTime = 0 with initial state i_state drawn from iDistribution
        accTime = 0
        i_state = int(np.random.choice(n_states, 1, p=iDistribution)) #initial state
        #print("i = " + encodeState[i_state])

        #initialize the sequence with the initial state
        sequence[seqIter] = '0.' + encodeState[i_state]

        for k in range(0, n_events - 1):

            q_ii = -genMatrix[i_state][i_state]
            if q_ii == 0:
                break
            else:
                #exponential holding time
                holdTime = float(np.random.exponential(q_ii, 1))
                accTime = accTime + holdTime
                #next state drawn from the i_th row of the stochastic matrix
                i_state = int(np.random.choice(n_states, 1, p=jumpMatrix[i_state][:])) 
                #update the sequence with the time elapsed holdTime and the current state i_state
                sequence[seqIter] = sequence[seqIter] + ',' + str(round(holdTime,2)) + '.' + str(encodeState[i_state])

        return sequence
    
    
    with Pool(4) as pool:
        sequenceList = pool.map(create_sequence, range(0, n_sequences))

    #initialize the dataframe structure that will support the sequences
    dataFrame = pd.DataFrame.from_dict(ChainMap(*sequenceList), orient='index', columns = ['aux_encode'])
    #print(dataFrame)

    return dataFrame

def compute_jump_matrix(genMatrix):
#Receives the generator matrix and generates the stochastic (jump) matrix

    size = len(genMatrix)

    jumpMatrix = np.zeros((size, size))

    for i in range(0, size):
        q_ii = -genMatrix[i][i]
        if q_ii == 0:
            #Absorbing state
            jumpMatrix[i][:] = 0
            jumpMatrix[i][i] = 1
        else:
            jumpMatrix[i][:] = genMatrix[i][:]/q_ii
            jumpMatrix[i][i] = 0

    return jumpMatrix

def save_sequence(sequence, sequenceList):
    sequenceList.append(sequence)