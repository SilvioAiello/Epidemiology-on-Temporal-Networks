"""
Functions contained in this script allow:
    1) to spread an epidemic in SI mode over a temporal network
    2) to measure virulence of each node

There are two epidemic-state: Susceptible (0) and Infected (1).
Nodes virulence is determined by:
    * making each node, at a time, the index case, and 
    * computing min or avg time infection takes to spread to a certain % of the whole network
Each propagation is a stochastic process, so several performances should be performed, with same initial condition.
Since tempnet structures are themself stochastic, the whole set of iterations should be repeated over several ones.

Functions work in Pyhon3, and may require the following libraries (so, check if they are installed):
    * numpy, used for its data structures and anaylisis, and to get random functions 
    * pickle, used to store, in an efficient way, the complex information generated
    * scipy.integrate, just for function "quad", usate to integrate

[If you want to get some plots, you should also use: matplotlib.pyplot]

From this module you can extract the following functions:

    
"""
#a temporal network of n nodes is described by t adiacency matrices (n x n), which have the sctructure of numpy array. 
#desease state is a dictionary that maps each node to 0 (susceptible) or 1 (infected).

import numpy as np #Used for its random functions and data structures
import matplotlib.pyplot as plt #Used for graphic belluries
import pickle


# TABLE OF CONTENTS: after initial arrangements, two sections develop and perform SI propagations, and two sections develop and perform centrality measures.
# Final section makes comparisons and prints result. Note that centrality measures are indipendent of the desease, and may be performed before its propagation.
# 1) Set of parameters and network import
# 2) Containers and functions for epidemic
# 3) Epidemic propagation
# 4) Containers and functions for centrality measures
# 5) Centrality measures
# 6) Comparisons and results print

###                     PARAMETERS OF THE SYSTEM & NETWORKS IMPORT         ###
fig_count = 0 #several figures may be printed


def network_load(N,T,start,k=1,isDAR=True,P=1):
    name = str()
    if isDAR:
        name = "Networks/N"+str(N)+"_T"+str(T)+"_DAR"+str(P)+"_"+start+"/realization"+str(k)+"/network.txt"
    else:
        name = "Networks/N"+str(N)+"_T"+str(T)+"_TGRG_"+start+"/realization"+str(k)+"/network.txt" 
    with open(name, 'rb') as f:
        return pickle.load(f)

#TODO: VERIFICA CHE QUESTI DUE SONO OK (PER ORA SONO COMMENTATI PERCHE' DEVO VERIFICARE ALTRE COSE)
        #(secondo me non va perche' anche se è salvato con quel titolo, c'è qualcosa di strano)
#temporal_dar = network_load(100,100,'alphaeqs_xieqs',k=1,isDAR=True,P=1)
#temporal_fitn= load_network(100,100,isDAR=True,alleq,k=1)

#TODO: determinare il destino di questo vecchio modo di caricare
#temporal_dar = np.loadtxt('Networks\N100_T100_DAR1_alphaeqs_xieqs\realization1'+'%i_N%i_wholenetwork_T%i.txt' %(P,N,T))
#temporal_dar = temporal_dar.reshape((T,N,N)) #salvatolo come unica colonna, confermo che vuole tempo righe colonne
#temporal_fitn= np.loadtxt('Examples/FITN'+'_N%i_wholenetwork_T%i.txt' %(N,T))
#temporal_fitn= temporal_fitn.reshape((T,N,N))


###                      CONTAINERS AND FUNCTIONS FOR EPIDEMIC           ###

# CONTAINERS #
#Initialization of dictionaries that keep track of node states evolution. 
#This will be the syntax: label_[indexcase][repetition][instant]. If you set these 3 params, you get a dict of nodes state for that index, that repetition, at that time.
label_dar = []
label_fitn= []
    
#These lists save the a number of time steps needed to infect 60% of the network, for each node, for each iteration (they will be a list of N lists with K items):
score_dar = []
score_fitn =[]
#In order to perform a direct comparison with structural centralities, average virulence score is computed for each node,
#And virulence ranking will be stored in a list of ordered nodes, from those with lowest to those with highest average time to infect.
avg_dar = []
avg_fitn = []


# EPIDEMIC FUNCTIONS BUILDING #
#The idea of SI model is that each I node, at t-1, can make infect one of its S neighbours, at t-1, with a certain probability, whose rate is beta.
#So, there will be usefuls some simple functions that build the neighbourhood of a node, or find S/I nodes from a list at a given time, or compute links' temporal duration.

def neighbourhood(adiacency,node):
    """Extracts the neighbourhood reachable by a node a given time. 
    
    Parameters
    ----------
    adiacency : np.array
        Adiacency matrix (dar, tgrg...), N*N-array of values 0 or 1, expressing existence state of a link between two nodes at a certain time
        As an example, it may be a temporal_dar[t]
    node: int
        Node of interest
    
    Returns
    -------
    neigh: set
        Set of neighbours linked to the input node
    """
    neigh = {i for i in range(len(adiacency)) if adiacency[node,i]==1}
    return neigh

def onlyzeros(nodes_set,states_dict):
    """Extracts the susceptible ones from a list of nodes.
    
    Parameters
    ----------
    nodes_set: tuple/list/set
        Ensemble of nodes (not necessairly all network's nodes)
    states_dict: dict
        Dictionary mapping state of required nodes to a 0 or 1 value.
        It can have more keys than nodes_set: the extra ones won't be scanned.
        As an example, state_dict may be a label_dar[t]
    
    Returns
    -------
    selected: set
        Set of susceptible nodes
    """
    selected = {node for node in nodes_set if states_dict[node]==0} #0 means susceptible
    return selected

def contact_lasting(adiacency,state,t,infected_node,susceptible_node):
    #UN MODO PER TESTARLA E' VEDERE CHE AD OGNI ISTANTE SUCCESSIVO E' DIVERSO DAL PRECEDENTE, PER LA STESSA COPPIA
    #This function computes the duration of a contact (in number of temporal steps) for a couple of nodes I-S, as long as I is infected (otherwise, it couldn't propagate the epidemic)
    #This is accomplished by checking backwards the existence of the link and the state of the I node, increasing the value of a counter variable until these conditions are satisfied
    if state[t][infected_node] != 1:
        print("Error: infected node %i is not infected at the latest instant %i" %(infected_node,t))
        raise AssertionError
    assert state[t][susceptible_node] == 0, "Error: susceptible node is not susceptible at the latest instant"
    assert adiacency[t,infected_node,susceptible_node] == 1, "Error: nodes couple is not linked at latest instant"
    counter = 0
    for instant in range(t+1):
        if (adiacency[t-instant,infected_node,susceptible_node] == 1 and state[t-instant][infected_node] != state[t-instant][susceptible_node]):
            counter +=1
        else:
            break
    return counter #this should be included in [0,t]

#OBSOLETE:
def infect_extraction(probab): #Remember: the outcome of this function is stochastic
    #This function extracts a random number from Uniform Distribution, and compares it to the given probability
    #If the random number is lower, function return True, to signal that contagion can occur
    x = np.random.uniform(0,1) #soglia di probabilità da superare
    if probab > x:
        return (True)
    else:
        return (False)

def onlyones(nodes_set,states_dict):
    """Extracts the infected ones from a list of nodes.
    
    Parameters
    ----------
    nodes_set: tuple/list/set
        Ensemble of nodes (not necessairly all network's nodes)
    states_dict: dict
        Dictionary mapping state of required nodes to a 0 or 1 value.
        It can have more keys than nodes_set: the extra ones won't be scanned
    
    Returns
    -------
    selected: set
        Set of infected nodes
    """
    selected = {node for node in nodes_set if states_dict[node]==1} #1 means infected
    return selected
#FINE OBSOLETE

#%%
    #le prime tre righe erano sopra contact lasting
#Epidemic spread follows Chen approach: time of infection (in unit steps) follows a Poissonian distribution, normalized to return beta for 1 step, integrated within link duration.
#(note: beta is the probability rate of contagion [1/s], but also the actual probability of contagion after 1 unit time: infact, P(1) = beta*1 u.t. = beta [dimensionless]).
#Chen's algorithm make a contagion happen by performing an extraction from Unif(0,1): if this number is lower than the Poisson integral, contagion takes place at that time.

#This function performs the actual propagation, from 0 to T, given an index case.
#It finds the S neighboroughs of all I nodes at a given time, and makes them infect at t+1, with probability beta
#This function will be evoked for each time step but the first and the last one
def propagation(tempnet,index_case,probabilities): #Remember: the outcome of this function is stochastic
    '''
    Produces the evolution of disease states over a temporal network
    
    Parameters
    ----------
    network: np.array
        T*N*N (the functions extracts by itself T and N)
    index_case: int
        Index-case node (defining the initial state)
        
    Returns
    -------
    states_sequence: dict
        It's a dictonary of dictionaries: first key selects time, second one selects node.
        So, states_sequence[0] is a dictionary, in which each key stands for the node-state at that time
        states_sequence[0][i] is the state at the time for i-th node.

    '''
    def set_infected(node,t):
        for instant in range(t,T):
            states_sequence[instant][node] = 1
        return 
    
    #La funzione produce il dizionario iniziale, e poi un dizionario temporaneo per ogni step, 
    #che fa corrispondere alla valore del key i del principale
    T = tempnet.shape[0] #shape return tuple (T,N,N)
    N = tempnet.shape[1] #shape return tuple (T,N,N)
    
    
    #Outout initialization:
    states_sequence = dict()
    for t in range(T):
        states_sequence[t] = dict.fromkeys(range(N),0)
    set_infected(index_case,0)
    
    #Sets initialization
    susceptibles = onlyzeros(range(N),states_sequence[0]) #"targets"
    infecteds = {index_case} #they will be decisive to change target state
    
    #First, it is inialized the state-dict at t+1, equal to the previous one, so who is infected stays infected
    for t in range(1,T): #TEMPI RICHIESTI IN: NEIGHBOUR, CONTACT LASTING, SET INFECT, E VEDI CHI E' -1 E CHI NO
        print("Determino lo stato all'istante " +str(t)+" in cui gli infetti sono " + str(infecteds))
        for s in susceptibles.copy(): #copy avoids rising an error when the iteration set changes
            infectneighbourhood = neighbourhood(tempnet[t-1],s).intersection(infecteds) #takes the infected neigbs. of that node
            for i in infectneighbourhood.copy(): #then, for each infected neighbour, perform the extraction, according to link duration
                    if probabilities[contact_lasting(tempnet,states_sequence,t-1,i,s)]>np.random.uniform(0,1):
                        print("Nuova infezione! " + str(i)+ " Ha infettato " + str(s))
                        set_infected(s,t) #if successful, change the state of the node, at next t
                        susceptibles.remove(s)
                        infecteds.add(s)
                        break
    return(states_sequence) #in questo momento non servono nè onlyones nè infect extraction

def infected_counter(infected):
    #This function counts the number of infected nodes in a set:
    counter = 0
    for i in range(len(infected)):
        if infected[i]==1:
            counter+=1
    return counter

#This function, provided of the whole states-evolution, returns the time step at wich the disease has reached a given fraction of the network, using infected_counter.
#If that fraction has not been reached, it return the total time of propagation:
def time_score(scores,fraction):
    """
    Returns the time each node spent to infect a fraction of network, for one whole iteration
    
    The output is the first time step when the fraction of infected nodes is bigger than the given fraction.
    If this never happens, function just returns the last time step of network evolution.
    
    Parameters
    ----------
    scores: dict of dicts
        Sequence of T dictionaries of the states of all N nodes (T and N are extracted bt function iteself)
    fraction: float
        Real number, from 0 to 1, representing network fraction
        
    Returns
    -------
    time_spent: int
        Score for that iteration
    """
    #Scores sarebbe il label dar da 0 a T, e lui deve trovare in quale label_dar[t] la soglia è stata superata
    T,N = len(scores),len(scores[0])
    assert fraction > 0, "Error, only values between 0 and 1 are allowed"
    assert fraction < 1, "Error, only values between 0 and 1 are allowed"
    #asserire che sum deve essere sempre <= N?
    time_spent = T-1 #initialized as the final temporal step
    for t in range(T):
        if infected_counter(scores[t])>=fraction*N:
            time_spent = t
            break
    return time_spent #ha senso restituire il total time? magari inf, magari errore, vediamo