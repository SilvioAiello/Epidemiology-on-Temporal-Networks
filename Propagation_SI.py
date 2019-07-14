"""
Functions contained in this script allow:
    1) to spread an epidemic in SI mode over a temporal network
    2) to measure virulence of each node
    
Functions work in Pyhon3, and may require the following libraries (so, check if they are installed):
    * numpy, used for its data structures and anaylisis, and to get random functions 
    * pickle, used to store, in an efficient way, the complex information generated
[If you want to get some plots, you may use matplotlib.pyplot, for plots belluries]

From this module you can extract the following functions:
#TODO: UPDATE

For further understandings on how this script operates, check file "howto.md"
For further theoretical understandings, check file "explanation.md"
"""
import numpy as np #Used for its random functions and data structures
import pickle

#%%
def network_load(N,T,start,k=1,isDAR=True,P=1):
    """ Loads a previously generated temporal network using pickle (so it must be installed), if it's present
    
    Parameters
    ----------
    N: int
        number of nodes
    T: int
        temporal duration
    start: string
        identificative name of the network
    k: int (default = 1)
        iteration of network realization
    isDAR: bool (default = True)
        defines whether to search a DAR(P) or TGRG
    P: int (default = 1)
        defines DAR order
    
    Returns
    -------
    pickle.load(f): np.array
        If path and file exist, the TNN-tempnet is returned
    """
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

#%%
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
#%%
def propagation(tempnet,index_case,probabilities): #Remember: the outcome of this function is stochastic
    '''
    Produces the evolution of disease states over a temporal network
    
    Parameters
    ----------
    tempnet: np.array
        T*N*N (the functions extracts by itself T and N)
    index_case: int
        Index-case node (defining the initial state)
    probabilities: dict
        Dictionary with T keys, from 0 to T-1, expressing probability of infection for each possible duration of a contact
        
    Returns
    -------
    states_sequence: dict
        It's a dictonary of dictionaries: first key selects time, second one selects node.
        So, states_sequence[0] is a dictionary, in which each key stands for the node-state at that time
        states_sequence[0][i] is the state at the time for i-th node.

    '''
    T = tempnet.shape[0] #shape of array returns (T,N,N)-tuple
    N = tempnet.shape[1] #shape of array returns (T,N,N)-tuple
    
    def set_infected(node,t): #once one is infected,it stays infeced
        for instant in range(t,T):
            states_sequence[instant][node] = 1
        return     
    
    #Output initialization:
    states_sequence = dict()
    for t in range(T):
        states_sequence[t] = dict.fromkeys(range(N),0)
    set_infected(index_case,0)
    
    #Sets initialization
    susceptibles = onlyzeros(range(N),states_sequence[0]) #"targets"
    infecteds = {index_case} #they will be decisive to change target state
    
    for t in range(1,T): #TEMPI RICHIESTI IN: NEIGHBOUR, CONTACT LASTING, SET INFECT, E VEDI CHI E' -1 E CHI NO
        print("DETERMINO STATO ALL'ISTANTE " +str(t))
        for s in susceptibles.copy(): #copy avoids rising an error when the iteration set changes
            infectneighbourhood = neighbourhood(tempnet[t-1],s).intersection(infecteds)
            for i in infectneighbourhood.copy(): 
                if probabilities[contact_lasting(tempnet,states_sequence,t-1,i,s)]>np.random.uniform(0,1): #rand extraction
                        set_infected(s,t) #if successful, change the state of the node, at next t
                        susceptibles.remove(s)
                        infecteds.add(s)
                        break
            else:
                continue # only executed if the inner loop did NOT break
            break  # only executed if the inner loop DID break
    return(states_sequence)

#%%
def infected_counter(set_of_nodes):
    """
    Counts the number of infected nodes in a set
    """
    counter = 0
    for i in range(len(set_of_nodes)):
        if set_of_nodes[i]==1:
            counter+=1
    return counter
def time_score(scores,fraction):
    """
    Returns the time each node spent to infect a fraction of network, for one whole iteration.
    
    The output is the first time step when the fraction of infected nodes is bigger than the given fraction.
    If this never happens, function just returns the last time step of network evolution.
    Output is computed using external function "infected_counter".
    
    Parameters
    ----------
    scores: dict of dicts
        Sequence of T dictionaries of the states of all N nodes (T and N are extracted by function iteself)
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
    return time_spent 
#TODO: ha senso restituire il total time? magari inf, magari errore, vediamo