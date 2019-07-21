"""
Functions in this script work in Pyhon3, may require numpy (v1.16) and allow to:
    1) spread an epidemic, in SI mode, over a temporal network
    2) measure virulence of each node

From this module you can extract the following functions:
    * neighbourhood, onlyzeros, contact_lasting, poisson_probability
    * propagation
    * infected_counter, time_score

For further understandings on how this script operates, check file "howto.md".
For further theoretical understandings, check file "explanation.md".
"""
import numpy as np #Used for its random functions and data structures
import Test_suite
#%% EASING FUNCTIONS
def neighbourhood(adjacency,node):
    """Extracts the neighbourhood reachable by a node a given time.
    So, it doesn't care about network's directness.
    
    Parameters
    ----------
    adjacency : np.array
        Adiacency matrix (dar, tgrg...), N*N-array of values 0 or 1, expressing existence state of a link between two nodes at a certain time
        As an example, it may be a temporal_dar[t]
    node: int
        Node of interest
    
    Returns
    -------
    neigh: set
        Set of neighbours linked to the input node
    """
    #ASSERTS
    Test_suite.assert_ndarray(adjacency,2)
    Test_suite.assert_square(adjacency)
    Test_suite.assert_nulldiagonal(adjacency)
    
    assert isinstance(node, int)
    assert node <= len(adjacency), "Error: node not present"
     
    #FUNCTION
    neigh = {i for i in range(len(adjacency)) if adjacency[node,i]==1}
    return neigh

def onlyzeros(nodes_set,states_dict):
    """Extracts the susceptible ones from a list of nodes.
    
    Parameters
    ----------
    nodes_set: tuple/list/set
        Ensemble of nodes (not necessairly all network's nodes)
    states_dict: N-dict
        Dict mapping state of required nodes to a 0 or 1 value.
        It can have more keys than nodes_set: the extra ones won't be scanned.
        As an example, state_dict may be a label_dar[t]
    
    Returns
    -------
    selected: set
        Set of susceptible nodes
    """
    assert type(states_dict) == dict, "states_dict is not a dictionary"
    
    selected = {node for node in nodes_set if states_dict[node]==0} #0 means susceptible
    return selected

def contact_lasting(tempnet,states_sequence,t,infected_node,susceptible_node):
    """
    Computes the number of consectuive temporal steps a couple of nodes I-S are linked, as long as they hold these states.
    
    Parameters
    ----------
    tempnet: np.array
        TNN temporal network
    states_sequence: dict
        TN dict, states evolution at each time step
    t: int
        time step from which going backwards in time
    infected_node: int
        
    susceptible_node: int
        
    Returns
    -------
    counter: int
        SI contact duration; if 0, nodes are not linked, if t, they have been linked since network is born
    """
    #ASSERTS
    assert isinstance(t, int)
    assert isinstance(infected_node, int)
    assert isinstance(susceptible_node, int)
    assert isinstance(states_sequence, dict), "states_sequence is not a dictionary"
    
    T = tempnet.shape[0]
    N = tempnet.shape[1]
    assert t <= T, "t is bigger than the total network duration"
    assert infected_node<N, "infected_node doesn't exist"
    assert susceptible_node<N, "susceptible_node doesn't exist"
    
    Test_suite.assert_ndarray(tempnet,3) #ask it to be a square array of the first dimension of the network
    [Test_suite.assert_square(tempnet[i]) for i in range(t)] #check square for each step
    [Test_suite.assert_nulldiagonal(tempnet[i]) for i in range(t)] #check null diagonal for each step
    
    assert infected_node != susceptible_node, "Nodes are the same"
    assert states_sequence[t][infected_node] == 1, "Infected node is not infected at time %i"%(infected_node,t)
    assert states_sequence[t][susceptible_node] == 0, "Node %i is not susceptible at time %i"%(susceptible_node,t)
    assert tempnet[t,infected_node,susceptible_node] == 1, "Nodes %i and %i are not linked at time %i"%(infected_node,susceptible_node,t)
    
    #FUNCTION
    counter = 0
    for instant in range(t+1):
        if (tempnet[t-instant,infected_node,susceptible_node] == 1 and states_sequence[t-instant][infected_node] != states_sequence[t-instant][susceptible_node]):
            counter +=1
        else:
            break
    return counter

def poisson_probability(t,beta):
    """
    This function reproduces the Poisson PDF, whose average depends on beta.
    Its integral is the probability of having and infection for a I-S contact lasting t.
    
    Parameters
    ----------
    t: int
        time duration of a contact
    beta: float
        infection rate
        
    Returns
    -------
    A float, representing PDF value for that t, with given beta.
    """
    assert beta <1, "Error: beta must be <1"
    assert beta >=0, "Error: beta must be >=0"
    lamda = -np.log(1-beta) # Chen, Benzi use 60
    return(lamda*np.exp(-lamda*t))
    
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
    
    #ASSERTS
    assert isinstance(index_case, int)
    assert index_case<N, "index_case doesn't exist"
    
    Test_suite.assert_ndarray(tempnet,3) #ask it to be a square array of the first dimension of the network
    [Test_suite.assert_square(tempnet[t]) for t in range(T)] #check square for each step
    [Test_suite.assert_nulldiagonal(tempnet[t]) for t in range(T)] #check null diagonal for each step
    
    assert isinstance(probabilities, dict), "Probabilities is not a dictionary"
    assert len(probabilities) == T, "Probabilities doesn't meet tempnet duration"
    
    #FUNCTION
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
    
    for t in range(1,T):
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

#%% EPIDEMIC SCORE COMPUTING FUNCTIONS
def infected_counter(set_of_nodes):
    """
    Counts the number of infected nodes in a states-set at a certain time
    """
    counter = 0
    for i in range(len(set_of_nodes)):
        if set_of_nodes[i]==1:
            counter+=1
    return counter
def time_score(scores_evolution,fraction):
    """
    Returns the time each node spent to infect a fraction of network, for one whole iteration.
    
    The output is the first time step when the fraction of infected nodes is bigger than the given fraction.
    If this never happens, function just returns the last time step of network evolution.
    Output is computed using external function "infected_counter".
    
    Parameters
    ----------
    scores: TN-dict
        Sequence of T dictionaries of the states of all N nodes (T and N are extracted by function iteself)
    fraction: float
        Real number, from 0 to 1, representing network fraction
        
    Returns
    -------
    time_spent: int
        Score for that iteration
    """
    T = len(scores_evolution)
    N = len(scores_evolution[0])
    
    #ASSERTS
    assert fraction > 0, "Error, only values between 0 and 1 are allowed"
    assert fraction < 1, "Error, only values between 0 and 1 are allowed"
    assert isinstance(scores_evolution, dict), "scores_evolution is not a dictionary"
    
    #FUNCTION
    time_spent = T-1 #initialized as the final temporal step
    for t in range(T):
        if infected_counter(scores_evolution[t])>=fraction*N:
            time_spent = t
            break
    return time_spent