""" 
The functions contained in this module allow:
    1) to generate and update a network according to the dar(p) or trgr laws of evolution.
    (if you don't know what they are, check the documentation)
    2) to perform analysis of network's structure, such as degree evolution...
    3)... and centrality measures (Communicability, AD, BD)

A network is described through the temporal sequence of its adiacency matrices.
Each update is a stochastic process, so several performances should be performed.

It requires the following libraries (so, check if they are installed):
    * numpy, used for its data structures and anaylisis, and to get random functions 
    * matplotlib.pyplot, used for graphics plot and belluries
    * pickle, used to store, in an efficient way, the complex information generated
    * networkx, just used to plot small networks

From this module you can extract the following functions:
    * network_generation_dar, that generates a DAR(P) network in form of np.array
    * network_generation_dar, that generates a TGRG network in form of np.array
    
    * degree_node, degree_mean_t, degree_sequence, that return degree of a node, 
    mean degree of all nodes at a time step, mean degree evolution over time
    #TODO: CENTRALITA
    
    * make_name, that generate the path where storing nodes and results,
    according to the syntax illustrated in documentation
    * network save, that saves a generated network through pickle
    * plot save, if you want to save something in automatized way, 
    with automatically generated name

ONLY UNDERICTED ARE DELT WITH, AT THE MOMENT. DOVREI METTERE "SE UND, FAI COSI', ECC"
"""

import numpy as np 
import networkx as nx 
import matplotlib.pyplot as plt 
import pickle

### EVOLUTIONS

#First, an initial state is created, i.e. the adiacency matrix for T=0.
#In future development, this state will be the same for both evolutions.
#Then, it's initialized a 3D tensor, N*N*T, i.e. a sequence of adiacency matrices evolving in time.
#Each of the T adiacences matrices (A) will be updated according to the values of the previous one.
#Adiacences are kept simmetric with 0-diagonal.
#This scheme will be the same for both evolutions, just the law of updating will be different.

#In DAR, there's a probability alpha_ij of having A_ij(t) = A_ij(t-1)
#Otherwise, the new state will be extacted by a "coin toss", with probability xi_ij of getting "1"
#The evolution is performed just for the upper triangular adiancency, and than the matrix is simmetrized

def assert_probability(structure):
    assert (structure<= 1).all(),"Error: at least one element in a probability matrix is >1, so it's not a probability"
    assert (structure>0).all(),"Error: at least one element in probability matrix is <0, so it's not a probability"

def assert_natural(number):
    assert type(number) == int, "Error: you should provide an integer"
    assert number>0, "Error: you should prove a positive number"

def network_generation_dar(alpha,xi,N=100,T=100, P=1, directed = False):
    """
    Generates a DAR(P) network in form of np.array
    
    Parameters
    ----------
    alpha: np.array 
        N*N matrix expressing probability of "persistence" of previous state of each link
    xi: np.array 
        N*N matrix expressing "tossing coin" probability for each link state
    N: int (default = 100)
        Natural number expressing network dimension
    T: int (default = 100)
        Natural number expressing duration of temporal evolution
    P: int (default = 1)
        Natural number expressing order of Discrete AutoRegressive model
    directed: Boolean (default = False)
        Determines symmetry of the network, at each time step
        
    
    Returns
    -------
    temporal_network: np.array
        T*(N*N) array, expressing adiacencies for each time step 
        (so, it's entries are 0 or 1, and it can be or not symmetryc)
    
    """
    # Asserts #
    assert len(alpha) == N, "Error: you provided an alpha matrix whose dimension is different from that of the network"
    assert len(xi) == N, "Error: you provided an alpha matrix whose dimension is different from that of the network"    
    assert_probability(alpha)
    assert_probability(xi)
    
    assert_natural(N)
    assert_natural(T)
    assert_natural(P)
    assert N<1000, "Error: for computational ease, this functuion accepts only networks with dimension < 1000"
    assert T<1000, "Note: for computational ease, this functuion accepts only networks with duration < 1000"
    assert P < T, "Error: you're trying to generate a DAR(P) model where P is higher or equal to evolution's lasting"
    
    assert type(directed) == bool, "Error: only bool type is allowed for variable directed"
    
    #TODO: SPIEGA ALTROVE CHE L'IDEA E' DI LAVORARE SEMPRE E COMUNQUE PRIMA CON LA TRIANGOLARE SUPERIORE E POI, A SECONDA, AGGIORNARE L'INFERIORE
        #SE SERVE, AGGIUNGI QUESTO, SENNO' ELIMINA:
    #First comprehension says: for each row i, taking only the upper triangle (j>i), choose the value 
    #from the same cell of the T-P matrix, with probability alpha, or from a coin tossing (whose probability of getting 1 is xi).
    #The second comprehension just simmetrizes the matrix
    
    #Initial adiacency, as tossing simple coin (0.5); if undirected, it's made symmetric
    initial_state = [[np.random.choice((0,1),p=(0.5,1-0.5)) if j>i else 0 for j in range(N)] for i in range(N)] #upper triangle first (j>i)
    if directed == False: #setting lower triangle as directed value
        initial_state = [[initial_state[j][i] if j<i else initial_state[i][j] for j in range(N)] for i in range(N)]
    else:
        initial_state = [[np.random.choice((0,1),p=(0.5,1-0.5)) if j<i else initial_state[i][j] for j in range(N)] for i in range(N)]
    
    #Temporal network initilialization
    temporal_network = np.zeros((T,N,N))
    temporal_network[0] = initial_state
    
    #Network evolution
    for t in range(1,T):
        #UPPER TRIANGLE (j>i):
        temporal_network[t] = [[np.random.choice((temporal_network[t-P,i,j],np.random.choice((1,0), p=(xi[i,j],1-xi[i,j]))),
                        p=(alpha[i,j],1-alpha[i,j])) if j>i else 0 for j in range(N)] for i in range(N)] #DAR(P) generation
        #LOWER TRIANGLE (j<i)
        if directed == False:    
            temporal_network[t] = [[temporal_network[t,j,i] if j<i else temporal_network[t,i,j] for j in range(N)] for i in range(N)]
            #copy upper triangle, if undirected
        else:
            temporal_network[t] = [[np.random.choice((temporal_network[t-P,i,j],np.random.choice((1,0), p=(xi[i,j],1-xi[i,j]))),
                        p=(alpha[i,j],1-alpha[i,j])) if j<i else temporal_network[t,i,j] for j in range(N)] for i in range(N)]
            #follow the same rule as upper, if directed
    return temporal_network

def network_generation_tgrg(phi0,phi1,sigma,N=100,T=100):
    """
    Generates a DAR(P) network in form of np.array
    
    Parameters
    ----------
    bla: int
        bla bla
    
    Returns
    -------
    bla: int
        bla
    
    """
    #phi0, phi1, sigma (DEV ST!) sono vettori N, per es
    #phi0 = 0.27*np.ones(N) #intercepts, typical of each node
    #phi1 = 0.18*np.ones(N) #slopes, typical of each node
    
    theta = np.zeros((N,T)) #fitness matrix creation
    #So, theta is inizialized and then the random walk is performed
    theta[:,0] = np.random.normal(0,sigma,size=N)
    for t in range(1,T): #costruzione dei random walks
        for i in range(N):
            theta[i,t] = phi0[i] + phi1[i]*theta[i,t-1] + np.random.normal(0,sigma)
    
    #Adiacency computation
    #Now everything is set to let the adiacencies evolve.
    
    #First, an initial state is created, i.e. the adiacency matrix for T=0.
    #In future development, this state will be the same for both evolutions.
    #Then, it's initialized a 3D tensor, N*N*T, i.e. a sequence of adiacency matrices evolving in time.
    #Each of the T adiacences matrices (A) will be updated according to the values of the previous one.
    #Adiacences are kept simmetric with 0-diagonal.
    #This scheme will be the same for both evolutions, just the law of updating will be different.
    
    #In FITNESS Model, the probability of having a link between two nodes is a sort of logistic of the sum of their fitnesses.
    #Being a probability, it obviously takes values from 0 to 1: exp(theta1+theta2)/(1+exp(theta1+theta2)).
    #So, for each temporal step, the algorithm takes the nodes i and j, computes the probability and chooses whether rising a link or not.
    #The evolution is performed just for the upper triangular adiancencies, and than the matrix is simmetrized
    
    temporal_network = np.zeros((T,N,N)) #Creation of Adiacency tensor
    for t in range(T):
        for i in range(N):
            for j in range(N):
                if j>i: #working on the upper triangoular matrix
                    prob = np.exp(theta[i,t]+theta[j,t])/(1+np.exp(theta[i,t]+theta[j,t])) #probability for these nodes at these time
                    temporal_network[t,i,j] = np.random.choice((1,0),p=(prob,1-prob))
                    if temporal_network[t,i,j] == 1: #symmetrization
                        temporal_network[t,j,i] = 1
    
    return temporal_network, theta
    

# NETWORK ANALYSYS #
def degree_node(network,node,out = True):
    """
    Generates a DAR(P) network in form of np.array
    
    Parameters
    ----------
    bla: int
        bla bla
    
    Returns
    -------
    bla: int
        bla
    
    """
    assert node <= len(network), "Error: node not present"
    #you give an adiacency[t], and a node: this function computes node degree
    #if out = true, sums the line (link that start from node)
    #if out = false, sums the column (link that reach that node)
    #if graph is undirected, out = in
    if out:
        return sum(network[node])
    else:
        return sum(network[:,node])

def degree_mean_t(network,out = True):
    """
    Generates a DAR(P) network in form of np.array
    
    Parameters
    ----------
    bla: int
        bla bla
    
    Returns
    -------
    bla: int
        bla
    
    """
    #this function is not "pure", since it relies on another function
    #network must be a N*N np.array
    #you give an adiacency for a certain t: it computes mean degree at that t, using degree_node
    #depending in out value, it returns out degree or in degree. If undirected, it's the same
    degrees = []
    for node in range(len(network)):
        degrees.append(degree_node(network,node,out))
    return(np.average(degrees)) #it returns a number

def degree_mean_sequence(network,T,out = True):
    """
    Generates a DAR(P) network in form of np.array
    
    Parameters
    ----------
    bla: int
        bla bla
    
    Returns
    -------
    bla: int
        bla
    
    """
    #here you give the whole set of adiacencies
    d_seq = []
    for t in range(T):
        d_seq.append(degree_mean_t(network[t],out))
    return d_seq #it return a list
    

# NETWORK SAVE #

def make_name(N,T,isDAR, start,P=1):
    """ isDAR = True -> add "DAR(P); False -> add "TGRG"
    start should be like a str like "alphaeqsxieqs"
    """
    if isDAR:
        return "Networks/N"+str(N)+"_T"+str(T)+"_DAR"+str(P)+"_"+start
    if isDAR == False:
        return "Networks/N"+str(N)+"_T"+str(T)+"_TGRG_"+start

def network_save(network, N,T,isDAR, start,P=1,k=1):
    #Il motivo per cui lo metto qua è che questa funzione genera l'intero percorso dove storo tutto
    name = make_name(N,T,isDAR,start,P)+"/realization"+str(k)+"/network.txt"
    with open(name, 'wb') as handle:
        return pickle.dump(network,handle)

# PLOT FUNCTIONS #
def networkplot(graph,figcount,savefig=False, figname = None, k = None, t = None):
    """Requires an adiacency[t] and the figure count
    If you want to save the image, you have to say savefig=True, and pass a 2-ple and the value of k and t.
    Tuple figname contains the values needed by make_name: N,T,isDAR, start,P (P MUST BE SPECIFIED HERE)"""
    figcount+=1
    plt.figure(figcount)
    plt.title("")
    nxgraph = nx.from_numpy_matrix(graph) #convert np-matrix in a Network, readable by nx library:
    nx.draw(nxgraph, pos = nx.drawing.layout.circular_layout(nxgraph), with_labels = True) #draw a circular plot and show the number of each node
    if savefig and figname and k and t:
        plt.savefig(make_name(figname)+"Results"+'_realization%i_t=%i.pdf' %(k,t))
    return figcount+1, plt.show()


### TESTS:
def test_symmetry(temporal_network,T=100):
    for t in range(T):
        assert (temporal_network[t] == temporal_network[t].T).any, "Error: network at t = %i is not symmetric" %t

def test_nulldiagonal(temporal_network,T=100):
    for t in range(T):
        assert sum(np.diag(temporal_network[t])) == 0, "Error: network at t = %i has not-0 diagonal" %t

#Save the network for further use
#np.savetxt(start_name+'%i_N%i_wholenetwork_T%i.txt' %(P,N,T), temporal_network.reshape(T*N*N,1))
#To import:
#new_data = np.loadtxt('start_name+'%i__N%i_wholenetwork_T%i.txt' %(P,N,T))
#new_data = new_data.reshape((T,N,N))