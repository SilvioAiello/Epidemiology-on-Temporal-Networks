""" 
Functions in this script work in Pyhon3, may require numpy (v1.16) and allow to:
    1) generate and update a temporal network, according to the dar(p) or trgr laws of evolution.
    2) perform structural analysis, such as degree evolution and centrality measures (BC/RC, AD, BD)

From this module you can extract the following functions:
    * network_generation_dar, network_generation_tgrg
    * degree_node, degree_mean_t, degree_sequencee
    * communicability, broadcast_ranking, receive_ranking
    #TODO: ADD CENTRALITIES AND UPDATE THIS LIST

For further understandings on how this script operates, check file "howto.md".
For further theoretical understandings, check file "explanation.md".
"""
import numpy as np
import Test_suit
#%%
def network_generation_dar(alpha,xi,P=1, T=100, directed = False):
    """
    Generates a DAR(P) network in form of np.array
    
    Parameters
    ----------
    alpha: np.array 
        N*N matrix expressing probability of "persistence" of previous state of each link; null-diagonal only
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
        (so, it's entries are 0 or 1, and it can be or not symmetryc; null-diagonal only)
    """
    # PRELIMINARY ASSERTS, TO UNSURE FUNCTION WORKS PROPERLY
    Test_suit.assert_ndarray(alpha,2)
    Test_suit.assert_square(alpha)
    Test_suit.assert_probability(alpha)
    
    Test_suit.assert_ndarray(xi,2)
    Test_suit.assert_square(xi)
    Test_suit.assert_probability(xi)

    N = len(alpha) #if everything is ok, get info about number of nodes
    assert N<1000, "Error: for computational ease, this functuion accepts only networks with dimension < 1000"
    
    Test_suit.assert_natural(T)
    Test_suit.assert_natural(P)
    assert T<1000, "Note: for computational ease, this functuion accepts only networks with duration < 1000"
    assert P < T, "Error: you're trying to generate a DAR(P) model where P is higher or equal to evolution's lasting"
    
    assert type(directed) == bool, "Error: only bool type is allowed for variable directed"
    
    #EVOLUTION
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

#%%
def network_generation_tgrg(phi0,phi1,sigma,T=100, directed = False):
    """
    Generates a TGRG in form of np.array
    
    Parameters
    ----------
    phi0: np.array 
        Vector of length N expressing intercept for each of the N nodes
    phi1: np.array 
        Vector of length N expressing angoular coefficient for each of the N nodes
    epsilon: np.array 
        Vector of length N expressing standard deviation for the gaussian stochastic term for each of the N nodes
    N: int (default = 100)
        Natural number expressing network dimension
    T: int (default = 100)
        Natural number expressing duration of temporal evolution
    directed: Boolean (default = False)
        Determines symmetry of the network, at each time step
    
    Returns
    -------
    temporal_network: np.array
        T*(N*N) tensor, expressing adiacencies for each time step 
        (so, it's entries are 0 or 1, and it can be or not symmetryc)
    theta: np.array
        N*T matrix, expressing evolution of fitnesses for each node
    """
    #phi0 = 0.27*np.ones(N) #intercepts, typical of each node (values from Mazzarisi-Lillo paper)
    #phi1 = 0.18*np.ones(N) #slopes, typical of each node
    
    # PRELIMINARY ASSERTS, TO UNSURE FUNCTION WORKS PROPERLY 
    Test_suit.assert_ndarray(phi0,1)
    Test_suit.assert_ndarray(phi1,1)
    Test_suit.assert_ndarray(sigma,1)
    
    N = len(phi0) #if everything is ok, get info about number of nodes
    assert N<1000, "Error: for computational ease, this functuion accepts only networks with dimension < 1000"
    
    Test_suit.assert_natural(T)
    assert T<1000, "Note: for computational ease, this functuion accepts only networks with duration < 1000"
    
    assert type(directed) == bool, "Error: only bool type is allowed for variable directed"
    
    #FITNESSES EVOLUTION
    theta = np.zeros((N,T)) #definition
    theta[:,0] = phi0
    for t in range(1,T):
        theta[:,t] = phi0 + phi1*theta[:,t-1] + np.random.normal(0,sigma,size=N) #evolution for each node and each time

    #ADIACENCIES COMPUTATION
    temporal_network = np.zeros((T,N,N)) #tensor definition
    for t in range(T):
        for i in range(N):
            for j in range(N):
                #TGRG FOR UPPER TRIANGLE
                temporal_network[t] = [[np.random.choice((1,0), p=(np.exp(theta[i,t]+theta[j,t])/(1+np.exp(theta[i,t]+theta[j,t])),
                                1-np.exp(theta[i,t]+theta[j,t])/(1+np.exp(theta[i,t]+theta[j,t])))) if j>i else 0 for j in range(N)] for i in range(N)]
                #LOWER TRIANGLE (j<i)
                if directed == False:    
                    temporal_network[t] = [[temporal_network[t,j,i] if j<i else temporal_network[t,i,j] for j in range(N)] for i in range(N)]
                    #copy upper triangle, if undirected
                else:
                    temporal_network[t] = [[np.random.choice((1,0), p=(np.exp(theta[i,t]+theta[j,t])/(1+np.exp(theta[i,t]+theta[j,t])),
                                    1-np.exp(theta[i,t]+theta[j,t])/(1+np.exp(theta[i,t]+theta[j,t])))) if j<i else temporal_network[t,i,j] for j in range(N)] for i in range(N)]
    return temporal_network, theta
    
#%% NETWORK ANALYSYS
def degree_node(network,node,out = True):
    """
    Returns out- or in-going degree of a node at a certain time
    
    Parameters
    ----------
    network: np.array
        N*N adiacency (so, for a selected time; null diagonal)
    node: int
        index of the node of interest
    out: bool
        Expresses interest in out or in-going degree; if undirected graph, both give the same result
    
    Returns
    -------
    degree: int
    
    """
    #ASSERTS
    Test_suit.assert_ndarray(network,2)
    Test_suit.assert_square(network)
    Test_suit.assert_nulldiagonal(network)
    
    Test_suit.assert_natural(node)
    assert node <= len(network), "Error: node not present"
     
    #FUNCTION
    if out:
        return sum(network[node])
    else:
        return sum(network[:,node])

def degree_mean_t(network,out = True):
    """
    Returns out- or in-going mean degree of a network at a certain time.
    
    Parameters
    ----------
    network: np.array
        N*N adiacency (so, for a selected time; null diagonal)
    out: bool
        Expresses interest in out or in-going degree; if undirected graph, both give the same result
    
    Returns
    -------
    degree: float
        Mean degree at that time
    
    """
    #ASSERTS
    Test_suit.assert_ndarray(network,2) 
    Test_suit.assert_square(network)
    Test_suit.assert_nulldiagonal(network)
    
    #FUNCTION
    degrees = []
    for node in range(len(network)):
        if out:
            degrees.append(sum(network[node]))
        else:
            degrees.append(sum(network[:,node]))
    return(np.average(degrees)) #it returns a float

def degree_mean_sequence(network,T, initial = 0, out = True):
    """
    Returns out- or in-going mean degree of a whole sequence of networks
    
    It makes use of function np.average, so check if numpy is available first
    
    Parameters
    ----------
    network: np.array
        T*N*N adiacency (so, an adiacency for each time step; null diagonal)
    T: int
        Natural number expressing last instant to check
    initial: int (default = 0)
        Natural number expressing first instant to check
    out: bool
        Expresses interest in out or in-going degree; if undirected graph, both give the same result
    
    Returns
    -------
    degree: float
        Mean degree
    
    """
    #ASSERTS
    Test_suit.assert_natural(T)
    assert T>initial, "Error: something wrong in initial-final time step"
    
    Test_suit.assert_ndarray(network,3) #ask it to be a square array of the first dimension of the network
    [Test_suit.assert_square(network[t]) for t in range(len(network))] #check square for each step
    [Test_suit.assert_nulldiagonal(network[t]) for t in range(initial,T)] #check null diagonal for each step
    
    #FUNCTION
    d_seq = []
    for t in range(T):
        degrees = []
        for node in range(len(network)):
            if out:
                degrees.append(sum(network[node]))
            else:
                degrees.append(sum(network[:,node]))
        d_seq.append(np.average(degrees))
    return d_seq #it return a list

#%% CENTRALITY
def communicability(temporal): 
    """
    Return Communicability matrix of a tempnetwork, as defined by Grindrod, and max spectral radius.
    
    It uses several numpy functions
    
    Parameters
    ----------
    temporal: np.array
        T*N*N adiacencies (so, an adiacency for each time step; null diagonal)
    
    Returns
    -------
    maxradius: float
        Reciprocal of maximum eigenvalue of the whole set of adiacencies
    Q: np.array
        N*N matrix (N being number of nodes), expressing "how well information can be passed from i to j"
    """
    #At the moment, this function takes as default, as coefficient, a quarter of the inverse of max spectral radius
    
    #ASSERTS
    Test_suit.assert_ndarray(temporal,3) #ask it to be a square array of the first dimension of the network
    [Test_suit.assert_square(temporal[t]) for t in range(len(temporal))] #check square for each step
    [Test_suit.assert_nulldiagonal(temporal[t]) for t in range(len(temporal))] #check null diagonal for each step
    
    #FUNCTION
    T = temporal.shape[0]
    N = temporal.shape[1]
    #Find max spectral radius:
    spec = []
    for t in range(T):
        spec.append(np.real(max(np.linalg.eigvals(temporal[t])))) #find eigenval with max real part for each adiacency
    maxradius = 1/max(spec) #reciprocal of the maximum eigenvalue
    #Communicability builing:
    Q = np.identity(N)/np.linalg.norm(np.identity(N)) #initialization (and normalization)
    for t in range(T):
        inv = np.linalg.inv(np.identity(N)-0.25*maxradius*temporal[t]) #inverse factor, which has to be multiplicated to the previous Q
        Q = np.matmul(Q,inv)/np.linalg.norm(np.matmul(Q,inv)) #updating and normalizing of Q
    return(maxradius,Q) 

def broadcast_ranking(Q):
    """
    Provided with a Communicabiltiy, return a list of nodes sorted by Broadcast ranking, and the scores
    
    It uses several numpy functions
    
    Parameters
    ----------
    Q: np.array
        N*N communicability (according to Grindrod definition and previously generated here)
    
    Returns
    -------
    lines_sum: np.array
        N-list of floats representing Broacast Centrality scores of nodes
    rank: np.array
        N list of nodes, sorted by ranking (rank[0] is the most central node)
    """
    #FUNCTION    
    lines_sum = np.sum(Q, axis = 1) #Broadcast -> sum over lines:
    rank = np.flip(np.argsort(lines_sum)) #argsort -> nodes ordered by increasing score; flip -> decreasing
    return(lines_sum,rank)
    
def receive_ranking(Q):
    """
    Provided with a Communicabiltiy, returns a list of nodes sorted by Receiving ranking, and the scores
    
    It uses several numpy functions
    
    Parameters
    ----------
    Q: np.array
        N*N communicability (according to Grindrod definition and previously generated here)
    
    Returns
    -------
    lines_sum: np.array
        N-list of floats representing Receiveing Centrality scores of nodes
    rank: np.array
        N list of nodes, sorted by ranking (rank[0] is the most central node)
    """
    #FUNCTION
    lines_sum = np.sum(Q, axis = 0) #Broadcast -> sum over columns:
    rank = np.flip(np.argsort(lines_sum)) #argsort -> increasing score; flip -> decreasing
    return(lines_sum,rank)