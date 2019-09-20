""" 
Functions in this script work in Pyhon3, may require numpy (v1.16) and allow to:
    1) generate and update a temporal network, according to the dar(p) or trgr laws of evolution.
    2) perform structural analysis, such as degree evolution and centrality measures (BC/RC, AD, BD)

From this module you can extract the following functions:
    * network_generation_dar, network_generation_tgrg
    * degree_node, degree_mean_t, degree_sequencee
    * communicability, broadcast_ranking, receive_ranking
    #TODO: ADD CENTRALITIES AND UPDATE THIS LIST

For further understandings on how this script operates, check file "docs/howto.md".
For further theoretical understandings, check file "docs/explanation.md".
"""
import numpy as np
import Assertions_suite
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
        T*(N*N) array, expressing adjacencies for each time step; this is a stochastic output
    
    Examples
    --------
    
        >>> network_generation_dar(0.6*np.ones((3,3)),0.5*np.ones((3,3)),P=1, T=3, directed = False)
        array([[[0., 1., 0.],
        [1., 0., 1.],
        [0., 1., 0.]],
        [[0., 1., 0.],
        [1., 0., 1.],
        [0., 1., 0.]],
        [[0., 0., 0.],
        [0., 0., 0.],
        [0., 0., 0.]]])
    
    """
    # PRELIMINARY ASSERTS, TO UNSURE FUNCTION WORKS PROPERLY
    assert Assertions_suite.check_is_ndarray(alpha,2)
    assert Assertions_suite.check_is_square(alpha)
    assert Assertions_suite.check_is_probability(alpha)
    
    assert Assertions_suite.check_is_ndarray(xi,2)
    assert Assertions_suite.check_is_square(xi)
    assert Assertions_suite.check_is_probability(xi)

    N = len(alpha) #if everything is ok, get info about number of nodes
    assert N<1000, "Error: for computational ease, this functuion accepts only networks with dimension < 1000"
    
    assert Assertions_suite.check_is_natural(T)
    assert Assertions_suite.check_is_natural(P)
    assert T<1000, "Note: for computational ease, this functuion accepts only networks with duration < 1000"
    assert P < T, "Error: you're trying to generate a DAR(P) model where P is higher or equal to evolution's lasting"
    
    assert type(directed) == bool, "Error: only bool type is allowed for variable directed"
    
    #EVOLUTION
    #Initial adjacency, as tossing simple coin (0.5); if undirected, it's made symmetric
    initial_state = [[np.random.choice((1,0),p=(xi[i,j],1-xi[i,j])) if j>i else 0 for j in range(N)] for i in range(N)] #upper triangle first (j>i)
    if directed == False: #setting lower triangle as directed value
        initial_state = [[initial_state[j][i] if j<i else initial_state[i][j] for j in range(N)] for i in range(N)]
    else:
        initial_state = [[np.random.choice((1,0),p=(xi[i,j],1-xi[i,j])) if j<i else initial_state[i][j] for j in range(N)] for i in range(N)]
    
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
    sigma: np.array 
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
        T*(N*N) tensor, expressing adjacencies for each time step 
        (so, it's entries are 0 or 1, and it can be or not symmetryc); this is a stochastic output
    theta: np.array
        N*T matrix, expressing evolution of fitnesses for each node; this is a stochastic output
    
    Examples
    --------
    
        >>> network_generation_tgrg(0.6*np.ones(3),0.4*np.ones(3),0.1*np.ones(3),T=3, directed = False)
        ( array([[[0., 1., 0.],
         [1., 0., 1.],
         [0., 1., 0.]],
        [[0., 1., 1.],
         [1., 0., 1.],
         [1., 1., 0.]],
        [[0., 1., 1.],
         [1., 0., 1.],
         [1., 1., 0.]]]), 
        array([[0.6       , 0.99664444, 1.09950828],
        [0.6       , 0.9285223 , 1.00334695],
        [0.6       , 0.73700156, 0.85323848]]) )
    
    """
    # PRELIMINARY ASSERTS, TO UNSURE FUNCTION WORKS PROPERLY 
    assert Assertions_suite.check_is_ndarray(phi0,1)
    assert Assertions_suite.check_is_ndarray(phi1,1)
    assert Assertions_suite.check_is_ndarray(sigma,1)
    
    N = len(phi0) #if everything is ok, get info about number of nodes
    assert N<1000, "Error: for computational ease, this functuion accepts only networks with dimension < 1000"
    
    assert Assertions_suite.check_is_natural(T)
    assert T<1000, "Note: for computational ease, this functuion accepts only networks with duration < 1000"
    
    assert type(directed) == bool, "Error: only bool type is allowed for variable directed"
    
    #FITNESSES EVOLUTION
    theta = np.zeros((N,T)) #definition
    theta[:,0] = phi0
    for t in range(1,T):
        theta[:,t] = phi0 + phi1*theta[:,t-1] + np.random.normal(0,sigma,size=N) #evolution for each node and each time

    #ADJACENCIES COMPUTATION
    temporal_network = np.zeros((T,N,N)) #tensor definition
    for t in range(T):
        prob = np.array([[np.exp(theta[i,t]+theta[j,t])/(1+np.exp(theta[i,t]+theta[j,t])) for j in range(N)] for i in range(N)])
        #TGRG FOR UPPER TRIANGLE
        temporal_network[t] = np.array([[np.random.choice((1,0), p=(prob[i,j],1-prob[i,j])) if j>i else 0 for j in range(N)] for i in range(N)])
        #LOWER TRIANGLE (j<i)
        if directed == False:    
            temporal_network[t] = temporal_network[t] + temporal_network[t].T
        else:
            temporal_network[t] = temporal_network[t] + np.array([[np.random.choice((1,0), p=(prob[i,j],1-prob[i,j])) if j<i else 0 for j in range(N)] for i in range(N)])
    return temporal_network, theta
    
#%% NETWORK ANALYSYS
def degree_node(network,node,out = True):
    """
    Returns out- or in-going degree of a node at a certain time
    
    Parameters
    ----------
    network: np.array
        N*N adjacency (so, for a selected time; null diagonal)
    node: int
        index of the node of interest
    out: bool
        Expresses interest in out or in-going degree; if undirected graph, both give the same result
    
    Returns
    -------
    degree: int
    
    Examples
    --------
    
        >>> degree_node(np.array([[0,0,1],[1,0,1],[1,0,0]]),0)
        1
        
        >>> degree_node(np.array([[0,0,1],[1,0,1],[1,0,0]]),1)
        2
        
    """
    #ASSERTS
    assert Assertions_suite.check_is_ndarray(network,2)
    assert Assertions_suite.check_is_square(network)
    assert Assertions_suite.check_is_nulldiagonal(network)
    
    assert type(node) == int, "Error: index for nodes must be an integer"
    assert node <= len(network), "Error: node not present"
     
    #FUNCTION
    if out:
        return sum(network[node])
    else:
        return sum(network[:,node])
    
def degree_mean(tempnet, out = True):
    """
    Returns out- or in-going mean degree of a whole temporal network, even of duration 1 (i.e. a single adjacency).
    It returns a list of mean degrees over time.
    If you're interested in just one instant, just select the proper output entry.
    
    Parameters
    ----------
    tempnet: np.array
        TNN temporal network
    out: bool
        Expresses interest in out or in-going degree; if undirected graph, both give the same result
    
    Returns
    -------
    d_seq: list
        list of floats, one per temporal step, expressing mean degree at that time
    
    Examples
    --------
    
        >>> degree_mean(np.array([[0,0,1],[1,0,1],[1,0,0]]))
        [1.3333333333333333]
        
        >>> degree_mean(np.array([[[0,0,1],[1,0,1],[1,0,0]],
        [[0,1,1],[1,0,1],[1,0,0]]]))
        [1.3333333333333333, 1.6666666666666667]
        
    """
    assert isinstance(tempnet,np.ndarray)
    
    N = tempnet.shape[1] #this is true anyway
    #Infer temporal duration and perform assertions
    if len(tempnet.shape) == 3:
        T = tempnet.shape[0]
        for t in range(T):
            assert Assertions_suite.check_is_square(tempnet[t]) #check square for each step
            assert Assertions_suite.check_is_nulldiagonal(tempnet[t]) #check null diagonal for each step
        network = tempnet #same new for 3- or 2-dimensional arrays
    elif len(tempnet.shape) == 2:
        T = 1
        assert Assertions_suite.check_is_square(tempnet) #check square for each step
        assert Assertions_suite.check_is_nulldiagonal(tempnet) #check null diagonal for each step
        network = np.zeros((T,N,N))
        network[0] = tempnet #same new for 3- or 2-dimensional arrays
    else:
        raise AssertionError ("You must provide an adjacency or a sequence of adjacencies")
        
    #FUNCTION
    d_seq = []
    for t in range(T):
        degrees = []
        for node in range(N):
            if out:
                degrees.append(sum(network[t][node]))
            else:
                degrees.append(sum(network[t][:,node]))
        d_seq.append(np.average(degrees))
    return d_seq #output is always a list

#%% CENTRALITY
def communicability(temporal): 
    """
    Return Communicability matrix of a tempnetwork, as defined by Grindrod, and max spectral radius.
    
    It uses several numpy functions
    
    Parameters
    ----------
    temporal: np.array
        T*N*N adjacencies (so, an adjacency for each time step; null diagonal)
    
    Returns
    -------
    rec_maxradius: float
        Reciprocal of maximum eigenvalue of the whole set of adjacencies
    Q: np.array
        N*N matrix (N being number of nodes), expressing "how well information can be passed from i to j"
    
    Examples
    --------
    
        >>> communicability(np.array([[[0,0,1],[1,0,1],[1,0,0]],
        [[0,1,1],[1,0,1],[1,0,0]]]))
        ( 0.6180339887498951, 
        array([[0.54377529, 0.0840179 , 0.17483766],
        [0.20185156, 0.52294101, 0.20185156],
        [0.16411784, 0.0253576 , 0.53305547]]))
    """
    #At the moment, this function takes as default, as coefficient, a quarter of the inverse of max spectral radius
    
    T = temporal.shape[0]
    N = temporal.shape[1]
    
    #ASSERTS
    assert Assertions_suite.check_is_ndarray(temporal,3) #ask it to be a square array of the first dimension of the network
    for t in range(T):
        assert Assertions_suite.check_is_square(temporal[t])
        assert Assertions_suite.check_is_nulldiagonal(temporal[t]) #check null diagonal for each step
        
    #FUNCTION
    T = temporal.shape[0]
    N = temporal.shape[1]
    #Find max spectral radius:
    spec = []
    for t in range(T):
        spec.append(np.real(max(np.linalg.eigvals(temporal[t])))) #find eigenval with max real part for each adjacency
    rec_maxradius = 1/max(spec) #reciprocal of the maximum eigenvalue
    #Communicability builing:
    Q = np.identity(N)/np.linalg.norm(np.identity(N)) #initialization (and normalization)
    for t in range(T):
        inv = np.linalg.inv(np.identity(N)-0.25*rec_maxradius*temporal[t]) #new factor for that time step
        Q = np.matmul(Q,inv)/np.linalg.norm(np.matmul(Q,inv)) #Q updating and normalizing
    return(rec_maxradius,Q)

def communicability_onelink(temporal): 
    """
    Return Communicability matrix of a tempnetwork, counting only one link per snapshot, and max spectral radius.
    
    It uses several numpy functions
    
    Parameters
    ----------
    temporal: np.array
        T*N*N adjacencies (so, an adjacency for each time step; null diagonal)
    
    Returns
    -------
    rec_maxradius: float
        Reciprocal of maximum eigenvalue of the whole set of adjacencies
    Q: np.array
        N*N matrix (N being number of nodes), expressing "how well information can be passed from i to j"
    
    Examples
    --------
    
        >>> communicability(np.array([[[0,0,1],[1,0,1],[1,0,0]],
        [[0,1,1],[1,0,1],[1,0,0]]]))
        ( 0.6180339887498951, 
        array([[0.54377529, 0.0840179 , 0.17483766],
        [0.20185156, 0.52294101, 0.20185156],
        [0.16411784, 0.0253576 , 0.53305547]]))
    """
    #At the moment, this function takes as default, as coefficient, a quarter of the inverse of max spectral radius
    
    T = temporal.shape[0]
    N = temporal.shape[1]
    
    #ASSERTS
    assert Assertions_suite.check_is_ndarray(temporal,3) #ask it to be a square array of the first dimension of the network
    for t in range(T):
        assert Assertions_suite.check_is_square(temporal[t])
        assert Assertions_suite.check_is_nulldiagonal(temporal[t]) #check null diagonal for each step
        
    #FUNCTION
    T = temporal.shape[0]
    N = temporal.shape[1]
    #Find max spectral radius:
    spec = []
    for t in range(T):
        spec.append(np.real(max(np.linalg.eigvals(temporal[t])))) #find eigenval with max real part for each adjacency
    inv_maxradius = 1/max(spec) #reciprocal of the maximum eigenvalue
    #Communicability builing:
    Q = np.identity(N)/np.linalg.norm(np.identity(N)) #initialization (and normalization)
    for t in range(T):
        matr = np.identity(N)+0.25*inv_maxradius*temporal[t] #new factor for that time step
        Q = np.matmul(Q,matr)/np.linalg.norm(np.matmul(Q,matr)) #Q updating and normalizing
    return(inv_maxradius,Q) 

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
    
    Examples
    --------
    
        >>> broadcast_ranking(np.array([[0.54377529, 0.0840179 , 0.17483766],
        [0.20185156, 0.52294101, 0.20185156],
        [0.16411784, 0.0253576 , 0.53305547]]))
        (array([0.80263085, 0.92664413, 0.72253091]), 
        array([1, 0, 2], 
        dtype=int64))
    
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
    
    Examples
    --------
    
        >>> receive_ranking(np.array([[0.54377529, 0.0840179 , 0.17483766],
        [0.20185156, 0.52294101, 0.20185156],
        [0.16411784, 0.0253576 , 0.53305547]]))
        (array([0.90974469, 0.63231651, 0.90974469]), 
        array([2, 0, 1], 
        dtype=int64))
    """
    #FUNCTION
    lines_sum = np.sum(Q, axis = 0) #Broadcast -> sum over columns:
    rank = np.flip(np.argsort(lines_sum)) #argsort -> increasing score; flip -> decreasing
    return(lines_sum,rank)