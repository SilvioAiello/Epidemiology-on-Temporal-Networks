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
import scipy.stats
import matplotlib.pyplot as plt
import os
import Assertions_suite
#%%
def input_sampling(samples,N, isDAR):
    """
    This sample generates inputs for the network generation functions, if these input have to be sampled from an external distribution.
    
    Parameters
    ----------
    samples: array
        Matrix, or vector (it depends on the observation)
    N: int
        nodes of the RESULTING distribution you want to get
    isDAR: bool
    
    Returns
    -------
    sampled: np.array
        If it's a DAR network, the output will be a N*N matrix, otherwise a N-vector
    """
    norm = samples.shape[0]*samples.shape[1]
    resh = samples.reshape(norm,1)
    empiric_values = np.linspace(resh.min(),resh.max(),N)
    probabilities = plt.hist(resh, bins = N)[0]/norm 
    assert abs(sum(probabilities)-1) <=0.001
    if isDAR:
        size = (N,N)
    else:
        size = N
    sampled = np.random.choice(empiric_values, size=size,p=probabilities)
    return sampled

def beta_distribution(dataset,shape):
    """
    Returns random extractions by beta distribution sampled from an empiric dataset
    """
    mu = np.mean(dataset)
    sigma = np.std(dataset)
    A = ((1-mu)/(sigma**2) - 1/mu)*(mu**2)
    B = A*(1/mu - 1)
    output = scipy.stats.beta.rvs(A,B, size = shape)
    return output

def hist_plots(fig_count,data,N,isDAR,string_title,saving_path):
    """
    Sequence of recurring operations to plot an histogram.
    
    Parameters
    ----------
    fig_count: int
        It is also an output
    data: N or N*N array 
        data you want to perform an histogram over
    size: int or tuple
    """
    assert N == data.shape[0]
    
    plt.figure(fig_count)
    fig_count +=1
    if isDAR:
        plt.hist(data.reshape(N*N,1))
    else:
        plt.hist(data)
    plt.xlabel("Values")
    plt.ylabel("Occurences")
    plt.grid()
    plt.title(string_title)
    os.makedirs(os.path.dirname(saving_path), exist_ok=True)
    plt.savefig(saving_path)
    return fig_count
    

def network_generation_dar(alpha,xi,P=1, T=100, directed = False):
    """
    Generates a DAR(P) network in form of np.array
    
    Parameters
    ----------
    alpha: np.array 
        N*N matrix expressing probability of "persistence" of previous state of each link; null-diagonal only
    xi: np.array 
        N*N matrix expressing "tossing coin" probability for each link state
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
    #assert N<1000, "Error: for computational ease, this functuion accepts only networks with dimension < 1000"
    
    assert Assertions_suite.check_is_natural(T)
    assert Assertions_suite.check_is_natural(P)
    assert T<1000, "Note: for computational ease, this functuion accepts only networks with duration < 1000"
    assert P < T, "Error: you're trying to generate a DAR(P) model where P is higher or equal to evolution's lasting"
    
    assert type(directed) == bool, "Error: only bool type is allowed for variable directed"
    
    #EVOLUTION
    #Initializazion (with just "coin tossing")
    initial_state = [[np.random.choice((1,0),p=(xi[i,j],1-xi[i,j])) if j>i else 0 for j in range(N)] for i in range(N)] #upper triangle first (j>i)
    if directed == False: #setting lower triangle as directed value
        initial_state = [[initial_state[j][i] if j<i else initial_state[i][j] for j in range(N)] for i in range(N)]
    else:
        initial_state = [[np.random.choice((1,0),p=(xi[i,j],1-xi[i,j])) if j<i else initial_state[i][j] for j in range(N)] for i in range(N)]
    
    temporal_network = np.zeros((T,N,N))
    temporal_network[0] = initial_state
    
    #Network evolution
    for t in range(1,T):
        #UPPER TRIANGLE (j>i):
        temporal_network[t] = [[np.random.choice((temporal_network[t-P,i,j],np.random.choice((1,0), p=(xi[i,j],1-xi[i,j]))),
                        p=(alpha[i,j],1-alpha[i,j])) if j>i else 0 for j in range(N)] for i in range(N)] #DAR(P) generation
        #LOWER TRIANGLE (j<i)
        if directed == False:    
            temporal_network[t] = temporal_network[t] + temporal_network[t].T  #copy upper triangle, if undirected
        else:
            temporal_network[t] = temporal_network[t] + np.array([[np.random.choice((temporal_network[t-P,i,j],np.random.choice((1,0), p=(xi[i,j],1-xi[i,j]))),
                        p=(alpha[i,j],1-alpha[i,j])) if j<i else 0 for j in range(N)] for i in range(N)]) #first choice has as 1° argument a 2-tuple
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
    
    def logistic(x,y):
        return np.exp(x+y)/(1+np.exp(x+y))
    
    #ADJACENCIES COMPUTATION
    temporal_network = np.zeros((T,N,N)) #tensor definition
    for t in range(T):
        prob = np.array([[logistic(theta[i,t],theta[j,t]) if logistic(theta[i,t],theta[j,t])< np.exp(705) else 1 for j in range(N)] for i in range(N)])
        #TGRG FOR UPPER TRIANGLE
        temporal_network[t] = np.array([[np.random.choice((1,0), p=(prob[i,j],1-prob[i,j])) if j>i else 0 for j in range(N)] for i in range(N)])
        #LOWER TRIANGLE (j<i)
        if directed == False:    
            temporal_network[t] = temporal_network[t] + temporal_network[t].T
        else:
            temporal_network[t] = temporal_network[t] + np.array([[np.random.choice((1,0), p=(prob[i,j],1-prob[i,j])) if j<i else 0 for j in range(N)] for i in range(N)])
    return temporal_network, theta

#%% CENTRALITY
def recipr_max_eigen(temporal):
    T = temporal.shape[0]
    spec = []
    for t in range(T):
        spec.append(np.real(max(np.linalg.eigvals(temporal[t])))) #find eigenval with max real part for each adjacency
    a = 1/max(spec) #reciprocal of the maximum eigenvalue
    return 1/a    

def communicability(temporal, a, length_one=True): 
    """
    Return Communicability matrix of a tempnetwork, as defined by Grindrod, and max spectral radius.
    User must set the parameter of computation, a, and can choose if counting walks of all length or just of length 1.
    Lenght 1 is useful for undirected networks?
    
    It uses several numpy functions
    
    Parameters
    ----------
    temporal: np.array
        T*N*N adjacencies (so, an adjacency for each time step; null diagonal)
    
    a: float
        Parameter used in matrix computing, which must be within (0,1]. It might be a quarter of max spec rad
    
    length_one: bool
        If true, the version without matrix inverse is selected.
    
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
    assert a<=1
    assert a>0
        
    #FUNCTION
    T = temporal.shape[0]
    N = temporal.shape[1]
    
    #Communicability builing:
    Q = np.zeros((T,N,N))
    
    if length_one:
        Q[0] = np.identity(N)+a*temporal[0]
        for t in range(1,T):
            matr = np.identity(N)+a*temporal[t] #new factor for that time step
            Q[t] = np.matmul(Q[t-1],matr) #Q updating, no normalization
    else:
        Q[0] = np.linalg.inv(np.identity(N)+a*temporal[0])/np.linalg.norm(np.linalg.inv(np.identity(N)+a*temporal[0]))
        for t in range(1,T):
            matr = np.linalg.inv(np.identity(N)-a*temporal[t]) #new factor for that time step
            Q[t] = np.matmul(Q[t-1],matr)/np.linalg.norm(np.matmul(Q[t-1],matr)) #Q updating and normalizing
    return Q


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

def aggregate_degree(temporal, directed=False):
    """
    Provided with a temporal network, returns AD scores for all nodes.
    Aggregate matrix: sum of 1s, for each couple of nodes, over time.
    If network is directed, it returns 2 ranking vectors, out and in.
    
    It uses several numpy functions
    
    Parameters
    ----------
    temporal: np.array
        T*N*N adjacencies (so, an adjacency for each time step; null diagonal)
    
    Returns
    -------
    lines_sum1: np.array
        N-list of floats representing out-AD scores of nodes
    
    lines_sum2: np.array [if directed = True]
        N-list of floats representing in-AD scores of nodes
    """
    T = temporal.shape[0]
    N = temporal.shape[1]
    
    ad = np.zeros((N,N))
    for t in range(T):
        ad += temporal[t]
    
    if directed:
        out_ranking = np.zeros(N)
        in_ranking = np.zeros(N)
        for i in range(N):
            out_ranking[i] = sum(ad[i])
            in_ranking[i] = sum(ad.T[i])
        return out_ranking,in_ranking
    else:
        ranking = np.zeros(N)
        for i in range(N):
            ranking[i] = sum(ad[i])
        return ranking

def binarized_degree(temporal, directed=False):
    """
    Provided with a temporal network, returns BD scores for all nodes.
    Binarized matrix: 1, for each couple of nodes, if a link exists over time.
    If network is directed, it returns 2 ranking vectors, out and in.
    
    It uses several numpy functions
    
    Parameters
    ----------
    temporal: np.array
        T*N*N adjacencies (so, an adjacency for each time step; null diagonal)
    
    Returns
    -------
    lines_sum1: np.array
        N-list of floats representing out-BD scores of nodes
    
    lines_sum2: np.array [if directed = True]
        N-list of floats representing in-BD scores of nodes
    """
    N = temporal.shape[1]
    
    bd = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            if np.sum(temporal[:,i,j]) != 0:
                bd[i,j] = 1
    
    if directed:
        out_ranking = np.zeros(N)
        in_ranking = np.zeros(N)
        for i in range(N):
            out_ranking[i] = sum(bd[i])
            in_ranking[i] = sum(bd.T[i])
        return out_ranking,in_ranking
    else:
        ranking = np.zeros(N)
        for i in range(N):
            ranking[i] = sum(bd[i])
        return ranking

#%% OLD FUNCTIONS FOR NETWORK ANALYSIS
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