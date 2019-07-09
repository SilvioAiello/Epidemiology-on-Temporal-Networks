""" 
Functions contained in this script allow:
    1) to generate and update a network according to the dar(p) or trgr laws of evolution.
    (if you don't know what they are, check the documentation)
    2) to perform analysis of network's structure, such as degree evolution...
    3)... and centrality measures (Communicability, AD, BD)

A network is described through the temporal sequence of its adiacency matrices.
Each update is a stochastic process, so several performances should be performed.

It requires the following libraries (so, check if they are installed):
    * numpy, used for its data structures and anaylisis, and to get random functions 
    * pickle, used to store, in an efficient way, the complex information generated

[If you want to plot, you should use these libraries:
    * matplotlib.pyplot, used for graphics plot and belluries
    * networkx, just used to plot small networks]

From this module you can extract the following functions:
    * network_generation_dar, that generates a DAR(P) network in form of np.array
    * network_generation_dar, that generates a TGRG network in form of np.array
    
    * degree_node, degree_mean_t, degree_sequence, that return degree of a node, 
    mean degree of all nodes at a time step, mean degree evolution over time
    
    * communicability, broadcast_ranking, receive_ranking
    build communicability matrix and extract from it broacdcast and receive centrality for each node
    #TODO: ADD CENTRALITIES AND UPDATE THIS LIST
    
    * network save, that saves a generated network through pickle, in a path and
    with a name that follows syntax illustrated in documentation
    * plot save, if you want to save something in automatized way, 
    with automatically generated name
"""

import numpy as np 
import pickle

def assert_ndarray(matrix,dimensions):
     #check if matrix is a np.array, has a certain dimension
    assert isinstance(matrix,np.ndarray), "Error: matrix must be a numpy array"
    assert len(matrix.shape) == dimensions, "Error: matrix is not 2-dimensional"

def assert_square(matrix):
    #check if matrix is a square by comparing each side length with the first side
    for i in range(len(matrix.shape)):
        assert matrix.shape[i] == matrix.shape[0], "Error: matrix is not a square"

def assert_probability(matrix):
    #check that all elements are probabilities
    assert (matrix<= 1).all(),"Error: at least one element in a probability matrix is >1, so it's not a probability"
    assert (matrix>=0).all(),"Error: at least one element in probability matrix is <0, so it's not a probability"

def assert_natural(number):
    #check that a number is a positive integer
    assert isinstance(number, int), "Error: %f is not an integer, but it should be" %number
    assert number>0, "Error: %i is not positive, but it should be" %number

def assert_nulldiagonal(matrix):
    #check that a matrix has null diagonal
    assert sum(np.diag(matrix)) == 0, "Error: network has not-0 diagonal"
    
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
    assert_ndarray(alpha,2)
    assert_square(alpha)
    assert_probability(alpha)
    
    assert_ndarray(xi,2)
    assert_square(xi)
    assert_probability(xi)

    N = len(alpha) #if everything is ok, get info about number of nodes
    assert N<1000, "Error: for computational ease, this functuion accepts only networks with dimension < 1000"
    
    assert_natural(T)
    assert_natural(P)
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
    #phi0 = 0.27*np.ones(N) #intercepts, typical of each node
    #phi1 = 0.18*np.ones(N) #slopes, typical of each node
    
    #TODO: SPIEGARE DA QUALCHE PARTE CHE PRIMA FACCIAMO EVOLVERE COMPLETAMENTE LE FITNESS E POI PERFORMIAMO L'EVOLUZIONE
        #SE SERVE, AGGIUNGI QUESTO:
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
    
    # PRELIMINARY ASSERTS, TO UNSURE FUNCTION WORKS PROPERLY 
    assert_ndarray(phi0,1)
    assert_ndarray(phi1,1)
    assert_ndarray(sigma,1)
    
    N = len(phi0) #if everything is ok, get info about number of nodes
    assert N<1000, "Error: for computational ease, this functuion accepts only networks with dimension < 1000"
    
    assert_natural(T)
    assert T<1000, "Note: for computational ease, this functuion accepts only networks with duration < 1000"
    
    assert type(directed) == bool, "Error: only bool type is allowed for variable directed"
    
    #FITNESSES EVOLUTION
    theta = np.zeros((N,T)) #definition
    theta[:,0] = np.random.normal(0,sigma,size=N) #random-normal initialization
    for t in range(1,T):
        theta[:,t] = phi0 + phi1*theta[:,t-1] + np.random.normal(0,sigma,size=N) #evolution for each node and each time

    #ADIACENCIES COMPUTATION
    temporal_network = np.zeros((T,N,N)) #tensor definition
    for t in range(T):
        for i in range(N):
            for j in range(N):
                #TGRG FOR UPPER TRIANGLE
                temporal_network[t] = [[np.random.choice((1,0), 
                                p=(np.exp(theta[i,t]+theta[j,t])/(1+np.exp(theta[i,t]+theta[j,t])),1-np.exp(theta[i,t]+theta[j,t])/(1+np.exp(theta[i,t]+theta[j,t])))) if j>i else 0 for j in range(N)] for i in range(N)]
                #LOWER TRIANGLE (j<i)
                if directed == False:    
                    temporal_network[t] = [[temporal_network[t,j,i] if j<i else temporal_network[t,i,j] for j in range(N)] for i in range(N)]
                    #copy upper triangle, if undirected
                else:
                    [[np.random.choice((1,0), 
                                p=(np.exp(theta[i,t]+theta[j,t])/(1+np.exp(theta[i,t]+theta[j,t])),1-np.exp(theta[i,t]+theta[j,t])/(1+np.exp(theta[i,t]+theta[j,t])))) if j<i else temporal_network[t,i,j] for j in range(N)] for i in range(N)]
    return temporal_network, theta
    
#%%
# NETWORK ANALYSYS #
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
    assert_ndarray(network,2)
    assert_square(network)
    assert_nulldiagonal(network)
    
    assert_natural(node)
    assert node <= len(network), "Error: node not present"
    
    #TODO: SPIEGARE COSA FA, PUOI USARE:
    #if out = true, sums the line (link that start from node)
    #if out = false, sums the column (link that reach that node)
    #if graph is undirected, out = in
    
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
    assert_ndarray(network,2) 
    assert_square(network)
    assert_nulldiagonal(network)
    
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
    assert_natural(T)
    assert T>initial, "Error: something wrong in initial-final time step"
    
    assert_ndarray(network,3) #ask it to be a square array of the first dimension of the network
    [assert_square(network[t]) for t in range(len(network))] #check square for each step
    [assert_nulldiagonal(network[t]) for t in range(initial,T)] #check null diagonal for each step
    
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

#%%
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
    #TODO: DOCUMENTA COME FUNZIONA (ANCHE ALTROVE DEVI FARLO, DICENDO CHE TROVA DA SOLO LE LUNGHEZZE E PERCHE')
        #As known, to compute communicability one as to choose a coefficient that multiplicates adiacencies
        #At the moment, this function takes as default, as coefficient, a quarter of the inverse of max spectral radius
    
    #ASSERTS
    assert_ndarray(temporal,3) #ask it to be a square array of the first dimension of the network
    [assert_square(temporal[t]) for t in range(len(temporal))] #check square for each step
    [assert_nulldiagonal(temporal[t]) for t in range(len(temporal))] #check null diagonal for each step
    
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
    #TODO: Nello spiegone documentato puoi usare questo:
    #Next two functions compute broadcast/receiving centralities and node rankings
    #For centralities, they use function np.sum, where one can choose to sum of lines (BC) or columns (RC)
    #For rankings, they use np.argsort, whose input is a list and output is a list of the indices of the input, sorted according to their decreasing values
    #So, the first element of the output list has the highest rank
    
    lines_sum = np.sum(Q, axis = 1) #Broadcast -> sum over lines:
    rank = np.flip(np.argsort(lines_sum)) #argsort -> increasing score; flip -> decreasing
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
    lines_sum = np.sum(Q, axis = 0) #Broadcast -> sum over columns:
    rank = np.flip(np.argsort(lines_sum)) #argsort -> increasing score; flip -> decreasing
    return(lines_sum,rank)


# NETWORK SAVE #
def network_save(network, start,isDAR=True,P=1, k=1):
    """ Saves network using pickle (so it must be present) and giving it its proper name (with identification provided by user) and folder
    
    Parameters
    ----------
    network: np.array
        T*N*N (the functions extracts by itself T and N)
    start: str
        Name you choose for the function
    isDAR: bool (default = True)
        User is required to specify wheter network is DAR or TGRG
    P: int (default = 1)
        Natural number expressing the order of DAR. If network is TGRG, its value doesn't affect the result
    k: int (default = 1)
        Natural number expressing what iteration of the same network this is
    
    Returns
    -------
    name: network.txt
        A file txt with the network, in a folder that follows syntax presented in documentation (if path doesn't exist, it's created)
    """
    #ASSERTS
    assert_natural(P) #there's no need to perform other checks, since they have been already performed
    assert_natural(k)
    
    #FUNCTION
    T = network.shape[0]
    N = network.shape[1]
    
    name = str()
    if isDAR:
        name = "Networks/N"+str(N)+"_T"+str(T)+"_DAR"+str(P)+"_"+start+"/realization"+str(k)+"/network.txt"
    else:
        name = "Networks/N"+str(N)+"_T"+str(T)+"_TGRG_"+start+"/realization"+str(k)+"/network.txt"
    
    with open(name, 'wb') as handle:
        return pickle.dump(network,handle)