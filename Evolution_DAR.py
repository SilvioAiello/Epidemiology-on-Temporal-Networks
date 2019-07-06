#THIS MODULE ALLOWS TO UPDATE A NETWORK ACCORDING TO THE DAR(P) LAW OF EVOLUTION.
#IF YOU DON'T KNOW WHAT'S A DAR(P), CHECK THE DOCUMENTATION
#THE NETWORK IS DESCRIBED THROUGH IT'S ADIACENCY MATRIX
#EACH UPDATE IS A STOCHASTIC PROCESS, SO SEVERAL PERFORMANCES ARE NEEDED, IN ORDER TO MAKE STATISTICS
#ONLY UNDIRECTED GRAPHS WITH NO AUTO-LOOPS ARE DEALT WITH, AT THE MOMENT

import numpy as np #Used for its random functions and data structures
import networkx as nx #Used to extract information from networks, and plot them
import matplotlib.pyplot as plt #Used for graphic belluries
import pickle

### EVOLUTION

#First, an initial state is created, i.e. the adiacency matrix for T=0.
#In future development, this state will be the same for both evolutions.
#Then, it's initialized a 3D tensor, N*N*T, i.e. a sequence of adiacency matrices evolving in time.
#Each of the T adiacences matrices (A) will be updated according to the values of the previous one.
#Adiacences are kept simmetric with 0-diagonal.
#This scheme will be the same for both evolutions, just the law of updating will be different.

#In DAR, there's a probability alpha_ij of having A_ij(t) = A_ij(t-1)
#Otherwise, the new state will be extacted by a "coin toss", with probability xi_ij of getting "1"
#The evolution is performed just for the upper triangular adiancency, and than the matrix is simmetrized


def network_generation_dar(N,T,alpha,xi, P=1):
    #potresti aggiungere una variabile indicatrice tipo undirected = True
    #Genera un intero network di evoluzione
    #Initial state is generated randomly
    initial_state = [[np.random.choice((0,1),p=(0.5,1-0.5)) if j>i else 0 for j in range(N)] for i in range(N)]
    initial_state = [[initial_state[j][i] if j<i else initial_state[i][j] for j in range(N)] for i in range(N)]
    #Tensor generation, T=0 equal to initial_state
    temporal_network = np.zeros((T,N,N))
    temporal_network[0] = initial_state
    #Network update
    #First comprehension says: for each row i, taking only the upper triangle (j>i), choose the value 
    #from the same cell of the T-P matrix, with probability alpha, or from a coin tossing (whose probability of getting 1 is xi).
    #The second comprehension just simmetrizes the matrix
    for t in range(1,T):
        temporal_network[t] = [[np.random.choice((temporal_network[t-P,i,j],np.random.choice((1,0), p=(xi[i,j],1-xi[i,j]))),p=(alpha[i,j],1-alpha[i,j])) if j>i else 0 for j in range(N)] for i in range(N)]
        temporal_network[t] = [[temporal_network[t,j,i] if j<i else temporal_network[t,i,j] for j in range(N)] for i in range(N)]
    return temporal_network

# NETWORK ANALYSYS #
def degree_node(network,node,out = True):
    #you give an adiacency for a certain t, and a node: this function computes node degree
    #if out = true, sums the line (link that start from node)
    #if out = false, sums the column (link that reach that node)
    #if graph is undirected, out = in
    if out:
        return sum(network[node])
    else:
        return sum(network[:,node])

def degree_mean_t(network,out = True):
    #this function is not "pure", since it relies on another function
    #network must be a N*N np.array
    #you give an adiacency for a certain t: it computes mean degree at that t, using degree_node
    #depending in out value, it returns out degree or in degree. If undirected, it's the same
    degrees = []
    for node in range(len(network)):
        degrees.append(degree_node(network,node,out))
    return(np.average(degrees)) #it returns a number

def degree_mean_sequence(network,T,out = True):
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
        return pickle.dump()

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

### RESULTS PRINT AND SAVE:
start_name = 'Examples/DAR' #begin of the name of files that will be saved (txt,img...)
degree_evolution = [] #here the sequence of mean degrees will be held
    #Since it's not possible to directly handle the mean degree, here it is computed:
    degree_evolution.append(sum([grafo.degree(i) for i in range(N)])/N)
        #nx.draw(grafo, pos = nx.drawing.layout.circular_layout(grafo), with_labels = True)
    #plt.savefig(start_name+'%i_N%i_image_t=%i.pdf' %(P,N,t))
    #plt.show()
print(degree_evolution)

### TESTS:
def test_symmetry(temporal_network,T):
    for t in range(T):
        assert (temporal_network[t] == temporal_network[t].T).any, "Error: network at t = %i is not symmetric" %t

def test_nulldiagonal(temporal_network,T):
        assert sum(np.diag(temporal_network[t])) == 0, "Error: network at t = %i has not-0 diagonal" %t

#Save the network for further use
#np.savetxt(start_name+'%i_N%i_wholenetwork_T%i.txt' %(P,N,T), temporal_network.reshape(T*N*N,1))
#To import:
#new_data = np.loadtxt('start_name+'%i__N%i_wholenetwork_T%i.txt' %(P,N,T))
#new_data = new_data.reshape((T,N,N))