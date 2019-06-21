#THIS MODULE ALLOWS TO UPDATE A NETWORK ACCORDING TO THE DAR(P) LAW OF EVOLUTION.
#IF YOU DON'T KNOW WHAT'S A DAR(P), CHECK THE DOCUMENTATION
#THE NETWORK IS DESCRIBED THROUGH IT'S ADIACENCY MATRIX
#EACH UPDATE IS A STOCHASTIC PROCESS, SO SEVERAL PERFORMANCES ARE NEEDED, IN ORDER TO MAKE STATISTICS
#ONLY UNDIRECTED GRAPHS WITH NO AUTO-LOOPS ARE DEALT WITH, AT THE MOMENT

import numpy as np #Used for its random functions and data structures
import networkx as nx #Used to extract information from networks, and plot them
import matplotlib.pyplot as plt #Used for graphic belluries

### PARAMETERS OF THE SYSTEM

N = 10 #number of nodes of the network
T = 10 #number of steps of temporal evolution
K = 100 #number of repetitions 
P = 1 #DAR(P)
fig_count = 0 #several figures will be printed

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

#At the moment, alpha will be the same for all entries; alpha = 0.66 and xi = 0.5:
alpha = 0.66*np.ones((N,N))
xi = 0.5*np.ones((N,N))

#Initial state is generated randomly
#np.random.seed(1)
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

### RESULTS PRINT AND SAVE:
start_name = 'Examples/DAR' #begin of the name of files that will be saved (txt,img...)
degree_evolution = [] #here the sequence of mean degrees will be held
for t in range(T):
    plt.figure(fig_count)
    fig_count+=1
    #Function "from_numpy_matrix" converts a np matrix in a Network, readable by nx library:
    grafo = nx.from_numpy_matrix(temporal_network[t])
    #Since it's not possible to directly handle the mean degree, here it is computed:
    degree_evolution.append(sum([grafo.degree(i) for i in range(N)])/N)
    #Function draw is provided of two informations: draw a circular plot and show the label (i.e. the number) of each node
    nx.draw(grafo, pos = nx.drawing.layout.circular_layout(grafo), with_labels = True)
    #plt.savefig(start_name+'%i_N%i_image_t=%i.pdf' %(P,N,t))
    plt.show()
print(degree_evolution)

### TESTS:
for t in range(T):
    assert (temporal_network[t] == temporal_network[t].T).any, "Error: network at t = %i is not symmetric" %t
    assert sum(np.diag(temporal_network[t])) == 0, "Error: network at t = %i has not-0 diagonal" %t
#Save the network for further use
np.savetxt(start_name+'%i_N%i_wholenetwork_T%i.txt' %(P,N,T), temporal_network.reshape(T*N*N,1))
#To import:
#new_data = np.loadtxt('start_name+'%i__N%i_wholenetwork_T%i.txt' %(P,N,T))
#new_data = new_data.reshape((T,N,N))