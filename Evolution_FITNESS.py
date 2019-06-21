#THIS MODULE ALLOWS TO UPDATE A NETWORK ACCORDING TO THE FITNESS LAW OF EVOLUTION.
#IF YOU DON'T KNOW WHAT'S A FITNESS EVOLUTION, CHECK THE DOCUMENTATION
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
fig_count = 0 #several figures will be printed

### EVOLUTION

#Fitness computation
#First, the evolution of fitnesses, which is a kind of random walk, is separately performed.
#It is stored in a matrix, theta, whose entry ij is the value of i-node's fitness at t = j.
#The random walk has, for each step, an increment extracted by a gaussian with 0-mean and 1-variance
phi0 = 0.27*np.ones(N) #intercepts, typical of each node
phi1 = 0.18*np.ones(N) #slopes, typical of each node
theta = np.zeros((N,T)) #fitness matrix creation

#So, theta iS inizialized and then the random walk is performed
theta[:,0] = np.random.normal(0,0.1,size=N)
for t in range(1,T): #costruzione dei random walks
    for i in range(N):
        theta[i,t] = phi0[i] + phi1[i]*theta[i,t-1] + np.random.normal(0,0.1)

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

### RESULTS PRINT AND SAVE:
start_name = 'Examples/FITN' #begin of the name of files that will be saved (txt,img...)
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
    #plt.savefig(start_name+'_N%i_image_t=%i.pdf' %(N,t))
    plt.show()
print(degree_evolution)

### TESTS:
for t in range(T):
    assert (temporal_network[t] == temporal_network[t].T).any, "Error: network at t = %i is not symmetric" %t
    assert sum(np.diag(temporal_network[t])) == 0, "Error: network at t = %i has not-0 diagonal" %t

#Save the network for further use
np.savetxt(start_name+'_N%i_wholenetwork_T%i.txt' %(N,T), temporal_network.reshape(T,N*N))
#To import:
#new_data = np.loadtxt(start_name+'_N%i_wholenetwork_T%i.txt' %(N,T))
#new_data = new_data.reshape((T,N,N))