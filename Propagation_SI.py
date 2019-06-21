#THIS MODULE MAKES A DISEASE SPREADING, OVER PREVIOUSLY GENERATED TEMPORAL NETWORK (DAR AND FITN), IN SI MODE.
#IF YOU DON'T KNOW WHAT ARE FITNESS/DAR, OR THE SI PROPAGATION, CHECK THE DOCUMENTATION AND THE "EVOLUTION" PYTHON MODULES IN THIS FOLDER.
#THE NETWORK IS DESCRIBED THROUGH IT'S ADIACENCY MATRIX. ITS DESEASE STATE (SUSCEPTIBLE, INFECTED) IS A LABEL FOR EACH NODE.
#SUSCEPTIBLE -> 0, INFECTED -> 1
#EACH UPDATE OF THE DESEASE IS A STOCHASTIC PROCESS, SO SEVERAL PERFORMANCES ARE NEEDED, IN ORDER TO MAKE STATISTICS
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


### NETWORKS IMPORT
temporal_network_dar = np.loadtxt('Examples/DAR'+'%i__N%i_wholenetwork_T%i.txt' %(P,N,T))
temporal_network_dar = temporal_network_dar.reshape((T,N,N))
temporal_network_fitn= np.loadtxt('Examples/FITN'+'_N%i_wholenetwork_T%i.txt' %(N,T))
temporal_network_fitn= temporal_network_fitn.reshape((T,N,N))


### EPIDEMIC FUNCTION BUILDING
#The idea of SI model is that one infect node, at t-1, can make infect one of its neighboroughs at t-1, with probability beta
beta = 0.15
#The following matrices (T,N) will keep track of node-labels evolution in time, so the entry ij means "state at time i for node j":
label_dar = np.zeros((T,N))
label_fitn = np.zeros((T,N))
#They are both initialized randomly, according to beta
label_dar[0]=np.random.choice((0,1),p=(1-beta,beta), size=N)
label_fitn[0]=np.random.choice((0,1),p=(1-beta,beta), size=N)

#Since the spread will be performed two times, on two different network, a function can be useful:
def propagation(adiacency,state):
    #This function takes the neighboroughs of all infected nodes, and infects them with probability beta
    nextstate = state #next state vector is initialized as equal to the previous, so who is infected will stay infected
    graph = nx.from_numpy_matrix(adiacency) #make adiacency a nx graph, so finding neighboroughs is straightforward
    for i in range(N): #for each node...
        if state[i] == 1: #... if the node is infected...
            neighb = list(nx.neighbors(graph, i)) #...take its neighb. and put them in a list.
            for n in neighb: # Then, for each neighb...
                if state[n] == 0: #...if it is susceptible...
                    nextstate[n] = np.random.choice((0,1),p=(1-beta,beta)) #...infect him with probabilty beta
    return(nextstate)

### EPIDEMIOLOGY ON DAR NETWORK:
for t in range(T-1): #updating label state for both networks, providing the previous state
    label_dar[t+1] = propagation(temporal_network_dar[t],label_dar[t])
    label_fitn[t+1]=propagation(temporal_network_fitn[t],label_fitn[t])
 
### RESULTS PRINT
#To get some "readable" graphics, one can tell nx to show in red infected nodes, and in blue the others
#So, here is a functions that thakes the vector of state at time t and return a sequence of colours
def colorstate(state):
    colormap = []
    for node in range(N):
        if state[node] == 1:
            colormap.append("red")
        elif state[node]==0:
            colormap.append("blue")
    return colormap


#Plot DAR
for t in range(T):
    #What to plot:
    graph = nx.from_numpy_matrix(temporal_network_dar[t])
    colors = colorstate(label_dar[t])
    #Plotting settings
    plt.figure(fig_count)
    fig_count+=1
    plt.title(r"DAR(%i) SI-Epidemiology,$\beta$ = %.2f, N=%i,T=%i" %(P,beta,N,T))
    nx.draw(graph, pos = nx.drawing.layout.circular_layout(graph),node_color = colors, with_labels = True)
    plt.show()

#Plot FITN
for t in range(T):
    #What to plot:
    graph = nx.from_numpy_matrix(temporal_network_fitn[t])
    colors = colorstate(label_fitn[t])
    #Plotting settings
    plt.figure(fig_count)
    fig_count+=1
    plt.title(r"FITN SI-Epidemiology,$\beta$ = %.2f, N=%i,T=%i" %(beta,N,T))
    nx.draw(graph, pos = nx.drawing.layout.circular_layout(graph),node_color = colors, with_labels = True)
    plt.show()
    
### NETWORK ANALYSIS
    
start_name = 'Examples/SI_EPIDEMIC_' #begin of the name of files that will be saved (txt,img...)
#np.savetxt(start_name+'DAR1_N%i_T%i.txt' %(N,T), temporal_network_dar.reshape(T,N*N))
#np.savetxt(start_name+'_LABELS_DAR1_N%i_T%i.txt' %(N,T), label_dar) #saving labels in a separate file
#np.savetxt(start_name+'FITN_N%i_T%i.txt' %(N,T), temporal_network_fitn.reshape(T,N*N))
#np.savetxt(start_name+'_LABELS_FITN_N%i_T%i.txt' %(N,T), label_fitn) #saving labels in a separate file

#Communicability analysis
#Building of Communicability matrix, just by applying definition and using some useful np functions:
def communicability(temporal_network): 
    #As known, to compute communicability one as to choose a coefficient that multiplicates adiacencies
    #At the moment, this function takes as default, as coefficient, a quarter of the inverse of max spectral radius
    #Find max spectral radius:
    spec = []
    for t in range(T):
        spec.append(np.real(max(np.linalg.eigvals(temporal_network[t])))) #find max eigenvalue for each adiacency, taking only RE
    maxradius = 1/max(spec) #maximum is the reciprocal of the maximum eigenvalue
    #Communicability builing
    Q = np.identity(N) 
    for t in range(T):
        inv = np.linalg.inv(np.identity(N)-0.25*maxradius*temporal_network[t]) #inverse factor, which has to be multiplicated to the previous Q
        Q = np.matmul(Q,inv)/np.linalg.norm(np.matmul(Q,inv)) #updating and normalizing of Q
    return(maxradius,Q) #just for sake of completeness, also tha max spectral radius is returned

#Next two functions compute broadcast/receiving centralities and node rankings
#For centralities, they use function np.sum, where one can choose to sum of lines (BC) or columns (RC)
#For rankings, they use np.argsort, whose input is a list and output is a list of the indices of the input, sorted according to their decreasing values
#So, the first element of the output list has the highest rank
def broadcast_ranking(Q):
    #Broadcast = sum over lines:
    lines_sum = np.sum(Q, axis = 1) #this is a vector which reports the BC for each node
    #Using argsort (flipping the result) to get the ranks vector: rank_i means the i-th best node
    rank = np.flip(np.argsort(lines_sum)) #argsort returns the increasing order
    return(lines_sum,rank)
def receive_ranking(Q):
    #Receive = sum over columns:
    lines_sum = np.sum(Q, axis = 0) #this is a vector which reports the RC for each node
    #Using argsort (flipping the result) to get the ranks vector: rank_i means the i-th best node
    rank = np.flip(np.argsort(lines_sum)) #argsort returns the increasing order
    return(lines_sum,rank)

#Everything is set to perform analysis over the imported netowrks:
#DAR
spec_radius_dar, Q_dar = communicability(temporal_network_dar)
nodes_Bcentrality_dar, nodes_Brank_dar = broadcast_ranking(Q_dar) #scores, node rankings
nodes_Rcentrality_dar, nodes_Rrank_dar = receive_ranking(Q_dar) #cores, node rankings
print("Top 3-ranked Broadcast centrality, and their score [DAR]:")
for i in range(3):
    print(nodes_Brank_dar[i], nodes_Bcentrality_dar[nodes_Brank_dar[i]])

#FITN
spec_radius_fitn, Q_fitn = communicability(temporal_network_fitn)
nodes_Bcentrality_fitn, nodes_Brank_fitn = broadcast_ranking(Q_fitn)
nodes_Rcentrality_fitn, nodes_Rrank_fitn = receive_ranking(Q_fitn)
print("Top 3-ranked Receive centrality, and their score [FITN]:")
for i in range(3):
    print(nodes_Brank_fitn[i], nodes_Bcentrality_fitn[nodes_Brank_fitn[i]])