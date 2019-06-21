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
N = 100 #number of nodes of the network
T = 100 #number of steps of temporal evolution
K = 100 #number of repetitions 
P = 1 #DAR(P)
fig_count = 0 #several figures will be printed


### NETWORKS IMPORT
temporal_dar = np.loadtxt('Examples/DAR'+'%i_N%i_wholenetwork_T%i.txt' %(P,N,T))
temporal_dar = temporal_dar.reshape((T,N,N)) #salvatolo come unica colonna, confermo che vuole tempo righe colonne
temporal__fitn= np.loadtxt('Examples/FITN'+'_N%i_wholenetwork_T%i.txt' %(N,T))
temporal__fitn= temporal__fitn.reshape((T,N,N))


### EPIDEMIC FUNCTIONS BUILDING
#The idea of SI model is that one infected node, at t-1, can infect one of its neighboroughs with probability beta
beta = 0.15
#Since the spread will be performed two times, on two different network, a function can be useful:
#This function takes the neighboroughs of all infected nodes, and infects them with probability beta
#It also updates the number of infections caused by that node
def propagation(adiacency,state,killed):
    nextstate = np.zeros(N)
    for i in range(N):
        if state[i]==1:
            nextstate[i] = 1 #next state-vector is state be initially equal to the previous, so who is infected will stay infected
    graph = nx.convert_matrix.from_numpy_matrix(adiacency) #make adiacency a nx-graph, so finding neighboroughs is straightforward
    infected = [z for z in range(N) if state[z]==1]
    #print(infected)
    for i in infected: #for each infected node...
        neighb = list(nx.neighbors(graph, i)) #...take its neighboroughs and put them in a list.
        #print("Nodo infetto # %i, i suoi vicini sono " %i + str(neighb))
        for n in neighb: # Then, for each neighb...
            if nextstate[n] == 0: #...if it is (still!) susceptible...
                nextstate[n] = np.random.choice((0,1),p=(1-beta,beta)) #...infect him with probabilty beta...
                if nextstate[n] == 1: #...and if there is infection... 
                    #print("New infection! %i infected by %i" %(n,i))
                    killed[i]+=1 #...update the infecting score of the i-th node
    #print("Giro concluso, nuovo label state:")
    #print(nextstate)
    #print("")
    return(nextstate,killed)

#It may be useful a function that counts the number of infected nodes at a certain time:
def infected_counter(infected):
    assert len(infected)==N, "Error: array has not length N"
    assert sum(infected)<=N, "Error: there are more infections then the actual number of nodes"
    counter = 0
    for i in range(N):
        if infected[i]==1:
            counter+=1
    return counter

### EPIDEMIOLOGY ON DAR NETWORK:
#The following (T,N)matrices will keep track of node-labels evolution in time, so the entry ij means "state, at time i, of node j":
label_dar = np.zeros((T,N))
label_fitn = np.zeros((T,N))
#They are both initialized randomly, according to beta:
label_dar[0]=np.random.choice((0,1),p=(1-beta,beta), size=N)

label_fitn[0]=np.random.choice((0,1),p=(1-beta,beta), size=N)
#Infective score is also created now; it will be used later to make comparisons with structural measures:
infect_score_dar = np.zeros(N)
infect_score_fitn = np.zeros(N)
#Actual propagation:
for t in range(T-1): #updating label state for both networks, providing the previous state
#    print("t = %i" %t)
#    print("Stato infettivo:")
#    print(label_dar[t])
    label_dar[t+1],infect_score_dar =propagation(temporal_dar[t],label_dar[t],infect_score_dar)
    label_fitn[t+1],infect_score_fitn=propagation(temporal__fitn[t],label_fitn[t],infect_score_fitn)
#Normalization of infection scores, so each i-th entry is the percentage of the network infected by the i-th node:
infect_score_dar = 100*infect_score_dar/N
infect_score_fitn = 100*infect_score_fitn/N
#Using flip sort (L'HO SPIEGATA PIU AVANTI) we obtain a list of most infective nodes
rankings_infections_dar = np.flip(np.argsort(infect_score_dar))
rankings_infections_fitn = np.flip(np.argsort(infect_score_fitn))


### PRINTS
print("Epidemiological results, DAR:")
print("Infection rate at 0: %.f %%" %(100*infected_counter(label_dar[0])/N))
print("Infection rate at T/2: %.f %%" %(100*infected_counter(label_dar[int(T/2)])/N))
print("Infection rate at T: %.f %%" %(100*infected_counter(label_dar[T-1])/N))
print("Best 10 infecting node:")
print(rankings_infections_dar[0:10])
print("and their percentages:")
print(infect_score_dar[rankings_infections_dar[0:10]])
print("")
print("")
print("Epidemiological results, FITN:")
print("Infection rate at 0: %.f %%" %(100*infected_counter(label_fitn[0])/N))
print("Infection rate at T/2: %.f %%" %(100*infected_counter(label_fitn[int(T/2)])/N))
print("Infection rate at T: %.f %%" %(100*infected_counter(label_fitn[T-1])/N))
print("Best 10 infecting node:")
print(rankings_infections_fitn[0:10])
print("and their percentages:")
print(infect_score_fitn[rankings_infections_fitn[0:10]])
### PLOTS
##To get some "readable" graphics, one can tell nx to show in red infected nodes, and in blue the others
##So, here is a functions that thakes the vector of state at time t and return a sequence of colours
#def colorstate(state):
#    colormap = []
#    for node in range(N):
#        if state[node] == 1:
#            colormap.append('r')
#        elif state[node]==0:
#            colormap.append('b')
#    return colormap
##Plot DAR
#for t in range(T):
#    #What to plot:
#    graph = nx.convert_matrix.from_numpy_matrix(temporal_dar[t])
#    colors = colorstate(label_dar[t])
#    #Plotting settings
#    plt.figure(fig_count)
#    fig_count+=1
#    plt.title(r"DAR(%i) SI-Epidemiology,$\beta$ = %.2f, N=%i,T=%i" %(P,beta,N,T))
#    nx.draw(graph, pos = nx.drawing.layout.circular_layout(graph),node_color = colors, with_labels = True)
#    plt.show()
#Plot FITN
#for t in range(T):
#    #What to plot:
#    graph = nx.convert_matrix.from_numpy_matrix(temporal__fitn[t])
#    colors = colorstate(label_fitn[t])
#    #Plotting settings
#    plt.figure(fig_count)
#    fig_count+=1
#    plt.title(r"FITN SI-Epidemiology,$\beta$ = %.2f, N=%i,T=%i" %(beta,N,T))
#    nx.draw(graph, pos = nx.drawing.layout.circular_layout(graph),node_color = colors, with_labels = True)
#    plt.show()

###SAVINGS
#start_name = 'Examples/SI_EPIDEMIC_' #begin of the name of files that will be saved (txt,img...)
#np.savetxt(start_name+'DAR1_N%i_T%i.txt' %(N,T), temporal__dar.reshape(T,N*N))
#np.savetxt(start_name+'_LABELS_DAR1_N%i_T%i.txt' %(N,T), label_dar) #saving labels in a separate file
#np.savetxt(start_name+'FITN_N%i_T%i.txt' %(N,T), temporal__fitn.reshape(T,N*N))
#np.savetxt(start_name+'_LABELS_FITN_N%i_T%i.txt' %(N,T), label_fitn) #saving labels in a separate file


### CENTRALITY MEASURES
    
#Communicability analysis
#Building of Communicability matrix, just by applying definition and using some useful np functions:
def communicability(temporal_): 
    #As known, to compute communicability one as to choose a coefficient that multiplicates adiacencies
    #At the moment, this function takes as default, as coefficient, a quarter of the inverse of max spectral radius
    #Find max spectral radius:
    spec = []
    for t in range(T):
        spec.append(np.real(max(np.linalg.eigvals(temporal_[t])))) #find max eigenvalue for each adiacency, taking only RE
    maxradius = 1/max(spec) #maximum is the reciprocal of the maximum eigenvalue
    #Communicability builing
    Q = np.identity(N) 
    for t in range(T):
        inv = np.linalg.inv(np.identity(N)-0.25*maxradius*temporal_[t]) #inverse factor, which has to be multiplicated to the previous Q
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
spec_radius_dar, Q_dar = communicability(temporal_dar)
nodes_Bcentrality_dar, nodes_Brank_dar = broadcast_ranking(Q_dar) #scores, node rankings
nodes_Rcentrality_dar, nodes_Rrank_dar = receive_ranking(Q_dar) #cores, node rankings
print("Top 3-ranked Broadcast centrality, and their score [DAR]:")
for i in range(3):
    print(nodes_Brank_dar[i], nodes_Bcentrality_dar[nodes_Brank_dar[i]])

#FITN
spec_radius_fitn, Q_fitn = communicability(temporal__fitn)
nodes_Bcentrality_fitn, nodes_Brank_fitn = broadcast_ranking(Q_fitn)
nodes_Rcentrality_fitn, nodes_Rrank_fitn = receive_ranking(Q_fitn)
print("Top 3-ranked Receive centrality, and their score [FITN]:")
for i in range(3):
    print(nodes_Brank_fitn[i], nodes_Bcentrality_fitn[nodes_Brank_fitn[i]])