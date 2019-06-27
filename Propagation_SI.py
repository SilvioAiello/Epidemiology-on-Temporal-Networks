#THIS MODULE MAKES A DISEASE SPREADING, OVER PREVIOUSLY GENERATED TEMPORAL NETWORK (DAR AND FITN), IN SI MODE.
#IF YOU DON'T KNOW WHAT ARE FITNESS/DAR, OR THE SI PROPAGATION, CHECK THE DOCUMENTATION AND THE "EVOLUTION" PYTHON MODULES IN THIS FOLDER.
#THE NETWORK IS DESCRIBED THROUGH IT'S ADIACENCY MATRIX. ITS DESEASE STATE (SUSCEPTIBLE, INFECTED) IS A LABEL FOR EACH NODE.
#SUSCEPTIBLE -> 0, INFECTED -> 1
#EACH UPDATE OF THE DESEASE IS A STOCHASTIC PROCESS, SO SEVERAL PERFORMANCES ARE NEEDED, IN ORDER TO MAKE STATISTICS
#ONLY UNDIRECTED GRAPHS WITH NO AUTO-LOOPS ARE DEALT WITH, AT THE MOMENT

import numpy as np #Used for its random functions and data structures
import networkx as nx #Used to extract information from networks, and plot them
import matplotlib.pyplot as plt #Used for graphic belluries

#TODO: Spiegare

# TABLE OF CONTENTS: after initial arrangements, two sections develop and perform SI propagations, and two sections develop and perform centrality measures.
# Final section makes comparisons and prints result. Note that centrality measures are indipendent of the desease, and may be performed before its propagation.
# 1) Set of parameters and network import
# 2) Containers and functions for epidemic
# 3) Epidemic propagation
# 4) Containers and functions for centrality measures
# 5) Centrality measures
# 6) Comparisons and results print

###                     PARAMETERS OF THE SYSTEM & NETWORKS IMPORT         ###

#Here you can choose number of nodes, evolution length, number of iterations, DAR type, beta
N = 100 #number of nodes of the network
T = 100 #number of steps of temporal evolution
K = 5 #number of repetitions 
P = 1 #DAR(P)
beta = 0.015 #infection probability
fig_count = 0 #several figures will be printed


temporal_dar = np.loadtxt('Examples/DAR'+'%i_N%i_wholenetwork_T%i.txt' %(P,N,T))
temporal_dar = temporal_dar.reshape((T,N,N)) #salvatolo come unica colonna, confermo che vuole tempo righe colonne
temporal_fitn= np.loadtxt('Examples/FITN'+'_N%i_wholenetwork_T%i.txt' %(N,T))
temporal_fitn= temporal_fitn.reshape((T,N,N))


###                      CONTAINERS AND FUNCTIONS FOR EPIDEMIC           ###

# CONTAINERS #
#The following (K,T,N) matrices will keep track of node-state evolution in time, so the entry kij means "state of node j, at time i and for the k-th iteration":
#Quindi devi anticipare, da qualche parte sopra, che si itererà tutto K volte
label_dar = np.zeros((T,N)) #infective state for DAR-network
label_fitn = np.zeros((T,N))#infective state for FITN-network

#This array saves the number of time steps needed to infect 60% of the network, for each of K iterations, for each node:
score_dar = np.zeros((N,K))
score_fitn = np.zeros((N,K))
#The average result, for each node, is stored in a proper array
avg_score_dar = np.zeros(N)
avg_score_fitn = np.zeros(N)
#This one will rank the nodes according to their scores. So, he will be involved in the comparison with the centralities:
ranked_nodes_dar = np.zeros(N)
ranked_nodes_fitn = np.zeros(N)


# EPIDEMIC FUNCTIONS BUILDING #
#The idea of SI model is that one infected node, at t-1, can infect one of its neighboroughs with probability beta
#The infection, appens as in Chen -> SPIEGARE

from scipy.integrate import quad
def poisson_probability(t):
    #This function does.......
    lam = -np.log(1-beta) #applico la formula del Chen (pag 17), ma invece di 60 faccio direttamente 1
    return(lam*np.exp(-lam*t))
I = quad(poisson_probability,0,np.inf) #questo è per verificare che integrato all'inf fa 1

#THIS FUNCTION COMPUTES THE DURATION OF CONTACT (number of temporal steps) FOR A COUPLE OF NODES
def contact_lasting(adiacency,t,i,j):
    counter = 0
    for instant in range(t):
        if adiacency[t-instant,i,j] == 1:
            counter +=1
        else:
            break
    return counter

def infect_extraction(adiacency,t,i,j): #sarebbe bastato il solo t ma voglio far funzionare tutto veloce
    # Potrei farla più garbata definendo status, così ritorna un solo valore
    #This functions does.....
    x = np.random.uniform(0,1) #soglia di probabilità da superare
    if quad(poisson_probability,0,contact_lasting(adiacency,t,i,j))[0] > x:
        return True #se l'integrale è maggiore, allora dai l'ok
    else:
        return False

#Since the spread will be performed two times, on two different network, a function can be useful:
#This function takes the neighboroughs of all infected nodes, and infects them with probability beta
#It also updates the number of infections caused by that node
def propagation(adiacency,state,t): #quello che gli darai in pasto è t-1, e lui ti calcola lo stato a t
    nextstate = np.zeros(N)
    for i in range(N):
        if state[t,i]==1:
            nextstate[i] = 1 #next state-vector is state be initially equal to the previous, so who is infected will stay infected
    infected = [z for z in range(N) if state[t,z]==1]
    for i in infected: #for each infected node...
        neighb = [nodes for nodes in range(N) if adiacency[t,i,nodes] == 1] #...TAKE THE VICINI, VEDIAMO SE LO SONO...
        for n in neighb: # Then, for each neighb...
            if nextstate[n] == 0: #...if it is (still!) susceptible...
                if infect_extraction(adiacency,t,i,n) == True: #...and the extraction is successful...
                    nextstate[n] = 1 #...infect the node
                    #print("New infection! %i infected by %i" %(n,i))
    return(nextstate)

#This function counts the number of infected nodes at a certain time:
#Ma poi scusa, il numero infetti non è semplicemente sum(label_dar[t])?? Perché complicarsi la vita?
def infected_counter(infected):
    counter = 0
    for i in range(N):
        if infected[i]==1:
            counter+=1
    return counter

#Next function returns the time step at wich the disease has reached 60% of the network, if it exist:
def score_60(scores):
    #asserire che sum deve essere sempre <= N?
    time_60 = T-1 #initialized as the final temporal step
    for t in range(T):
        if sum(scores[t])>=0.6*N:
            time_60 = t
            break
    return time_60


###                         PROPAGATIONS                         ###
#Here the actual propagation occurs. First, the initial states are set for both networks, then propagation-function is called
#Rembember that, while for DAR there's only 1 initial infected, there will be more for FITN.

# DAR #
import time
for i in range(N):
    initial_state = np.zeros(N) 
    initial_state[i] = 1 #node i is the "index case" ("paziente zero"), i.e. the first infected node
    #Now, everything is ready to perform propagation, K times
    for k in range(K):
        start = time.time()
        label_dar[0] =initial_state
        for t in range(1,T): 
            label_dar[t] =propagation(temporal_dar,label_dar,t-1)
        score_dar[i,k] = score_60(label_dar) #score updating, dovresti mettere un "if == T-1, allora non ci è arrivato"
        print(time.time()-start)
    avg_score_dar[i] = np.mean(score_dar[i])


#Ranking computation:
# It's used argsort, which takes a vector and returns a vector of the input vector indices, in increasing order. 
# So, the first entry of the following output is the node who took less time to infect; the last entry is the worst node.
ranked_nodes_dar = np.argsort(avg_score_dar)


# FITN #
#Working progress


###                      CONTAINERS AND FUNCTIONS FOR CENTRALITIES      ###

# COMMUNICABILITY #
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
    Q = np.identity(N)/np.linalg.norm(np.identity(N))
    for t in range(T):
        inv = np.linalg.inv(np.identity(N)-0.25*maxradius*temporal_[t]) #inverse factor, which has to be multiplicated to the previous Q
        Q = np.matmul(Q,inv)/np.linalg.norm(np.matmul(Q,inv)) #updating and normalizing of Q
    return(maxradius,Q) #just for sake of completeness, also the max spectral radius is returned

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

###                      CENTRALITIES MEASURES                      ###
# DAR #
spec_radius_dar, Q_dar = communicability(temporal_dar)
nodes_Bcentrality_dar, nodes_Brank_dar = broadcast_ranking(Q_dar) #scores, node rankings
nodes_Rcentrality_dar, nodes_Rrank_dar = receive_ranking(Q_dar) #cores, node rankings

# FITN #
spec_radius_fitn, Q_fitn = communicability(temporal_fitn)
nodes_Bcentrality_fitn, nodes_Brank_fitn = broadcast_ranking(Q_fitn)
nodes_Rcentrality_fitn, nodes_Rrank_fitn = receive_ranking(Q_fitn)



###                     RESULTS PRINT                               ###
print("###   DAR NETWORK   ###")
print("Top 10- infective nodes, and their scores:")
for i in range(10):
    print(ranked_nodes_dar[i], avg_score_dar[ranked_nodes_dar[i]])
print("Top 10-ranked Broadcast centrality, and their scores:")
for i in range(10):
    print(nodes_Brank_dar[i], nodes_Bcentrality_dar[nodes_Brank_dar[i]])
print("Common nodes between infective and BC:")
#Function intesection shows the common top-ranking nodes, but lists have to be converted in sets
print(list(set(ranked_nodes_dar[0:9]).intersection(set(nodes_Brank_dar[0:9])))) 

print("")
print("")

#For FITN it will be the same


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