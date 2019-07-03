#THIS MODULE PREFORMS DISEASES SPREADING, OVER PREVIOUSLY GENERATED TEMPORAL NETWORKS (DAR AND FITN), IN SI MODE.
#IF YOU DON'T KNOW WHAT ARE FITNESS/DAR, OR THE SI PROPAGATION, CHECK THE DOCUMENTATION AND THE "EVOLUTION" PYTHON MODULES IN THIS FOLDER.

#THE MAIN PURPOSE IS TO EXTRACT AN INDEX OF CORRELATION BETWEEN NODES VIRULENCE AND CENTRALITY.
#NODE VIRULENCE IS DETERMINED BY MAKING EACH NODE THE INDEX CASE (PAZIENTE 0), AND COMPUTING THE TIME NEEDED TO INFECT A CERTAIN % OF THE NETWORK.
#EACH EPIDEMIC SPREAD IS A STOCHASTIC PROCESS, SO SEVERAL (K) SPREAD ITERATIONS ARE NEEDED, WITH THE SAME INITIAL CONDITION, OVER THE SAME NETWORK.
#TEMPORAL NETWORKS ARE THEMSELVES STOCHASTIC, SO THE WHOLE SET OF ITERATIONS HAS TO BE REPEATED OVER SEVERAL TEMPNETS.
#AFTER THIS, ONE CAN DETERMINE IF THERE'S AN ACTUAL CORRELATION BETWEEN VIRULENCE AND CENTRALITIES SCORES.

#A TEMPORAL NETWORK OF N NODES IS DESCRIBED BY T ADIACENCY MATRICES (N x N), WHICH HAVE THE SCTRUCTURE OF NUMPY ARRAY. 
#DESEASE STATE IS A DICTIONARY THAT MAPS EACH NODE TO 0 (SUSCEPTIBLE) OR 1 (INFECTED).

#ONLY UNDIRECTED GRAPHS WITH NO AUTO-LOOPS ARE DEALT WITH, AT THE MOMENT

import numpy as np #Used for its random functions and data structures
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
K = 50 #number of repetitions 
P = 1 #DAR(P)
beta = 0.1 #infection rate (probability of infecting a node within unit time)
fig_count = 0 #several figures may be printed


temporal_dar = np.loadtxt('Examples/DAR'+'%i_N%i_wholenetwork_T%i.txt' %(P,N,T))
temporal_dar = temporal_dar.reshape((T,N,N)) #salvatolo come unica colonna, confermo che vuole tempo righe colonne
temporal_fitn= np.loadtxt('Examples/FITN'+'_N%i_wholenetwork_T%i.txt' %(N,T))
temporal_fitn= temporal_fitn.reshape((T,N,N))


###                      CONTAINERS AND FUNCTIONS FOR EPIDEMIC           ###

# CONTAINERS #
#Initialization of dictionaries that keep track of node states evolution. 
#Each key, representing the instant, has as value a dictionary, that maps each node to the state at that instant.
label_dar = dict()
label_fitn= dict()
    
#These dicts save the average (over K iterations) number of time steps needed to infect 60% of the network, for each node:
score_dar = dict.fromkeys(range(N),0)
score_fitn =dict.fromkeys(range(N),0)

#These lists rank the nodes according to their scores. So, they will be involved in the comparisons with centralities:
sorted_nodes_dar = []
sorted_nodes_fitn = []


# EPIDEMIC FUNCTIONS BUILDING #
#The idea of SI model is that each I node, at t-1, can make infect one of its S neighbours, at t-1, with a certain probability, whose rate is beta.
#So, there will be usefuls some simple functions that build the neighbourhood of a node, or find S/I nodes from a list at a given time, or compute links' temporal duration.

def neighbourhood(adiacency,t,node):
    #This function takes and adiacency matrix (dar, tgrg...), a particoular time step and a given node
    #And returns a SET with all nodes that have a link with that node at that time
    neigh = {i for i in range(N) if adiacency[t,node,i]==1}
    return neigh

def onlyones(state,nodes_list,t):
    #This function takes a list of nodes (not necessarily the whole set of nodes), whole nodes-state dict and a given time
    #And return a SET with just infected nodes (i.e. those whose state is 1).
    selected = {node for node in nodes_list if state[t][node]==1}
    return selected

def onlyzeros(state,nodes_list,t):
    #Same as onlyones, but returns suceptible nodes.
    selected = {node for node in nodes_list if state[t][node]==0}
    return selected

def contact_lasting(adiacency,state,t,infected_node,susceptible_node):
    #This function computes the duration of a contact (in number of temporal steps) for a couple of nodes I-S, as long as I is infected (otherwise, it couldn't propagate the epidemic)
    #This is accomplished by checking backwards the existence of the link and the state of the I node, increasing the value of a counter variable until these conditions are satisfied
    counter = 0
    for instant in range(t+1):
        if (adiacency[t-instant,infected_node,susceptible_node] == 1 and state[t-instant][infected_node]==1):
            counter +=1
        else:
            break
    return counter #this should be included in [0,t]


#Epidemic spread follows Chen approach: time of infection (in unit steps) follows a Poissonian distribution, normalized to return beta for 1 step, integrated within link duration.
#(note: beta is the probability rate of contagion [1/s], but also the actual probability of contagion after 1 unit time: infact, P(1) = beta*1 u.t. = beta [dimensionless]).
#Chen's algorithm make a contagion happen by performing an extraction from Unif(0,1): if this number is lower than the Poisson integral, contagion takes place at that time.

from scipy.integrate import quad #this function performs integration
def poisson_probability(t):
    #This function reproduces the Poisson PDF, where the average is set in function of beta, who is a parameter of the disease
    lam = -np.log(1-beta) #See Chen (who uses 60)
    return(lam*np.exp(-lam*t))
#I = quad(poisson_probability,0,np.inf) #verify the correct normalization
#The same integrals will be performed many times, varying at most the total duration, whose maximum value can be T.
#So, they are computed here once for all, and the results stored in a dict, that maps each duration to the integral
probabilities = dict()
for t in range(T):
    probabilities[t] = quad(poisson_probability,0,t)[0]


def infect_extraction(probab): #Remember: the outcome of this function is stochastic
    #This function extracts a random number from Uniform Distribution, and compares it to the given probability
    #If the random number is lower, function return True, to signal that contagion can occur
    x = np.random.uniform(0,1) #soglia di probabilità da superare
    if probab > x:
        return (True)
    else:
        return (False)

#This function performs the actual propagation. It requires the adiacency and the states at t, and return the states at t+1
#It finds the S neighboroughs of all I nodes at a given time, and makes them infect at t+1, with probability beta
#This function will be evoked for each time step but the first and the last one
def propagation(adiacency,state,t): #Remember: the outcome of this function is stochastic
    #First, it is inialized the state-dict at t+1, equal to the previous one, so who is infected stays infected
    nextstate = dict.fromkeys(range(N),0)
    for i in range(N):
        if state[t][i]==1:
            nextstate[i] = 1

    #Then, it finds susceptible nodes...
    susc = onlyzeros(state,range(N),t) #un modo per migliorare ancora sarebbe crearla una tantum e rimuovere volta per volta
    #...and for each of them, it performs the extraction, evaluating the whole neighbourhood
    for s in susc:
        infectneighbourhood = onlyones(state,neighbourhood(adiacency,t,s),t) #takes the infected neigbs. of that node
        if len(infectneighbourhood): #if the set is empty, there can't be infection
            for i in infectneighbourhood:
                p = 0
                p += probabilities[contact_lasting(adiacency,state,t,i,s)] #compute the total probability due to neighbourhood
            p /= len(infectneighbourhood) #normalize
            if infect_extraction(p): #evoke the actual extraction
                nextstate[s] = 1 #if successful, change the state of the node, at next t
    return(nextstate)

#This function counts the number of infected nodes in a set:
def infected_counter(infected):
    counter = 0
    for i in range(N):
        if infected[i]==1:
            counter+=1
    return counter

#This function, provided of the whole states-evolution, returns the time step at wich the disease has reached a given fraction of the network, using infected_counter.
#If that fraction has not been reached, it return the total time of propagation:
def time_score(scores,fraction):
    #Scores sarebbe il label dar totale, e lui deve trovare in quale label_dar[t] la soglia è stata superata
    assert fraction > 0, "Error, only values between 0 and 1 are allowed"
    assert fraction < 1, "Error, only values between 0 and 1 are allowed"
    #asserire che sum deve essere sempre <= N?
    time_spent = T-1 #initialized as the final temporal step
    for t in range(T):
        if infected_counter(scores[t])>=fraction*N:
            time_spent = t
            break
    return time_spent #ha senso restituire il total time? magari inf, magari errore, vediamo


###                         PROPAGATIONS                         ###
#Here the propagation is performed. First, the initial states are set for both networks, then function is called
#Initial state is: all nodes Susceptible but one. Each node at a time is the index case, and for each there are K iterations

import time
start = time.time()

# DAR #
for i in range(N): #do it for each node
    print("Processing node %i" %i)
    initial_state = dict.fromkeys(range(N),0) #fromkeys vuole una tupla e un valore
    initial_state[i] = 1 #make node the index case
    #Now, everything is ready to perform propagation, K times
    for k in range(K): #repeat K time for the same network (same intial condition, different outcome since propagation is stochastic)
        label_dar[0] = initial_state #import initial state
        for t in range(1,T):
            label_dar[t] =propagation(temporal_dar,label_dar,t-1) #perform progation for each time step (but not the first and the last)
        score_dar[i] += time_score(label_dar,0.6) #score updating, dovresti mettere un "if == T-1, allora non ci è arrivato"
    score_dar[i] /= K #occhio al fatto di T-1 se non è raggiunta la percentuale!!!

#Ranking computation:
sorted_nodes_dar = sorted(score_dar.keys(), key= score_dar.get) #list of nodes, sorted by their score

print(time.time()-start)

plt.plot(range(T),[infected_counter(label_dar[t]) for t in range(T)]) #così, messo qui, te lo mostra solo per l'ultimo nodo del ciclo
plt.xlabel("Time step")
plt.ylabel("Percentage of infected nodes")
plt.title("Node 99, iteration %i" %k)
plt.show()
# FITN #
#Working progress
#%%
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
    rank = np.flip(np.argsort(lines_sum)) #argsort returns the list of nodes the increasing order (from the lowest centrality to the highset)
    return(lines_sum,rank)
def receive_ranking(Q):
    #Receive = sum over columns:
    lines_sum = np.sum(Q, axis = 0) #this is a vector which reports the RC for each node
    #Using argsort (flipping the result) to get the ranks vector: rank_i means the i-th best node
    rank = np.flip(np.argsort(lines_sum)) #argsort returns the list of nodes the increasing order (from the lowest centrality to the highset)
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
    print(sorted_nodes_dar[i], score_dar[sorted_nodes_dar[i]])
print("Top 10-ranked Broadcast centrality, and their scores:")
for i in range(10):
    print(nodes_Brank_dar[i], nodes_Bcentrality_dar[nodes_Brank_dar[i]])
print("Common nodes between infective and BC:")
#Function intesection shows the common top-ranking nodes, but lists have to be converted in sets
print(list(set(sorted_nodes_dar[0:9]).intersection(set(nodes_Brank_dar[0:9])))) 

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