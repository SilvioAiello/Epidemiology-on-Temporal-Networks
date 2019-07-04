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
beta = 0.005 #infection rate (probability of infecting a node within unit time)
fig_count = 0 #several figures may be printed


temporal_dar = np.loadtxt('Examples/DAR'+'%i_N%i_wholenetwork_T%i.txt' %(P,N,T))
temporal_dar = temporal_dar.reshape((T,N,N)) #salvatolo come unica colonna, confermo che vuole tempo righe colonne
temporal_fitn= np.loadtxt('Examples/FITN'+'_N%i_wholenetwork_T%i.txt' %(N,T))
temporal_fitn= temporal_fitn.reshape((T,N,N))


###                      CONTAINERS AND FUNCTIONS FOR EPIDEMIC           ###

# CONTAINERS #
#Initialization of dictionaries that keep track of node states evolution. 
#This will be the syntax: label_[indexcase][repetition][instant]. If you set these 3 params, you get a dict of nodes state for that index, that repetition, at that time.
label_dar = []
label_fitn= []
    
#These lists save the a number of time steps needed to infect 60% of the network, for each node, for each iteration (they will be a list of N lists with K items):
score_dar = []
score_fitn =[]
#In order to perform a direct comparison with structural centralities, average virulence score is computed for each node,
#And virulence ranking will be stored in a list of ordered nodes, from those with lowest to those with highest average time to infect.
avg_dar = []
avg_fitn = []


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
    #UN MODO PER TESTARLA E' VEDERE CHE AD OGNI ISTANTE SUCCESSIVO E' DIVERSO DAL PRECEDENTE, PER LA STESSA COPPIA
    #This function computes the duration of a contact (in number of temporal steps) for a couple of nodes I-S, as long as I is infected (otherwise, it couldn't propagate the epidemic)
    #This is accomplished by checking backwards the existence of the link and the state of the I node, increasing the value of a counter variable until these conditions are satisfied
    counter = 0
    for instant in range(t+1):
        if (adiacency[t-instant,infected_node,susceptible_node] == 1 and state[t-instant][infected_node]==1 and state[t-instant][susceptible_node]==0):
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
        for i in infectneighbourhood: #then, for each infected neighbour, perform the extraction, according to link duration
                if infect_extraction(probabilities[contact_lasting(adiacency,state,t,i,s)]):
                    nextstate[s] = 1 #if successful, change the state of the node, at next t
                if nextstate[s] == 1: 
                    break #in order to avoid to scan all the nodes, if state is already changed, break from the loop
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
    score_dar.append([])
    label_dar.append([]) #create the i-th entry of label dar, which has K lists of T dictionaries, and each of the K has the same dict for t0, i.e. "node i index case"
    initial_state = dict.fromkeys(range(N),0) #fromkeys requires a tuple of keys and one value
    initial_state[i] = 1 #make node i the index case
    #Now, everything is ready to perform propagation, K times
    for k in range(K): #repeat K time for the same network (same intial condition, different outcome since propagation is stochastic)
        label=dict() #label è un dict di dict che viene sovrascritto ad ogni k-esima iterazione
        label[0] = initial_state #import initial state
        for t in range(1,T):
            label[t] = propagation(temporal_dar,label,t-1) #perform progation for each time step (but not the first and the last)
        label_dar[i].append(label)
        score_dar[i].append(time_score(label,0.6)) #score updating, dovresti mettere un "if == T-1, allora non ci è arrivato"
    avg_dar.append(np.average(score_dar[i]))

#Ranking computation:
sorted_nodes_dar = np.argsort(avg_dar) #list of nodes, sorted by their score, from the best to the worst

print(time.time()-start)

# FITN #
#Working progress

#TEST
def new_infected(state2,state1):
    #This function takes two dictionaries, describing system state at two different time steps, and returns the new infecteds
    #It is used for testing
    new_infeced_nodes = []
    for i in range(N):
        if state2[i] != state1[i]:
            new_infeced_nodes.append(i)
    return new_infeced_nodes

def intersection(list1,list2):
    #Returns intersection of two lists. Is used in test to check that the new infected were neigbs of the infecteds
    temp = set(list2)  #use set to improve performance
    list3 = [value for value in list1 if value in temp] 
    return list3

#Actual test:
for indexcase in range(N):
    for rep in range(K):
        for t in range(1,T):
            for new in new_infected(label_dar[indexcase][rep][t],label_dar[indexcase][rep][t-1]):
                #se il nuovo infetto non aveva nemmeno un vicino infetto c'è un errore
                if len(intersection(neighbourhood(temporal_dar,t-1,new),onlyones(label_dar[indexcase][rep],range(N),t-1))) == 0: 
                    print("Test failed for indexcase = %i, repetion %i, at t = %i" %(indexcase,rep,t))
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
    print(sorted_nodes_dar[i], avg_dar[sorted_nodes_dar[i]])
print("Top 10-ranked Broadcast centrality, and their scores:")
for i in range(10):
    print(nodes_Brank_dar[i], nodes_Bcentrality_dar[nodes_Brank_dar[i]])
print("Common nodes between infective and BC:")
#Function intesection shows the common top-ranking nodes, but lists have to be converted in sets
print(list(set(sorted_nodes_dar[0:9]).intersection(set(nodes_Brank_dar[0:9])))) 
print("")
print("")

#For FITN it will be the same

### PLOT
plt.figure(fig_count)
plt.plot(range(T),[infected_counter(label_dar[99][k][t]) for t in range(T)]) #così, messo qui, te lo mostra solo per l'ultimo nodo del ciclo
plt.xlabel("Time step")
plt.ylabel("Percentage of infected nodes")
plt.title("Node 99, iteration %i, beta = %.3f" %(k,beta))
plt.show()
fig_count+=1

x = np.zeros(N)
y = np.zeros(N)
for i in range(N):
    x[i] = avg_dar[i]/max(avg_dar)
    y[i] = nodes_Bcentrality_dar[sorted_nodes_dar[i]]/max(nodes_Bcentrality_dar)

plt.figure(fig_count)
plt.scatter(x,y)
plt.xlabel("Virulenza")
plt.ylabel("Centralità")
plt.show()
fig_count+=1
###SAVINGS

#Since we're dealing with dict of dict, it's better to save the results using 

#start_name = 'Examples/SI_EPIDEMIC_' #begin of the name of files that will be saved (txt,img...)
#np.savetxt(start_name+'DAR1_N%i_T%i.txt' %(N,T), temporal__dar.reshape(T,N*N))
#np.savetxt(start_name+'_LABELS_DAR1_N%i_T%i.txt' %(N,T), label_dar) #saving labels in a separate file
#np.savetxt(start_name+'FITN_N%i_T%i.txt' %(N,T), temporal__fitn.reshape(T,N*N))
#np.savetxt(start_name+'_LABELS_FITN_N%i_T%i.txt' %(N,T), label_fitn) #saving labels in a separate file