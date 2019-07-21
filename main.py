"""
From this script, user can import and use functions from all the others, providing its own parameters.
Functions in this script work in Pyhon3, may require numpy (v1.16) and function "quad" from scipy.integrate (scipy v1.3).

In this module are defined the following functions:
    * poisson_probability

You can set your parameters in "USER ACTION" section.

For further understandings on how this script operates, check file "howto.md".
For further theoretical understandings, check file "explanation.md".
"""
import numpy as np
from scipy.integrate import quad #used in dictionary of probabilities

import Evolutions
import Propagation_SI
import Saves
#%%         USER ACTION
    #GENERAL PARAMETERS
N=10 #nodes
T=100 #duration
isDAR = True #type
isDIRECTED = False #network symmetry
    #NUMBER OF ITERATIONS
NET_REAL = 1 #tempnet
K = 5 #infection
    #ANALYSIS
net_name = 'CHOOSE NAME' #identification name
it_chosen = 1 #what network iteration to infect

    #DAR INPUTS
P = 1 #order
alpha = 0.6*np.ones((N,N)) #you can change values, but keep it (N,N)
xi = 0.5*np.ones((N,N)) #you can change value, but keep it (N,N)
    #TGRG INPUTS
phi0 = 0.5*np.ones(N) #change values, but keep it (N)
phi1 = 0.5*np.ones(N) #change values, but keep it (N)
epsilon=0.5*np.ones(N)#change values, but keep it (N)

    #EPIDEMIC PARAMETERS
beta = 0.005 #infection rate

#%% OUTPUTS GENERATION
    #preliminar assertions
assert NET_REAL >= 1, "NET_REAL should be >=1"
assert K > 1, "K should be >1"
assert it_chosen >=1, "it_chosen should be >=1"
assert it_chosen <=NET_REAL, "Realization not found"

    #TEMPNETS GENERATION
for k in range(1,NET_REAL+1):  #so first realization has index 1
    if isDAR: #use the proper functiond wheter user selected dar or tgrg
        temp = Evolutions.network_generation_dar(alpha,xi,P=P,T=T,directed=isDIRECTED) 
    else:
        temp = Evolutions.network_generation_tgrg(alpha,xi,P=P,T=T,directed=isDIRECTED) 
    Saves.network_save(temp,net_name, isDAR = isDAR, k=k, P=1)

    #Network to analyze
temporal_network = Saves.network_load(N=N,T=T,start=net_name, k=it_chosen)

    #SI PROPAGATION
    #Probabilities dict
probabilities = dict() #dict building
for t in range(T):
    probabilities[t] = quad(Propagation_SI.poisson_probability,0,t, args = beta)[0] #quad produces several outputs, integral is the first

    #Function evoking
label = [] 
for index_case in range(N):
    label.append([]) #create the i-th entry
    for iteration in range(K):
        label[index_case].append(Propagation_SI.propagation(temporal_network, index_case, probabilities))
        # OUTPUT TESTS
        assert label[index_case][iteration][0][index_case] == 1, "An index case appears to be uninfected"
        assert sum(label[index_case][iteration][0].values()) == 1, "There should be only 1 infect at the beginning"
assert [[label[index_case][iteration][0] == label[index_case][iteration-1][0] for iteration in range(1,K)] for index_case in range(N)], "Initial condition is not equal for all iterations" 
#asserire che sum deve essere sempre <= N? (usa .values())

    #Centrality measures
rec_spec_radius, Q = Evolutions.communicability(temporal_network)
nodes_Bcentrality, nodes_Brank = Evolutions.broadcast_ranking(Q) #scores, node rankings
nodes_Rcentrality, nodes_Rrank = Evolutions.receive_ranking(Q) #scores, node rankings

    #Virulence measures
virulence = [] #i-th entry is virulence of i-th node
for index_case in range(N):
    virulence.append([]) #create the i-th entry
    for iteration in range(K):
        virulence[index_case].append(Propagation_SI.time_score(label[index_case][iteration],0.6))
    virulence[index_case] = np.mean(virulence[index_case])
virulence_rank = np.argsort(virulence)

#Results print
print("Top B-Centrality nodes:")
print(nodes_Brank[0:9])
print("Top Virulence nodes:")
print(virulence_rank[0:9])
print("Common nodes")
print(set(nodes_Brank[0:9]).intersection(set(virulence_rank[0:9])))