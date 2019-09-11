"""
From this script, user can generate temporal networks and perform epidemic upon them.
System parameters must be set in "inputs.ini" file.
Functions in this script work in Pyhon3, may require numpy (v1.16) and function "quad" from scipy.integrate (scipy v1.3).

In this module are defined the following functions:
    * poisson_probability

For further understandings on how this script operates, check file "howto.md".
For further theoretical understandings, check file "explanation.md".
"""
import numpy as np
from scipy.integrate import quad #used in dictionary of probabilities

import Evolutions
import Propagation_SI
import Saves

import pickle
import matplotlib.pyplot as plt

#%%         CONFIGURATION PARAMETERS
net_name = "THESIS_TRY"
N = 100
T = 90
isDAR = True
P = 1
isDIRECTED = False
beta = 0.005 #infection rate

NET_REAL = 25
K = 100 #infection realizations

alpha_constant = 0.25
xi_constant = 0.5

phi0_constant = 0.1
phi1_constant = 0.35
sigma_constant = 0.05

#%% MATRICES BUILDING
    #DAR MATRICES
alpha = alpha_constant*np.ones((N,N))
xi = xi_constant*np.ones((N,N))
    #TGRG MATRICES
phi0 = phi0_constant*np.ones(N)
phi1 = phi1_constant*np.ones(N)
sigma= sigma_constant*np.ones(N)
#%%     OUTPUTS GENERATION
    #preliminar assertions
assert NET_REAL >= 1, "NET_REAL should be >=1"
assert K > 1, "K should be >1"

    #TEMPNETS GENERATION
for k in range(1,NET_REAL+1):  #so first realization has index 1
    if isDAR: #use the proper functiond wheter user selected dar or tgrg
        temporal_network = Evolutions.network_generation_dar(alpha,xi,P=P,T=T,directed=isDIRECTED) 
    else:
        temporal_network = Evolutions.network_generation_tgrg(phi0,phi1,sigma,T=T,directed=isDIRECTED) 
    Saves.network_save(temporal_network,net_name, isDAR = isDAR, isDIRECTED = isDIRECTED, k=k, P=1)
    
        #CENTRALITIES
    Q = Evolutions.communicability(temporal_network)[1] #as there's 1 tempnet, there's 1 Q per k
    nodes_Bcentrality = Evolutions.broadcast_ranking(Q)[0]
    nodes_Rcentrality = Evolutions.receive_ranking(Q)[0]
        #SAVES
    starting_name = str()
    if isDAR:
        if isDIRECTED:
            starting_name = "Networks/N"+str(N)+"_T"+str(T)+"_DIRECT"+"_DAR"+str(P)+"_"+net_name+"/realization"+str(k)+"/infections_beta"+str(beta)
        else:
            starting_name = "Networks/N"+str(N)+"_T"+str(T)+"_UNDIRECT"+"_DAR"+str(P)+"_"+net_name+"/realization"+str(k)+"/infections_beta"+str(beta)            
    else:
        if isDIRECTED:
            starting_name = "Networks/N"+str(N)+"_T"+str(T)+"_DIRECT"+"_TGRG_"+net_name+"/realization"+str(k)+"/infections_beta"+str(beta)
        else:
            starting_name = "Networks/N"+str(N)+"_T"+str(T)+"_UNDIRECT"+"_TGRG_"+net_name+"/realization"+str(k)+"/infections_beta"+str(beta)           

    with open(starting_name+"_BCENTR.pkl", 'wb') as handle:
        pickle.dump(nodes_Bcentrality,handle)
    with open(starting_name+"_RCENTR.pkl", 'wb') as handle:
        pickle.dump(nodes_Rcentrality,handle)
        
        #SI PROPAGATION
        #Probabilities dict
    probabilities = dict() #dict building
    for t in range(T):
        probabilities[t] = quad(Propagation_SI.poisson_probability,0,t, args = beta)[0] #quad produces several outputs, integral is the first
        #Function evoking
    label = []
    virulence=[]
    for index_case in range(N):
        label.append([]) #create the i-th entry
        for iteration in range(K):
            label[index_case].append(Propagation_SI.propagation(temporal_network, index_case, probabilities))
        virulence.append(np.mean([Propagation_SI.time_score(label[index_case][k],0.6) for k in range(K)]))
        #SAVES
    Saves.infection_save(label,N,T,beta, net_name, isDAR = isDAR, isDIRECTED=isDIRECTED, k=k, P=1)
    with open(starting_name+"_VIRULENCE.pkl", 'wb') as handle:
        pickle.dump(virulence,handle)

#%% ANALYSIS
nodes_Bcentrality = []
nodes_Rcentrality = []
virulence = []

for k in range(1,NET_REAL+1):  #so first realization has index 1
    starting_name = str()
    if isDAR:
        if isDIRECTED:
            starting_name = "Networks/N"+str(N)+"_T"+str(T)+"_DIRECT"+"_DAR"+str(P)+"_"+net_name+"/realization"+str(k)+"/infections_beta"+str(beta)
        else:
            starting_name = "Networks/N"+str(N)+"_T"+str(T)+"_UNDIRECT"+"_DAR"+str(P)+"_"+net_name+"/realization"+str(k)+"/infections_beta"+str(beta)            
    else:
        if isDIRECTED:
            starting_name = "Networks/N"+str(N)+"_T"+str(T)+"_DIRECT"+"_TGRG_"+net_name+"/realization"+str(k)+"/infections_beta"+str(beta)
        else:
            starting_name = "Networks/N"+str(N)+"_T"+str(T)+"_UNDIRECT"+"_TGRG_"+net_name+"/realization"+str(k)+"/infections_beta"+str(beta)           
    
    with open(starting_name+"_BCENTR.pkl", 'rb') as f:
        nodes_Bcentrality.append(pickle.load(f))
    with open(starting_name+"_RCENTR.pkl", 'rb') as f:
        nodes_Rcentrality.append(pickle.load(f))
    with open(starting_name+"_VIRULENCE.pkl", 'rb') as f:
        virulence.append(pickle.load(f))

nodes_B = []
[nodes_B.extend(el) for el in nodes_Bcentrality]
nodes_R = []
[nodes_R.extend(el) for el in nodes_Rcentrality]
vir = []
[vir.extend(el) for el in virulence]

plt.figure(1)
plt.scatter(nodes_B, vir)
plt.xlabel("Broadcast Centrality")
plt.ylabel('Epidemic score')
plt.title(r"Undirected DAR(1) network, $\alpha$ = %.2f, $\chi$ = %.2f; $\beta$ = %.3f; N=%i, T=%i; Netw iter = %i, Epid iter = %i" %(alpha_constant,xi_constant,beta,N,T,NET_REAL,K))

plt.figure(2)
plt.scatter(nodes_R, vir)
plt.xlabel("Receive Centrality")
plt.ylabel('Epidemic score')
plt.title(r"Undirected DAR(1) network, $\alpha$ = %.2f, $\chi$ = %.2f; $\beta$ = %.3f; N=%i, T=%i; Netw iter = %i, Epid iter = %i" %(alpha_constant,xi_constant,beta,N,T,NET_REAL,K))


import scipy.stats
print(scipy.stats.pearsonr(nodes_B, vir))
print(scipy.stats.pearsonr(nodes_R, vir))