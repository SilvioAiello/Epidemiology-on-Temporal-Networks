"""
From this script, user can perform structural and epidemiological measures upon networks.
System parameters must be set in "inputs.ini" file
Functions in this script work in Pyhon3, may require numpy (v1.16).

or further understandings on how this script operates, check file "howto.md".
For further theoretical understandings, check file "explanation.md".
"""
import numpy as np
import matplotlib.pyplot as plt

import configparser

import pickle
import time
start = time.time()
#%%         CONFIGURATION READING
config = configparser.ConfigParser()
config.read('run1.ini') #CHANGE THE NUMBER OF THE RUN TO PERFORM YOUR SIMULATION
N = config['simulation'].getint('N') #nodes
T = config['simulation'].getint('T') #duration
isDAR = config['simulation'].getboolean('isDAR') #type
isDIRECTED = config['simulation'].getboolean('isDIRECTED') #network symmetry

NET_REAL  = config['simulation'].getint('NET_REAL') #network realizations
K = config['simulation'].getint('K') #infection realizations

net_name = config['simulation']['net_name'] #identification name

alpha_constant = config['simulation'].getfloat('alpha_constant')
xi_constant = config['simulation'].getfloat('xi_constant')

P = config['simulation'].getint('P') #DAR order

beta = config['simulation'].getfloat('beta') #infection rate

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


##Results print
#print("")
#print("RESULTS FOR NETWORK: "+net_name)
#print("Top B-Centrality nodes:")
##print(nodes_Brank[0:9])
#print("Top Virulence nodes:")
##print(virulence_rank[0:9])
#print("Common nodes")
##print(set(nodes_Brank[0:9]).intersection(set(virulence_rank[0:9])))
#
#end = time.time()
#print(end-start)