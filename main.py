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
import os
import matplotlib.pyplot as plt

#%%         CONFIGURATION PARAMETERS
net_name = input("Specify a name for this simulation: ")

isSAMPLED = False #if True, sample from the empiric distribution
isDAR = True
P = 1
isDIRECTED = False

N = 100
T = 100

beta = 0.01 #infection rate

NET_REAL = 20
K = 100 #infection realizations
assert NET_REAL >= 1, "NET_REAL should be >=1"
assert K > 1, "K should be >1"

alpha_sigma = 0.05

eigen_fraction = 0.25 #fraction of max eigenvalue to use in communicability
multiple_infections= True #allow or not multiple infections per time step
infective_fraction = 0.5

file_name = Saves.make_basics(isDAR=isDAR,P=P,isDIRECTED=isDIRECTED,isSAMPLED=isSAMPLED) + "_N"+str(N)+"_T"+str(T)+"_"+net_name+"/input_parameters.txt"
os.makedirs(os.path.dirname(file_name), exist_ok=True)
with open(file_name, 'w') as f:
    f.write("NAME = " + net_name + "\n\n")
    f.write("sampled = " + str(isSAMPLED)+ "\n")
    f.write("dar = " + str(isDAR)+ "\n")
    f.write("P = %i\n" %P)
    f.write("directed = " + str(isDIRECTED)+ "\n\n")
    f.write("N = %i\n" %N)
    f.write("T = %i\n" %T)
    f.write("beta = %.2f \n\n" %beta)
    f.write("NET_REAL = %i\n" %NET_REAL)
    f.write("K = %i\n\n" %K)
    f.write("alpha_sigma = %.2f \n\n" %alpha_sigma)
    f.write("eigen_fraction = %.2f \n" %eigen_fraction)
    f.write("multiple_infections = " + str(multiple_infections) + "\n")
    f.write("infective_fraction = %.2f" %infective_fraction)

#%%                         INPUTS BUILDING
if isSAMPLED:
    #DAR MATRICES
    if isDAR:
        with open('Empiric_Data/alphaDAR1.pkl', 'rb') as f:
            alphaDAR1 = pickle.load(f)
        with open('Empiric_Data/chiDAR1.pkl', 'rb') as f:
            chiDAR1 = pickle.load(f)
        #ALPHAS
        resh = alphaDAR1.reshape(98*98,1)
        empiric_alphas = np.linspace(0,0.99,N)
        alphas_probabilities = plt.hist(resh, bins = N)[0]/(98*98) #entry [0] is a count, so it has to be normalized
        alpha = np.random.choice(empiric_alphas, size=(N,N),p=alphas_probabilities)
        #CHIS
        resh = chiDAR1.reshape(98*98,1)
        empiric_chis = np.linspace(0,0.99,N)
        chis_probabilities = plt.hist(resh, bins = N)[0]/(98*98) #entry [0] is a count, so it has to be normalized
        chi = np.random.choice(empiric_chis, size=(N,N),p=chis_probabilities)
        assert abs(sum(alphas_probabilities)-1) <=0.001
        assert abs(sum(chis_probabilities)-1) <= 0.001
    else:   
        #TGRG MATRICES
        with open('Empiric_Data/phi0TGRG.pkl', 'rb') as f:
            phi0TGRG = pickle.load(f)
        with open('Empiric_Data/phi1TGRG.pkl', 'rb') as f:
            phi1TGRG = pickle.load(f)
        with open('Empiric_Data/sigmaTGRG.pkl', 'rb') as f:
            sigmaTGRG = pickle.load(f)
        #PHI0S
        resh = phi0TGRG[0:98] #FOR DIRECTED CASE, ALL 196 SHOULD BE TOOK
        empiric_phi0s = np.linspace(resh.min(),resh.max(),N)
        phi0s_probabilities = plt.hist(resh, bins = N)[0]/98 #entry [0] is a count, so it has to be normalized
        phi0 = np.random.choice(empiric_phi0s, size=N,p=phi0s_probabilities)
        #PHI1S
        resh = phi1TGRG[0:98] #FOR DIRECTED CASE, ALL 196 SHOULD BE TOOK
        empiric_phi1s = np.linspace(resh.min(),resh.max(),N)
        phi1s_probabilities = plt.hist(resh, bins = N)[0]/98 #entry [0] is a count, so it has to be normalized
        phi1 = np.random.choice(empiric_phi1s, size=N,p=phi1s_probabilities)
        #SIGMAS
        resh = sigmaTGRG[0:98] #FOR DIRECTED CASE, ALL 196 SHOULD BE TOOK
        empiric_sigmas = np.linspace(resh.min(),resh.max(),N)
        sigmas_probabilities = plt.hist(resh, bins = N)[0]/98 #entry [0] is a count, so it has to be normalized
        sigma = np.random.choice(empiric_sigmas, size=N,p=sigmas_probabilities)
        assert abs(sum(phi0s_probabilities)-1) <=0.001
        assert abs(sum(phi1s_probabilities)-1) <= 0.001
        assert abs(sum(sigmas_probabilities)-1) <= 0.001
else:
    if isDAR:
        alpha_sigma = 0.09
        alpha = np.random.normal(0.5, alpha_sigma, size=(N,N))
        chi = np.random.normal(0.5, alpha_sigma, size=(N,N))
    else:
        phi0 = np.random.poisson(3, size = N)
        phi1 = np.random.uniform(-0.99,0.99, size = N)
        sigma = np.random.poisson(1, size = N)
#%%                         OUTPUTS GENERATION
import time
start = time.time()

#Following lists contain NET_REAL arrays; each of these lists has N entries, with the average scores for each node
nodes_Bcentrality = []
nodes_Rcentrality = []
nodes_AD = []
nodes_BD = []
virulence = []
    
for k in range(1,NET_REAL+1):  #so first realization has index 1; for each k one tempnet and one k are overwritten
    #TEMPNETS GENERATION AND SAVE
    if isDAR: #use the proper functiond wheter user selected dar or tgrg
        temporal_network = Evolutions.network_generation_dar(alpha,chi,P=P,T=T,directed=isDIRECTED) 
    else:
        temporal_network = Evolutions.network_generation_tgrg(phi0,phi1,sigma,T=T,directed=isDIRECTED)[0]
    Saves.network_save(temporal_network,net_name, isDAR = isDAR, isDIRECTED = isDIRECTED, isSAMPLED=isSAMPLED, k=k, P=P)
    
    #CENTRALITIES GENERATION AND SAVE
    inv_maxradius, Q = Evolutions.communicability(temporal_network, eigen_fraction=eigen_fraction, length_one= not isDIRECTED) #if is undirected, use length 1
    
    singleiter_nodes_Bcentrality = Evolutions.broadcast_ranking(Q)[0]
    Saves.analysis_save(singleiter_nodes_Bcentrality,"BCENTR", net_name, N,T,isDAR=isDAR,isDIRECTED=isDIRECTED,isSAMPLED=isSAMPLED, k=k,P=P)
    nodes_Bcentrality.append(singleiter_nodes_Bcentrality)
    
    singleiter_nodes_Rcentrality = Evolutions.receive_ranking(Q)[0]
    Saves.analysis_save(singleiter_nodes_Rcentrality, "RCENTR", net_name, N,T,isDAR=isDAR,isDIRECTED=isDIRECTED,isSAMPLED=isSAMPLED, k=k,P=P)
    nodes_Rcentrality.append(singleiter_nodes_Rcentrality)   
    
    singleiter_nodes_AD = Evolutions.aggregate_degree(temporal_network, directed=False) #if dir = True, there are 2 outputs
    Saves.analysis_save(singleiter_nodes_AD, "AGGDEG", net_name, N,T,isDAR=isDAR,isDIRECTED=isDIRECTED,isSAMPLED=isSAMPLED, k=k,P=P)
    nodes_AD.append(singleiter_nodes_AD)

    singleiter_nodes_BD = Evolutions.binarized_degree(temporal_network, directed=False) #if dir = True, there are 2 outputs
    Saves.analysis_save(singleiter_nodes_BD, "BINDEG", net_name, N,T,isDAR=isDAR,isDIRECTED=isDIRECTED,isSAMPLED=isSAMPLED, k=k,P=P)
    nodes_BD.append(singleiter_nodes_BD)
    
    #SI PROPAGATION
    probabilities = dict() #probabilities dict
    for t in range(T):
        probabilities[t] = quad(Propagation_SI.poisson_probability,0,t, args = beta)[0] #quad produces several outputs, integral is the first
    label = []
    singleiter_virulece = []
    for index_case in range(N):
        label.append([]) #create the i-th entry
        for iteration in range(K):
            label[index_case].append(Propagation_SI.propagation(temporal_network, index_case, probabilities, multiple_infections = multiple_infections))
        singleiter_virulece.append(np.mean([Propagation_SI.time_score(label[index_case][iteration],infective_fraction) for iteration in range(K)]))
    virulence.append(singleiter_virulece)
    #SAVES
    Saves.infection_save(label,N,T,beta, net_name, isDAR = isDAR, isDIRECTED=isDIRECTED,isSAMPLED=isSAMPLED, k=k, P=P)
    Saves.analysis_save(singleiter_virulece, "VIRULENCE"+str(beta), net_name, N,T,isDAR=isDAR,isDIRECTED=isDIRECTED,isSAMPLED=isSAMPLED, k=k,P=P)
print(time.time()-start)  
#Re-loading of single-iteration lists to re-build the whole ones is possible but not implemented here
#%% ANALYSIS
nodes_B = []
[nodes_B.extend(el) for el in nodes_Bcentrality]
nodes_R = []
[nodes_R.extend(el) for el in nodes_Rcentrality]
nod_A = []
[nod_A.extend(el) for el in nodes_AD]
nod_B = []
[nod_B.extend(el) for el in nodes_BD]
vir = []
[vir.extend(el) for el in virulence]

fig_count = 2
plt.figure(fig_count)
fig_count+=1
plt.scatter(nodes_B, vir)
plt.xlabel("Broadcast Centrality")
plt.ylabel('Epidemic score')
plt.title(r"Undirected DAR(1) network, $\beta$ = %.3f; N=%i, T=%i; Netw iter = %i, Epid iter = %i" %(beta,N,T,NET_REAL,K))

plt.figure(fig_count)
fig_count+=1
plt.scatter(nodes_R, vir)
plt.xlabel("Receive Centrality")
plt.ylabel('Epidemic score')
plt.title(r"Undirected DAR(1) network; $\beta$ = %.3f; N=%i, T=%i; Netw iter = %i, Epid iter = %i" %(beta,N,T,NET_REAL,K))

plt.figure(fig_count)
fig_count+=1
plt.scatter(nod_A, vir)
plt.xlabel("Aggregate Degree")
plt.ylabel('Epidemic score')
plt.title(r"Undirected DAR(1) network; $\beta$ = %.3f; N=%i, T=%i; Netw iter = %i, Epid iter = %i" %(beta,N,T,NET_REAL,K))

plt.figure(fig_count)
fig_count+=1
plt.scatter(nod_B, vir)
plt.xlabel("Binarized Degree")
plt.ylabel('Epidemic score')
plt.title(r"Undirected DAR(1) network; $\beta$ = %.3f; N=%i, T=%i; Netw iter = %i, Epid iter = %i" %(beta,N,T,NET_REAL,K))


import scipy.stats
print(scipy.stats.pearsonr(nodes_B, vir))
print(scipy.stats.pearsonr(nodes_R, vir))
print(scipy.stats.pearsonr(nod_A, vir))
print(scipy.stats.pearsonr(nod_B, vir))

print("Common nodes in first 10 positions, BCENTR vs VIR")
for k in range(1,NET_REAL+1):
    print(set(np.flip(np.argsort(nodes_Bcentrality[k-1]))[0:10]).intersection(set(np.argsort(virulence[k-1])[0:10])))
    #Highest BCENTR should meet lowest virulence, which is a score of how much time it takes. So virulence is not flipped