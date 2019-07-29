"""
From this script, user can perform structural and epidemiological measures upon networks.
System parameters must be set in "inputs.ini" file
Functions in this script work in Pyhon3, may require numpy (v1.16).

or further understandings on how this script operates, check file "howto.md".
For further theoretical understandings, check file "explanation.md".
"""
import numpy as np

import Evolutions
import Propagation_SI
import Saves

import configparser

#%%         CONFIGURATION READING
config = configparser.ConfigParser()
config.read('inputs.ini')
for section in config.sections(): #if you want to analyze few networks, just ???
    N = config[section].getint('N') #nodes
    T = config[section].getint('T') #duration
    isDAR = config[section].getboolean('isDAR') #type
    isDIRECTED = config[section].getboolean('isDIRECTED') #network symmetry
    
    NET_REAL  = config[section].getint('NET_REAL') #network realizations
    K = config[section].getint('K') #infection realizations
    
    net_name = config[section]['net_name'] #identification name
    
    P = config[section].getint('P') #DAR order
    
    beta = config[section].getfloat('beta') #infection rate
    
    #Imports
    for k in range(1,NET_REAL+1):  #so first realization has index 1
        temporal_network = Saves.network_load(N,T,net_name,isDAR=isDAR,k=k,P=P)
        label = Saves.infection_load(N,T,beta,net_name,isDAR=True,k=1, P=1)
        
    #%% MEASURES
        #Centrality
    rec_spec_radius, Q = Evolutions.communicability(temporal_network)
    nodes_Bcentrality, nodes_Brank = Evolutions.broadcast_ranking(Q) #scores, node rankings
    nodes_Rcentrality, nodes_Rrank = Evolutions.receive_ranking(Q) #scores, node rankings
    
        #Virulence
    virulence = [] #i-th entry is virulence of i-th node
    for index_case in range(N):
        virulence.append([]) #create the i-th entry
        for iteration in range(K):
            virulence[index_case].append(Propagation_SI.time_score(label[index_case][iteration],0.6))
        virulence[index_case] = np.mean(virulence[index_case])
    virulence_rank = np.argsort(virulence)
    
    #Results print
    print("")
    print("RESULTS FOR NETWORK: "+net_name)
    print("Top B-Centrality nodes:")
    print(nodes_Brank[0:9])
    print("Top Virulence nodes:")
    print(virulence_rank[0:9])
    print("Common nodes")
    print(set(nodes_Brank[0:9]).intersection(set(virulence_rank[0:9])))