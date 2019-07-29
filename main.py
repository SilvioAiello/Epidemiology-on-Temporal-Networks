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

import configparser
#%%         CONFIGURATION READING
config = configparser.ConfigParser()
config.read('inputs.ini')
for section in config.sections():
    N = config[section].getint('N') #nodes
    T = config[section].getint('T') #duration
    isDAR = config[section].getboolean('isDAR') #type
    isDIRECTED = config[section].getboolean('isDIRECTED') #network symmetry
    
    NET_REAL  = config[section].getint('NET_REAL') #network realizations
    K = config[section].getint('K') #infection realizations
    
    net_name = config[section]['net_name'] #identification name
    
    P = config[section].getint('P') #DAR order
    alpha_constant = config[section].getfloat('alpha_constant') #alpha matrix constant
    xi_constant = config[section].getfloat('xi_constant') #xi matrix constant
    
    phi0_constant = config[section].getfloat('phi0_constant') #phi0 vector constant
    phi1_constant = config[section].getfloat('phi1_constant') #phi1 vector constant
    epsilon_constant = config[section].getfloat('epsilon_constant') #epsilon vector constant

    #%% MATRICES BUILDING
        #DAR MATRICES
    alpha = alpha_constant*np.ones((N,N))
    xi = xi_constant*np.ones((N,N))
        #TGRG MATRICES
    phi0 = phi0_constant*np.ones(N)
    phi1 = phi1_constant*np.ones(N)
    epsilon=epsilon_constant*np.ones(N)
    
        #EPIDEMIC PARAMETERS
    beta = 0.005 #infection rate
    
    #%% OUTPUTS GENERATION
        #preliminar assertions
    assert NET_REAL >= 1, "NET_REAL should be >=1"
    assert K > 1, "K should be >1"
    
        #TEMPNETS GENERATION
    for k in range(1,NET_REAL+1):  #so first realization has index 1
        if isDAR: #use the proper functiond wheter user selected dar or tgrg
            temporal_network = Evolutions.network_generation_dar(alpha,xi,P=P,T=T,directed=isDIRECTED) 
        else:
            temporal_network = Evolutions.network_generation_tgrg(alpha,xi,P=P,T=T,directed=isDIRECTED) 
        Saves.network_save(temporal_network,net_name, isDAR = isDAR, k=k, P=1)
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
        Saves.infection_save(label,N,T,net_name, isDAR = isDAR, k=k, P=1)