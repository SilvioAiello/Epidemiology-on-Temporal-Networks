"""
From this script, user can generate temporal networks and perform epidemic upon them.
System parameters must be set in "inputs.ini" file.
Functions in this script work in Pyhon3, may require numpy (v1.16) and function "quad" from scipy.integrate (scipy v1.3).

For further understandings on how this script operates, check file "howto.md".
For further theoretical understandings, check file "explanation.md".
"""
import Saves
import Main_functions

import time
import numpy as np
#%%         CONFIGURATION PARAMETERS
import configparser
config = configparser.ConfigParser()
config.read('inputs.ini')
for section in config.sections():
    net_name = config[section]['net_name']
    
    N = config[section].getint('N')
    T = config[section].getint('T')
    beta = config[section].getfloat('beta')    
    
    isSAMPLED = config[section].getboolean('isSAMPLED') #directly sampled from empiric distribution
    isDAR = config[section].getboolean('isDAR')
    P = config[section].getint('P')
    isDIRECTED = config[section].getboolean('isDIRECTED')
    
    NET_REAL = config[section].getint('NET_REAL')
    K = config[section].getint('K')
    
    eigen_fraction = config[section].getfloat('eigen_fraction')
    multiple_infections = config[section].getboolean('multiple_infections')
    infective_fraction = config[section].getfloat('infective_fraction')
    
    directory_name = Saves.main_directory(beta, isDIRECTED,isDAR,P,isSAMPLED,N,T,net_name)
    
    #Tempo infezione al variare di alfa:
    #alfa = [0.1,0.25,0.5,0.75,0.99]
    alfa = 0.5
    chi = 0.05
    fig_count = 2
    q_param = beta


#%%                         INPUTS BUILDING AND SAVE     
    if isDAR:
        start = time.time()
        temporal_networks_RTNN, alpha, chi, avgdegree_evolution_RT,communicabilities_RTNN,nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN, fig_count = Main_functions.STRUCTURE_GEN_ANALYSIS(N,T, alfa, NET_REAL, P, q_param, isDIRECTED, isDAR,isSAMPLED,fig_count, directory_name)
        infection_matrices_RNITN,virulence_RNI,time_tobe_infected_RNIN,epidemic_size_evolution_RNIT,nodes_Bcentrality_RNT,mean_epidemic_RNT,correlations_RT = Main_functions.INFECTION_GEN_ANALYSIS(N,T, beta, NET_REAL, K, isDAR, infective_fraction,multiple_infections, temporal_networks_RTNN,communicabilities_RTNN, q_param, directory_name)
        end = time.time()
#    else:
#        start = time.time()
#        temporal_networks_RTNN, phi0,phi1,sigma, avgdegree_evolution_RT,communicabilities_RTNN,nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN, fig_count = Main_functions.STRUCTURE_GEN_ANALYSIS(N,T, alfa, NET_REAL, P, q_parameters, isDIRECTED, isDAR,isSAMPLED,fig_count, directory_name)
#        infection_matrices_RNITN,virulence_RNI,time_tobe_infected_RNIN,epidemic_size_evolution_RNIT,nodes_Bcentrality_QRNT,mean_epidemic_RNT,correlations_RT = Main_functions.INFECTION_GEN_ANALYSIS(N,T, beta, NET_REAL, K, isDAR, infective_fraction,multiple_infections, temporal_networks_RTNN,communicabilities_QRTNN, directory_name)
#        end = time.time()
    print(end-start)
    #%%                         ANALYSIS AND RESULTS
        #virulence_scores,timetobe_scores,avg_epid_ev,avg_correlations_evolution,correlationsBR,correlationsBV,correlationsRV,correlationsAV,correlationsBinV,correlationsRTim, fig_count = Main_functions.results(N,T, beta, K, alfa, isDAR,P, isDIRECTED, isSAMPLED, NET_REAL,net_name, multiple_infections,eigen_fraction, infective_fraction, fig_count, virulence_RNI, directory_name+'_alfa'+str(alfa)+'_beta'+str(beta), time_tobe_infected_RNIN,epidemic_size_evolution_RNIT, correlations_ABRT[alfa][beta],nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN,avgdegree_evolution_RT)
    virulence_scores,timetobe_scores,avg_epid_ev,avg_correlations_evolution, fig_count = Main_functions.results(N,T, beta, K, alfa, isDAR,P, isDIRECTED, isSAMPLED, NET_REAL,net_name, multiple_infections,eigen_fraction, infective_fraction, fig_count, virulence_RNI, directory_name, time_tobe_infected_RNIN,epidemic_size_evolution_RNIT, correlations_RT,q_param,nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN,avgdegree_evolution_RT)
    #virulence_scores,avg_epid_ev,avg_correlations_evolution, fig_count = Main_functions.results(N,T, beta, K, alfa, isDAR,P, isDIRECTED, isSAMPLED, NET_REAL,net_name, multiple_infections,eigen_fraction, infective_fraction, fig_count, virulence_RNI, directory_name, time_tobe_infected_RNIN,epidemic_size_evolution_RNIT, correlations_QRT,q_parameters,nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN,avgdegree_evolution_RT)
    
    #Per queste analisi ulteriori, lasciamo sempre che le funzioni normalmente agiscano fino a R, per livelli superiori agire da qua
    