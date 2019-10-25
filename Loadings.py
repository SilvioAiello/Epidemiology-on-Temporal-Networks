# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 08:31:03 2019

@author: Silvio
"""
#Ricorda di settare l'inputs.ini in modo da corrispondere con i dati su cui vuoi lavorare
import Saves
import Main_functions
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
    alfa_values = [0.2,0.4,0.6,0.8]
    beta_values = [0.1,0.2,0.3,0.5,0.6,0.9]
    chi = 0.05
    fig_count = 1

    Saves.main_directory(beta, isDIRECTED,isDAR,P,isSAMPLED,N,T,net_name)
    correlations_ABRT = {alfa:{beta:{} for beta in beta_values} for alfa in alfa_values}
    
    for alfa in alfa_values:
        for beta in beta_values:
            temporal_networks_RTNN = Saves.network_load(directory_name+'_alfa'+str(alfa)+'_beta'+str(beta),'networks')
            avgdegree_evolution_RT = Saves.network_load(directory_name+'_alfa'+str(alfa)+'_beta'+str(beta),'DEGREE_EVOLUTION')
            communicabilities_RTNN = Saves.network_load(directory_name+'_alfa'+str(alfa)+'_beta'+str(beta),'Communicabilities_evolutions')
            nodes_Bcentrality_RN = Saves.network_load(directory_name+'_alfa'+str(alfa)+'_beta'+str(beta),'BROADCAST_CENTRALITIES')
            nodes_Rcentrality_RN = Saves.network_load(directory_name+'_alfa'+str(alfa)+'_beta'+str(beta),'RECEIVER_CENTRALITIES')
            nodes_AD_RN = Saves.network_load(directory_name+'_alfa'+str(alfa)+'_beta'+str(beta),'AGGREGATE_DEGREES')
            nodes_BD_RN = Saves.network_load(directory_name+'_alfa'+str(alfa)+'_beta'+str(beta),'BINARIZED_DEGREES')
            infection_matrices_RNITN = Saves.network_load(directory_name+'_alfa'+str(alfa)+'_beta'+str(beta),'infective_state_matrices')
            virulence_RNI = Saves.network_load(directory_name+'_alfa'+str(alfa)+'_beta'+str(beta),'VIRULENCE_SCORES')
            time_tobe_infected_RNIN = Saves.network_load(directory_name+'_alfa'+str(alfa)+'_beta'+str(beta),'TIMEINFECTED_SCORES')
            epidemic_size_evolution_RNIT = Saves.network_load(directory_name+'_alfa'+str(alfa)+'_beta'+str(beta),'Epidemic_size_evolution')
            correlations_ABRT[alfa][beta] = Saves.network_load(directory_name+'_alfa'+str(alfa)+'_beta'+str(beta),'correlations')
            
            #virulence_scores,timetobe_scores,avg_epid_ev,avg_correlations_evolution,correlationsBR,correlationsBV,correlationsRV,correlationsAV,correlationsBinV,correlationsRTim, fig_count = Main_functions.results(N,T, beta, K, alfa, isDAR,P, isDIRECTED, isSAMPLED, NET_REAL,net_name, multiple_infections,eigen_fraction, infective_fraction, fig_count, virulence_RNI, directory_name+'_alfa'+str(alfa)+'_beta'+str(beta), time_tobe_infected_RNIN,epidemic_size_evolution_RNIT, correlations_ABRT[alfa][beta],nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN,avgdegree_evolution_RT)
            
            virulence_scores,avg_epid_ev,avg_correlations_evolution, fig_count = Main_functions.correlations_analysis(N,T,beta,NET_REAL,fig_count,virulence_RNI,correlations_ABRT[alfa][beta],epidemic_size_evolution_RNIT,alfa,beta, chi,directory_name)