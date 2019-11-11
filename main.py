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
    
    fig_count= 2
    
    #Tempo infezione al variare di alfa:
    #alfa = [0.1,0.25,0.5,0.75,0.99]
#%%                         INPUTS BUILDING AND SAVE     
#    if isDAR:
#        start = time.time()
#        temporal_networks_RTNN, alpha, chi, avgdegree_evolution_RT, \
#         communicabilities_RTNN,nodes_Bcentrality_RN,nodes_Rcentrality_RN, nodes_AD_RN,nodes_BD_RN, fig_count \
#           = Main_functions.STRUCTURE_GEN_ANALYSIS(N,T, alfa, NET_REAL, P, beta, isDIRECTED, isDAR,isSAMPLED,fig_count, directory_name)
#        
#        infection_matrices_RNITN,virulence_RNI,time_tobe_infected_RNIN,epidemic_size_evolution_RNIT, \
#              nodes_Bcentrality_RNT,mean_epidemic_RNT,correlations_RT \
#                = Main_functions.INFECTION_GEN_ANALYSIS(N,T, beta, NET_REAL, K, isDAR, infective_fraction,multiple_infections, temporal_networks_RTNN,communicabilities_RTNN, directory_name)
##    else:
##        start = time.time()
##        temporal_networks_RTNN, phi0,phi1,sigma, avgdegree_evolution_RT,communicabilities_RTNN,nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN, fig_count = Main_functions.STRUCTURE_GEN_ANALYSIS(N,T, alfa, NET_REAL, P, q_parameters, isDIRECTED, isDAR,isSAMPLED,fig_count, directory_name)
##        infection_matrices_RNITN,virulence_RNI,time_tobe_infected_RNIN,epidemic_size_evolution_RNIT,nodes_Bcentrality_QRNT,mean_epidemic_RNT,correlations_RT = Main_functions.INFECTION_GEN_ANALYSIS(N,T, beta, NET_REAL, K, isDAR, infective_fraction,multiple_infections, temporal_networks_RTNN,communicabilities_QRTNN, directory_name)
##        end = time.time()
    #%%                         ANALYSIS AND RESULTS
#        virulence_scores,timetobe_scores,avg_epid_ev,avg_correlations_evolution, fig_count \
#           = Main_functions.results(N,T, beta, K, alfa, isDAR,P, isDIRECTED, isSAMPLED, NET_REAL,net_name, multiple_infections,eigen_fraction, infective_fraction, fig_count, virulence_RNI, directory_name, time_tobe_infected_RNIN,epidemic_size_evolution_RNIT, correlations_RT,beta,nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN,avgdegree_evolution_RT)


import Evolutions
import Main_functions
import scipy.stats
import matplotlib.pyplot as plt
import os
file_name = directory_name+"/grapichs/"
os.makedirs(os.path.dirname(file_name), exist_ok=True)
#%% QPARAM DA DIVERSI BETA
#alfa = 0.5
#q_parameters = [0.3,0.35,0.37,0.4] #0.01, 0.05,0.1,0.15,0.5 per beta = 0.1
#correlations_QRT = {q_param:np.zeros((NET_REAL, T-1)) for q_param in q_parameters}
#nodes_Bcentrality_QRNT = {q_param:np.zeros((NET_REAL, N,T)) for q_param in q_parameters}
#fig_count+=1
#plt.figure(fig_count)
#
#temporal_networks_RTNN = Saves.network_load(directory_name,'networks')
#epidemic_size_evolution_RNIT = Saves.network_load(directory_name,'Epidemic_size_evolution')
#mean_epidemic_RNT = np.zeros((NET_REAL,N,T))
#for netw_realiz in range(NET_REAL):
#    for index_case in range(N):
#        mean_epidemic_RNT[netw_realiz][index_case] = np.average(epidemic_size_evolution_RNIT[netw_realiz][index_case], axis=0)
#
#for q_param in q_parameters:
#    for netw_realiz in range(NET_REAL):
#        nodes_Bcentrality_QRNT[q_param][netw_realiz]= np.array([Evolutions.broadcast_ranking(Main_functions.produce_communicRTNN(N,T,NET_REAL, isDIRECTED, q_param, directory_name, temporal_networks_RTNN)[0][netw_realiz][t])[0] for t in range(T)]).T
#        correlations_QRT[q_param][netw_realiz] = np.array([scipy.stats.pearsonr(mean_epidemic_RNT[netw_realiz].T[t+1],nodes_Bcentrality_QRNT[q_param][netw_realiz].T[t])[0] if scipy.stats.pearsonr(mean_epidemic_RNT[netw_realiz].T[t+1],nodes_Bcentrality_QRNT[q_param][netw_realiz].T[t])[0]==scipy.stats.pearsonr(mean_epidemic_RNT[netw_realiz].T[t+1],nodes_Bcentrality_QRNT[q_param][netw_realiz].T[t])[0] else 0 for t in range(T-1)])
#    #Saves.network_save(correlations_QRT, directory_name,'Qcorrelations')
#correlazione_media = {q_param:np.zeros(T-1) for q_param in q_parameters}
#for q_param in q_parameters:
#    correlazione_media[q_param] = np.average(correlations_QRT[q_param],axis=0)
#for q_param in q_parameters:
#    plt.plot(np.linspace(0,T-2,T-1),correlazione_media[q_param], label="q_param = "+str(q_param))
#plt.legend()
#plt.grid()
#plt.title(r"alfa = %.2f, chi = 0.9, beta = %.2f" %(alfa,beta))
#plt.xlabel("Temporal step")
#plt.ylabel(r"c(t)")
#plt.savefig(file_name+"qparam"+str(q_param)+"correlation_evolutions555.pdf")


#%% ALPHA VS TIME PER LE SIMULAZIONI GIA' FATTE #
fig_count+=1
plt.figure(fig_count)

ki = 0.5
alfa_parameters = [0.2,0.4,0.6,0.8]
beta_parameters = [0.1,0.2,0.3,0.5,0.6,0.9]
temporal_networks_ABRTNN = {alfa:{beta:Saves.network_load(Saves.main_directory(0.05, isDIRECTED,isDAR,P,isSAMPLED,N,T,net_name)+"_alfa"+str(alfa)+"_beta"+str(beta),'networks') for beta in beta_parameters} for alfa in alfa_parameters}
infection_matrices_ABRNITN = {alfa:{beta:Saves.network_load(Saves.main_directory(0.05, isDIRECTED,isDAR,P,isSAMPLED,N,T,net_name)+"_alfa"+str(alfa)+"_beta"+str(beta),'infective_state_matrices') for beta in beta_parameters} for alfa in alfa_parameters}
virulence_ABRNI = {alfa:{beta:Saves.network_load(Saves.main_directory(0.05, isDIRECTED,isDAR,P,isSAMPLED,N,T,net_name)+"_alfa"+str(alfa)+"_beta"+str(beta),'VIRULENCE_SCORES') for beta in beta_parameters} for alfa in alfa_parameters}

network_score = {alfa:{beta:0 for beta in beta_parameters} for alfa in alfa_parameters}
for alfa in alfa_parameters:
    for beta in beta_parameters:
        network_score[alfa][beta] = np.average([np.average([np.average(virulence_ABRNI[alfa][beta][net_r][index],axis=0) for index in range(N)],axis=0) for net_r in range(NET_REAL)],axis=0)

fig_count+=1
plt.figure(fig_count)
for beta in beta_parameters:
    plt.plot(alfa_parameters, [network_score[alfa][beta] for alfa in alfa_parameters], label = "beta = "+str(beta))
plt.legend()
plt.grid()
plt.title(r"Avg time (over iterations, index cases and realizations); chi = %.2f" %ki)
plt.xlabel(r"$\alpha$ values")
plt.ylabel(r"Average infective time")
plt.savefig(file_name+"_chi"+str(ki)+"timevsalfa.pdf")



#%% CHI PICCOLI
#alfa = 0.65
#beta = 0.4
#ki_parameters = [0.001, 0.01, 0.1,0.5,0.9]
#
#directory_name = Saves.main_directory(beta, isDIRECTED,isDAR,P,isSAMPLED,N,T,net_name)
#
#temporal_networks_CRTNN = {ki:[] for ki in ki_parameters}
#infection_matrices_CRNITN = {ki:[] for ki in ki_parameters}
#correlations_CRT = {ki:np.zeros((NET_REAL, T-1)) for ki in ki_parameters}
#
#for ki in ki_parameters:
#    temporal_networks_CRTNN[ki], alpha, chi, avgdegree_evolution_RT, \
#     communicabilities_RTNN,nodes_Bcentrality_RN,nodes_Rcentrality_RN, nodes_AD_RN,nodes_BD_RN, fig_count \
#       = Main_functions.STRUCTURE_GEN_ANALYSIS(N,T, alfa,ki, NET_REAL, P, beta, isDIRECTED, isDAR,isSAMPLED,fig_count, directory_name)
#    print("tempnet ok")
#       
#    infection_matrices_CRNITN[ki],virulence_RNI,time_tobe_infected_RNIN,epidemic_size_evolution_RNIT, \
#          nodes_Bcentrality_RNT,mean_epidemic_RNT,correlations_CRT[ki] \
#            = Main_functions.INFECTION_GEN_ANALYSIS(N,T, beta, NET_REAL, K, isDAR, infective_fraction,multiple_infections, temporal_networks_CRTNN[ki],communicabilities_RTNN, directory_name)
#
#fig_count+=1
#for ki in ki_parameters:
#    plt.plot(np.linspace(0,T-2,T-1), np.average(correlations_CRT[ki], axis = 0), label = "ki = "+str(ki))
#plt.legend()
#plt.grid()
#plt.title(r"<N>-BC correlations by varying chi parameter")
#plt.xlabel(r"$\alpha$ values")
#plt.ylabel(r"Correlation")
#plt.savefig(file_name+net_name+".pdf")


#%% SCORPORATI: N,BC,"integrazione", dalla simulazione già fatta

##PRIMA COSA: MODIFICA L'INPUTS.INI PER ANDARE A RIPESCARE, E ASSICURATI CHE SOPRA TUTTO SIA COMMENTATO
##Importo i due oggetti da cui calcoli le correlazioni
#alfa = 0.6
#beta =  0.1
#directory_name += "_alfa"+str(alfa)+"_beta"+str(beta)
#a_epidemic_size_evolution_RNIT= Saves.network_load(directory_name,'Epidemic_size_evolution')
#a_communicabilities_RTNN= Saves.network_load(directory_name,'Communicabilities_evolutions')
#a_temporal_network_RTNN = Saves.network_load(directory_name,'networks')
##Ma devo inportare anche l'evoluzione degli stati, a fungere da interrogazione alla natura:
#a_infection_matrices_RNITN = Saves.network_load(directory_name,'infective_state_matrices')
#
##Controllare a cipolla perché non è detto che abbia scritto cose sensate
#n_T = np.average([np.average([np.average(a_epidemic_size_evolution_RNIT[net_real][index_case], axis = 0) for index_case in range(N)], axis=0) for net_real in range(NET_REAL)], axis=0)
#bc_T= np.average([np.average([Evolutions.broadcast_ranking(a_communicabilities_RTNN[netw_realiz][t])[0] for t in range(T)],axis=1) for netw_realiz in range(NET_REAL)], axis = 0)
#
#x = np.zeros((NET_REAL,N,N,T)) #per ogni realizzazione, probab infetto nel tempo quando parte dall'index case
#for netw_real in range(NET_REAL):
#    for index_case in range(N):
#        x[netw_real][index_case][index_case] = np.ones(T)
#
#for netw_real in range(NET_REAL):
#    for index_case in range(N):
#        for t in range(1,T):
#            for n in [x for x in range(N) if x != index_case]:
#                x[netw_real][index_case][n][t] = x[netw_real][index_case][n][t-1] + beta*(1-x[netw_real][index_case][n][t-1])*sum([a_temporal_network_RTNN[netw_real][t-1][k][n]*x[netw_real][index_case][k][t-1] for k in range(N)])/N
##->aggiungere il fatto che il <N> è la somma sui nodi (questo per ogni index case)
#x_plot = np.average([np.average([sum(x[net_real][i]) for i in range(N)],axis=0) for net_real in range(NET_REAL)],axis=0)
#
#fig_count+=1
#plt.figure(fig_count+1050)
#TT = 25
#plt.scatter(np.linspace(0,TT-1,TT),n_T[0:TT], label="<N>(t)")
##plt.scatter(np.linspace(0,TT-1,TT),bc_T[0:TT], label="BC(t)")
#plt.plot(np.linspace(0,TT-1,TT),x_plot[0:TT], label="Somma probabilità")
#plt.legend()
#plt.grid()
#plt.title(r"Temporal evolution of average quantities")
#plt.xlabel(r"time")
#plt.ylabel(r"Average values")


