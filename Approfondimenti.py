# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 11:41:03 2019

@author: Silvio
"""
import numpy as np
import Evolutions
import Saves
import Main_functions
import time
import matplotlib.pyplot as plt
import os
import scipy.stats

#%% 0.3 BC con parametri diversi da beta
#N = 50; T=35; P =1
#beta = 0.4
#alfa=0.6
#ki = 0.95
#
#NET_REAL = 50; #K = 30
#isDAR=True;isDIRECTED = False; isSAMPLED = True
#infective_fraction = 1; multiple_infections = True; fig_count=1; final_fig = 1000
#net_name= "INFSOLOBETA 0.3 Q con parametri diversi da beta"
#directory_name = Saves.main_directory(beta, isDIRECTED,isDAR,P,isSAMPLED,N,T,net_name)
#
#q_parameters = [0.1,0.4,0.7] #0.01, 0.05,0.1,0.15,0.5 per beta = 0.1
# Guarda che esiste la funzione in MainFunc "produce communicability"!

#temporal_networks_RTNN, alpha, chi, avgdegree_evolution_RT,communicabilities_RTNN,\
# nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN, fig_count \
#  = Main_functions.STRUCTURE_GEN_ANALYSIS(N,T, alfa,ki, NET_REAL, P, q_param, isDIRECTED, isDAR,isSAMPLED,fig_count, directory_name, optionals=False)
#nodes_Bcentrality_QRNT= {qparam:np.zeros((NET_REAL,T,N)) for qparam in q_parameters}
#BC_Q = {qparam:{} for qparam in q_parameters}
#fig_count+=100
#plt.figure(fig_count)
#for q_param in q_parameters:
#    for netw_realiz in range(NET_REAL):
#        nodes_Bcentrality_QRNT[q_param][netw_realiz] = np.array([Evolutions.broadcast_ranking(Evolutions.communicability(temporal_networks_RTNN[netw_realiz], q_param, length_one=True)[t])[0] for t in range(T)])
#    BC_Q[q_param] = np.average(nodes_Bcentrality_QRNT[q_param].T[0],axis=1)
#    plt.plot(np.linspace(0,T-1,T),np.log(BC_Q[q_param]), label = "qparam = "+str(q_param))
#plt.legend()
#plt.grid()
#plt.title(r"log(BC) Centrality comparisons; $\alpha$=%.2f, $\beta$=%.2f,$\chi$=%.2f, N=%i,T=%i,#REAL=%i" %(alfa,beta,ki,N,T,NET_REAL))

##%% 1. Correlation <N>-BC by varying beta and chi
#start = time.time()
##Le uniche variazioni qui sono da apportare a seconda che tu scelga Newman (statico, alfa=0), quindi cambi i nomi, o Williams (viceversa viceversa)
#N = 50; T=55; P =1
#
#NET_REAL = 5; K = 20
#
#isDAR=True;isDIRECTED = False; isSAMPLED = True
#infective_fraction = 1; multiple_infections = True; fig_count=1; final_fig = 1000
#net_name="INFSOLOBETA 1. Correlation N-BC varying beta chi"
#directory_name = Saves.main_directory(0.00, isDIRECTED,isDAR,P,isSAMPLED,N,T,net_name)
#
#alfa = 0.7
#beta_parameters=[0.1,0.5,0.9]
#ki_parameters = [0.1,0.5,0.9]
#
#temporal_networks_BCRTNN = {beta:{ki:[] for ki in ki_parameters} for beta in beta_parameters}
#infection_matrices_BCRNITN = {beta:{ki:[] for ki in ki_parameters} for beta in beta_parameters}
#correlations_BCRT = {beta:{ki:np.zeros((NET_REAL, T-1)) for ki in ki_parameters} for beta in beta_parameters}
#
#for beta in beta_parameters:
#    for ki in ki_parameters:
        #temporal_networks_RTNN, alpha, chi, avgdegree_evolution_RT,communicabilities_RTNN,\
        # nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN, fig_count \
        #  = Main_functions.STRUCTURE_GEN_ANALYSIS(N,T, alfa,ki, NET_REAL, P, beta, isDIRECTED, isDAR,isSAMPLED,fig_count, directory_name, optionals=False)
#        print("tempnet ok")  
#        infection_matrices_BCRNITN[beta][ki],virulence_RNI,time_tobe_infected_RNIN,epidemic_size_evolution_RNIT, \
#              nodes_Bcentrality_RNT,mean_epidemic_RNT,correlations_BCRT[beta][ki] \
#                = Main_functions.INFECTION_GEN_ANALYSIS(N,T, beta, NET_REAL, K, isDAR, infective_fraction,multiple_infections, temporal_networks_BCRTNN[beta][ki],communicabilities_RTNN, directory_name+"_betaki"+str(beta)+str(ki))
#        print("tempnet ok")
#
#fig_count+=2000
#plt.figure(fig_count)
#plt.subplot(131)
#for ki in ki_parameters:
#    plt.plot(np.linspace(0,T-2,T-1), np.average(correlations_BCRT[0.1][ki], axis = 0), label = "ki = "+str(ki))
#plt.legend()
#plt.grid()
#plt.xlabel('Temporal step')
#plt.ylabel('Correlation')
#
#plt.subplot(132)
#for ki in ki_parameters:
#    plt.plot(np.linspace(0,T-2,T-1), np.average(correlations_BCRT[0.5][ki], axis = 0), label = "ki = "+str(ki))
#plt.xlabel(r"Temporal step")
#plt.ylabel(r"Correlation")
#plt.legend()
#plt.grid()
#
#plt.subplot(133)
#for ki in ki_parameters:
#    plt.plot(np.linspace(0,T-2,T-1), np.average(correlations_BCRT[0.9][ki], axis = 0), label = "ki = "+str(ki))
#plt.xlabel(r"Temporal step")
#plt.ylabel(r"Correlation")
#plt.legend()
#plt.grid()
#
#plt.suptitle(r"<N>-BC correlations, $\alpha$=%.2f, $\beta$=0.1/0.5/0.9, N=%i,T=%i,#REAL=%i,K=%i,inf_frac=%.2f;" %(alfa,N,T,NET_REAL,K,infective_fraction))
#
#
#file_name = directory_name+"/grapichs/beta"+str(beta)
#os.makedirs(os.path.dirname(file_name), exist_ok=True)
#plt.savefig(file_name+net_name+".pdf")
#fig_count+=1
#
##Se funonzia, ripetere per chi
#plt.figure(fig_count)
#plt.subplot(131)
#for beta in beta_parameters:
#    plt.plot(np.linspace(0,T-2,T-1), np.average(correlations_BCRT[beta][0.1], axis = 0), label = "beta = "+str(beta))
#plt.legend()
#plt.grid()
#plt.xlabel('Temporal step')
#plt.ylabel('Correlation')
#
#plt.subplot(132)
#for beta in beta_parameters:
#    plt.plot(np.linspace(0,T-2,T-1), np.average(correlations_BCRT[beta][0.5], axis = 0), label = "beta = "+str(beta))
#plt.xlabel(r"Temporal step")
#plt.ylabel(r"Correlation")
#plt.legend()
#plt.grid()
#
#plt.subplot(133)
#for beta in beta_parameters:
#    plt.plot(np.linspace(0,T-2,T-1), np.average(correlations_BCRT[beta][0.9], axis = 0), label = "beta = "+str(beta))
#plt.xlabel(r"Temporal step")
#plt.ylabel(r"Correlation")
#plt.legend()
#plt.grid()
#
#plt.suptitle(r"<N>-BC correlations, $\alpha$=%.2f, $\chi$=0.1/0.5/0.9, N=%i,T=%i,#REAL=%i,K=%i,inf_frac=%.2f;" %(alfa, N,T,NET_REAL,K,infective_fraction))
#os.makedirs(os.path.dirname(file_name), exist_ok=True)
#plt.savefig(file_name+net_name+".pdf")


#%% 2. Correlation between alfa and time of infection
start = time.time()

N = 5; T=15; P =1

ki = 0.33
alfa_values = [0.1,0.35,0.6,0.85] #quando tutto funzionerà, aumentare i valori
beta_values = [0.1,0.3]

NET_REAL = 10; K = 20

isDAR=True;isDIRECTED = False; isSAMPLED = True
infective_fraction = 2*1/N; multiple_infections = True; fig_count=1; final_fig = 1000

net_name="INFSOLOBETA 2. Correlation between alfa and time of infection - TUTTIPUNTI"

directory_name = Saves.main_directory(0.00, isDIRECTED,isDAR,P,isSAMPLED,N,T,net_name)

temporal_networks_ARTNN = {alfa:{beta:{} for beta in beta_values} for alfa in alfa_values}
infection_matrices_ARNITN = {alfa:{beta:{} for beta in beta_values} for alfa in alfa_values}
virulence_ARNI = {alfa:{beta:{} for beta in beta_values} for alfa in alfa_values}

communicabilities_ARTNN= {alfa:{beta:{} for beta in beta_values} for alfa in alfa_values}
epidemic_size_evolution_ARNIT= {alfa:{beta:{} for beta in beta_values} for alfa in alfa_values}
infection_states_ARNITN= {alfa:{beta:{} for beta in beta_values} for alfa in alfa_values}
time_tobe_infected_ARNIN= {alfa:{beta:{} for beta in beta_values} for alfa in alfa_values}

for alfa in alfa_values:
    for beta in beta_values:
        temporal_networks_ARTNN[alfa][beta], alpha, chi, avgdegree_evolution_RT,communicabilities_ARTNN[alfa][beta],\
         nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN, fig_count \
          = Main_functions.STRUCTURE_GEN_ANALYSIS(N,T, alfa,ki, NET_REAL, P, beta, isDIRECTED, isDAR,isSAMPLED,fig_count, directory_name, optionals=True)
        print(time.time()-start)
        infection_states_ARNITN[alfa][beta],virulence_ARNI[alfa][beta],time_tobe_infected_ARNIN[alfa][beta],epidemic_size_evolution_ARNIT[alfa][beta],nodes_Bcentrality_RNT,mean_epidemic_RNT,correlations_RT \
          =Main_functions.INFECTION_GEN_ANALYSIS(N,T, beta, NET_REAL, K, isDAR, infective_fraction,multiple_infections, temporal_networks_ARTNN[alfa][beta],communicabilities_ARTNN[alfa][beta], directory_name+str(beta))
        print(time.time()-start)

network_score = {alfa:{} for alfa in alfa_values}
for alfa in alfa_values:
    for beta in beta_values:
        network_score[alfa][beta] = np.average([np.average([np.average(virulence_ARNI[alfa][beta][net_r][index],axis=0) for index in range(N)],axis=0) for net_r in range(NET_REAL)],axis=0)

Saves.network_save(network_score,directory_name,'network_score')


beta_colors = {0.1:"blue", 0.3:"green", 0.5:"orange"}
fig_count+=1
plt.figure(fig_count)
for beta in beta_values:
    plt.scatter(alfa_values, [network_score[alfa][beta] for alfa in alfa_values], color=beta_colors[beta], label = "beta = "+str(beta))
plt.legend()
plt.grid()
plt.suptitle(r'Avg (over NetR,EpiR, Index) $\tau(\alpha)$, chi=%.3f, N=%i, T=%i,' %(ki, N,T,))
plt.title(r'#REAL=%i, K=%i, inf_frac=%.2f; dots = simulation results' %(NET_REAL,K,infective_fraction))
plt.xlabel(r"$\alpha$ values")
plt.ylabel(r"Average infective time")
file_name = directory_name+"/grapichs/"
os.makedirs(os.path.dirname(file_name), exist_ok=True)
plt.savefig(file_name+"_chi"+str(ki)+"timevsalfa.pdf")

##Provo a calcolare l'aspettazione teorica
#tau = {alfa:{} for alfa in alfa_values}
#for alfa in alfa_values:
#    T_matrix = np.zeros((2**P,2**P)); TL_matrix = np.zeros((2**P,2**P))
#    T_matrix[0][0] = 1-(1-alfa)*ki; T_matrix[1][0] = (1-alfa)*(1-ki)
#    TL_matrix[0][0] = 1-(1-alfa)*ki; TL_matrix[1][0] = (1-alfa)*(1-ki)
#    T_matrix[0][1] = (1-alfa)*ki; T_matrix[1][1] = alfa+(1-alfa)*ki
#    tau_vector = [sum(np.linalg.inv(np.identity(2**P) - (1-beta)*T_matrix - beta*TL_matrix)[i]*(2**(-P))) for i in range(2**P)]
#    tau[alfa] = tau_vector[0]*(1-ki) + tau_vector[1]*ki

#Provo a calcolare l'aspettazione teorica 2
#Questo strumento, anche finché non è provato che torni con le simulazioni, è prezioso per capire quanto far durare il network
tau_2 = {alfa:{beta:{} for beta in beta_values} for alfa in alfa_values}
plt.figure(fig_count)
for beta in beta_values:
    for alfa in alfa_values:
        T_matrix_2 = np.zeros((2**P,2**P)); TL_matrix_2 = np.zeros((2**P,2**P))
        T_matrix_2[0][0] = 1-(1-alfa)*ki; T_matrix_2[1][0] = (1-alfa)*(1-ki)
        TL_matrix_2[0][0] = 1-(1-alfa)*ki; TL_matrix_2[1][0] = (1-alfa)*(1-ki)
        T_matrix_2[0][1] = (1-beta)*(1-alfa)*ki; T_matrix_2[1][1] = (1-beta)*(alfa+(1-alfa)*ki)
        tau_vector_2 = [sum(np.linalg.inv(np.identity(2**P) - T_matrix_2)[i]*(2**(-P))) for i in range(2**P)]
        tau_2[alfa][beta] = tau_vector_2[0]*(1-ki) + tau_vector_2[1]*ki
    plt.plot(alfa_values, [tau_2[alfa][beta] for alfa in alfa_values], color=beta_colors[beta])

#%% 3. Scorporati e Reproduction of Newman and Williams results
#start = time.time()
##Le uniche variazioni qui sono da apportare a seconda che tu scelga Newman (statico, alfa=0), quindi cambi i nomi, o Williams (viceversa viceversa)
#N = 100; T=20; P =1
#beta = 0.99
#
##alfa = 0; alpha = alfa*np.ones(shape = (N,N)) #Dynamic network
#
#alfa = 1; alpha = alfa*(np.ones(shape = (N,N)) - np.diag([1 for i in range(N)])) #Static network
#ki_values = [0.1]
#
#NET_REAL = 10; K = 2
#isDAR=True;isDIRECTED = False; isSAMPLED = True
#infective_fraction = 1; multiple_infections = True; fig_count=1; final_fig = 1000
#
#net_name="3. Reproduction of Newman and Williams results - Riprova"
#
#directory_name = Saves.main_directory(beta, isDIRECTED,isDAR,P,isSAMPLED,N,T,net_name)
#
#temporal_networks_CRTNN = {ki:{} for ki in ki_values}
#infection_states_CRNITN = {ki:{} for ki in ki_values}
#virulence_CRNI= {ki:{} for ki in ki_values};time_tobe_infected_CRNIN= {ki:{} for ki in ki_values};epidemic_size_evolution_CRNIT= {ki:{} for ki in ki_values};nodes_Bcentrality_CRNT= {ki:{} for ki in ki_values};
#mean_epidemic_CRNT= {ki:{} for ki in ki_values};correlations_CRT= {ki:{} for ki in ki_values}; communicabilities_CRTNN = {ki:{} for ki in ki_values}
#for ki in ki_values:
##    if ki == 0.1:
##        T = 20
##    else:
##        T=100
#    chi = ki*(np.ones(shape = (N,N)) - np.diag([1 for i in range(N)]))
#    temporal_networks_RTNN, alpha, chi, avgdegree_evolution_RT,communicabilities_RTNN,\
#     nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN, fig_count \
#      = Main_functions.STRUCTURE_GEN_ANALYSIS(N,T, alfa,ki, NET_REAL, P, beta, isDIRECTED, isDAR,isSAMPLED,fig_count, directory_name, optionals=False)
#    step1=time.time(); print(step1-start)
#    #Propagation
#    infection_states_CRNITN[ki],virulence_CRNI[ki],time_tobe_infected_CRNIN[ki],epidemic_size_evolution_CRNIT[ki],nodes_Bcentrality_CRNT[ki],mean_epidemic_CRNT[ki],correlations_CRT[ki] \
#     = Main_functions.INFECTION_GEN_ANALYSIS(N,T, beta, NET_REAL, K, isDAR, infective_fraction,multiple_infections, temporal_networks_CRTNN[ki],communicabilities_CRTNN[ki], directory_name)
#    print(time.time()-step1)
#    print("ok")
#
#start = time.time()
#
#x_CRNNT = {ki:np.zeros((NET_REAL,N,N,T)) for ki in ki_values} #per ogni realizzazione, probab infetto nel tempo quando parte dall'index case
#
#for ki in ki_values:
##    if ki == 0.1:
##        T = 20
##    else:
##        T=100
#    for netw_real in range(NET_REAL):
#        for index_case in range(N):
#            x_CRNNT[ki][netw_real][index_case][index_case] = np.ones(T)
#for ki in ki_values:
##    if ki == 0.1:
##        T = 20
##    else:
##        T=100
#    for netw_real in range(NET_REAL):
#        for index_case in range(N):
#            for t in range(1,T):
#                for n in [x for x in range(N) if x != index_case]:
#                    numero = len({infettatore for infettatore in range(N) if temporal_networks_CRTNN[ki][netw_real][t-1][infettatore][n] ==1 and x_CRNNT[ki][netw_real][index_case][infettatore][t-1] != 0})
#                    if numero == 0:
#                        x_CRNNT[ki][netw_real][index_case][n][t] = x_CRNNT[ki][netw_real][index_case][n][t-1]
#                    else:
#                        x_CRNNT[ki][netw_real][index_case][n][t] = (x_CRNNT[ki][netw_real][index_case][n][t-1] + beta*(1-x_CRNNT[ki][netw_real][index_case][n][t-1])*sum([temporal_networks_CRTNN[ki][netw_real][t-1][k][n]*x_CRNNT[ki][netw_real][index_case][k][t-1] for k in range(N)])/numero)
#print(time.time()-start)
#
#bc_T= np.average([np.average([Evolutions.broadcast_ranking(communicabilities_CRTNN[ki][netw_realiz][t])[0] for t in range(T)],axis=1) for netw_realiz in range(NET_REAL)], axis = 0)
#
#os.makedirs(os.path.dirname(directory_name+"/grapichs/"), exist_ok=True)
#for ki in ki_values:
##    if ki == 0.1:
##        T = 20
##    else:
##        T=100
#    plt.figure(final_fig)
#    plt.plot(np.linspace(0,T-1,T),np.average([np.average([np.average(epidemic_size_evolution_CRNIT[ki][net_real][index_case], axis = 0) for index_case in range(N)], axis=0) for net_real in range(NET_REAL)], axis=0)/N, label="<N>(t)")
#    plt.plot(np.linspace(0,T-1,T),np.average([np.average([sum(x_CRNNT[ki][net_real][index_case]) for index_case in range(N)],axis=0)/N for net_real in range(NET_REAL)],axis = 0), label="first-order", color = "orange")
#    final_fig+=1
#    plt.legend()
#    plt.grid()
#    plt.xlabel("Temporal step")
#    plt.ylabel("Fraction of infected nodes")
#    plt.suptitle(r'chi=%.3f, N=%i, T=%i, beta = %.2f, alpha = %.2f' %(ki, N,T,beta,alfa))
#    plt.title(r'#REAL=%i, K=%i, inf_frac=%.2f' %(NET_REAL,K,infective_fraction))
#    plt.show()
#    plt.savefig(directory_name+"/grapichs/"+"3. TempNET"+str(ki)+".png")
##    plt.figure(final_fig)
##    plt.scatter(np.linspace(0,T-1,T),bc_T[0:T], label="BC(t)")
##    plt.grid()
##    plt.title("chi=%.3f, N=%i, T=%i, beta=%.2f, alfa=%.2f, #REAL=%i, K=%i, Williams, BC(t)" %(ki, N,T,beta,alfa,NET_REAL,K))
##    plt.show()
##    plt.savefig(directory_name+"/grapichs/"+"WilliamsFigure_chi"+str(ki)+"communicability"+".png")
##    final_fig+=1
#
#Saves.network_save(x_CRNNT,directory_name,"prob_newman_williams")
    