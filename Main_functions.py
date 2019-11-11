# -*- coding: utf-8 -*-
"""
Created on Sun Oct 20 18:37:31 2019

@author: Silvio
"""
import pickle
import numpy as np
import Evolutions
import Propagation_SI
import Assertions_suite
import Saves
import os
from scipy.integrate import quad #used in dictionary of probabilities
import scipy.stats
import matplotlib.pyplot as plt



def STRUCTURE_GEN_ANALYSIS(N,T, alfa,ki, NET_REAL, P, q_param, isDIRECTED, isDAR,isSAMPLED,fig_count, directory_name, optionals=True):
    """
    Returns temporal_networks_RTNN, alpha, chi, avgdegree_evolution_RT,
    communicabilities_RTNN,nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN, fig_count.

    """
    fig_count+=1
    alpha,chi,phi0,phi1,sigma,temporal_networks_RTNN,fig_count = produce_inpunt_tempnetRTNN(N,T,alfa,ki,NET_REAL,isDAR,P,isSAMPLED,isDIRECTED,directory_name,fig_count)
    
    #CENTRALITIES GENERATION AND SAVE
    communicabilities_RTNN, nodes_Bcentrality_RN, nodes_Rcentrality_RN= produce_communicRTNN(N, T,NET_REAL, isDIRECTED, q_param, directory_name, temporal_networks_RTNN)
    if optionals ==True:
        avgdegree_evolution_RT = [Evolutions.degree_mean(temporal_networks_RTNN[netw_realiz]) for netw_realiz in range(NET_REAL)]  
        nodes_AD_RN = [Evolutions.aggregate_degree(temporal_networks_RTNN[netw_realiz], directed=isDIRECTED)for netw_realiz in range(NET_REAL)]  #if dir = True, there are 2 outputs
        nodes_BD_RN = [Evolutions.binarized_degree(temporal_networks_RTNN[netw_realiz], directed=isDIRECTED)for netw_realiz in range(NET_REAL)]  #if dir = True, there are 2 outputs
    else:
        avgdegree_evolution_RT=[];nodes_AD_RN=[];nodes_BD_RN=[]
    #Saves.network_save(thetas_RNT,directory_name,'TGRG_thetas')
    Saves.network_save(temporal_networks_RTNN, directory_name, 'networks')
    Saves.network_save(avgdegree_evolution_RT, directory_name,'DEGREE_EVOLUTION')
    Saves.network_save(nodes_AD_RN, directory_name,'AGGREGATE_DEGREES')
    Saves.network_save(nodes_BD_RN, directory_name,'BINARIZED_DEGREES')
        
    if isDAR:
        return temporal_networks_RTNN, alpha, chi, avgdegree_evolution_RT,communicabilities_RTNN,nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN, fig_count
    else:
        return 0
        #return temporal_networks_RTNN, phi0,phi1,sigma, avgdegree_evolution_RT,communicabilities_RTNN,nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN, fig_count



def INFECTION_GEN_ANALYSIS(N,T, beta, NET_REAL, K, isDAR, infective_fraction,multiple_infections, temporal_networks_RTNN,communicabilities_RTNN, directory_name):
    #R -> NET_REAL, I -> infection iteration, N -> N, T->T
    infection_matrices_RNITN = []
    virulence_RNI = [] #scores computed when each node is the index case
    epidemic_size_evolution_RNIT = []
    time_tobe_infected_RNIN = [] #scores computed when each node is not the index case

    nodes_Bcentrality_RNT = np.zeros((NET_REAL,N,T))
    correlations_RT = np.zeros((NET_REAL,T-1))
    mean_epidemic_RNT = np.zeros((NET_REAL,N,T))

    for netw_realiz in range(NET_REAL):            
        #SI PROPAGATION
#        probabilities = dict() #probabilities dict
#        for t in range(T):
#            probabilities[t] = quad(Propagation_SI.poisson_probability,0,t, args = beta)[0] #quad produces several outputs, integral is the first
        
        infection_matrix = []
        iteration_virulences = []
        iterations_timetobe = []
        epidemic_size_evolution = []
        for index_case in range(N):
            infection_matrix.append([]) #create the i-th entry
            epidemic_size_evolution.append([])
            for iteration in range(K):
                infection_matrix[index_case].append(Propagation_SI.propagation_solobeta(temporal_networks_RTNN[netw_realiz], index_case, beta, multiple_infections = multiple_infections))
                epidemic_size_evolution[index_case].append([sum(infection_matrix[index_case][iteration][t]) for t in range(T)])
            iteration_virulences.append([Propagation_SI.time_score_v2(infection_matrix[index_case][iteration],infective_fraction) for iteration in range(K)])
            iterations_timetobe.append([list(Propagation_SI.when_is_infected(infection_matrix[index_case][iteration], index_case)) for iteration in range(K)])
        #Propagation last appends after loops
        infection_matrices_RNITN.append(infection_matrix)
        epidemic_size_evolution_RNIT.append(epidemic_size_evolution)
        virulence_RNI.append(iteration_virulences)
        time_tobe_infected_RNIN.append(iterations_timetobe)

        #Correlations <N>-BC
        for index_case in range(N):
            mean_epidemic_RNT[netw_realiz][index_case] = np.average(epidemic_size_evolution_RNIT[netw_realiz][index_case], axis=0)
            nodes_Bcentrality_RNT[netw_realiz][index_case]= np.array([Evolutions.broadcast_ranking(communicabilities_RTNN[netw_realiz][t])[0][index_case] for t in range(T)])
        correlations_RT[netw_realiz] = np.array([scipy.stats.pearsonr(mean_epidemic_RNT[netw_realiz].T[t+1],nodes_Bcentrality_RNT[netw_realiz].T[t])[0] if scipy.stats.pearsonr(mean_epidemic_RNT[netw_realiz].T[t+1],nodes_Bcentrality_RNT[netw_realiz].T[t])[0]==scipy.stats.pearsonr(mean_epidemic_RNT[netw_realiz].T[t+1],nodes_Bcentrality_RNT[netw_realiz].T[t])[0] else 0 for t in range(T-1)])
        
    #SAVES
    Saves.network_save(temporal_networks_RTNN,directory_name,'networks')
    Saves.network_save(infection_matrices_RNITN,directory_name,'infective_state_matrices')
    Saves.network_save(virulence_RNI, directory_name,'VIRULENCE_SCORES')
    Saves.network_save(time_tobe_infected_RNIN, directory_name,'TIMEINFECTED_SCORES')
    Saves.network_save(epidemic_size_evolution_RNIT, directory_name,'Epidemic_size_evolution')
    Saves.network_save(correlations_RT, directory_name,'correlations')
    
    if isDAR:
        return infection_matrices_RNITN,virulence_RNI,time_tobe_infected_RNIN,epidemic_size_evolution_RNIT,nodes_Bcentrality_RNT,mean_epidemic_RNT,correlations_RT
    else:
        return infection_matrices_RNITN,virulence_RNI,time_tobe_infected_RNIN,epidemic_size_evolution_RNIT,nodes_Bcentrality_RNT,mean_epidemic_RNT,correlations_RT
    

def results(N,T, beta, K, alfa, isDAR,P, isDIRECTED, isSAMPLED, NET_REAL,net_name, multiple_infections,eigen_fraction, infective_fraction, fig_count, virulence_RNI, directory_name, time_tobe_infected_RNIN,epidemic_size_evolution_RNIT, correlations_RT,q_param,nodes_Bcentrality_RN,nodes_Rcentrality_RN,nodes_AD_RN,nodes_BD_RN,avgdegree_evolution_RT):
    fig_count+=50
    virulence_scores = [list(1/np.average(virulence_RNI[realization], axis=1)) for realization in range(NET_REAL)] #axis 1 -> averages over epidemic Iterations
    timetobe_scores = [1/np.average([np.average(time_tobe_infected_RNIN[realization][index], axis = 0) for index in range(N)],axis=0) for realization in range(NET_REAL)]
    avg_epid_ev = [np.average([np.average(epidemic_size_evolution_RNIT[realization][index],axis=0) for index in range(N)],axis=0) for realization in range(NET_REAL)]
    avg_correlations_evolution = np.average(correlations_RT,axis=0)

    file_name = directory_name+"/grapichs/"
    os.makedirs(os.path.dirname(file_name), exist_ok=True)
    
#    if isDAR == False:
#        thetas = Saves.network_load(directory_name,'TGRG_thetas')[NET_REAL-1]
#        fig_count+=1
#        plt.figure(fig_count)
#        plt.errorbar(np.linspace(0,T-1,T),np.average(thetas,axis=0),yerr = np.std(thetas,axis=0))
#        plt.grid()
#        plt.title("Evolution of avg (over nodes) fitness with error bars (stdev), last realization")
#        plt.xlabel("Temporal step")
#        plt.ylabel(r"$<\theta>(t)$")
#        plt.savefig(directory_name+"/inputs/TGRG_fitness_evolution_lastrealizat.pdf")
    
    fig_count+=1
    plt.figure(fig_count)
    fig_count+=1
    for netw_realiz in range(NET_REAL):
        plt.plot(np.linspace(0,T-2,T-1),correlations_RT[netw_realiz], label="netw_realiz = "+str(netw_realiz))
    plt.legend()
    plt.grid()
    plt.title(r"q_param = %.2f, beta = 0.1" %q_param)
    plt.xlabel("Temporal step")
    plt.ylabel(r"c(t)")
    plt.savefig(file_name+"qparam"+str(q_param)+"correlation_evolutions.pdf")
    
#    f, axes = plt.subplots(2,2)
#    for realization in range(NET_REAL):
#        axes[0, 0].scatter(nodes_Bcentrality_RN[realization], virulence_scores[realization], c = "blue")
#        axes[0, 1].scatter(nodes_Rcentrality_RN[realization], virulence_scores[realization], c = "blue")
#        axes[1, 0].scatter(nodes_AD_RN[realization], virulence_scores[realization], c = "blue")
#        axes[1, 1].scatter(nodes_BD_RN[realization], virulence_scores[realization], c = "blue")
#    axes[0, 0].set(xlabel="Broadcast Centrality", ylabel="Epidemic score")
#    axes[0, 1].set(xlabel="Receiver Centrality", ylabel="Epidemic score")
#    axes[1, 0].set(xlabel="Aggregate Degree", ylabel="Epidemic score")
#    axes[1, 1].set(xlabel="Binarized Degree", ylabel="Epidemic score")
#    f.suptitle(r"Direct=%s, DAR=%s; $\beta$=%.3f; N=%i, T=%i; Net_itr=%i, Epi_itr=%i" %(str(isDIRECTED),str(isDAR),beta,N,T,NET_REAL,K))
#    plt.savefig(file_name+"Virulece_correlations.pdf")
#    
#    f, axes = plt.subplots(1,2)
#    for realization in range(NET_REAL):
#        axes[0].scatter(nodes_Rcentrality_RN[realization],nodes_Bcentrality_RN[realization], c = "blue")
#        axes[1].scatter(nodes_Rcentrality_RN[realization],timetobe_scores[realization], c = "blue")
#    axes[0].set(xlabel="Receiver Centrality", ylabel="Broadcast Centrality")
#    axes[1].set(xlabel="Receiver Centrality", ylabel="Time to be infected")
#    f.suptitle(r"Direct=%s, DAR=%s; $\beta$=%.3f; N=%i, T=%i; Net_itr=%i, Epi_itr=%i" %(str(isDIRECTED),str(isDAR),beta,N,T,NET_REAL,K))
#    plt.savefig(file_name+"Receiver_correlations.pdf")
#    
#    f, (ax1, ax2) = plt.subplots(1,2)
#    ax1.plot(np.linspace(0,T-1,T), np.average(avgdegree_evolution_RT,axis=0)/N, label="Average degree(t)")
#    ax1.plot(np.linspace(0,T-1,T),min(min(avgdegree_evolution_RT))*np.ones(T),linestyle='--', label="Total min degree")
#    ax1.plot(np.linspace(0,T-1,T),max(max(avgdegree_evolution_RT))*np.ones(T),linestyle='--', label="Total max degree")
#    ax2.plot(np.linspace(0,T-1,T), np.average(avg_epid_ev, axis=0)/N)
#    ax1.legend()
#    ax1.grid()
#    ax1.set(xlabel="Temporal step", ylabel="Networks' avg percent-degree (over realizations)") 
#    ax2.grid()
#    ax2.set(xlabel="Temporal step", ylabel="Average percentual epidemic size (over realizations and iters)")
#    f.suptitle(r"Direct=%s, DAR=%s; $\beta$=%.3f; N=%i, T=%i; Net_itr=%i, Epi_itr=%i" %(str(isDIRECTED),str(isDAR),beta,N,T,NET_REAL,K))
#    plt.savefig(file_name+"Degree_and_virulence_evolutions.pdf")  
#    
#    correlationsBR = [scipy.stats.pearsonr(nodes_Bcentrality_RN[realization], nodes_Rcentrality_RN[realization])[0] for realization in range(NET_REAL)]
#    correlationsBV = [scipy.stats.pearsonr(nodes_Bcentrality_RN[realization], virulence_scores[realization])[0] for realization in range(NET_REAL)]
#    correlationsRV = [scipy.stats.pearsonr(nodes_Rcentrality_RN[realization], virulence_scores[realization])[0] for realization in range(NET_REAL)]
#    correlationsAV = [scipy.stats.pearsonr(nodes_AD_RN[realization], virulence_scores[realization])[0] for realization in range(NET_REAL)]
#    correlationsBinV = [scipy.stats.pearsonr(nodes_BD_RN[realization], virulence_scores[realization])[0] for realization in range(NET_REAL)]
#    correlationsRTim = [scipy.stats.pearsonr(nodes_Rcentrality_RN[realization], timetobe_scores[realization])[0] for realization in range(NET_REAL)]
#    file_name = directory_name+"/resume.txt"
#    os.makedirs(os.path.dirname(file_name), exist_ok=True)
#    with open(file_name, 'w') as f:
#        f.write("### INPUTS ###\n\n")
#        f.write("NAME = " + net_name + "\n\n")
#        f.write("sampled = " + str(isSAMPLED)+ "\n")
#        f.write("dar = " + str(isDAR)+ "\n")
#        f.write("P = %i\n" %P)
#        f.write("directed = " + str(isDIRECTED)+ "\n\n")
#        f.write("N = %i\n" %N)
#        f.write("T = %i\n" %T)
#        f.write("beta = %.2f \n\n" %beta)
#        f.write("NET_REAL = %i\n" %NET_REAL)
#        f.write("K = %i\n\n" %K)
#        f.write("eigen_fraction = %.2f \n" %eigen_fraction)
#        f.write("multiple_infections = " + str(multiple_infections) + "\n")
#        f.write("infective_fraction = %.2f" %infective_fraction)
#        print("\n\n### OUTPUTS ###\n", file=f)
#        print("BC - RC average correlation:", file=f)
#        print("%.5f" %np.average(correlationsBR), file=f)
#        print("BC - vir average correlation:", file=f)
#        print("%.5f" %np.average(correlationsBV), file=f)
#        print("RC - vir average correlation:", file=f)
#        print("%.5f" %np.average(correlationsRV), file=f)
#        print("AD - vir average correlation:", file=f)
#        print("%.5f" %np.average(correlationsAV), file=f)
#        print("BD - vir average correlation:", file=f)
#        print("%.5f" %np.average(correlationsBinV), file=f)
#        print("RC - tim average correlation:", file=f)
#        print("%.5f" %np.average(correlationsRTim), file=f)
#        
#        print("Common nodes in first 10 positions, BCENTR vs VIR, for each network realization", file=f)
#        for k in range(NET_REAL):
#            print(set(np.flip(np.argsort(nodes_Bcentrality_RN[k]))[0:10]).intersection(set(np.flip(np.argsort(virulence_scores[k]))[0:10])), file=f)
#            #Highest BCENTR should meet lowest virulence, which is a score of how much time it takes. So virulence is not flipped
#        
#        print("Common nodes in first 10 positions, RCENTR vs TIM, for each network realization", file=f)
#        for k in range(NET_REAL):
#            print(set(np.flip(np.argsort(nodes_Rcentrality_RN[k]))[0:10]).intersection(set(np.flip(np.argsort(timetobe_scores[k]))[0:10])), file=f)
#            #Highest BCENTR should meet lowest virulence, which is a score of how much time it takes. So virulence is not flipped
#        
#        print("BC - RC correlations:", file=f)
#        print(correlationsBR, file=f)
#        print("BC - vir correlations:", file=f)
#        print(correlationsBV, file=f)
#        print("RC - vir correlations:", file=f)
#        print(correlationsRV, file=f)
#        print("AD - vir correlations:", file=f)
#        print(correlationsAV, file=f)
#        print("BD - vir correlations:", file=f)
#        print(correlationsBV, file=f)
#        print("RC - tim correlations:", file=f)
#        print(correlationsRTim, file=f)
  
    return virulence_scores,timetobe_scores,avg_epid_ev,avg_correlations_evolution, fig_count
#   return virulence_scores,timetobe_scores,avg_epid_ev,avg_correlations_evolution,correlationsBR,correlationsBV,correlationsRV,correlationsAV,correlationsBinV,correlationsRTim, fig_count

def produce_communicRTNN(N,T,NET_REAL, isDIRECTED, q_param, directory_name, temporal_networks_RTNN):
    communicabilities_RTNN = np.zeros((NET_REAL,T,N,N))
    nodes_Bcentrality_RN = []
    nodes_Rcentrality_RN = []
    #communicabilities_QRTNN = {q_param:[[] for n in range(NET_REAL)] for q_param in q_parameters}
    for netw_realiz in range(NET_REAL):
        communicabilities_RTNN[netw_realiz] = Evolutions.communicability(temporal_networks_RTNN[netw_realiz], a = q_param, length_one= not isDIRECTED) #if undirected, use length 1; Evolutions.recipr_max_eigen
    
    nodes_Bcentrality_RN.append(Evolutions.broadcast_ranking(communicabilities_RTNN[netw_realiz][T-1])[0])
    nodes_Rcentrality_RN.append(Evolutions.receive_ranking(communicabilities_RTNN[netw_realiz][T-1])[0])
    
    Saves.network_save(communicabilities_RTNN, directory_name,'Communicabilities_evolutions')
    Saves.network_save(nodes_Bcentrality_RN, directory_name,'BROADCAST_CENTRALITIES')
    Saves.network_save(nodes_Rcentrality_RN, directory_name,'RECEIVER_CENTRALITIES')

    return communicabilities_RTNN, nodes_Bcentrality_RN, nodes_Rcentrality_RN

def produce_inpunt_tempnetRTNN(N,T,alfa,ki,NET_REAL,isDAR,P,isSAMPLED,isDIRECTED,directory_name,fig_count):
    if isDAR:   #DAR NETWORK GENERATION#
        phi0=[];phi1=[];sigma=[] #generate, as empty, the other inputs
        #SAMPLED DATA IMPORT
        with open('Empiric_Data/alphaDAR1.pkl', 'rb') as f:
            alphaDAR1 = pickle.load(f)
        with open('Empiric_Data/chiDAR1.pkl', 'rb') as f:
            chiDAR1 = pickle.load(f)
        
        #PARAMETERS GENERATION
        if isSAMPLED: #extract from empiric distribution
#            alpha = Evolutions.input_sampling(alphaDAR1,N=N,isDAR = isDAR)
#            chi = Evolutions.input_sampling(chiDAR1,N=N,isDAR = isDAR)
            
            #Mi serve un'alpha omogeneo:
            alpha = alfa*(np.ones(shape = (N,N)) - np.diag([1 for i in range(N)]))
            chi = ki*(np.ones(shape = (N,N)) - np.diag([1 for i in range(N)]))
            
        else: #generate from analytic beta distribution
            alpha = Evolutions.beta_distribution(alphaDAR1,(N,N))
            chi = Evolutions.beta_distribution(chiDAR1,(N,N))
        
        #CHECKS AND SAVES
        Assertions_suite.check_is_probability(alpha)
        Assertions_suite.check_is_probability(chi)  
        #Figures
        #alpha
#        if isSAMPLED:
#            fig_count = Evolutions.hist_plots(fig_count,alpha,N,isDAR,"alpha values generated by empiric distribution",directory_name+"/inputs/histALPHA.pdf")
#            fig_count = Evolutions.hist_plots(fig_count,chi,N,isDAR,"chi values generated by empiric distribution",directory_name+"/inputs/histCHI.pdf")
#        else:
#            fig_count = Evolutions.hist_plots(fig_count,alpha,N,isDAR,"alpha values generated by beta distribution",directory_name+"/inputs/histALPHA.pdf")
#            fig_count = Evolutions.hist_plots(fig_count,chi,N,isDAR,"chi values generated by empiric distribution",directory_name+"/inputs/histCHI.pdf")
        #Pickles
        os.makedirs(os.path.dirname(directory_name+"/inputs/alpha.pkl"), exist_ok=True)
        with open(directory_name+"/inputs/alpha.pkl", 'wb') as handle: #wb = write binary
            pickle.dump(alpha,handle)
        os.makedirs(os.path.dirname(directory_name+"/inputs/chi.pkl"), exist_ok=True)
        with open(directory_name+"/inputs/chi.pkl", 'wb') as handle: #wb = write binary
            pickle.dump(chi,handle)
    else:
        alpha=[];chi=[]
                                #TGRG NETWORK GENERATION#
        #SAMPLED DATA IMPORT
        with open('Empiric_Data/phi0TGRG.pkl', 'rb') as f:
            phi0TGRG = pickle.load(f)
        with open('Empiric_Data/phi1TGRG.pkl', 'rb') as f:
            phi1TGRG = pickle.load(f)
        with open('Empiric_Data/sigmaTGRG.pkl', 'rb') as f:
            sigmaTGRG = pickle.load(f)
#        if isDIRECTED == False: #if it's undirected, only first 98 observations are took
#            phi0TGRG = phi0TGRG[1:98]
#            phi1TGRG = phi1TGRG[1:98]
#            sigmaTGRG = sigmaTGRG[1:98]

        #PARAMETERS GENERATION
        if isSAMPLED:
            phi0 = Evolutions.input_sampling(phi0TGRG,N=N,isDAR = isDAR)
            phi1 = Evolutions.input_sampling(phi1TGRG,N=N,isDAR = isDAR)
            for i in range(len(phi1)):
                if phi1[i]>=1:
                    phi1[i] = 0.99
            sigma = Evolutions.input_sampling(sigmaTGRG,N=N,isDAR = isDAR)
        else: #extract from analytic distribution
            phi0 = np.random.poisson(np.average(phi0TGRG), size = N)
            phi1 = np.random.uniform(-0.99,0.99, size = N)
            sigma = np.random.poisson(np.average(sigmaTGRG), size = N)
        
        #CHECKS AND SAVES
        assert abs(min(phi1))>=0
        assert abs(max(phi1))<=1
#        #Figures
#        #alpha
#        if isSAMPLED:
#            fig_count = Evolutions.hist_plots(fig_count,phi0,N,isDAR,r"$\phi_0$ values generated by empiric distribution",directory_name+"/inputs/histPHI0.pdf")
#            fig_count = Evolutions.hist_plots(fig_count,phi1,N,isDAR,r"$\phi_1$ values generated by empiric distribution",directory_name+"/inputs/histPHI1.pdf")
#            fig_count = Evolutions.hist_plots(fig_count,sigma,N,isDAR,r"$\sigma$ values generated by empiric distribution",directory_name+"/inputs/histSIGMA.pdf")
#        else:
#            fig_count = Evolutions.hist_plots(fig_count,phi0,N,isDAR,r"$\phi_0$ values generated by analytic distribution",directory_name+"/inputs/histPHI0.pdf")
#            fig_count = Evolutions.hist_plots(fig_count,phi1,N,isDAR,r"$\phi_1$ values generated by analytic distribution",directory_name+"/inputs/histPHI1.pdf")
#            fig_count = Evolutions.hist_plots(fig_count,sigma,N,isDAR,r"$\sigma$ values generated by analytic distribution",directory_name+"/inputs/histSIGMA.pdf")
        #Pickles        
        os.makedirs(os.path.dirname(directory_name+"/inputs/phi0.pkl"), exist_ok=True)
        with open(directory_name+"/inputs/phi0.pkl", 'wb') as handle: #wb = write binary
            pickle.dump(phi0,handle)
        os.makedirs(os.path.dirname(directory_name+"/inputs/phi1.pkl"), exist_ok=True)
        with open(directory_name+"/inputs/phi1.pkl", 'wb') as handle: #wb = write binary
            pickle.dump(phi1,handle)
        os.makedirs(os.path.dirname(directory_name+"/inputs/sigma.pkl"), exist_ok=True)
        with open(directory_name+"/inputs/sigma.pkl", 'wb') as handle: #wb = write binary
            pickle.dump(sigma,handle)

    #TEMPNETS GENERATION
    temporal_networks_RTNN = []
    for netw_realiz in range(NET_REAL):
        if isDAR: #use the proper functiond wheter user selected dar or tgrg
            temporal_network = Evolutions.network_generation_dar(alpha,chi,P=P,T=T,directed=isDIRECTED) 
        else:
            temporal_network, thetas = Evolutions.network_generation_tgrg(phi0,phi1,sigma,T=T,directed=isDIRECTED)
        temporal_networks_RTNN.append(temporal_network)
    return alpha,chi,phi0,phi1,sigma,temporal_networks_RTNN,fig_count