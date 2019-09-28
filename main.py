"""
From this script, user can generate temporal networks and perform epidemic upon them.
System parameters must be set in "inputs.ini" file.
Functions in this script work in Pyhon3, may require numpy (v1.16) and function "quad" from scipy.integrate (scipy v1.3).

For further understandings on how this script operates, check file "howto.md".
For further theoretical understandings, check file "explanation.md".
"""
import numpy as np
import scipy.stats
from scipy.integrate import quad #used in dictionary of probabilities

import Evolutions
import Propagation_SI
import Saves

import pickle
import os

import matplotlib.pyplot as plt
fig_count = 2

#%%         CONFIGURATION PARAMETERS
import configparser
config = configparser.ConfigParser()
config.read('inputs.ini')
for section in config.sections():
    net_name = config[section]['net_name']
    
    N = config[section].getint('N')
    T = config[section].getint('T')
    beta = config[section].getfloat('beta')    
    
    isSAMPLED = config[section].getboolean('isSAMPLED')
    isDAR = config[section].getboolean('isDAR')
    P = config[section].getint('P')
    isDIRECTED = config[section].getboolean('isDIRECTED')
    
    NET_REAL = config[section].getint('NET_REAL')
    K = config[section].getint('K')
    
    eigen_fraction = config[section].getfloat('eigen_fraction')
    multiple_infections = config[section].getboolean('multiple_infections')
    infective_fraction = config[section].getfloat('infective_fraction')
#%%                         INPUTS BUILDING AND SAVE
    directory_name = Saves.main_directory(beta, isDIRECTED,isDAR,P,isSAMPLED,N,T,net_name)
    
    if isDAR:
        #SAMPLED DATA WILL BE USED ANYWAY
        with open('Empiric_Data/alphaDAR1.pkl', 'rb') as f:
            alphaDAR1 = pickle.load(f)
        with open('Empiric_Data/chiDAR1.pkl', 'rb') as f:
            chiDAR1 = pickle.load(f)
            
        if isSAMPLED: #extract from empiric distribution
            alpha = Evolutions.input_sampling(alphaDAR1,N=N,isDAR = isDAR)
            chi = Evolutions.input_sampling(chiDAR1,N=N,isDAR = isDAR)
        else: #generate from analytic distribution
            
    #        mu_alpha = np.average(alphaDAR1)
    #        sigma_alpha= np.std(alphaDAR1)
    #        def lognorm(x,mu,sigma): #my attempt
    #            return np.exp(-((np.log(x)-mu)**2)/(2*sigma**2))/x
    #        factor = quad(lognorm,0,1, args=(mu_alpha, sigma_alpha))[0]
    #        prob = [quad(lognorm,x-0.005,x+0.005, args=(mu_alpha, sigma_alpha))[0]/factor for x in np.linspace(0.01,0.99,100)]
    #        
    #        plt.plot(np.linspace(0.01,0.99,100),prob)
    #        plt.xlabel(r"$\alpha_{ij}$ possible values")
    #        plt.ylabel("Probability")
    #        plt.title(r"$\mu =$ %.3f, $\sigma =$ %.3f"%(mu_alpha,sigma_alpha))
    #        plt.show()
            
            #truncated gaussians
            lower = 0
            upper = 1
            N=1000
            mu = np.mean(alphaDAR1)
            sigma = np.std(alphaDAR1)
            alpha = scipy.stats.truncnorm.rvs((lower-mu)/sigma,(upper-mu)/sigma,loc=mu,scale=sigma,size=(N,N))
            plt.figure(1)
            plt.hist(alpha.reshape(N*N,1))
            plt.title(r"Generazione di valori per $\alpha$ da Gaussiana troncata, $\mu$ e $\sigma$ uguali alla distr $\alpha$ empiriche, num = %i$^2$" %N)
            plt.xlabel(r"$\alpha_{ij}$")
            plt.ylabel("Occorrenze")
            
            mu = np.mean(chiDAR1)
            sigma = np.std(chiDAR1)
            chi = scipy.stats.truncnorm.rvs((lower-mu)/sigma,(upper-mu)/sigma,loc=mu,scale=sigma,size=(N,N))
            plt.figure(2)
            plt.hist(chi.reshape(N*N,1))
            plt.title(r"Generazione di valori per $\chi$ da Gaussiana troncata, $\mu$ e $\sigma$ uguali alla distr $\chi$ empiriche, num = %i$^2$" %N)
            plt.xlabel(r"$\chi_{ij}$")
            plt.ylabel("Occorrenze")
        
        os.makedirs(os.path.dirname(directory_name+"/alpha.pkl"), exist_ok=True)
        with open(directory_name+"/alpha.pkl", 'wb') as handle: #wb = write binary
            pickle.dump(alpha,handle)
        os.makedirs(os.path.dirname(directory_name+"/chi.pkl"), exist_ok=True)
        with open(directory_name+"/chi.pkl", 'wb') as handle: #wb = write binary
            pickle.dump(chi,handle)
    else:
        #SAMPLED DATA WILL BE USED ANYWAY
        with open('Empiric_Data/phi0TGRG.pkl', 'rb') as f:
            phi0TGRG = pickle.load(f)
        with open('Empiric_Data/phi1TGRG.pkl', 'rb') as f:
            phi1TGRG = pickle.load(f)
        with open('Empiric_Data/sigmaTGRG.pkl', 'rb') as f:
            sigmaTGRG = pickle.load(f)
        if isDIRECTED == False: #if it's undirected, only first 98 observations are took
            phi0TGRG = phi0TGRG[1:98]
            phi1TGRG = phi0TGRG[1:98]
            sigmaTGRG = sigmaTGRG[1:98]
        
        if isSAMPLED:
            #PHI0S
            phi0 = Evolutions.input_sampling(phi0TGRG,N=N,isDAR = isDAR)
            #PHI1S
            phi1 = Evolutions.input_sampling(phi1TGRG,N=N,isDAR = isDAR)
            #SIGMAS
            sigma = Evolutions.input_sampling(sigmaTGRG,N=N,isDAR = isDAR)
        else: #extract from analytic distribution
            phi0 = np.random.poisson(3, size = N)
            phi1 = np.random.uniform(-0.99,0.99, size = N)
            sigma = np.random.poisson(1, size = N)
        
        os.makedirs(os.path.dirname(directory_name+"/phi0.pkl"), exist_ok=True)
        with open(directory_name+"/phi0.pkl", 'wb') as handle: #wb = write binary
            pickle.dump(phi0,handle)
        os.makedirs(os.path.dirname(directory_name+"/phi1.pkl"), exist_ok=True)
        with open(directory_name+"/phi1.pkl", 'wb') as handle: #wb = write binary
            pickle.dump(phi1,handle)
        os.makedirs(os.path.dirname(directory_name+"/sigma.pkl"), exist_ok=True)
        with open(directory_name+"/sigma.pkl", 'wb') as handle: #wb = write binary
            pickle.dump(sigma,handle)
    
#%%                         OUTPUTS GENERATION
    import time
    start = time.time()
    
    #Following lists contain NET_REAL arrays; each of these lists has N entries, with the average scores for each node
    nodes_Bcentrality = []
    nodes_Rcentrality = []
    nodes_AD = []
    nodes_BD = []
    
    virulence = [] #scores computed when each node is the index case
    time_tobe_infected = [] #scores computed when each node is not the index case
        
    for k in range(1,NET_REAL+1):  #so first realization has index 1; for each k one tempnet and one k are overwritten
        #TEMPNETS GENERATION AND SAVE
        if isDAR: #use the proper functiond wheter user selected dar or tgrg
            temporal_network = Evolutions.network_generation_dar(alpha,chi,P=P,T=T,directed=isDIRECTED) 
        else:
            temporal_network = Evolutions.network_generation_tgrg(phi0,phi1,sigma,T=T,directed=isDIRECTED)[0]
        Saves.network_save(temporal_network,directory_name,k,'network')
        
        #CENTRALITIES GENERATION AND SAVE
        inv_maxradius, Q = Evolutions.communicability(temporal_network, eigen_fraction=eigen_fraction, length_one= not isDIRECTED) #if is undirected, use length 1
        
        singleiter_nodes_Bcentrality = Evolutions.broadcast_ranking(Q)[0]
        Saves.network_save(singleiter_nodes_Bcentrality,directory_name,k,'BCENTR')
        nodes_Bcentrality.append(singleiter_nodes_Bcentrality)
        
        singleiter_nodes_Rcentrality = Evolutions.receive_ranking(Q)[0]
        Saves.network_save(singleiter_nodes_Rcentrality, directory_name,k,'RCENTR')
        nodes_Rcentrality.append(singleiter_nodes_Rcentrality)   
        
        singleiter_nodes_AD = Evolutions.aggregate_degree(temporal_network, directed=False) #if dir = True, there are 2 outputs
        Saves.network_save(singleiter_nodes_AD, directory_name,k,'AGGDEG')
        nodes_AD.append(singleiter_nodes_AD)
    
        singleiter_nodes_BD = Evolutions.binarized_degree(temporal_network, directed=False) #if dir = True, there are 2 outputs
        Saves.network_save(singleiter_nodes_BD, directory_name,k,'BINDEG')
        nodes_BD.append(singleiter_nodes_BD)
        
        #SI PROPAGATION
        probabilities = dict() #probabilities dict
        for t in range(T):
            probabilities[t] = quad(Propagation_SI.poisson_probability,0,t, args = beta)[0] #quad produces several outputs, integral is the first
        label = []
        singleiter_virulece = []
        singleiter_time_tobe_infected = []
        for index_case in range(N):
            label.append([]) #create the i-th entry
            for iteration in range(K):
                label[index_case].append(Propagation_SI.propagation(temporal_network, index_case, probabilities, multiple_infections = multiple_infections))
            singleiter_virulece.append(np.mean([Propagation_SI.time_score(label[index_case][iteration],infective_fraction) for iteration in range(K)]))
            singleiter_time_tobe_infected.append(np.mean([Propagation_SI.when_is_infected(label[index_case][iteration], index_case) for iteration in range(K)],axis=0))
        virulence.append(singleiter_virulece)
        time_tobe_infected.append(singleiter_time_tobe_infected)
        #SAVES
        Saves.network_save(label,directory_name,k,'infections')
        Saves.network_save(virulence, directory_name,k,'VIRULENCE')
        Saves.network_save(time_tobe_infected, directory_name,k,'TIMEINFECTED')
    end = time.time()
    print(end-start)  
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
    vir = [1/i for i in vir] #so lowest times have correctly greatest score
    time_tobe_infected = np.mean(time_tobe_infected,axis=1) #averaged over all the iterations
    tim = []
    [tim.extend(el) for el in time_tobe_infected] 
    tim = [1/i for i in tim] #so lowest times have correctly greatest score
    
    
#%% RESULTS OUTPUT
    file_name = directory_name+"/grapichs/"
    os.makedirs(os.path.dirname(file_name), exist_ok=True)
    
    plt.figure(fig_count)
    fig_count+=1
    plt.scatter(nodes_B, vir)
    plt.xlabel("Broadcast Centrality")
    plt.ylabel('Epidemic score')
    plt.title(r"Direct = %s, DAR = %s; $\beta$ = %.3f; N=%i, T=%i; Netw iter = %i, Epid iter = %i" %(str(isDIRECTED),str(isDAR),beta,N,T,NET_REAL,K))
    plt.savefig(file_name+"BC.pdf")
    
    plt.figure(fig_count)
    fig_count+=1
    plt.scatter(nodes_R, vir)
    plt.xlabel("Receive Centrality")
    plt.ylabel('Epidemic score')
    plt.title(r"Direct = %s, DAR = %s; $\beta$ = %.3f; N=%i, T=%i; Netw iter = %i, Epid iter = %i" %(str(isDIRECTED),str(isDAR),beta,N,T,NET_REAL,K))
    plt.savefig(file_name+"RC.pdf")
    
    plt.figure(fig_count)
    fig_count+=1
    plt.scatter(nod_A, vir)
    plt.xlabel("Aggregate Degree")
    plt.ylabel('Epidemic score')
    plt.title(r"Direct = %s, DAR = %s; $\beta$ = %.3f; N=%i, T=%i; Netw iter = %i, Epid iter = %i" %(str(isDIRECTED),str(isDAR),beta,N,T,NET_REAL,K))
    plt.savefig(file_name+"AD.pdf")
    
    plt.figure(fig_count)
    fig_count+=1
    plt.scatter(nod_B, vir)
    plt.xlabel("Binarized Degree")
    plt.ylabel('Epidemic score')
    plt.title(r"Direct = %s, DAR = %s; $\beta$ = %.3f; N=%i, T=%i; Netw iter = %i, Epid iter = %i" %(str(isDIRECTED),str(isDAR),beta,N,T,NET_REAL,K))
    plt.savefig(file_name+"BD.pdf")
    
    plt.figure(fig_count)
    fig_count+=2
    plt.scatter(nodes_R, tim)
    plt.xlabel("Receive Centrality")
    plt.ylabel('Time to be infected')
    plt.title(r"Direct = %s, DAR = %s; $\beta$ = %.3f; N=%i, T=%i; Netw iter = %i, Epid iter = %i" %(str(isDIRECTED),str(isDAR),beta,N,T,NET_REAL,K))
    plt.savefig(file_name+"TIM.pdf")
    
    file_name = directory_name+"/resume.txt"
    os.makedirs(os.path.dirname(file_name), exist_ok=True)
    with open(file_name, 'w') as f:
        f.write("### INPUTS ###\n\n")
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
        f.write("eigen_fraction = %.2f \n" %eigen_fraction)
        f.write("multiple_infections = " + str(multiple_infections) + "\n")
        f.write("infective_fraction = %.2f" %infective_fraction)
        print("\n\n### OUTPUTS ###\n", file=f)
        print("BC - vir correlation:", file=f)
        print(scipy.stats.pearsonr(nodes_B, vir), file=f)
        print("RC - vir correlation:", file=f)
        print(scipy.stats.pearsonr(nodes_R, vir), file=f)
        print("AD - vir correlation:", file=f)
        print(scipy.stats.pearsonr(nod_A, vir), file=f)
        print("BD - vir correlation:", file=f)
        print(scipy.stats.pearsonr(nod_B, vir), file=f)
        print("RC - tim correlation:", file=f)
        print(scipy.stats.pearsonr(nodes_R, tim), file=f)
        
        print("Common nodes in first 10 positions, BCENTR vs VIR, for each network realization", file=f)
        for k in range(1,NET_REAL+1):
            print(set(np.flip(np.argsort(nodes_Bcentrality[k-1]))[0:10]).intersection(set(np.argsort(virulence[k-1])[0:10])), file=f)
            #Highest BCENTR should meet lowest virulence, which is a score of how much time it takes. So virulence is not flipped
        
        print("Common nodes in first 10 positions, RCENTR vs TIM, for each network realization", file=f)
        for k in range(1,NET_REAL+1):
            print(set(np.flip(np.argsort(nodes_Rcentrality[k-1]))[0:10]).intersection(set(np.argsort(time_tobe_infected[k-1])[0:10])), file=f)
            #Highest BCENTR should meet lowest virulence, which is a score of how much time it takes. So virulence is not flipped