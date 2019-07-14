import numpy as np

import matplotlib.pyplot as plt
import networkx as nx

import Evolutions
import Propagation_SI
#%% PLOT FUNCTIONS, MANEGGIARE CON CURA#

def networkplot(graph,figcount,savefig=False, figname = None, k = None, t = None):
    """Requires an adiacency[t] and the figure count
    If you want to save the image, you have to say savefig=True, and pass a 2-ple and the value of k and t.
    Tuple figname contains the values needed by make_name: N,T,isDAR, start,P (P MUST BE SPECIFIED HERE)"""
    figcount+=1
    plt.figure(figcount)
    plt.title("")
    nxgraph = nx.from_numpy_matrix(graph) #convert np-matrix in a Network, readable by nx library:
    nx.draw(nxgraph, pos = nx.drawing.layout.circular_layout(nxgraph), with_labels = True) #draw a circular plot and show the number of each node
    if savefig and figname and k and t:
        plt.savefig(figname+"Results"+'_realization%i_t=%i.pdf' %(k,t))
    return figcount+1, plt.show()

#Save the network for further use
#np.savetxt(start_name+'%i_N%i_wholenetwork_T%i.txt' %(P,N,T), temporal_network.reshape(T*N*N,1))
#To import:
#new_data = np.loadtxt('start_name+'%i__N%i_wholenetwork_T%i.txt' %(P,N,T))
#new_data = new_data.reshape((T,N,N))


#%%             USER ACTION        ###
#Tempnet generation
N=10 #nodes
T=100 #duration
alpha = 0.6*np.ones((N,N)) #give the shape you want but respect the rule
xi = 0.5*np.ones((N,N)) #give the shape you want but respect the rule
temporal_dar = Evolutions.network_generation_dar(alpha,xi,T=T,directed=False)
#andrebbe usata la funzione di salvataggio

#TODO: troppe parole ci sono
K = 5 #to generate more tempnet is up to you; but iterating propagation for nodes is "compulsory": scores MUST be averages
beta = 0.005 #infection rate (probability of infecting a node within unit time)

#CREANDO DIZIONARIO#
from scipy.integrate import quad #this function performs integration
def poisson_probability(t):
    #This function reproduces the Poisson PDF, where the average is set in function of beta, who is a parameter of the disease
    lam = -np.log(1-beta) #See Chen (who uses 60)
    return(lam*np.exp(-lam*t))
#I = quad(poisson_probability,0,np.inf) #verify the correct normalization
#The same integrals will be performed many times, varying at most the total duration, whose maximum value can be T.
#So, they are computed here once for all, and the results stored in a dict, that maps each duration to the integral
probabilities = dict()
for t in range(T):
    probabilities[t] = quad(poisson_probability,0,t)[0]
# FINE DINIZIONARIO

#You also must generate, before, a dict of probabilities, which is totally defined once beta (and T) is set
label_dar = [] #list of N lists (one per index case), each of which is a list of K (one per iteration) sequences of T dictionaries
#(forse anche sti commenti lunghi li puoi mettere direttamente in documentazione)
#label dar 0: tutto riferito allo 0 nodo come index, label_dar[0,3]: 3 iterazione
for index_case in range(N):
    print("Processing node %i" %index_case)
    label_dar.append([]) #create the i-th entry
    for iteration in range(K):
        label_dar[index_case].append(Propagation_SI.propagation(temporal_dar, index_case, probabilities))
        #TEST#
        assert label_dar[index_case][iteration][0][index_case] == 1, "L'index case non sembra esserlo"
        assert sum(label_dar[index_case][iteration][0].values()) == 1, "Ci dovrebbe essere solo un infetto all'inizio"

assert [[label_dar[index_case][iteration][0] == label_dar[index_case][iteration-1][0] for iteration in range(1,K)] for index_case in range(N)], "Error: some initial condition is not equal for all iterations" 

#%%                      CENTRALITIES MEASURES                      ###
# DAR #
spec_radius_dar, Q_dar = Evolutions.communicability(temporal_dar)
nodes_Bcentrality_dar, nodes_Brank_dar = Evolutions.broadcast_ranking(Q_dar) #scores, node rankings
nodes_Rcentrality_dar, nodes_Rrank_dar = Evolutions.receive_ranking(Q_dar) #cores, node rankings

## FITN #
#spec_radius_fitn, Q_fitn = Evolution_DAR.communicability(temporal_fitn)
#nodes_Bcentrality_fitn, nodes_Brank_fitn = Evolutions.broadcast_ranking(Q_fitn)
#nodes_Rcentrality_fitn, nodes_Rrank_fitn = Evolutions.receive_ranking(Q_fitn)

#%%                     VIRULENCE MISURES                            ###
#DAR #
score_dar = [] #list of N lists (one per index case), each of which is a list of K (one per iteration) floats (score of that iteration)
for index_case in range(N):
    score_dar.append([]) #create the i-th entry
    for iteration in range(K):
        score_dar[index_case].append(Propagation_SI.time_score(label_dar[index_case][iteration],0.6))