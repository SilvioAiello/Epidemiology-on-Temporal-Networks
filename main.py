"""
This script manages all the others, allowing user to actually create and save networks and epidemics, and to perform analysis.

Functions work in Pyhon3, and may require the following libraries (so, check if they are installed):
    * numpy, used for its data structures and anaylisis, and to get random functions 
    * pickle, used to store, in an efficient way, the complex information generated
    * function "quad" from scipy.integrate
[If you want to get some plots, you may use matplotlib.pyplot, for plots belluries, and networkx, to plot small networks]

#TODO: Stucture?

For further understandings on how this script operates, check file "howto.md"
For further theoretical understandings, check file "explanation.md"
"""

import numpy as np
from scipy.integrate import quad #used in dictionary of probabilities

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

#%%             USER ACTION        ###
# PARAMETERS SETTING #
N=10 #nodes
T=100 #duration
K = 5 #infection iterations
beta = 0.005 #infection rate
    #dar inputs
alpha = 0.6*np.ones((N,N)) #you can change values, but keep it (N,N)
xi = 0.5*np.ones((N,N)) #you can change value, but keep it (N,N)
    #tgrg inputs
phi0 = 0.5*np.ones(N)
phi1 = 0.5*np.ones(N)
epsilon=0.5*np.ones(N)

    #Probabilities dict, according to beta (same for DAR and TGRG)
def poisson_probability(t): #function definition
    """
    This function reproduces the Poisson PDF, whose average depends on beta.
    Its integral is the probability of having and infection for a I-S contact lasting t.
    
    Parameter: t (int), representing time duration of a contact.
    Returns: a float, representing PDF for that t.
    """
    lam = -np.log(1-beta) # Chen, Benzi use 60
    return(lam*np.exp(-lam*t))
probabilities = dict() #dict building
for t in range(T):
    probabilities[t] = quad(poisson_probability,0,t)[0]

# TEMPNETS GENERATION AND SAVE
temporal_dar = Evolutions.network_generation_dar(alpha,xi,T=T,directed=False) #tempnet generation
#andrebbe usata la funzione di salvataggio

# SI PROPAGATION
label_dar = [] 
#TODO: AGGIUNGERE TGRG OVUNQUE
for index_case in range(N):
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

## TGRG #
#spec_radius_tgrg, Q_tgrg = Evolution_DAR.communicability(temporal_tgrg)
#nodes_Bcentrality_tgrg, nodes_Brank_tgrg = Evolutions.broadcast_ranking(Q_tgrg)
#nodes_Rcentrality_tgrg, nodes_Rrank_tgrg = Evolutions.receive_ranking(Q_tgrg)

#%%                     VIRULENCE MEASURES                            ###
#DAR #
score_dar = [] #list of N lists (one per index case), each of which is a list of K (one per iteration) floats (score of that iteration)
for index_case in range(N):
    score_dar.append([]) #create the i-th entry
    for iteration in range(K):
        score_dar[index_case].append(Propagation_SI.time_score(label_dar[index_case][iteration],0.6))