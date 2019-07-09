import numpy as np

import matplotlib.pyplot as plt
import networkx as nx

import Evolutions

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
alpha = 0.6*np.ones((10,10)) #give the shape you want but respect the rule
xi = 0.5*np.ones((10,10)) #give the shape you want but respect the rule
temporal_dar = Evolutions.network_generation_dar(alpha,xi)


#%%                      CENTRALITIES MEASURES                      ###
# DAR #
spec_radius_dar, Q_dar = Evolutions.communicability(temporal_dar)
nodes_Bcentrality_dar, nodes_Brank_dar = Evolutions.broadcast_ranking(Q_dar) #scores, node rankings
nodes_Rcentrality_dar, nodes_Rrank_dar = Evolutions.receive_ranking(Q_dar) #cores, node rankings

## FITN #
#spec_radius_fitn, Q_fitn = Evolution_DAR.communicability(temporal_fitn)
#nodes_Bcentrality_fitn, nodes_Brank_fitn = Evolutions.broadcast_ranking(Q_fitn)
#nodes_Rcentrality_fitn, nodes_Rrank_fitn = Evolutions.receive_ranking(Q_fitn)