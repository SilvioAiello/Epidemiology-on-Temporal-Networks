"""
Functions in this script work in Pyhon3, require os, pickle, and allow to save results.

In this module are defined the following functions:
    * network_save and network_load
    * plot save

For further understandings on how this script operates, check file "docs/howto.md"
"""
import pickle
import os

import Test_suite

def network_save(network, start,isDAR=True,k=1, P=1):
    """ Saves network using pickle (so it must be installed) and giving it its proper name (with an identification provided by user) and folder
    
    Parameters
    ----------
    network: np.array
        T*N*N tensor (the functions extracts by itself T and N)
    start: str
        Name you choose for the function
    isDAR: bool (default = True)
        User is required to specify wheter network is DAR or TGRG
    P: int (default = 1)
        DAR order. If network is TGRG, its value doesn't affect the result
    k: int (default = 1)
        Iteration of the same network this is
    
    Returns
    -------
    /PATH/network.pkl
        If PATH didn't exist, it's created
    """
    #ASSERTS
    Test_suite.assert_natural(P) #there's no need to perform other checks, since they have been already performed
    Test_suite.assert_natural(k)
    assert type(start) == str

    #FUNCTION
    T = network.shape[0]
    N = network.shape[1]
    
    name = str()
    if isDAR:
        name = "Networks/N"+str(N)+"_T"+str(T)+"_DAR"+str(P)+"_"+start+"/realization"+str(k)+"/network.pkl"
    else:
        name = "Networks/N"+str(N)+"_T"+str(T)+"_TGRG_"+start+"/realization"+str(k)+"/network.pkl"
    os.makedirs(os.path.dirname(name), exist_ok=True)
    with open(name, 'wb') as handle: #wb = write binary
        pickle.dump(network,handle)

def network_load(N,T,start,isDAR=True,k=1,P=1):
    """ Loads a previously generated temporal network using pickle (so it must be installed), if it's present
    
    Parameters
    ----------
    N: int
        number of nodes
    T: int
        temporal duration
    start: string
        identificative name of the network
    k: int (default = 1)
        iteration of network realization
    isDAR: bool (default = True)
        whether to search a DAR(P) or TGRG
    P: int (default = 1)
        DAR order
    
    Returns
    -------
    pickle.load(f): np.array
        TNN-tempnet is returned
    """
    name = str()
    if isDAR:
        name = "Networks/N"+str(N)+"_T"+str(T)+"_DAR"+str(P)+"_"+start+"/realization"+str(k)+"/network.pkl"
    else:
        name = "Networks/N"+str(N)+"_T"+str(T)+"_TGRG_"+start+"/realization"+str(k)+"/network.pkl" 
    with open(name, 'rb') as f:
        return pickle.load(f)
#temporal_dar = network_load(100,100,'alphaeqs_xieqs',k=1,isDAR=True,P=1)
#temporal_fitn= load_network(100,100,isDAR=True,alleq,k=1)