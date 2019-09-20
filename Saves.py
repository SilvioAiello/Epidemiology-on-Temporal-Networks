"""
Functions in this script work in Pyhon3, require os, pickle, and allow to save results.

In this module are defined the following functions:
    * network_save and network_load
    * plot save

For further understandings on how this script operates, check file "docs/howto.md"
"""
import pickle
import os

import Assertions_suite

def make_basics(isDAR=True,P=1,isDIRECTED=False,isSAMPLED=True):
    name = str()
    if isDAR:
        name = "Networks/DAR" + str(P)
    else:
        name = "Networks/TGRG"
    
    if isSAMPLED:
        name+="_SAMPLED"
    else:
        name+="_UNSAMPLED"
    
    if isDIRECTED:
        name+="_DIRECTED"
    else:
        name+="_DIRECTED"
    return name

def network_save(network, start,isDAR=True,isDIRECTED=False,isSAMPLED=True,k=1, P=1):
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
    assert Assertions_suite.check_is_natural(P) #there's no need to perform other checks, since they have been already performed
    assert Assertions_suite.check_is_natural(k)
    assert type(start) == str
    
    #FUNCTION
    T = network.shape[0]
    N = network.shape[1]
    
    name = make_basics(isDAR=isDAR,P=P,isDIRECTED=isDIRECTED,isSAMPLED=isSAMPLED) + "_N"+str(N)+"_T"+str(T)+"_"+start+"/realization"+str(k)+"/network.pkl"

    os.makedirs(os.path.dirname(name), exist_ok=True)
    with open(name, 'wb') as handle: #wb = write binary
        pickle.dump(network,handle)

def network_load(N,T,start,isDAR=True,isDIRECTED=False,isSAMPLED=True,k=1,P=1):
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

    name = make_basics(isDAR=isDAR,P=P,isDIRECTED=isDIRECTED,isSAMPLED=isSAMPLED) + "_N"+str(N)+"_T"+str(T)+"_"+start+"/realization"+str(k)+"/network.pkl"

    with open(name, 'rb') as f:
        return pickle.load(f)

def analysis_save(centr,centr_name, start, N,T,isDAR=True,isDIRECTED=False,isSAMPLED=True,k=1,P=1):
    """
    Here N and T need to be specified
    
    Returns
    -------
    Saves an array of length N, containing the scores for each node
    """
    assert type(centr_name) == str
    
    name = make_basics(isDAR=isDAR,P=P,isDIRECTED=isDIRECTED,isSAMPLED=isSAMPLED) + "_N"+str(N)+"_T"+str(T)+"_"+start+"/realization"+str(k)+"/"+centr_name+".pkl"
    
    with open(name, 'wb') as handle: #wb = write binary
        pickle.dump(centr,handle)

def analysis_load(centr, centr_name, start, N,T,isDAR=True,isDIRECTED=False,isSAMPLED=True,k=1,P=1):
    """
    Here N and T need to be specified
    
    Returns
    -------
    Saves an array of length N, containing the scores for each node
    """
    assert type(centr_name) == str
    
    name = make_basics(isDAR=isDAR,P=P,isDIRECTED=isDIRECTED,isSAMPLED=isSAMPLED) + "_N"+str(N)+"_T"+str(T)+"_"+start+"/realization"+str(k)+"/"+centr_name+".pkl"
    
    with open(name, 'wb') as handle: #wb = write binary
        return pickle.load(handle)

def infection_save(label, N,T, beta, start,isDAR=True,isDIRECTED=False,isSAMPLED=True,k=1, P=1):
    """ Saves an infection data structure, defined by its different iterations, using pickle (so it must be installed).
    A proper name and folder is provided for this strucure.
    
    Parameters
    ----------
    label: NKTN data structure
        epidemic spread iterations
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
    /PATH/infections.pkl
        If PATH didn't exist, it's created
    """
    
    name = make_basics(isDAR=isDAR,P=P,isDIRECTED=isDIRECTED,isSAMPLED=isSAMPLED) + "_N"+str(N)+"_T"+str(T)+"_"+start+"/realization"+str(k)+"/infections_beta"+str(beta)+".pkl"

    os.makedirs(os.path.dirname(name), exist_ok=True)
    with open(name, 'wb') as handle: #wb = write binary
        pickle.dump(label,handle)

def infection_load(N,T, beta, start,isDAR=True,isDIRECTED=False,isSAMPLED=True,k=1, P=1):
    """ Loads a previously generated infection data structure, defined by its different iterations, using pickle (so it must be installed).
    
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
    /PATH/infections.pkl
        If PATH didn't exist, it's created
    """
    name = make_basics(isDAR=isDAR,P=P,isDIRECTED=isDIRECTED,isSAMPLED=isSAMPLED) + "_N"+str(N)+"_T"+str(T)+"_"+start+"/realization"+str(k)+"/infections_beta"+str(beta)+".pkl"
    os.makedirs(os.path.dirname(name), exist_ok=True)
    with open(name, 'rb') as f:
        return pickle.load(f)