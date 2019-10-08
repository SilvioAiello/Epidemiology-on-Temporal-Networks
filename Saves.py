"""
Functions in this script work in Pyhon3, require os, pickle, and allow to save results.

In this module are defined the following functions:
    * network_save and network_load
    * plot save

For further understandings on how this script operates, check file "docs/howto.md"
"""
import pickle
import os

def main_directory(beta, isDIRECTED,isDAR,P,isSAMPLED,N,T,title):
    assert type(beta) == float
    assert type(isDIRECTED) == bool
    assert type(isDAR) == bool
    assert type(P) == int
    assert type(isSAMPLED) == bool
    assert type(N) == int
    assert type(T) == int
    assert type(title) == str
    
    name = "Networks/beta"+str(beta)
    if isDAR:
        name = name + "_DAR" + str(P)
    else:
        name+= "_TGRG"
    
    if isSAMPLED:
        name+="_SAMPLED"
    else:
        name+="_UNSAMPLED"
    
    if isDIRECTED:
        name+="_DIRECTED"
    else:
        name+="_UNDIRECTED"
    
    name = name + "_N"+str(N)+ "_T"+str(T) + "_" + title
    
    return name

def network_save(network, directory_name, k, data_structure_name):
    """ Saves data structures using pickle. Be sure to provide the proper name:
        
    
    Parameters
    ----------
    network: np.array or list
        The data structure you want to save
    directory_name: string
        PATH of simulation directory (es. Networks/betaX_DAR_DIRECTED ecc)
    k: int 
        Number of network's realization
    data_structure_name: str
        One of these: alpha/chi/phi0/1/sigma; 
        AGGDEG/BINDEG,BCENTR/RCENTR, 
        TIMEINRECTED/VIRULECE; 
        infections/network
    
    Returns
    -------
    /directory_name/realizationk/data_structure_name.pkl
        If PATH doesn't exist, it's created
    """
    name = directory_name + "/realization"+str(k)+"/"+data_structure_name+".pkl"
    os.makedirs(os.path.dirname(name), exist_ok=True)
    with open(name, 'wb') as handle: #wb = write binary
        pickle.dump(network,handle)

def network_load(directory_name, k,data_structure):
    """ 
    Loads a previously generated temporal network using pickle (so it must be installed), if it's present.
    Network should be saved in "Networks/directory_name/realizationk/network.pkl
    
    Parameters
    ----------
    network: np.array or list
        The data structure you want to save
    directory_name: string
        PATH of simulation directory (es. Networks/betaX_DAR_DIRECTED ecc)
    k: int 
        Number of network's realization
    data_structure_name: str
        One of these: alpha/chi/phi0/1/sigma; 
        AGGDEG/BINDEG,BCENTR/RCENTR, 
        TIMEINRECTED/VIRULECE; 
        infections/network
    
    Returns
    -------
    pickle.load(f): np.array
        Data structure is returned
    """
    name = directory_name + "/realization"+str(k)+"/"+data_structure+".pkl"
    with open(name, 'rb') as f:
        return pickle.load(f)

