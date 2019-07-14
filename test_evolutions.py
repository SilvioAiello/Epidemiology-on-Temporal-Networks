"""
Functions in this script provide:
    1) tests to verify that functions in "Evolutions.py" and "Propagation_SI.py" run correctly.
    2) assertions functions to use at the beginning of those functions, to be sure their inputs are correct.

Functions work in Pyhon3, and may require the following libraries (so, check if they are installed):
    * numpy, used for its data structures and anaylisis, and to get random functions   

From this module you can extract the following functions:
#TODO: UPDATE

For further understandings on how this script operates, check file "howto.md"      
"""
import numpy as np

import Evolutions
import Propagation_SI
#%% FUNCTIONS DEFINITIONS
def check_symmetry(temporal_network):
    """
    Verifies if a matrix is simmetric, by appling definition using python built-in T module
    
    Remember that simmetry is defined only for squared matrices.
    """
    assert (temporal_network == temporal_network.T).all(), "Error: network is not symmetric"

def check_adiacency(network):
    """
    Verifies that the provided network only has 0 or 1 values
    
    This is accomplished by using numpy function "extract", which returns an array containing all values with a provided property.
    This function is applied twice: once to get an array of non-0 values, and once for non-1 values.
    If this array has 0 length, matrice is correctly and adiacency.
    """
    #network MUST be np.array
    #each value of network must be 0 or 1
    assert len(np.extract(np.extract(network!=0, network) !=1, network)) ==0, "Error: there is at least one non-0 or 1 value in adiacency"

def structural_suite(network,nodes_number,duration,symmetry = True):
    """
    Checks some ubiquitary properties of temporal networks in this work.
    
    If they are verified, it doesn't mean that networks is produced correctly, but just that it's structure is how it was supposed to.
    So, this just a preliminary test.
    
    These are checked properties:
        * 3dimensionality
        * proper duration and number of nodes
        * adiacencies actually have only 0 or 1 networks (this is not the case of having multiple links between 2 same nodes at the same time)
        * squareness and null-diagonal for adiacencies (using 2 assertion functions defined in "Evolutions" module)
        * symmetry (if required)
    
    Parameters:
    -----------
    network: np.array
        T*N*N sequence of matrices
    nodes_number: int
        supposed number of nodes
    duration: int
        supposed duration (i.e. number of adiacencies)
    """
    Evolutions.assert_ndarray(network,3) #3d np.array
    assert len(network) == duration, "Error: output network has not the right duration"
    assert len(network[0]) == nodes_number, "Error: number of nodes doesn't match with one of adiecency length"
    check_adiacency(network) #only 0 or 1 values for adiacencies
    [Evolutions.assert_square(network[t]) for t in range(duration)]  #each adiacency is a square of proper size
    [Evolutions.assert_nulldiagonal(network[t]) for t in range(duration)] #each adiacency has null-diagonal
    if symmetry:
        [check_symmetry(network[t]) for t in range(duration)]

#%% DAR TESTS
def test_DARgeneration1():
    #3 nodes, duration 10, a,xi->zeros
    a_input = np.zeros((3,3))
    xi_input = np.zeros((3,3))
    tested1 = Evolutions.network_generation_dar(a_input,xi_input, T = 10,directed = False)
    structural_suite(tested1,nodes_number=3,duration=10, symmetry = True) #structural test
    for t in range(1,10):
        assert (tested1[t] == np.zeros((3,3))).all()

def test_DARgeneration2():  
    #16 nodes, duration 20, a->zeros and xi->ones
    a_input = np.zeros((16,16))
    xi_input = np.ones((16,16))
    tested2 = Evolutions.network_generation_dar(a_input,xi_input, T = 20,directed = False)
    structural_suite(tested2,nodes_number=16,duration=20, symmetry = True) #structural test
    for t in range(1,20):
        assert (tested2[t] == (np.ones((16,16))-np.identity(16))).all()

def test_DARgeneration3():    
    #30 nodes, duration 5, a zeros and xi ones
    a_input = np.ones((30,30))
    xi_input = np.zeros((30,30))
    tested3 = Evolutions.network_generation_dar(a_input,xi_input, T = 5,directed = False)
    structural_suite(tested3,nodes_number=30,duration=5, symmetry = True) #structural test
    for t in range(1,5):
        assert (tested3[t] == tested3[0]).all()

#%% TGRG TESTS
def test_TGRGgeneration1():
    #20 nodes, duration 10, low phi0
    input_0 = 100*np.ones(20)
    input_1 = np.zeros(20)
    input_e = np.zeros(20)
    tested1 = Evolutions.network_generation_tgrg(input_0,input_1,input_e, T= 10, directed = False)[0] #tgrg return temporalnet and thetas
    structural_suite(tested1,nodes_number=20,duration=10, symmetry = True) #structural test
    for t in range(10):
        assert (tested1[t] == (np.ones(20)-np.identity(20))).all()
        
def test_TGRGgeneration2():
    #10 nodes, duration 15, low phi0
    input_0 = -100*np.ones(10)
    input_1 = np.zeros(10)
    input_e = np.zeros(10)
    tested2 = Evolutions.network_generation_tgrg(input_0,input_1,input_e, T= 15, directed = False)[0] #tgrg return temporalnet and thetas
    structural_suite(tested2,nodes_number=10,duration=15, symmetry = True) #structural test
    for t in range(15):
        assert (tested2[t] == np.zeros(10)).all()

#TODO : i test per la Katz possono essere fatti ritornando al caso statico
#TODO: test poisson probability using #I = quad(poisson_probability,0,np.inf) #verify the correct normalization