"""
Functions in this script provide:
    1) tests to verify that functions in "Evolutions.py" and "Propagation_SI.py" run correctly.
    2) assertions functions to use at the beginning of those functions, to be sure their inputs are correct.

Functions work in Pyhon3, and may require the following libraries (so, check if they are installed):
    * numpy, used for its data structures and anaylisis, and to get random functions   

From this module you can extract the following functions:
ASSERTIONS: adiacency, ndarray, nulldiagonal, probability, square, symmetry; natural.
STRUCTURAL:
    Structural suite
    Tests for DAR generation
    Tests for TGRG generation
    #TODO: TESTS FOR MEASURES
#TODO: AS TEACHER SAID, EXPLAIN HOW TO GET THE SAME RESULTS I GOT

For further understandings on how this script operates, check file "howto.md"      
"""
import numpy as np

import Evolutions
import Propagation_SI
#%% ASSERTIONS
# Matrices
def assert_adiacency(network):
    """
    Verifies that the provided network only has 0 or 1 values, so it's an adiacency with max 1 link per couple of nodes.
    
    np.extract returns an array containing input values with a provided property.
    If searching for non-(0,1) values returns an empty list, matrix is correctly and adiacency of desidered type.
    
    Parameters: np.array, of any dimension
    """
    #extact returns an array containing all values with a provided property
    assert len(np.extract(np.extract(network!=0, network) !=1, network)) ==0, "Error: there is at least one non-0 or non-1 value in adiacency"

def assert_ndarray(network,dimensions):
    """
    Check if a data structure is a np.array, with proper dimension expressed by an integer
    
    Parameters: np.array and an integer for dimension
    """
    assert isinstance(network,np.ndarray), "Error: matrix must be a numpy array"
    assert len(network.shape) == dimensions, "Error: matrix has not the proper dimension"

def assert_nulldiagonal(matrix):
    """
    Check that matrix has null diagonal
    """
    assert sum(np.diag(matrix)) == 0, "Error: network has not-0 diagonal"

def assert_probability(network):
    """
    Verify input matrix is a probability matrix: it's value must be whitin 0 and 1
    """
    #check that all elements are probabilities
    assert (network<= 1).all(),"Error: at least one element in a probability matrix is >1, so it's not a probability"
    assert (network>=0).all(),"Error: at least one element in probability matrix is <0, so it's not a probability"

def assert_square(network):
    """
    Checks if a data structure is square-like, i.e. each of its sides has the same length (as the first one)
    """
    #check if matrix is a square by comparing each side length with the first side
    for i in range(len(network.shape)):
        assert network.shape[i] == network.shape[0], "Error: matrix is not a square"

def assert_symmetry(temporal_network):
    """
    Verifies if a matrix is simmetric, by appling definition using python built-in T module
    
    Remember that simmetry is defined only for squared matrices.
    """
    assert (temporal_network == temporal_network.T).all(), "Error: network is not symmetric"

# Parameters
def assert_natural(number):
    """
    Check that input is a positive integer
    """
    assert isinstance(number, int), "Error: %f is not an integer, but it should be" %number
    assert number>0, "Error: %i is not positive, but it should be" %number

#%% STRUCTURAL PROPERTIES
def structural_suite(network,nodes_number,duration,symmetry = True):
    """
    Checks some ubiquitary properties of temporal networks in this work.
    
    Parameters:
    -----------
    network: np.array
        T*N*N sequence of matrices
    nodes_number: int
        supposed number of nodes
    duration: int
        supposed duration (i.e. number of adiacencies)
    symmetry: bool
    """
    assert_ndarray(network,3) #3d np.array
    assert len(network) == duration, "Error: output network has not the right duration"
    assert len(network[0]) == nodes_number, "Error: number of nodes doesn't match with one of adiecency length"
    assert_adiacency(network) #only 0 or 1 values for adiacencies
    [assert_square(network[t]) for t in range(duration)]  #each adiacency is a square of proper size
    [assert_nulldiagonal(network[t]) for t in range(duration)] #each adiacency has null-diagonal
    if symmetry:
        [assert_symmetry(network[t]) for t in range(duration)]

#%% DAR TESTS
def test_DARgeneration1():
    #3 nodes, duration 10, a,xi->zeros; check "docs/howto.md" for further clarifications
    a_input = np.zeros((3,3))
    xi_input = np.zeros((3,3))
    tested1 = Evolutions.network_generation_dar(a_input,xi_input, T = 10,directed = False)
    structural_suite(tested1,nodes_number=3,duration=10, symmetry = True) #structural test
    for t in range(1,10):
        assert (tested1[t] == np.zeros((3,3))).all()

def test_DARgeneration2():  
    #16 nodes, duration 20, a->zeros and xi->ones; check "docs/howto.md" for further clarifications
    a_input = np.zeros((16,16))
    xi_input = np.ones((16,16))
    tested2 = Evolutions.network_generation_dar(a_input,xi_input, T = 20,directed = False)
    structural_suite(tested2,nodes_number=16,duration=20, symmetry = True) #structural test
    for t in range(1,20):
        assert (tested2[t] == (np.ones((16,16))-np.identity(16))).all()

def test_DARgeneration3():    
    #30 nodes, duration 5, a zeros and xi ones; check "docs/howto.md" for further clarifications
    a_input = np.ones((30,30))
    xi_input = np.zeros((30,30))
    tested3 = Evolutions.network_generation_dar(a_input,xi_input, T = 5,directed = False)
    structural_suite(tested3,nodes_number=30,duration=5, symmetry = True) #structural test
    for t in range(1,5):
        assert (tested3[t] == tested3[0]).all()

#%% TGRG TESTS
def test_TGRGgeneration1():
    #20 nodes, duration 10, low phi0; check "docs/howto.md" for further clarifications
    input_0 = 100*np.ones(20)
    input_1 = np.zeros(20)
    input_e = np.zeros(20)
    tested1 = Evolutions.network_generation_tgrg(input_0,input_1,input_e, T= 10, directed = False)[0] #tgrg return temporalnet and thetas
    structural_suite(tested1,nodes_number=20,duration=10, symmetry = True) #structural test
    for t in range(10):
        assert (tested1[t] == (np.ones(20)-np.identity(20))).all()
        
def test_TGRGgeneration2():
    #10 nodes, duration 15, low phi0; check "docs/howto.md" for further clarifications
    input_0 = -100*np.ones(10)
    input_1 = np.zeros(10)
    input_e = np.zeros(10)
    tested2 = Evolutions.network_generation_tgrg(input_0,input_1,input_e, T= 15, directed = False)[0] #tgrg return temporalnet and thetas
    structural_suite(tested2,nodes_number=10,duration=15, symmetry = True) #structural test
    for t in range(15):
        assert (tested2[t] == np.zeros(10)).all()

#TODO : i test per la Katz possono essere fatti ritornando al caso statico
#TODO: test poisson probability using #I = quad(poisson_probability,0,np.inf) #verify the correct normalization
#TODO: UN MODO PER TESTARE contact_lasting E' VEDERE CHE AD OGNI ISTANTE SUCCESSIVO E' DIVERSO DAL PRECEDENTE, PER LA STESSA COPPIA
        #sicuramente il suo output non può essere maggiore del t di input