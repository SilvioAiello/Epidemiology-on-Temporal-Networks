"""
Functions in this script provide:
    1) assertions functions to use at the beginning of those functions, to be sure their inputs are correct.
    2) tests to verify that functions in "Evolutions.py" and "Propagation_SI.py" run correctly.
    
Functions in this script work in Pyhon3, may require numpy (v1.16) and function "quad" from scipy.integrate (scipy v1.3).

From this module you can extract the following assertion-functions: 
    adiacency, ndarray, nulldiagonal, probability, square, symmetry; natural.

#TODO: AS TEACHER SAID, EXPLAIN HOW TO GET THE SAME RESULTS I GOT

For further understandings on how this script operates, check file "docs/tests.md"      
"""
import numpy as np
from scipy.integrate import quad 

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

#%% DAR GENERATION 
def test_Evolutions_DARgeneration1():
    #3 nodes, duration 10, a,xi->zeros; check "docs/howto.md" for further clarifications
    a_input = np.zeros((3,3))
    xi_input = np.zeros((3,3))
    tested1 = Evolutions.network_generation_dar(a_input,xi_input, T = 10,directed = False)
    structural_suite(tested1,nodes_number=3,duration=10, symmetry = True) #structural test
    for t in range(1,10):
        assert (tested1[t] == np.zeros((3,3))).all()

def test_Evolutions_DARgeneration2():  
    #16 nodes, duration 20, a->zeros and xi->ones; check "docs/howto.md" for further clarifications
    a_input = np.zeros((16,16))
    xi_input = np.ones((16,16))
    tested2 = Evolutions.network_generation_dar(a_input,xi_input, T = 20,directed = False)
    structural_suite(tested2,nodes_number=16,duration=20, symmetry = True) #structural test
    for t in range(1,20):
        assert (tested2[t] == (np.ones((16,16))-np.identity(16))).all()

def test_Evolutions_DARgeneration3():    
    #30 nodes, duration 5, a zeros and xi ones; check "docs/howto.md" for further clarifications
    a_input = np.ones((30,30))
    xi_input = np.zeros((30,30))
    tested3 = Evolutions.network_generation_dar(a_input,xi_input, T = 5,directed = False)
    structural_suite(tested3,nodes_number=30,duration=5, symmetry = True) #structural test
    for t in range(1,5):
        assert (tested3[t] == tested3[0]).all()

#%% TGRG GENERATION 
def test_Evolutions_TGRGgeneration1():
    #20 nodes, duration 10, low phi0; check "docs/howto.md" for further clarifications
    input_0 = 100*np.ones(20)
    input_1 = np.zeros(20)
    input_e = np.zeros(20)
    tested1 = Evolutions.network_generation_tgrg(input_0,input_1,input_e, T= 10, directed = False)[0] #tgrg return temporalnet and thetas
    structural_suite(tested1,nodes_number=20,duration=10, symmetry = True) #structural test
    for t in range(10):
        assert (tested1[t] == (np.ones(20)-np.identity(20))).all()
        
def test_Evolutions_TGRGgeneration2():
    #10 nodes, duration 15, low phi0; check "docs/howto.md" for further clarifications
    input_0 = -100*np.ones(10)
    input_1 = np.zeros(10)
    input_e = np.zeros(10)
    tested2 = Evolutions.network_generation_tgrg(input_0,input_1,input_e, T= 15, directed = False)[0] #tgrg return temporalnet and thetas
    structural_suite(tested2,nodes_number=10,duration=15, symmetry = True) #structural test
    for t in range(15):
        assert (tested2[t] == np.zeros(10)).all()

#%% DEGREE AND COMMUNICABILITY 
def test_Evolutions_degree_node():
    for i in range(100):
        assert Evolutions.degree_node(np.zeros((100,100)),i) == 0
        assert Evolutions.degree_node(np.ones((100,100))-np.identity(100),i) == 99
    
    test_network = np.zeros((2,3,3))
    #T = 0, out and in degrees, that show to differ in nodes 0 and 1
    test_network[0,0,1] = 1
    assert Evolutions.degree_node(test_network[0],0,out=True) == 1
    assert Evolutions.degree_node(test_network[0],1,out=True) == 0
    assert Evolutions.degree_node(test_network[0],2,out=True) == 0
    assert Evolutions.degree_node(test_network[0],0,out=False) == 0
    assert Evolutions.degree_node(test_network[0],1,out=False) == 1
    assert Evolutions.degree_node(test_network[0],2,out=False) == 0
    #T=1, out and in degree
    test_network[1,0,1] = 1
    test_network[1,1,0] = 1
    test_network[1,2,0] = 1
    assert Evolutions.degree_node(test_network[1],0,out=True) == 1
    assert Evolutions.degree_node(test_network[1],1,out=True) == 1
    assert Evolutions.degree_node(test_network[1],2,out=True) == 1
    assert Evolutions.degree_node(test_network[1],0,out=False) == 2
    assert Evolutions.degree_node(test_network[1],1,out=False) == 1
    assert Evolutions.degree_node(test_network[1],2,out=False) == 0    

def test_Evolutions_degree_mean():
    assert Evolutions.degree_mean(np.zeros((100,100))) == [0]
    assert Evolutions.degree_mean(np.ones((100,100))-np.identity(100)) == [99]
    
    test_network = np.zeros((2,3,3))
    test_network[0,0,1] = 1
    test_network[1,1,0] = 1
    assert Evolutions.degree_mean(test_network, out=True) == [1/3, 1/3]
    assert Evolutions.degree_mean(test_network, out=False) == [1/3,1/3]
    #more links
    test_network[0,1,2] = 1
    assert Evolutions.degree_mean(test_network, out=True)[0] == [2/3]
    assert Evolutions.degree_mean(test_network, out=False)[0] == [2/3]
    test_network[1,0,1] = 1
    test_network[1,2,0] = 1
    assert Evolutions.degree_mean(test_network, out=True)[1] == [1]
    assert Evolutions.degree_mean(test_network, out=False)[1] == [1]

def test_Evolutions_communicability():
    #3 nodes and 3 time step, all nodes linked without autoloops
    #input building:
    temporal_test = np.ones((3,3,3)) - np.identity(3) #input for tested function
    adj = temporal_test[0] #input for manual building
    
    Q = np.identity(3)
    inv = np.linalg.inv(np.identity(3)-0.125*adj)
    for t in range(3):
        Q = np.matmul(Q,inv)/np.linalg.norm(np.matmul(Q,inv))
    assert (Q - Evolutions.communicability(temporal_test)[1] <= 1e-10).all()

#%% Propagation easing functions
def test_PropagationSI_neighbourhood():
    N = 10
    network_test = np.ones((N,N))-np.identity(N)
    for i in range(N):
        assert Propagation_SI.neighbourhood(np.zeros((N,N)),i) == set() #null network
        assert len(Propagation_SI.neighbourhood(network_test,i)) == N-1 #all-linked network
    
    for i in range(1,N):
        network_test[i,0] = 0
        assert len(Propagation_SI.neighbourhood(network_test,0)) == N-1-i
        assert len(Propagation_SI.neighbourhood(network_test,i)) == N-1 #confirm unaffection for reached node

def test_PropagationSI_contactlasting1():
    T = 50
    N = 20
    #tempnet initialization
    network_test = np.zeros((T,N,N))
    for t in range(T):
        for i in range(N):
            for j in range(N):
                if j != i: #null diagonal
                    network_test[t,i,j] = np.random.choice((0,1)) #random links
    
    #states initialization, only 0 infected, always
    states_sequence = dict()
    for t in range(T):
        states_sequence[t] = dict.fromkeys(range(N),0)
        states_sequence[t][0] = 1
    
    for t in range(T-2):
        for i in range(1,N):
            if network_test[t,0,i] == 1:
                contact_now = Propagation_SI.contact_lasting(network_test,states_sequence, t,0,i)
                assert contact_now <= t+1, "Contact lasts %i, which is bigger than %i"%(contact_now,t) #t+1 since python starts count from 0
                assert  contact_now >= 1
                if network_test[t+1,0,i] == 1:
                    assert contact_now < Propagation_SI.contact_lasting(network_test,states_sequence, t+1,0,i)
                
def test_PropagationSI_contactlasting2():
    T = 75
    N = 10
    #tempnet initialization
    network_test = np.zeros((T,N,N))
    for t in range(T):
        network_test[t,0,1] = 1 #0 and 1 always linked

    #states initialization, only 0 infected, always
    states_sequence = dict()
    for t in range(T):
        states_sequence[t] = dict.fromkeys(range(N),0)
        states_sequence[t][0] = 1
        
    #TEST2
    for t in range(T-1):
        assert Propagation_SI.contact_lasting(network_test,states_sequence, t,0,1) == t+1
    
    #TEST3
    network_test[10,0,1] =0 #break contact only at a certain instant
    assert Propagation_SI.contact_lasting(network_test,states_sequence, 74,0,1) == 64
    
    #TEST4
    network_test[10,0,1] =0 #turn back
    states_sequence[10][1] = 1 #infected only at a certain instant
    assert Propagation_SI.contact_lasting(network_test,states_sequence, 74,0,1) == 64

def test_PropagationSI_poissonprobability():
    beta = np.random.random() #for any beta it must work
    integral = quad(Propagation_SI.poisson_probability, 0, np.inf, args = beta)
    assert integral[0]-1 <= integral[1]
    
    
#%% PROPAGATION_SI
def test_PropagationSI1():
    #Inputs
    T = 25
    N = 10
    tempnet = np.ones((T,N,N))
    for t in range(T):
        tempnet[t] -= np.identity(N)
    
    #Tests
    #only index case stays infected
    beta = 1e-10
    probabilities = dict()
    for t in range(T):
        probabilities[t] = quad(Propagation_SI.poisson_probability,0,t, args = beta)[0]
    
    for index_case in range(N):
        final_state = dict.fromkeys(range(N),0)
        final_state[index_case] = 1
        assert Propagation_SI.propagation(tempnet,index_case,probabilities)[T-1] == final_state
    
    #all infected
    beta = 2/3
    probabilities = dict()
    for t in range(T):
        probabilities[t] = quad(Propagation_SI.poisson_probability,0,t, args = beta)[0]
    
    final_state = dict.fromkeys(range(N),1)
    for index_case in range(N):
        assert Propagation_SI.propagation(tempnet,index_case,probabilities)[T-1] == final_state

    #all infected but one node is isolated
    for t in range(T):
        tempnet[t,:,N-1] = np.zeros(N)
    final_state = dict.fromkeys(range(N),1)
    final_state[N-1] = 0
    for index_case in range(N-1):
        assert Propagation_SI.propagation(tempnet,index_case,probabilities)[T-1] == final_state

def test_PropagationSI2():
    K = 3
    T = 10
    N = 20
    temp = Evolutions.network_generation_dar(0.6*(np.ones((N,N))-np.identity(N)),0.6*(np.ones((N,N))-np.identity(N)),P=1,T=T,directed=False) 
    
    beta = 0.4
    probabilities = dict() #dict building
    for t in range(T):
        probabilities[t] = quad(Propagation_SI.poisson_probability,0,t, args = beta)[0] #quad produces several outputs, integral is the first
        
    label = [] 
    for index_case in range(N):
        label.append([]) #create the i-th entry
        for iteration in range(K):
            label[index_case].append(Propagation_SI.propagation(temp, index_case, probabilities))
    # OUTPUT TESTS
    assert label[index_case][iteration][0][index_case] == 1, "An index case appears to be uninfected"
    assert sum(label[index_case][iteration][0].values()) == 1, "There should be only 1 infect at the beginning"
    assert [[label[index_case][iteration][0] == label[index_case][iteration-1][0] for iteration in range(1,K)] for index_case in range(N)], "Initial condition is not equal for all iterations" 
#asserire che sum deve essere sempre <= N? (usa .values())