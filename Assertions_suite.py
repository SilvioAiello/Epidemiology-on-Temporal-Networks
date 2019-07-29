"""
Functions in this script provide "truth-statements" to assertions at the beginning of various functions among the project.
    
They work in Pyhon3 and may require numpy (v1.16) .

You can extract the following functions: 
    adiacency, ndarray, nulldiagonal, probability, square, symmetry; natural.

For further understandings on how this script operates, check file "docs/tests.md"      
"""
import numpy as np

# Matrices
def check_is_adiacency(network):
    """
    Verifies that the provided network only has 0 or 1 values, so it's an adiacency with max 1 link per couple of nodes.
    
    np.extract returns an array containing input values with a provided property.
    If searching for non-(0,1) values returns an empty list, matrix is correctly and adiacency of desidered type.
    
    Parameters: np.array, of any dimension
    """
    #extact returns an array containing all values with a provided property
    if len(np.extract(np.extract(network!=0, network) !=1, network)) ==0:
        return True
    else:
        return False #there is at least one non-0 or non-1 value in adiacency

def check_is_ndarray(network,dimensions):
    """
    Check if a data structure is a np.array, with proper dimension expressed by an integer
    
    Parameters: np.array and an integer for dimension
    """
    if isinstance(network,np.ndarray) and len(network.shape) == dimensions:
        return True
    else:
        return False #matrix is not a np array or has not the proper dimension

def check_is_nulldiagonal(matrix):
    """
    Check that matrix has null diagonal
    """
    if sum(np.diag(matrix)) == 0:
        return True
    else:
        return False

def check_is_probability(network):
    """
    Verify input matrix is a probability matrix: it's value must be whitin 0 and 1
    """
    #check that all elements are probabilities
    if (network<= 1).all() and (network>=0).all():
        return True
    else:
        return False

def check_is_square(network):
    """
    Checks if a data structure is square-like, i.e. each of its sides has the same length (as the first one).
    Checking is accomplished by a target variable that takes values "False" if just one time provides something unexpected.
    """
    label = True
    for i in range(len(network.shape)):
        label == label and network.shape[i] == network.shape[0]
    return label

def check_is_symmetry(temporal_network):
    """
    Verifies if a matrix is simmetric, by appling definition using python built-in T module
    
    Remember that simmetry is defined only for squared matrices.
    """
    if (temporal_network == temporal_network.T).all():
        return True 
    else:
        return False #network is not symmetric

# Parameters
def check_is_natural(number):
    """
    Check that input is a positive integer
    """
    if number>0 and isinstance(number, int):
        return True
    else:
        return False