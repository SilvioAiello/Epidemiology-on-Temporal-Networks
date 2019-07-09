In this file it is explained how functions of each module work.

# Evolutions module
It's made use of function from numpy and pickle libraries.

## network_generation_dar and network_generation_tgrg
These functions generate a DAR(P) and a TGRG network, according to definion you can find in Explanation (MAGARI POI LINKA).
The compulsory parameters they require are those defining the respective networks: alpha and xi matrices for the first, phi0,phi1,sigma for the latter. Default parameters the user can change are the number of past step to consider (only for DAR model, whose default is 1), time steps of the simulation (default is 100), and whether generate a directed or undirected graph.
The result is np.array with shape T\*N\*N, i.e. an ordered sequence of T square(N) matrices, each of wich is the adiacency at that time step.

Before performing the simulation, they check a couple of assertions: 
* matrices/vectors provided must be np.arrays with squadre or vectorial shape; 
* all these structures must share the same length;
* alpha and xi must be probabilities, i.e., all their values should be included within 0 and 1; 
* if all these are verified, function extract the information about the number of nodes (contained in matrices/vector "length") and verify it is < 1000 (it's improbable to work with more than 1000 nodes, so  if it's the case, it may be due to human error); 
* the provided T and P (for DAR(P)) must be natural numbers, and also T must be < 1000;
* P needs to be < T "for definition"

Then, evolutions begin; to better understand how they work, keep in mind that:
* in order to guarantee the possibility to generate a directed or undirected graph, they both first generate an upper triangoular matrix, and for the lower, depending on what the user chooses, they just copy the result or keep on extracting values;
* they both make use of the function "np.random.choice", that performs and extractions between (here) 2 numbers according to a (given) probability for each number; probability conservation (i.e. sum of them giving 1) is guaranteed by initial assertions (for DAR) or by mathematical properties of distribuion (for TGRG).
This is what the function actually do:
* In DAR, firstly, the initial adiacency is extracted (with 50% probability of having or not a link for each initial couple) and the whole T\*N\*N tensor is generated, copying in the the initial state. So, everything is set to continue evolution. Both initial state and evolutions are handled with two nested comprehensions and for cicles (that distinguish triangles by setting i> or <j), the inner one containing "choice". For the initial states, choice is between 0 and 1 with probability 0.5 for each, while for following moments, they are took from the previous P steps with probability alpha, or extracted by a tossing coin (whose "inner" probability is xi), just as DAR(P) says.
* In TGRG, firstly, all fitnesses evolutions are extracted, for each node, initializing them according to the same random-normal that disciplines their evolutions. Then, T\*N\*N tensor is generated and, making use of nested for-cicles, and distinguishing triangles by setting i> or <j, it is computed the probability of link for each time step and couple of nodes, and it is put in "choice".

DAR gives only back the network, TGRG also the matrix of fitnesses.

## degree functions
These functions extract various forms of degree from a network, if it's of the type of those produced by previous functions (so T\*N\*N np.array): degree for a single node, mean degree at certain time, sequence of degrees. They just require the adiacency (or their sequence, with the initial and final time step), the single node if necessary, and wheter is required out-going degree or in-going degree (defining it is beyond the purpose of this text).

Before operating, they make their arguments underpass a couple of assertions: network has to have the already mentioned form, and null diagonal; if node index is required, it has to be a natural number that doesn't exceed network length; if temporal steps of analysis are required, they have to be integers >= 0 and not exceed network evolution duration.

Note that each function doesn't rely on the other, so they are "pure".

* degree_node return out or in-going degree of the provided node, just by summing (respectively) column of line of node in provided network, and returns it as a natural number;
* degree_mean_t does the same but for each node of a certain adiacency, returning the mean degree at that time (real number);
* degree_mean_sequence does the same but for each adiacency from an initial to a final time step, return a list of degree evolution.

## communicability and rankings
These functions allow to determine the Communicability matrix for a provided temporal network, perform Broadcast and Receive Centrality and rank nodes accordingly.
The only assertions they perform are about shape (sequence of squadre matrices) and diagonal nullity.

* communicability extracts the number of nodes and evolution duration by the shape of the temporal network, then computes for each adiacency the maximum spectral radius (each of which makes use of np.linalg.eigvals) and selects the maximum from these. A quarter of the highest one is the factor applied to each adiacency in Communicability updating; communicability is initialized as a normalized (np.linalg.norm) identity matrix, and multiplied, at each time step, by the inverse (np.linalg.inv) of the difference between identity (np.identity) and the the weighted adiacency; all matrices are normalized step by step, and matrices multiplications are performed by np.matmul;
* broadcast_ranking requires such a communicability matrix, and returns a list with nodes BCs by summing over lines; then it generates a list, lines_sum, containing nodes indeces ordered by their scores, the firsts being those with highest centralities, using np.argsort (that return indeces but from the worst to the best) and np.flip (that invertes lists); note that "lines_sum" doesn't follow such a order;
* broadcast_ranking requires such a communicability matrix, and returns a list with nodes RCs by summing over columns; then it generates a list, lines_sum, containing nodes indeces ordered by their scores, the firsts being those with highest centralities, using np.argsort (that return indeces but from the worst to the best) and np.flip (that invertes lists); note that "lines_sum" doesn't follow such a order.

## network_save
This function saves, in an automatized way, a temporal network, using pickle (so it binarizes) and following the foldering/naming rule expressed in Documentation (so, distinguishing networks according to their parameters and realizations), allowing user to add a particoular identification name to the file.


# Propagation_SI
Working progress

# Propagation_LTM
Further developments
