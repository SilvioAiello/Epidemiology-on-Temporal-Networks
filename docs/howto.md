In this file it is explained how functions of each script work, what parameters require and why, and what outputs they produce.
If you're interested in just one function or script, use table of content to directly find it.

# Table of contents (scripts and their functions)
* [Main.py](#main)
* [Evolutions.py](#evolutions)
  * [Network generation functions](#network_generation_dar-and-network_generation_tgrg)
  * [Degree functions](#degree-functions)
  * [Communicability and rankings](#communicability-and-rankings)
  * [Network save](#network_save)
* [Propagation SI.py](#propagation_si)
  * [Network load](#network_load)
  * [Easing simulation functions](#easing-simulation)
  * [Propagation](#propagation)
  * [Epidemic scores functions](#epidemic-scores)
* [Propagation LTM](#propagation_ltm)
* [Tests.py](#test_suite)
  * [Assertions](#assertion-functions)
  * [Evolution tests](#evolution-functions-tests)
  * [Measures tests](#measures-functions-tests)

# Main
        What destiny for network plot?
This script manages all the others, allowing user to actually create and save networks and epidemics, and to perform analysis.
Modifing it, you can set the input parameters you prefer, and perform simulations according to them:
* N, T are the number of nodes and duration of temporal network
* K is the number of epidemic simulations:
* beta is the infection rate, that defines epidemic virulence; 
* alpha and xi are defined, for DAR networks, in explanation, as well as, phi0,phi1,epsilon for TGRG.
Può servirti usare sta roba: to generate more tempnet is up to you; but iterating propagation for nodes is "compulsory", since scores MUST be averages;
* **Poisson_Probability** function, that should not be modified, reproduces the Poisson PDF (Probability Density Function), whose average depends on beta, and whose definite integral is the probability of having infection at a certain time. Since we're dealing with discrete finite time, there's a finite number (T) of integrals/probabilities that can be computed, since a contact may last only from 0 to T. So, same integrals will be computed repeatedly, and it makes sense to compute them once for all, transfering results in a dictionary, whose keys are the durations and name is just **probabilities**. This may save much computational time (also check [Explanation](/explanation.md)).
* Epidemic scores are stored in label_dar/label_fitn data structures, which share this hierarchy: label_\[index_case]\[iteration]\[time_step]\[node]; a label_ is a list of N dictionaries (one for each node set as index case), containing a list of the results K iterations (one for each same epidemic simulation), which in turn are dictionary of dictionaries, describing evolution of nodes states (as you will see in propagation function)
* Centrality and virulence scores are saved in namesake lists.
* Everything is saved by...

              Quale destino per i salvataggi?

# Evolutions
It allows to generate and analyze temporal networks, making use of functions from numpy and pickle libraries, so be sure to have them installed. As already mentioned, linking-state between two nodes is described by zeros (no link at that time) and ones (link at that time). 

### network_generation_dar and network_generation_tgrg
These functions generate DAR(P) and TGRG networks, according to definion you can find in [Explanation](/explanation.md).
Temporal networks are described through the temporal sequence of their adiacency matrices, so this is what you get.
Each update is a stochastic process, so it's difficult to get twice the same results, and you should call these functions multiple time to get some statistically relevant results.
Let's have a look at what parameters these function require:

    network_generation_dar(alpha,xi,P=1, T=100, directed = False)

    network_generation_tgrg(phi0,phi1,sigma,T=100, directed = False)

The compulsory parameters user must provide depend on networks' nature: alpha and xi matrices for DAR(P), phi0,phi1,sigma vectors for TGRG. 
Some parameters are set by default, but user can change them: time steps of the simulation (default is 100), whether to generate a directed or undirected graph, and, only for DAR(P), the number P of past step to consider.
The output is np.array with shape T\*N\*N, i.e. an ordered sequence of T square(N) matrices, each of wich is the adiacency at that time step.

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
* In TGRG, firstly, all fitnesses evolutions are extracted (to see how this works, check "expalanation.md"), for each node, after initializing them to be equal to the respective "phi0". Then, T\*N\*N tensor is generated and, making use of nested for-cicles, and distinguishing triangles by setting i> or <j, it is computed the probability of link for each time step and couple of nodes, and it is put in "choice".

DAR gives only back the network, TGRG also the matrix of fitnesses.

### degree functions
These functions extract various forms of degree from a network, if it's of the type of those produced by previous functions (so T\*N\*N np.array): degree for a single node, mean degree at certain time, sequence of degrees. They just require the adiacency (or their sequence, with the initial and final time step), the single node if necessary, and wheter is required out-going degree or in-going degree (defining it is beyond the purpose of this text).

Before operating, they make their arguments underpass a couple of assertions: network has to have the already mentioned form, and null diagonal; if node index is required, it has to be a natural number that doesn't exceed network length; if temporal steps of analysis are required, they have to be integers >= 0 and not exceed network evolution duration.

Note that each function doesn't rely on the other, so they are "pure".

* **degree_node** return out or in-going degree of the provided node, just by summing (respectively) column of line of node in provided network, and returns it as a natural number;
* **degree_mean_t** does the same but for each node of a certain adiacency, returning the mean degree at that time (real number);
* **degree_mean_sequence** does the same but for each adiacency from an initial to a final time step, return a list of degree evolution.

### communicability and rankings
These functions allow to determine the Communicability matrix for a provided temporal network, perform Broadcast and Receive Centrality and rank nodes accordingly.
The only assertions they perform are about shape (sequence of squadre matrices) and diagonal nullity.

* **communicability** extracts the number of nodes and evolution duration by the shape of the temporal network, then computes for each adiacency the maximum spectral radius (each of which makes use of np.linalg.eigvals) and selects the maximum from these. A quarter of the highest one is the factor applied to each adiacency in Communicability updating (as in Chen paper), but in further developments it could be different, or it could be given the possibility to choose among different ones; communicability is initialized as a normalized (np.linalg.norm) identity matrix, and multiplied, at each time step, by the inverse (np.linalg.inv) of the difference between identity (np.identity) and the the weighted adiacency; all matrices are normalized step by step, and matrices multiplications are performed by np.matmul;
* **broadcast_ranking/receive_ranking** require such a communicability matrix, and returns a list with nodes BCs/RCs by summing over lines/columns; then it generates a list, lines_sum, containing nodes indeces ordered by their scores, the firsts being those with highest centralities, using np.argsort (that return indeces but from the worst to the best) and np.flip (that invertes lists); note that "lines_sum" doesn't follow such a order;

### network_save
This function saves, in an automatized way, a temporal network, using pickle (so it binarizes) and following the foldering/naming rule expressed in Documentation (so, distinguishing networks according to their parameters and realizations), allowing user to add a particoular identification name to the file.

# Propagation_SI
It allows to generate and analyze epidemics on temporal networks, making use of functions from numpy and pickle libraries, so be sure to have them installed. As already mentioned, epidemic state for a node is described by zeros (susceptible) and ones (infected). 
Each update is a stochastic process, so it's difficult to get twice the same results, and you should call these functions multiple time to get some statistically relevant results.

### network_load
       This function is still undergoing tests

### easing simulation
These functions allow* propagation function* to perform well and in an efficient way. They solve simple tasks analyzing adjacencies and node epidemic states, which is very useful since the idea of SI model is that each I-node can infect its *Susceptible neighbours*, at a certain time, according to a probability that depends on rate beta and *contact lasting*.
* **neighbourhood**: generates a set of neighbours of a node at a certain time, by checking the provided adiacency; set is create through a comprehension;
* **onlyzeros**: extracts the subset of susceptible nodes at that time, from a general set, list, etc, by checking the provided states dict; set is create through a comprehension;
* **contact_lasting**: computes the number of time steps an S and a I node were connected, starting from a provided time. This is accomplished by checking backwards, time by time, the existence of the link and the state of the I node, increasing the value of a counter variable as long these conditions are satisfied; so, the loop interrupts if it has gone so back that the link wasn't  yet present, or the I-node wasn't yet infected. Clearly, the higher the output is, the higher the probability of having an infection (the value is checked in the already mentioned probabilities dictionary).

### propagation
This function takes a temporal network, an index case and the probabilities dictionary, and performs an iteration of a whole epidemic propagation in SI mode, from the first to last snapshot of the tempnet; the result is a dictionary of dictionaries, describing nodes' states for each time step.

First operations are, as usual, definitions: T and N are extracted from the tempnet, and then an internal function, *set_infected*, is created with the purpose of setting input node states equals to 1 from a given time to the end (once infected, always infected); it will be used on every new infected; then, output is initialized by setting every node susceptible at every time, with the expection of the index case who is always infected (set_infected since t=0). Finally, infected and susceptibles nodes sets are initialized: they will be updated every time a new infection occurs. 
Now, following Chen approach, epidemic is updated time by time by processing each susceptible node (at the previous time step), finding its infective neighbourhood (at the previous time step), and performing the "infection extraction" for each of them: if infection occurs, infective state is set "1" at present time step, sets are updated and the program jumps directly to next susceptible. Infection extraction compares the integral of Poisson distribution with a random uniform (from 0 to 1) extraction, performed by *np.random.uniform* function. This is iterated until the end of network evolution.

### epidemic scores
They are computed making use of two functions:
* **infected_counter**: takes the states dictionary at a certain time, and counts of many nodes are infected;
* **time_score**: uses infected_counter for each time step of the propagation, and saves the first step when the fraction of infected nodes was higher than a given value. If this never happens, it returns the total duration of the propagation.

# Propagation_LTM
This section will be deepened in further developments.

# Test_suite
This script contains a suite of tests veryfing that every function in other scripts works properly, and some assertion functions, useful to these latter functions to check the inputs they are provided of. Most complex tests obviously deal with temporal networks and epidemics.

### assertion functions
This section is quite self-explaining: its function impose inputs to be an np.array of a certain dimension, or to have a square shape, or to be probability matrices (so accepting only values from 0 to 1), or to be a natural number (integer and >=0), or to be a matrix with 0-diagonal. Most of these functions are imported by Evolutions and Propagations scripts.

### evolution functions tests
This section contains both structural tests (checking the output temporal network having the right mathematical properties) and tests of the actual evolution (since it's a stochastic process, one cannot make assertions about the exact outcome values, aside from some limit-cases; these ones are found and tested).
* **structural_suite**: since some structural tests will be repetead multiple times in both evolutions, they are collected in this function, which performs some of the previous asserts and, in the end, checks that the output temporal network has the right parameters (number of nodes and duration) and mathematical properties, like being a succession of adiacencies, which in turn are square matrices with null diagonal and, if case, symmetrical (Explanation.md if you need to better understand these lines). This function is recalled in any DAR and TGRG test. If it is verified, it doesn't mean that networks are produced correctly, but just that their structure is how it was supposed to. So, this is just a preliminary test.
* **Evolution/dynamic tests**: as said, these tests check, besides the structural suite, some limit-case inputs, the only ones one can be sure of the outputs; since multiple combinations are possible, multiple tests of the same kind are performed, changing time by time some parameters like number of nodes and duration, just for sake of completenes.
  * DAR(1): we can be sure that, if matrix alpha = all zeros, all following states are determined by performing a random extraction (ruled by xi), while if alpha = all ones: all following states are equals to the first (total persistence); moreover, if matrix xi = all zeros, there's no way of getting state "1" for any link, if an extraction occurs, and vice versa for xi = all ones. So, these limit cases are tested: **alpha and xi = all zeros**, expecting a sequence of ONLY ZEROS adiacencies; **alpha = all zeros, xi = all ones** expecting a sequence of ONLY ONES (but null diagonal, of couse) adiacencies; **alpha = all ones, caringless of xi** expecting a sequence of adiacencies equal to the initial one.
  * TGRG: if "sigma" vector is set to 0, one removes randomness in fitnesses evolution (but not in link generation in following times); anyway, giving to each entry of phi0 and phi1 (or just phi0) very high (in module) values, one can jump into certainty domain. So, these limits will be checked, taking for guaranted that sigmas are 0: **very high values (100) for all phi0, and just 0 phi1**, expecting a sequence of ONLY ONES (but null diagonal) adiacencies; **very low values (-100) for all phi0, and just 0 phi1**, expecting a sequence of ONLY ZEROS adiacencies.

### measures functions tests
NEXT STEP
