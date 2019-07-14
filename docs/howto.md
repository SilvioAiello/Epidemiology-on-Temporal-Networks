In this file it is explained how functions of each script work, what parameters require and why, and what outputs they produce.
If you're interested in just one function or script, use table of content to directly find it.

# Table of contents (scripts and their functions)
* [Main.py](#main-script)
* [Evolutions.py](#evolutions-script)
  * [Network generation functions](#network_generation_dar-and-network_generation_tgrg)
  * [Degree functions](#degree-functions)
  * [Communicability and rankings](#communicability-and-rankings)
  * [Network save](#network_save)
* [Propagation SI.py](#propagation_si-script)
  * [Network load](#network-load)
  * [Easing simulation functions](#easing-simulation)
  * [Propagation](#propagation)
  * [Epidemic scores functions](#epidemic-scores)
* [Propagation LTM](#propagation_ltm)
* [Tests.py](#tests)
  * [Assertions](#assertions)

# Main script
Può servirti usare sta roba:
K -> #to generate more tempnet is up to you; but iterating propagation for nodes is "compulsory": scores MUST be averages
label_dar = #list of N lists (one per index case), each of which is a list of K (one per iteration) sequences of T dictionaries
#(forse anche sti commenti lunghi li puoi mettere direttamente in documentazione)
#label dar 0: tutto riferito allo 0 nodo come index, label_dar[0,3]: 3 iterazione

# Evolutions script
It allows to generate and analyze temporal networks, making use of functions from numpy and pickle libraries, so be sure to have them installed. As already mentioned, linking-state between two nodes is described by zeros (no link at that time) and ones (link at that time). 

## network_generation_dar and network_generation_tgrg
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

* communicability extracts the number of nodes and evolution duration by the shape of the temporal network, then computes for each adiacency the maximum spectral radius (each of which makes use of np.linalg.eigvals) and selects the maximum from these. A quarter of the highest one is the factor applied to each adiacency in Communicability updating (as in Chen paper), but in further developments it could be different, or it could be given the possibility to choose among different ones; communicability is initialized as a normalized (np.linalg.norm) identity matrix, and multiplied, at each time step, by the inverse (np.linalg.inv) of the difference between identity (np.identity) and the the weighted adiacency; all matrices are normalized step by step, and matrices multiplications are performed by np.matmul;
* broadcast_ranking requires such a communicability matrix, and returns a list with nodes BCs by summing over lines; then it generates a list, lines_sum, containing nodes indeces ordered by their scores, the firsts being those with highest centralities, using np.argsort (that return indeces but from the worst to the best) and np.flip (that invertes lists); note that "lines_sum" doesn't follow such a order;
* broadcast_ranking requires such a communicability matrix, and returns a list with nodes RCs by summing over columns; then it generates a list, lines_sum, containing nodes indeces ordered by their scores, the firsts being those with highest centralities, using np.argsort (that return indeces but from the worst to the best) and np.flip (that invertes lists); note that "lines_sum" doesn't follow such a order.

## network_save
This function saves, in an automatized way, a temporal network, using pickle (so it binarizes) and following the foldering/naming rule expressed in Documentation (so, distinguishing networks according to their parameters and realizations), allowing user to add a particoular identification name to the file.

# Propagation_SI script
It allows to generate and analyze epidemics on temporal networks, making use of functions from numpy and pickle libraries, so be sure to have them installed. As already mentioned, epidemic state for a node is described by zeros (susceptible) and ones (infected). 
Each update is a stochastic process, so it's difficult to get twice the same results, and you should call these functions multiple time to get some statistically relevant results.

## Network load

## Easing simulation
bla
* Poisson Probability: reproduces the Poisson PDF (Probability Density Function), whose average depends on beta; probability of having infection at a certain time is just the integral in time of this function, from 0 to linking duration of I-S nodes. So, since we're dealing with discrete finite time, there's a finite number of integrals/probabilities that can be computed, T, since a contact may last from 0 to T. So, same integrals will be computed repeatedly, and it makes sense to compute them once for all, and put results in a dictionary, whose keys are the durations. This may save much computational time.
Most of these ideas are due to [Chen, Benzi paper](https://pdfs.semanticscholar.org/0cd5/46424d279a5a41f4cff3e863c1e0416b067f.pdf)

## propagation
A parte spiegare tutto, sottolinea che questa funzione fa evolvere un'epidemia dal primo all'ultimo istante del temporal network che gli fornisci, ma per dare uno score ai nodi dovresti farne evolvere tante sullo stesso temporal network, cambiando volta per volta l'index case. E poi, per ragioni di tesi, ripetere tutto quanto anche su più network.
Puoi usare sta roba:
Nodes virulence is determined by:
    * making each node, at a time, the index case, and 
    * computing min or avg time infection takes to spread to a certain % of the whole network
The idea of SI model is that each I node, at t-1, can make infect one of its S neighbours, at t-1, with a certain probability, whose rate is beta.
So, there will be usefuls some simple functions that build the neighbourhood of a node, or find S/I nodes from a list at a given time, or compute links' temporal duration.

le prime tre righe erano sopra contact lasting
Epidemic spread follows Chen approach: time of infection (in unit steps) follows a Poissonian distribution, normalized to return beta for 1 step, integrated within link duration.
(note: beta is the probability rate of contagion [1/s], but also the actual probability of contagion after 1 unit time: infact, P(1) = beta*1 u.t. = beta [dimensionless]).
Chen's algorithm make a contagion happen by performing an extraction from Unif(0,1): if this number is lower than the Poisson integral, contagion takes place at that time.

This function performs the actual propagation, from 0 to T, given an index case.
It finds the S neighboroughs of all I nodes at a given time, and makes them infect at t+1, with probability beta
This function will be evoked for each time step but the first and the last one

## Epidemic Scores

# Propagation_LTM
Further developments

# Tests.py
This script contains a suite of tests veryfing that every function in other scripts works properly, and some assertion functions, useful to these latter functions to check the inputs they are provided of. Most complex tests obviously deal with temporal networks and epidemics.
* Assertions section is quite self-explaining: these function impose inputs to be an np.array of a certain dimension, or to have a square shape, or to be probability matrices (so accepting only values from 0 to 1), or to be a natural number (integer and >=0), or to be a matrix with 0-diagonal. Most of these functions are imported by Evolutions and Propagations scripts.
* Tests section contains both structural tests (checking the output temporal network having the right mathematical properties) and tests of the actual evolution (since it's a stochastic process, one cannot make assertions about the exact outcome values, aside from some limit-cases; these ones are found and tested).

## Evolutions test functions
* structural_suite: since some structural tests will be repetead multiple times in both evolutions, they are put in this function, which collects some of the previous asserts and, in the end, checks that the output temporal network has the right parameters (number of nodes and duration) and mathematical properties, like being a succession of adiacencies, square metrices with null diagonal and, if case, symmetric (Explanation.md if you need to better understand these lines). This function is recalled in any DAR and TGRG test.
* Evolution/dynamic tests: as said, these tests check, besides the structural suite, some limit-case inputs, the only ones one can be sure of the outputs; since multiple combinations are possible, multiple tests of the same kind are performed, changing time by time some parameters like number of nodes and duration, just for sake of completenes.
  * DAR(1): we can be sure that, if matrix alpha = all zeros, all following states are determined by performing a random extraction (ruled by xi), while if alpha = all ones: all following states are equals to the first (total persistence); moreover, if matrix xi = all zeros, there's no way of getting state "1" for any link, if an extraction occurs, and vice versa for xi = all ones. So, these limit cases are tested: **alpha and xi = all zeros**, expecting a sequence of ONLY ZEROS adiacencies; **alpha = all zeros, xi = all ones** expecting a sequence of ONLY ONES (but null diagonal, of couse) adiacencies; **alpha = all ones, caringless of xi** expecting a sequence of adiacencies equal to the initial one.
  * TGRG: if "sigma" vector is set to 0, one removes randomness in fitnesses evolution (but not in link generation in following times); anyway, giving to each entry of phi0 and phi1 (or just phi0) very high (in module) values, one can jump into certainty domain. So, these limits will be checked, taking for guaranted that sigmas are 0: **very high values (100) for all phi0, and just 0 phi1**, expecting a sequence of ONLY ONES (but null diagonal) adiacencies; **very low values (-100) for all phi0, and just 0 phi1**, expecting a sequence of ONLY ZEROS adiacencies.
