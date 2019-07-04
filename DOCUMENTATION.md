EPIDEMIOLOGY ON TEMPORAL NETWORKS

This repository, at this moment, is structured as follows:
Directory "Networks", strucutured as follows:
	A folder for each kind of network, defined by its type (DAR,TGRG), parameters (eg. p, alphas and xis for DAR(P)), number of nodes and temporal duration. Its name resumes this parameters. Each folder contains:
		A folder for each realization of that network, named "realizationX". Each folder contains:
			1 file named "network.txt", that stores the whole set of adiacency matrices defining the realization;
			A file, one for each value of beta (SI), named "betaX_labelTYPE_K%.txt containing the evolutions of states of node in time, for each index case, for each (of the total K) iterations.
			A file, one for each LTM. This will be delt with in the future.
Directory "Results", where the final results (images, maybe something else) will go.
2 python "Evolution" modules, each of which simulates the temporal evolution of the network (i.e. of links between nodes) alone.
1 python "Propagation" module, which simulate the SI spreading of a desease
1 test module
2 md files: this, and the README.
The gitignore file