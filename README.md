# Epidemiology on Temporal Networks
The purpose of this project is to simulate epidemics on temporal networks, try to correlate it to netwroks' structure, and apply results to real data, took from eMID interbank network.

TODO: add table of contents and references

A network is "temporal" when the existence of a link between two nodes depends on time: distribution of links evolves, while the number of nodes is held fixed.
Several laws of evolution exist; in this study we deal with **DAR(P)** and **TGRG**. Also, two ways of epidemic spreading will be analyzed: **SI** and **LTM**.

The analysis will become step by step more complex by violating symmetries such as network underictness, homogeneity of node-parameters; so, the scripts are written from the beginning in a form as general as possible, although many parameters may present the same value, at the moment. 

# Table of Contents
* [Project purposes](#project-purposes)
* [References](#references)
* [Strucutal evolution](#strucutal-evolution)
  * [DAR(P)](#darp)
  * [TGRG](#tgrg)
* [Epidemic Diffusion](#epidemic-diffusion)
  * [SI](#si)
  * [LTM](#ltm)
* [Documentation](https://github.com/SilvioAiello/Epidemiology-on-Temporal-Networks/blob/master/DOCUMENTATION.md)

# Project purposes
Epidemic's diffusion should relate to laws of evolution: understanding wheter, and how much, may be achieved by trying to quantify, wheter exists, correlation between *temporal centraility measures*, and *epidemic performances* of nodes. If a node is "central", we expect him to be epidemiologically virulent.

*Centrality measures* state a way to recognize nodes more connected to others, and each measure gives a particoular meaning to this sentence. There is an extended literature about centralities in static networks rather then in temporal, where there are changes and rising in complexity in computing, due to dependence on link existence. Here, are taken into account these measures: **Broadcast** and **Receive Centrality** (BC, RC), took from a matrix named Communicability [Grindrod, Communicability across evolving networks], **Aggregate Degree**, **Binarized Degree**.

*Virulence* of a node can be measured in various ways. For example, it can be expressed as the **average**, or **minimum**, number of time steps necessary to infect a certain percentage of node population, provived that that node is the "index case". So, for a set realization of a temporal network, one should simulate several realization of the same epidemic, defined by having that node as only initial infected, and extract the minimum of average time that led to an infection of, for example, 60% of nodes. This has to be iterated over all nodes, getting a score for each one. 

Once this is done, it is possible to comprare these scores to nodes centralities, extracting an index of correlation. To be sure that the informations one get are not just due to a fortuitous realtization, the whole process shoud be iterated over different temporal realtizations of the same network.

# References
1. [Mazzarisi, Lillo et. al., A dynamic network model with persistent links and node-specific latent variables, with an application to the interbank market](https://arxiv.org/pdf/1801.00185.pdf)
2. [Grindrod et. al., Communicability Across Evolving Networks](http://centaur.reading.ac.uk/19357/1/Coomunicability_accepted.pdf)
3. [Chen et. al., Dynamic communicability and epidemic spread: a case study on an empirical dynamic contact network](https://pdfs.semanticscholar.org/0cd5/46424d279a5a41f4cff3e863c1e0416b067f.pdf)

# Strucutal evolution
State of links is described by an adiacency matrix A(t) for each time step t. If a link between nodes i and j exists, A_ij = 1; otherwise, A_ij = 0. If more links are allowed between the same nodes (this not being the case here), A_ij can have value 2,3,... .
If we don't care about which from the two nodes "starts" the link, the graph is said to be _undirected_, and A_ij = A_ji, so the adiacency is symmetrical; otherwise, it is _directed_.
If auto-loops, i.e. link with self, are not contemplated, A_ii = 0 for each i.
Both evolutions under examination have a stochastic nature, so, to infer some properties of the networks, one should consider more realization of the graph, i.e. more evolutions with the same defining-parameters.

## DAR(P)
DAR(P) stands for Discrete Auto Regressive of order P. It means that the state of a link (on or off), at time t, depends on the states at the P previous time steps. Let's say P=1: at time t, the link between nodes i and j has probability alpha_ij of persisting, i.e. of being the same as at t-1, and probability 1-alpha_ij of being extracted by a "tossing coin", with probability xi_ij of getting "on" (so, the probability of having k links at a certain time is a binomial). 
If P>1, alpha_ij is just decomposed into a probability for each of the P previous states.
For a DAR(1), there are N^2 parameters to set, N being the number of nodes.

## TGRG
TGRG stands for Temporally Generalized Random Graphs. The existence of a link between two nodes depends on the sum of a quantity defined for each noded, named Fitness, theta(t). Fitness evolves in time accordind to this law:

> theta_i(t) = phi0_i + phi1_i*\theta_i(t-1) + epsilon_i(t)

So, at each time, its value depends to the one at the previous time step, according a factor phi1, to constant term phi0 that is added time by time, e to a stochastic term, epsilon, which is extracted from a distribution with mean = 0 and average = sigma^2_i (tpical of each node). 
Now, the existence of a link at a certain time is a stochastic extraction, whose probability relies on a kind sigmoid of the sum of the two nodes' fitnesses:

> P(A_ij(t) = 1) = exp(theta_i(t)+theta_j(t))/(1+exp(theta_i(t)+theta_j(t)))

So, the sum is high, there's almost certainty of having a link, and viceversa.
For a TGRG, there are 3N parameters to set.

# Epidemic Diffusion
When dealing with epidemiology, for each node is defined a state: epidemic maps each time step to the state of each node at that time. For SI and LTM epidemiologies, there are two possible state: infected ("1"), and susceptible ("0", meaning that that node can be infected).
Both diffusions are stochastic processes, so, in both case, one should perform more than one diffusion, in order to get some useful information.
Temporal varying structure of networks leads to different results than for static networks 

## SI
This propagation is not specific for network theory, and it's effective in many system of various nature.
It's defined by just one parameter, beta, expressing the probability of infection per unit time, i.e. the probability that one infected node infects a susceptible after one time step (in a system with discrete time).
According to [Chen et. at.](https://pdfs.semanticscholar.org/0cd5/46424d279a5a41f4cff3e863c1e0416b067f.pdf), this is an indipendent and memory-less Poisson process, where the mean (lambda\*) is set by setting beta: 

> lambda\* | beta\*(1 time unit) = P(infection in 1 time unit) = integral(from 0 to time unit) lambda\* \* exp(-lambda\*t).

In general, probability at time t is given by integrating from 0 to t.

## LTM
This section will be deepened in further developments.

