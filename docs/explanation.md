# Introduction
This text better clarifies the mathematical background of the project. If you already know what a temporal network, DAR(P) and TGRG evolutions, SI and LTM epidemiology are, you can skip reading; otherwise, it is suggested to give a look.

A **static network** (or graph)  can be defined by a set of **nodes** and **links** between nodes, where nodes represent "entities", and links a certain interaction between them. A static network can be described by an **adiacency matrix**, whose entry "ij" states the number of links between node i and j (if it can be >=2, we talk about **multi-graph**), and may not be equal to entry "ji", stating that the system is sensible to which nodes starts the interaction (**directed graph**).
A network is **temporal** when the existence of a link between two nodes depends on time: distribution of links evolves, while the number of nodes is held fixed, and we call a temporal network and _ordered_ sequence of static networks (described by an ordered sequence of adiacencies). 
The *arrow of time* provides such a structure of an inherent asymmetry, even if the single graphs are undirected, that emerges in information trasmission between nodes. Epidemiology deals with information/disease transmission in different systems, and gets an effective contribution by modelling systems as temporal networks.

Several laws of evolution exist, each of wich may effect differently an epidemic spread; in this study we deal with **DAR(P)** and **TGRG** evolutions. Also, two ways of epidemic spreads will be analyzed: **SI** and **LTM**.

As already mentioned in [Read Me](https://github.com/SilvioAiello/Epidemiology-on-Temporal-Networks/blob/master/README.md), we want to understand if there is a correlation between nodes' centrality and capability in broadcasting/receiving informations, and, if it is the case, to quantify it, or, at least, determine what centrality measures better catch it.

# Table of contents
* [Introduction](#introduction)
* [Strucutal evolution](#strucutal-evolution)
  * [DAR(P)](#darp)
  * [TGRG](#tgrg)
  * [Centrality measures](#centrality-measures)
* [Epidemic Diffusion](#epidemic-diffusion)
  * [SI](#si)
  * [LTM](#ltm)
  * [Epidemic scores](#epidemic-scores)
 

# Strucutal evolution
State of links is described by an adiacency matrix A(t) for each time step t. If a link between nodes i and j exists at time ![equation](https://latex.codecogs.com/gif.latex?t_{k}), then ![equation](https://latex.codecogs.com/png.latex?A_{ij}(t_{k})&space;=&space;1); otherwise, ![equation](https://latex.codecogs.com/png.latex?A_{ij}(t_{k})&space;=&space;0). 

If more links are allowed between the same nodes (this not being the case here), ![equation](https://latex.codecogs.com/png.latex?A_{ij}(t_{k})) can have value 2,3,... . If the graph is _undirected_, ![equation](https://latex.codecogs.com/png.latex?A_{ij}(t)&space;=&space;A_{ji}(t)), so the adiacency is symmetrical. If auto-loops, i.e. link with self, are not contemplated, ![equation](https://latex.codecogs.com/png.latex?A_{ii}(t)&space;=&space;0) ![equation](https://latex.codecogs.com/png.latex?\forall&space;t,i).

Both evolutions under examination have a stochastic nature: to infer some properties of the networks, one should consider more realization of the graph, i.e. more evolutions with the same defining-parameters.

## DAR(P)
DAR(P) stands for Discrete Auto Regressive of order P. It means that the state of a link (on or off), at time t, depends on the states at the P previous time steps. Let's say P=1: at time t, the link between nodes i and j has probability ![equation](https://latex.codecogs.com/png.latex?\alpha_{ij}) of **persisting**, i.e. of being the same as at t-1, and probability ![equation](https://latex.codecogs.com/png.latex?1-\alpha_{ij}) of being extracted by a "**tossing coin**", with probability ![equation](https://latex.codecogs.com/png.latex?\xi_{ij}) of getting "on" (so, the probability of having k links between i and j, after N tries, is a binomial of ![equation](https://latex.codecogs.com/png.latex?\xi)). 
If P>1, alpha_ij is just decomposed into a probability for each of the P previous states.
For a DAR(1), there are N^2 parameters to set, N being the number of nodes.

## TGRG
TGRG stands for Temporally Generalized Random Graphs. The existence of a link between two nodes depends on a nodes' property named Fitness, theta(t). Fitness evolves in time accordind to this law:

![equation](https://latex.codecogs.com/png.latex?\theta_i(t)&space;=&space;\phi_{0,i}&space;&plus;&space;\phi_{1,i}\theta_i(t-1)&space;&plus;&space;\epsilon_i(t))

So, at each time, its value depends to the one at the previous time step, according a factor phi1, to constant additive term phi0, and to a stochastic term, epsilon, which is extracted from a Normal(0, sigma) distribution (with sigma being typical of each node). 
Now, the existence of a link at each temporal step ![equation](https://latex.codecogs.com/png.latex?t_k) is a stochastic extraction, whose probability relies on a kind sigmoid function of the sum of the two nodes' fitnesses:

![equation](https://latex.codecogs.com/png.latex?P(A_{ij}(t_k)&space;=&space;1)&space;=&space;\frac{\exp(\theta_i(t_k)&plus;\theta_j(t_k))}{1&plus;\exp(\theta_i(t_k)&plus;\theta_j(t_k))})

So, the sum is high, there's almost certainty of having a link, and viceversa.
For a TGRG, there are 3N parameters to set.

## Centrality Measures
Express how a node is "central". Better: they define a way to quantify centrality of nodes, impling that it is not unique. Measures for *temporal networks* are dicrect generalization of static ones (if T->0, they return the same results), which are mostly based on _walks_, since information does not necessarly flow across geodesics, and walk counting is more resilient to missing or spourious edges, (cfr. Grindrod, a walk is a set of links that lead from a node to another, allowing to pass to same nodes multiple times; paths doesn't allow this). Non-commutativity of matrices product, in these measures, is how we account the arrow of time.

Here are some useful to know static measures:
* **Degree**: it can be in or out-going; if graph is undirected, they are the same.
* **Katz Centrality** [from Grindod et al.](http://centaur.reading.ac.uk/19357/1/Coomunicability_accepted.pdf): To quantify the propensity for a node i to communicate, or interact, with another node j, we may count how many walks of length w = 1,2,3,... there are from i to j. We may then combine these counts into a single, cumulative total over all w. Allowing for the fact that shorter walks are generally more important it makes sense to scale the counts according to the walk length. A particularly attractive choice is to scale walks of length w by a factor a^w, where a is a suitably chosen scalar. A basic identity from graph theory shows that the k-th power
of the adjacency matrix has i, j element that counts the **number of walks of length w from node i to node j**. This leads us to the expansion ![equation](https://latex.codecogs.com/png.latex?1&space;&plus;&space;aA&space;&plus;&space;a^2A^2&space;&plus;...), that converges to resolvent function 
![equation](https://latex.codecogs.com/png.latex?(I-aA)^{-1}),
if *a* is less than the reciprocal of the maximum spectral radius of all temporal adjacencies. This last expression is still a matrix, summarizing how well information can pass from node i to node j, so summing over the n-th row gives an index of well that nodes broadcasts information to the others and is called Katz Centrality, while summing over the n-th column expresses easiness to get information.

Dynamic measures used in this project are:
* **Broadcast/Receive Centrality**: are a generalization of Katz Centrality, taking into account that A have different values in time, so the final matrix has as many factor as many adjacencies, and is called **Communicability**: 
![equation](https://latex.codecogs.com/png.latex?Q&space;=&space;(I-aA^{[0]})^{-1}(I-aA^{[1]})^{-1}...(I-aA^{[T-1]})^{-1})
One can choose to normalize this matrix by diving each factor by its norm (it doesn't which one).
From this matrix, as in Katz, one can infer nodes' capability to broadcast (**BC**) information, by summing over lines, or to receive it (**RC**), by summing over columns,
* **Aggregated Degree**, ranks nodes according to the number of distinct contacts weighted by the duration of contact time,
in the sense that the longer the contact time between a pair of nodes, the higher their corresponding AD values.
* **Binarized Degree**, where most of temporal information is lost, just counts the number of times a contact has existed between nodes i and j, so is a unweighted version of AD.

# Epidemic Diffusion
When dealing with epidemiology, for each node is defined a state: epidemic maps each time step to the state of each node at that time. For SI and LTM epidemiologies, there are two possible state: infected ("1"), and susceptible ("0", meaning that that node can be infected).
Both diffusions are stochastic processes, so, in both case, one should perform more than one diffusion, in order to get some useful information.
Temporal varying structure of networks leads to different results than for static networks.

## SI
This propagation is not specific for complex networks theory, and it's effective in many systems of various nature.
It's defined by just one parameter, beta, expressing the probability of infection per unit time, i.e. the probability that one infected node infects a susceptible after one time step (in a system with discrete time).
According to [Chen et. at.](https://pdfs.semanticscholar.org/0cd5/46424d279a5a41f4cff3e863c1e0416b067f.pdf), this is an indipendent and memory-less Poisson process, whose mean (lambda\*) is set by setting beta and correlating it to the probability per unit time T (1 second, 1 minute, 1 day...):

![equation](https://latex.codecogs.com/png.latex?\lambda^*&space;:) P(infection in T) ![equation](https://latex.codecogs.com/png.latex?\equiv&space;\beta&space;T) = ![equation](https://latex.codecogs.com/png.latex?\int_0^T&space;\lambda^*\exp(-\lambda^*t))

Here we choose T = 1 time step. So, beta now is the probability rate of contagion (dimensionally 1/time), but also the actual probability of contagion after 1 unit time: infact, P(1) = beta*1 u.t. = beta (dimensionless); this bounds its values between 0 and 1.

Probability of infection, for a contact lasting t, is given by integrating Poisson function from 0 to t; normalization guarantees this probability to be 1 for t approaching infinite. When dealing with discrete finite time, there's a finite number (T) of integrals / probabilities that can be obtained, since a contact may last only from 0 to T. So, since same integrals will be computed repeatedly during simulations, it makes sense to do it once for all for each time step. This may save much computational time.


## LTM
This section will be deepened in further developments.

## Epidemic scores
Node virulence is determined by:
    * making a node an epidemic index case, so we can be sure that infection time depends on its, and neighbours', "centrality" and 
    * computing minimum or avaverageg time (over a range of K iterations, with the same initial condition) infection takes to reach to a certain percentage of the whole network.
Choosing minimum or average may lead to different scores; at the moment, only average is considered.
