# Introduction
This text provide the user of some theoretical knowledge necessary to better understand the mathematical background of the project. If you already know what a temporal network, DAR(P) and TGRG evolutions, SI and LTM epidemiology are, you can skip reading; otherwise, it is scrictly suggested to give a look.

A static network (or graph)  can be defined by a set of nodes and link between nodes, where nodes represent "entities", and links a certain interaction. A static network can be described by an adiacency matrix, whose entry "ij" states the number of links between node i and j (if it can be >=2, we talk about multi-graph), and may not be equal to entry "ji", stating that the system is sensible to which nodes starts the interaction (directed graph).
A network is "temporal" when the existence of a link between two nodes depends on time: distribution of links evolves, while the number of nodes is held fixed, and we call a temporal network and _ordered_ sequence of static networks (described by an ordered sequence of adiacencies). 
The arrow of time provides such a structure of an inherent asymmetry, even if the single graphs are undirected, that emerges in information trasmission between nodes. Epidemiology deals with information, or disease, transmission, and gets effective contribution my modelling systems as temporal networks.

Several laws of evolution exist, each of wich may effect differently an epidemic spread; in this study we deal with **DAR(P)** and **TGRG** evolutions. Also, two ways of epidemic spreads will be analyzed: **SI** and **LTM**.

As already mentioned in ReadMe, we want to to understand if there is a correlation between nodes' centrality and capability in broadcasting/receiving informations, and, if it is the case, to quantify it, or, at least, determine what centrality measures better catch it.

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

If more links are allowed between the same nodes (this not being the case here), ![equation](https://latex.codecogs.com/png.latex?A_{ij}(t_{k})) can have value 2,3,... .

If we don't care about which from the two nodes "starts" the link, the graph is said to be _undirected_, and ![equation](https://latex.codecogs.com/png.latex?A_{ij}(t)&space;=&space;A_{ji}(t)), so the adiacency is symmetrical; otherwise, it is _directed_.
If auto-loops, i.e. link with self, are not contemplated, ![equation](https://latex.codecogs.com/png.latex?A_{ii}(t)&space;=&space;0) ![equation](https://latex.codecogs.com/png.latex?\forall&space;t,i).

Both evolutions under examination have a stochastic nature, so, to infer some properties of the networks, one should consider more realization of the graph, i.e. more evolutions with the same defining-parameters.

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
(just basic scattered concepts)
Degree: it can be in or out-going; if graph is undirected, they are the same
These measures are dicrect generalization of static ones (if T->0, they return the same results), which are mostly based on _walks_ (information does not necessarly flow across geodesics, and walk counting is more resilient to missing or spourious edges, cfr. Grindrod); so we rely on temporal graphs,with Grindrod noticing that 1) repeated times are allowed: t1 < t2 = t3 < t4, i.e. in the same instant there can be more than 1 progression, 2) times are not required to be consecutive (if t2 > t1+1, times between t1 and t2 are "ignored", just like you don't care of what's appening in that interval). A the moment, times are not repeated and are consecutive.
Non commutativity of matrices product is how we account the arrow of time

Working progress (Communicability, AD, BD)

# Epidemic Diffusion
When dealing with epidemiology, for each node is defined a state: epidemic maps each time step to the state of each node at that time. For SI and LTM epidemiologies, there are two possible state: infected ("1"), and susceptible ("0", meaning that that node can be infected).
Both diffusions are stochastic processes, so, in both case, one should perform more than one diffusion, in order to get some useful information.
Temporal varying structure of networks leads to different results than for static networks 

## SI
This propagation is not specific for network theory, and it's effective in many systems of various nature.
It's defined by just one parameter, beta, expressing the probability of infection per unit time, i.e. the probability that one infected node infects a susceptible after one time step (in a system with discrete time).
According to [Chen et. at.](https://pdfs.semanticscholar.org/0cd5/46424d279a5a41f4cff3e863c1e0416b067f.pdf), this is an indipendent and memory-less Poisson process, whose mean (lambda\*) is set by setting beta and correlating it to the probability per unit time T (1 second, 1 minute, 1 day...):

![equation](https://latex.codecogs.com/png.latex?\lambda^*&space;:) P(infection in T) = ![equation](https://latex.codecogs.com/png.latex?\int_0^T&space;\lambda^*\exp(-\lambda^*t))

Probability of infection, for a contact lasting t, is given by integrating this function from 0 to t.
SISTEMA UN PO' STA ROBA
time of infection (in unit steps) follows a Poissonian distribution, normalized to return beta for 1 step, integrated within link duration. (note: beta is the probability rate of contagion [1/s], but also the actual probability of contagion after 1 unit time: infact, P(1) = beta*1 u.t. = beta [dimensionless]).
Most of these ideas are due to [Chen, Benzi paper](https://pdfs.semanticscholar.org/0cd5/46424d279a5a41f4cff3e863c1e0416b067f.pdf)

## LTM
This section will be deepened in further developments.

## Epidemic scores
Node virulence is determined by:
    * making a node an epidemic index case, and 
    * computing minimum or avaverageg time (over a range of K iterations, with the same initial condition) the infection takes to reach to a certain percentage of the whole network.
Choosing minimum or average may lead to different scores; at the moment, only average is considered.

