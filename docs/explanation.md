# Introduction
This text provide the user of some theoretical knowledge necessary to better understand the mathematical background of this project. If you already know what a temporal network is, what are DAR(P) and TGRG evolutions, what's SI and LTM epidemiology, you can skip reading; otherwise, it is scrictly suggested.

A network is "temporal" when the existence of a link between two nodes depends on time: distribution of links evolves, while the number of nodes is held fixed.
Several laws of evolution exist; in this study we deal with **DAR(P)** and **TGRG**. Also, two ways of epidemic spreading will be analyzed: **SI** and **LTM**.

# Table of contents
* [Introduction](#introduction)
* [Strucutal evolution](#strucutal-evolution)
  * [DAR(P)](#darp)
  * [TGRG](#tgrg)
* [Epidemic Diffusion](#epidemic-diffusion)
  * [SI](#si)
  * [LTM](#ltm)

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
