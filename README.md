# Epidemiology on Temporal Networks
The purpose of this project is to simulate epidemics on temporal networks, try to correlate it to netwroks' structure, and apply results to real data, took from eMID interbank network.

The analysis will become step by step more complex by violating symmetries such as network underictness, homogeneity of node-parameters; so, the scripts are written from the beginning in a form as general as possible, although many parameters may present the same value, at the moment. 

To use the scripts how this project as they are, you can just open file "main.py", modify the input structures and make it run to get results, that are saved in folder "Networks"

# Table of Contents
* [Project purposes](#project-purposes)
* [References](#references)
* [docs/Explanation (theoretical concepts)](https://github.com/SilvioAiello/Epidemiology-on-Temporal-Networks/blob/modules-integration/docs/explanation.md)
* [docs/How to (how code works)](https://github.com/SilvioAiello/Epidemiology-on-Temporal-Networks/blob/modules-integration/docs/howto.md)

# Project purposes
Epidemic's diffusion should relate to laws of evolution: understanding wheter, and how much, may be achieved by trying to quantify, wheter exists, correlation between *temporal centraility measures*, and *epidemic performances* of nodes. If a node is "central", we expect him to be epidemiologically virulent.

*Centrality measures* state a way to recognize nodes more connected to others, and each measure gives a particoular meaning to this sentence. There is an extended literature about centralities in static networks rather then in temporal, where there are changes and rising in complexity in computing, due to dependence on link existence. Here, are taken into account these measures: **Broadcast** and **Receive Centrality** (BC, RC), took from a matrix named Communicability [Grindrod, Communicability across evolving networks], **Aggregate Degree**, **Binarized Degree**.

*Virulence* of a node can be measured in various ways. For example, it can be expressed as the **average**, or **minimum**, number of time steps necessary to infect a certain percentage of node population, provived that that node is the "index case". So, for a set realization of a temporal network, one should simulate several realization of the same epidemic, defined by having that node as only initial infected, and extract the minimum of average time that led to an infection of, for example, 60% of nodes. This has to be iterated over all nodes, getting a score for each one. 

Once this is done, it is possible to comprare these scores to nodes centralities, extracting an index of correlation. To be sure that the informations one get are not just due to a fortuitous realtization, the whole process shoud be iterated over different temporal realtizations of the same network.

# References
1. [Mazzarisi, Lillo et. al., A dynamic network model with persistent links and node-specific latent variables, with an application to the interbank market](https://arxiv.org/pdf/1801.00185.pdf)
2. [Grindrod et. al., Communicability Across Evolving Networks](http://centaur.reading.ac.uk/19357/1/Coomunicability_accepted.pdf)
3. [Chen et. al., Dynamic communicability and epidemic spread: a case study on an empirical dynamic contact network](https://pdfs.semanticscholar.org/0cd5/46424d279a5a41f4cff3e863c1e0416b067f.pdf)



