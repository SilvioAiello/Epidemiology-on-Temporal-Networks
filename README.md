# Epidemiology on Temporal Networks

# Table of Contents
* [Repository scheme](#repository-scheme)
* [Project purposes](#project-purposes)
* [References](#references)
* [docs/Explanation (theoretical concepts)](/docs/explanation.md)
* [docs/How to (how code works)](/docs/howto.md)

# Repository scheme
This repository is structured as follows:
* **Main folder**: 
  * .gitgnore, \\_\_pycache\_\_\
  * README.md
  * Evolutions.py: suite of functions that generate and analyze temporal networks (and save results);
  * main.py: allows to smartly manage the whole set of functions and provide them of the desidered inputs;
  * Propagation_SI.py: suite of functions that perform SI propagations upon networks (generated by Evolutions or already available to the user), measure virulence scores and save results;
  * Saves.py: suite of functions allowing to save data structures;
  * Test_Suite.py,: all functions are verified to work correctly.
* **\docs**: complete documentation files
  * explanation.md: deepening of theorical concepts this project deals with;
  * howto.md: deepening of how each single script works, and how to use it properly
* **\Networks**: in this folder networks and analysis results are stored, as **pickle** (native data serialization module in Python, check [here](https://docs.python.org/3/library/pickle.html) for further explanations) files, with this "rule of thumb": one folder per network (defined by its input parameters), whose name keeps track of *number of nodes*, *duration*, *type of evolution* (DAR/TGRG) and an *identification name*; into each one, one can find 
  1) as many folders as iterations of network evolution (*/realization/*), with inside network (*/network.pkl*) and epidemic propagation (with a name that records beta value and propagation iteration) files; 
  2) a folder (*/Results/*) where epidemic and structural centrality measures are stored.
 
So, the final output will be something like: Networks/N_T_TYPE_name/realizationX/beta_k.pkl or /network.pkl.

**To make your personal use of this project, you can just open file "main.py", modify the input structures and let it run to get results. Scripts should not be modified; you should only modify values in "USER ACTION" section of main script.** Once main has run, you can keep on operating from your Python IDE; for example, to call a function from "Evolutions", you may just type Evolutions.FUNCTNAME, etc.

# Project purposes
The purpose of this project is to simulate epidemics on temporal networks, try to correlate it to netwroks' structure, and apply results to real data, took from eMID interbank network. The analysis will be made more complex, and realistic, by violating symmetries such as network underictness and homogeneity of node-parameters. 

Epidemic's diffusion should relate to laws of evolution: understanding wheter, and how much, may be achieved by trying to quantify, wheter exists, correlation between *temporal centraility measures*, and *epidemic performances* of nodes. If a node is "central", we expect him to be epidemiologically virulent.

*Centrality measures* state a way to recognize nodes more connected to others, and each measure gives a particoular meaning to this sentence. There is an extended literature about centralities in static networks rather then in temporal, where there are changes and rising in complexity in computing, due to dependence on link existence. Here, are taken into account these measures: **Broadcast** and **Receive Centrality** (BC, RC), took from a matrix named Communicability [Grindrod, Communicability across evolving networks], **Aggregate Degree**, **Binarized Degree**.

*Virulence* of a node can be measured in various ways. For example, it can be expressed as the **average**, or **minimum**, number of time steps necessary to infect a certain percentage of node population, provived that that node is the "index case". So, for a set realization of a temporal network, one should simulate several realization of the same epidemic, defined by having that node as only initial infected, and extract the minimum of average time that led to an infection of, for example, 60% of nodes. This has to be iterated over all nodes, getting a score for each one. 

Once this is done, it is possible to comprare these scores to nodes centralities, extracting an index of correlation. To be sure that the informations one get are not just due to a fortuitous realtization, the whole process shoud be iterated over different temporal realtizations of the same network.

# References
1. [Mazzarisi, Lillo et. al., A dynamic network model with persistent links and node-specific latent variables, with an application to the interbank market](https://arxiv.org/pdf/1801.00185.pdf)
2. [Grindrod et. al., Communicability Across Evolving Networks](http://centaur.reading.ac.uk/19357/1/Coomunicability_accepted.pdf)
3. [Chen et. al., Dynamic communicability and epidemic spread: a case study on an empirical dynamic contact network](https://pdfs.semanticscholar.org/0cd5/46424d279a5a41f4cff3e863c1e0416b067f.pdf)
