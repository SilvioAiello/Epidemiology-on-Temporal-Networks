# Epidemiology on Temporal Networks
The purpose of this project is to simulate epidemics on temporal networks, try to correlate it to netwroks' structure, and apply results to real data, took from eMID interbank network.

A network is "temporal" when the existence of a link between two nodes depends on time: distribution of links evolves, while the number of nodes is held fixed.
Several laws of evolution exist; in this study we deal with **DAR(p)** and **TGRG**. DEFINE. Also, two ways of epidemic spreading will be analyzed: **SI** and **LTM**. DEFINE.

The law of network's evolution influences epidemic's diffusion. Understanding wheter, and how much, is one of the mainly coped with topics of this study, and can be achieved by trying to quantify, whether exists, the correlation between *temporal centraility measures* (DEFINE), and *epidemic performances* of nodes (DEFINE).

The analysis will become step by step more complex, by violating symmetries such as network underictness, homogeneity of node-parameters; so, the scripts will be written from the beginning in a form as general as possible. 