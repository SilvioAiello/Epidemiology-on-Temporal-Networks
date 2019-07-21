This file provides further explanations about tests perfromed in Test Suite.
        Functions in this script provide:
        1) assertions functions to use at the beginning of those functions, to be sure their inputs are correct.
        2) tests to verify that functions in "Evolutions.py" and "Propagation_SI.py" run correctly.

# Table of contents
* [Assertions](#assertion-functions)
* [Evolution tests](#evolution-functions-tests)
* [Measures tests](#measures-functions-tests)

# Assertion functions
This project makes deep use of assertion as preliminary checks of **input** parameters for functions. It's known that assertions may be easily muted, so they don't guarantee high security standards, but for our purposes they're enough. 
Since some checks will require more than one assertions and more then one performance, they are collected in some functions, located here and mostly imported by **Evolutions** and **Propagations** scripts. These functions impose inputs to be an np.array of a certain dimension, or to have a square shape, or to be probability matrices (so accepting only values from 0 to 1), or to be a natural number (integer and >=0), or to be a matrix with 0-diagonal.

One larger function, used only in this script for **output** datas, is **Structural_suite**: since some structural tests will be repetead multiple times in DAR/TGRG tests, they are collected here, to ensure a temporal network has the right parameters (number of nodes and duration) and **mathematical properties**. If they are passed, it doesn't mean that output networks are produced correctly, but just that their structure is how it was supposed to.

# Evolution.py tests
Since it's a stochastic process, one cannot make assertions about the exact outcome values of network evolution functions, aside from some limit-cases, that are found and checked. The other functions are easier to prove right.

* **generations**: structural suite and limit cases are checked; since multiple limit combinations are possible, multiple tests of the same kind are performed, changing time by time some parameters like number of nodes and duration, just for sake of completenes.
  * DAR(1): we can be sure that, if matrix alpha = all zeros, all following states are determined by performing a random extraction (ruled by xi), while if alpha = all ones: all following states are equals to the first (total persistence); moreover, if matrix xi = all zeros, there's no way of getting state "1" for any link, if an extraction occurs, and vice versa for xi = all ones. So, these limit cases are tested: **alpha and xi = all zeros**, expecting a sequence of ONLY ZEROS adjacencies; **alpha = all zeros, xi = all ones** expecting a sequence of ONLY ONES (but null diagonal, of couse) adjacencies; **alpha = all ones, caringless of xi** expecting a sequence of adjacencies equal to the initial one.
  * TGRG: if "sigma" vector is set to 0, one removes randomness in fitnesses evolution (but not in link generation in following times); anyway, giving to each entry of phi0 and phi1 (or just phi0) very high (in module) values, one can jump into certainty domain. So, these limits will be checked, taking for guaranted that sigmas are 0: **very high values (100) for all phi0, and just 0 phi1**, expecting a sequence of ONLY ONES (but null diagonal) adjacencies; **very low values (-100) for all phi0, and just 0 phi1**, expecting a sequence of ONLY ZEROS adjacencies.

* **degrees**: functions are tested by comparing results one gets by applaying definition (sum over lines/columns), over some known networks, for both out and in degree

* **communicability**: test is performed by comparing a "manual" computation of a communicability matrix, and the one performed by the namesake function. Temporal network is a 3-3-3 sequence of equal adjacences, with all nodes linked and no auto-loops; max eigenvalue is known being 2, so the coefficient multipling each adjacency (check [Centrality Measure section](https://github.com/SilvioAiello/Epidemiology-on-Temporal-Networks/blob/master/docs/explanation.md#centrality-measures)) is 0.125.
Since the two computations follow slightly processes, some differences may raise at last decimals, so we don't state results to be equal, but just to differ within a very small number as 10^(-10), which is several orders of magnitude greater then the average deviation between all communicability values.

# Propagation_SI.py tests
Also in this case, stochastic functions are tested through their limit cases.

* **neighbourhood**: in this project, a node neighbourhood is the set of nodes that can be infected by another one. Still forbiding auto-loops, "empty" and "full" networks are verified to have all nodes with empty (0) and full (N-1) neighbourhoods. One node, from the full network, is also selected and deprived, time by time, of all outgoing links, verifying its neighbours (i.e. infectable nodes) to decrease, while all the others are unaffected. This proves that the way this functions has been built, avoids user of caring about out- or in- going links, since functions deals with it by itself.

* **contact_lasting** this fundamental function is tested in two ways:
1) by generating a random temporal network, and verifying that temporal durations are smaller than the inquired temporal step, greater than 1 if exist, and always smaller than duration at next step, if it still exists (this may sound obvious, but it's important to ensure temporal solidity of measures;
2) by generating a network with only one couple of nodes linked forever and with only one of them infected; duration is tested to be right step by step, then two "accidents" take place: links disappear at a certain time, and duration is checked again and expected to be smaller, then it is restoted but, at the same time (and only for that moment), susceptible node is made infected, and result is checked to be smaller, in the same way.

* **propagation_SI**:
