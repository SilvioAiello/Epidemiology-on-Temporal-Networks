This file provides further explanations about test perfromed in Test Suite.
        Functions in this script provide:
        1) assertions functions to use at the beginning of those functions, to be sure their inputs are correct.
        2) tests to verify that functions in "Evolutions.py" and "Propagation_SI.py" run correctly.

# Table of contents
* [Assertions](#assertion-functions)
* [Evolution tests](#evolution-functions-tests)
* [Measures tests](#measures-functions-tests)

### assertion functions
This section is quite self-explaining: its function impose inputs to be an np.array of a certain dimension, or to have a square shape, or to be probability matrices (so accepting only values from 0 to 1), or to be a natural number (integer and >=0), or to be a matrix with 0-diagonal. Most of these functions are imported by Evolutions and Propagations scripts.

### evolution functions tests
This section contains both structural tests (checking the output temporal network having the right mathematical properties) and tests of the actual evolution (since it's a stochastic process, one cannot make assertions about the exact outcome values, aside from some limit-cases; these ones are found and tested).
* **structural_suite**: since some structural tests will be repetead multiple times in both evolutions, they are collected in this function, which performs some of the previous asserts and, in the end, checks that the output temporal network has the right parameters (number of nodes and duration) and mathematical properties, like being a succession of adjacencies, which in turn are square matrices with null diagonal and, if case, symmetrical (Explanation.md if you need to better understand these lines). This function is recalled in any DAR and TGRG test. If it is verified, it doesn't mean that networks are produced correctly, but just that their structure is how it was supposed to. So, this is just a preliminary test.
* **Evolution/dynamic tests**: as said, these tests check, besides the structural suite, some limit-case inputs, the only ones one can be sure of the outputs; since multiple combinations are possible, multiple tests of the same kind are performed, changing time by time some parameters like number of nodes and duration, just for sake of completenes.
  * DAR(1): we can be sure that, if matrix alpha = all zeros, all following states are determined by performing a random extraction (ruled by xi), while if alpha = all ones: all following states are equals to the first (total persistence); moreover, if matrix xi = all zeros, there's no way of getting state "1" for any link, if an extraction occurs, and vice versa for xi = all ones. So, these limit cases are tested: **alpha and xi = all zeros**, expecting a sequence of ONLY ZEROS adjacencies; **alpha = all zeros, xi = all ones** expecting a sequence of ONLY ONES (but null diagonal, of couse) adjacencies; **alpha = all ones, caringless of xi** expecting a sequence of adjacencies equal to the initial one.
  * TGRG: if "sigma" vector is set to 0, one removes randomness in fitnesses evolution (but not in link generation in following times); anyway, giving to each entry of phi0 and phi1 (or just phi0) very high (in module) values, one can jump into certainty domain. So, these limits will be checked, taking for guaranted that sigmas are 0: **very high values (100) for all phi0, and just 0 phi1**, expecting a sequence of ONLY ONES (but null diagonal) adjacencies; **very low values (-100) for all phi0, and just 0 phi1**, expecting a sequence of ONLY ZEROS adjacencies.

### communicability test
Test is performed by comparing a "manual" computation of a communicability matrix, and the one performed by the namesake function.
Temporal network is a 3-3-3 sequence of equal adjacences, with all nodes linked and no auto-loops; max eigenvalue is known being 2, so the coefficient multipling each adjacency (check [Centrality Measure section](https://github.com/SilvioAiello/Epidemiology-on-Temporal-Networks/blob/master/docs/explanation.md#centrality-measures)) is 0.125.
Since the two computations follow slightly processes, roundoff may give some troubles in comparing the resulting matrices, so we don't state them to be equal, but just require the difference of each entry to be less the 10^(-10), which is several order of magnitude greater then the mean deviation between all entries.

### measures functions tests
NEXT STEP
