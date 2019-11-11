# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 18:41:08 2019

@author: Silvio
"""
import numpy as np
import Propagation_SI
from scipy.integrate import quad
import matplotlib.pyplot as plt

T = 50
N = 50
beta = 0.015

probabilities = dict() #probabilities dict
for t in range(T):
    probabilities[t] = quad(Propagation_SI.poisson_probability,0,t, args = beta)[0] #quad produces several outputs, integral is the first


temporal_network= np.random.choice((0,1),p=(0.5,0.5),size = (T,N,N))

for t in range(T):
    for i in range(N):
        temporal_network[t,i,i] = 0
    temporal_network[t] = np.triu(temporal_network[t]) + np.triu(temporal_network[t]).T
plt.figure(1)
plt.hist(temporal_network.reshape(T*N*N,1))

index_case = 0
labels, nodes_stats = Propagation_SI.propagation_senzaset(temporal_network, index_case, probabilities, multiple_infections = True)

inftime = np.zeros(N)
for node in range(N):
    inftime[node] = Propagation_SI.when_infected(labels,node,T)

plt.figure(2)
plt.hist(inftime)

#a = [nodes_stats[n]['cont last'] for n in range(1,N)]
#b = [nodes_stats[n]['inf time'] for n in range(1,N)]
#plt.figure(1)
#plt.hist(a)
#plt.figure(2)
#plt.hist(b)