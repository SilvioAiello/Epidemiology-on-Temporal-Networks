# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 16:54:54 2019

@author: Silvio
"""
#%% GENERATING FILES
import scipy.io

mat = scipy.io.loadmat('parametri_stimati_EMID.mat')

alphaDAR1 = mat['alphaDAR1']
chiDAR1 = mat['marginalProbabilityDAR1']
phi0TGRG = mat['phi0TGRG']
phi1TGRG = mat['phi1TGRG']
sigmaTGRG = mat['sigmaTGRG']

import pickle

with open('alphaDAR1.pkl', 'wb') as handle: #wb = write binary
    pickle.dump(alphaDAR1,handle)
with open('chiDAR1.pkl', 'wb') as handle: #wb = write binary
    pickle.dump(chiDAR1,handle)
with open('phi0TGRG.pkl', 'wb') as handle: #wb = write binary
    pickle.dump(phi0TGRG,handle)
with open('phi1TGRG.pkl', 'wb') as handle: #wb = write binary
    pickle.dump(phi1TGRG,handle)
with open('sigmaTGRG.pkl', 'wb') as handle: #wb = write binary
    pickle.dump(sigmaTGRG,handle)

import matplotlib.pyplot as plt
fig_count = 1
#ALPHAS PLOT
plt.figure(fig_count)
fig_count +=1
plt.hist(alphaDAR1.reshape(98*98,1))
plt.xlabel(r"$\alpha_{ij} campionati$")
plt.ylabel("Occorrenze")
plt.title(r"$\alpha$ empirici")
plt.savefig("istogrammaALPHA.pdf")
#CHI PLOT
plt.figure(fig_count)
fig_count +=1
plt.hist(chiDAR1.reshape(98*98,1))
plt.xlabel(r"$\chi_{ij} campionati$")
plt.ylabel("Occorrenze")
plt.title(r"$\chi$ empirici")
plt.savefig("istogrammaCHI.pdf")
#PHI0 PLOT
plt.figure(fig_count)
fig_count +=1
plt.hist(phi0TGRG)
plt.xlabel(r"$\phi_{0,i} campionati$")
plt.ylabel("Occorrenze")
plt.title(r"$\phi_0$ empirici")
plt.savefig("istogrammaPHI0.pdf")
#PHI1 PLOT
plt.figure(fig_count)
fig_count +=1
plt.hist(phi1TGRG)
plt.xlabel(r"$\phi_{1,i} campionati$")
plt.ylabel("Occorrenze")
plt.title(r"$\phi_1$ empirici")
plt.savefig("istogrammaPHI1.pdf")
#SIGMA PLOT
plt.figure(fig_count)
fig_count +=1
plt.hist(sigmaTGRG)
plt.xlabel(r"$\sigma_{0,i} campionati$")
plt.ylabel("Occorrenze")
plt.title(r"$\sigma$ empirici")
plt.savefig("istogrammaSIGMA.pdf")