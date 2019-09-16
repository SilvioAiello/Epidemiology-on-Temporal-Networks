# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 16:54:54 2019

@author: Silvio
"""
#%% GENERATING FILES
import scipy.io

mat = scipy.io.loadmat('C:/Users/Silvio/Desktop/parametri_stimati_EMID.mat')

alphaDAR1 = mat['alphaDAR1']
chiDAR1 = mat['marginalProbabilityDAR1']
phi0TGRG = mat['phi0TGRG']
phi1TGRG = mat['phi1TGRG']
sigmaTGRG = mat['sigmaTGRG']

import pickle

with open('C:/Users/Silvio/Documents/GitHub/Epidemiology-on-Temporal-Networks/Empiric_Data/alphaDAR1.pkl', 'wb') as handle: #wb = write binary
    pickle.dump(alphaDAR1,handle)
with open('C:/Users/Silvio/Documents/GitHub/Epidemiology-on-Temporal-Networks/Empiric_Data/chiDAR1.pkl', 'wb') as handle: #wb = write binary
    pickle.dump(chiDAR1,handle)
with open('C:/Users/Silvio/Documents/GitHub/Epidemiology-on-Temporal-Networks/Empiric_Data/phi0TGRG.pkl', 'wb') as handle: #wb = write binary
    pickle.dump(phi0TGRG,handle)
with open('C:/Users/Silvio/Documents/GitHub/Epidemiology-on-Temporal-Networks/Empiric_Data/phi1TGRG.pkl', 'wb') as handle: #wb = write binary
    pickle.dump(phi1TGRG,handle)
with open('C:/Users/Silvio/Documents/GitHub/Epidemiology-on-Temporal-Networks/Empiric_Data/sigmaTGRG.pkl', 'wb') as handle: #wb = write binary
    pickle.dump(sigmaTGRG,handle)