#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 16:07:01 2019

@author: Craig Chisholm
"""

import numpy as np
import RydbergFunctions as Rydberg

atom = '87Rb'

nnc6 = np.arange(35,83)

c6S1_2 = np.zeros(np.shape(nnc6))
c6D3_2 = np.zeros(np.shape(nnc6))
c6D5_2 = np.zeros(np.shape(nnc6))

for kk in range(0,len(nnc6.tolist())):
    c6S1_2[kk] = Rydberg.BlockadeShift(atom,nnc6[kk],0,0.5,0.5)[3]
    c6D3_2[kk] = Rydberg.BlockadeShift(atom,nnc6[kk],2,1.5,1.5)[3]
    c6D5_2[kk] = Rydberg.BlockadeShift(atom,nnc6[kk],2,2.5,2.5)[3]

np.savetxt('Data/'+atom+'C6dataS1_2.txt', np.array([nnc6, c6S1_2]).T,header='nnc6 c6S1_2')
np.savetxt('Data/'+atom+'C6dataD3_2.txt', np.array([nnc6, c6D3_2]).T,header='nnc6 c6D3_2')
np.savetxt('Data/'+atom+'C6datad5_2.txt', np.array([nnc6, c6D5_2]).T,header='nnc6 c6D5_2')