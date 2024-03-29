#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 16:27:43 2019

@author: Craig Chisholm
"""

#This script is just the blockade shift part of Demonstrations.py

import matplotlib.pyplot as plt
import numpy as np
import RydbergFunctions as Rydberg
from fractions import Fraction
from matplotlib import cm

plt.close("all")

#Input information
atom = '87Rb'
nn = 50
ll = 2
jj = 2.5
mj = 2.5

RRSI, theta, blockadeshiftGHzmesh, C_6val = Rydberg.BlockadeShift(atom,nn,ll,jj,mj)

stringlookup = 'SPD'
ratJ = Fraction(jj)
numer = ratJ.numerator
denom = ratJ.denominator
ratmJ = Fraction(mj)
numer1 = ratmJ.numerator
denom1 = ratmJ.denominator

plt.figure()
plt.pcolormesh(np.multiply(RRSI,1e6),theta,np.multiply(blockadeshiftGHzmesh,1e3),cmap=cm.inferno,shading='auto')
plt.xlabel('$R\, (\mu\mathrm{m})$')
plt.ylabel(r'$\theta\, (\mathrm{radians})$')
plt.title(f'Calculated Rydberg blockade shift for\n $|{nn}{stringlookup[ll]}_{{{numer}/{denom}}}, m_j = {numer1}/{denom1}\\rangle$ state of {atom}')
plt.colorbar()
plt.show()