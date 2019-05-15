#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 17:12:32 2019

@author: Craig Chisholm
"""

#Script for plotting radial wavefunctions

import matplotlib.pyplot as plt
import numpy as np
import RydbergFunctions as Rydberg
from SIunits import *

plt.close("all")

#Input information
atom = '87Rb'
nn = 5
ll = 0
jj = 0.5

#Get some constants (by adding more files with different data this code can be made to work for different atoms)
if (atom=='87Rb'):
    from Rb87Numbers import *

#Calculate wave function
normY_sol, rr = Rydberg.numerovfunc(atom,nn,ll,jj)

#Rescale for plotting
plotscale = np.sqrt(rr)
probamp = np.power(np.multiply(normY_sol,plotscale),2)

tString = atom + ' radial wavefunction n = ' + str(nn) + ', l = ' + str(ll) + ', j = ' + str(jj)

plt.figure()
plt.plot(plotscale,normY_sol)
if (nn>20):
    plt.xlim([np.sqrt(alpha_c**(1/3)),plotscale.max()])
plt.xlabel('$(r/a_0)^{1/2}$')
plt.ylabel('$r^{1/2}R(r)\, (a_0^{-1})$')
plt.title(tString)
plt.show()

plt.figure()
plt.plot(plotscale,probamp)
if (nn>20):
    plt.xlim([alpha_c**(1/3),plotscale.max()])
plt.xlabel('$(r/a_0)^{1/2}$')
plt.ylabel('$|r^{1/2}R(r)|.^2\, (a_0^{-1})$')
plt.title(tString)
plt.show()