#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 19:56:18 2019

@author: Craig Chisholm
"""

import matplotlib.pyplot as plt
import numpy as np
import RydbergFunctions as Rydberg
from SIunits import *
from fractions import Fraction

plt.close("all")

#Input information
atom = '87Rb'
nn = 50
ll = 0
jj = 0.5
mj = 0.5

#Get some constants (by adding more files with different data this code can be made to work for different atoms)
alpha_c = Rydberg.GetAtomParams(atom,nn,ll,jj)[3]

#Calculate wave function
normY_sol, rr = Rydberg.numerovfunc(atom,nn,ll,jj)

#Rescale for plotting
plotscale = np.sqrt(rr)
probamp = np.power(np.multiply(normY_sol,plotscale),2)

tString = atom + ' radial wavefunction |n,l,j> = |' + str(nn) + ',' + str(ll) + ',' + str(jj) + '>'

plt.figure()
plt.plot(plotscale,normY_sol)
if (nn>20):
    plt.xlim([np.sqrt(alpha_c**(1/3)),plotscale.max()])
plt.xlabel('$(r/a_0)^{1/2}$')
plt.ylabel('$r^{1/2}R(r)\, (a_0^{-1})$')
plt.title(tString)
plt.show()

plt.figure()
plt.plot(rr,probamp)
if (nn>20):
    plt.xlim([alpha_c**(1/3),rr.max()])
plt.xlabel('$r/a_0$')
plt.ylabel('$|rR(r)|.^2\, (a_0^{-1})$')
plt.title(tString)
plt.show()

#Compare groundstate to self-consistent field method
plotscale2, sqrtrR = np.loadtxt('Data/SCFcalc.txt', unpack=True)
normY_solg, rrg = Rydberg.numerovfunc('87Rb',5,0,0.5)

plt.figure()
plt.scatter(plotscale2,sqrtrR,label='Callaway')
plt.plot(np.sqrt(rrg),np.multiply(normY_solg,-1),label='Numerov')
plt.xlabel('$(r/a_0)^{1/2}$')
plt.ylabel('$r^{1/2}R(r)\, (a_0^{-1})$')
plt.legend()
plt.title('Rubidium Groundstate Radial Wavefunction')

#Lifetimes
sexpt, sexpterr = np.loadtxt('Data/LifetimesS.txt', unpack=True)
pexpt, pexpterr = np.loadtxt('Data/LifetimesP.txt', unpack=True)
dexpt, dexpterr = np.loadtxt('Data/LifetimesD.txt', unpack=True)

#S_1/2
nns = np.arange(28,46,1)
calcs = np.zeros(np.shape(nns))
itr = 0
while (itr<len(nns.tolist())):
    calcs[itr] = Rydberg.Radiative_Lifetimes(atom,nns[itr],0,0.5)*1e6
    itr+=1

#P_3/2
nnp = np.arange(34,45,1)
calcp = np.zeros(np.shape(nnp))
itr = 0
while (itr<len(nnp.tolist())):
    calcp[itr] = Rydberg.Radiative_Lifetimes(atom,nnp[itr],1,1.5)*1e6
    itr+=1

#D_5/2
nnd = np.arange(29,45,1)
calcd = np.zeros(np.shape(nnd))
itr = 0
while (itr<len(nnd.tolist())):
    calcd[itr] = Rydberg.Radiative_Lifetimes(atom,nnd[itr],2,2.5)*1e6
    itr+=1

plt.figure()
plt.scatter(nns,calcs,marker='x',label='Calculated')
plt.errorbar(nns,sexpt,yerr=sexpterr,fmt='o',label='Experimental')
plt.title('Plot of calculated and experimental lifetimes for $nS_{1/2}$ states')
plt.ylabel(r'$\tau\, (\mu\mathrm{s})$')
plt.xlabel('$n$')
plt.legend()
plt.show()

plt.figure()
plt.scatter(nnp,calcp,marker='x',label='Calculated')
plt.errorbar(nnp,pexpt,yerr=pexpterr,fmt='o',label='Experimental')
plt.title('Plot of calculated and experimental lifetimes for $nP_{3/2}$ states')
plt.ylabel(r'$\tau\, (\mu\mathrm{s})$')
plt.xlabel('$n$')
plt.legend()
plt.show()

plt.figure()
plt.scatter(nnd,calcd,marker='x',label='Calculated')
plt.errorbar(nnd,dexpt,yerr=dexpterr,fmt='o',label='Experimental')
plt.title('Plot of calculated and experimental lifetimes for $nD_{5/2}$ states')
plt.ylabel(r'$\tau\, (\mu\mathrm{s})$')
plt.xlabel('$n$')
plt.legend()
plt.show()

#Blockade shift

RRSI, theta, blockadeshiftGHzmesh, C_6val = Rydberg.BlockadeShift(atom,nn,ll,jj,mj)

stringlookup = 'SPD'
ratJ = Fraction(jj)
numer = ratJ.numerator
denom = ratJ.denominator
ratmJ = Fraction(mj)
numer1 = ratmJ.numerator
denom1 = ratmJ.denominator

plt.figure()
plt.pcolormesh(np.multiply(RRSI,1e6),theta,np.multiply(blockadeshiftGHzmesh,1e3))
plt.xlabel('$R\, (\mu\mathrm{m})$')
plt.ylabel(r'$\theta\, (\mathrm{radians})$')
plt.title('Calculated Rydberg blockade shift for |' + str(nn) + stringlookup[ll] + '_{' + str(numer) + '/' + str(denom) + '}, m_j = ' + str(numer1) + '/' + str(denom1) + '> state of Rb')
plt.colorbar()
plt.show()