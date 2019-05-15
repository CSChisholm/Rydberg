#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 17:19:48 2019

@author: Craig Chisholm
"""

#Functions for Rydberg library

import numpy as np
from SIunits import *

def numerovfunc(atom,nn,ll,jj):
    '''Function for getting solution to the radial Schrodinger equation using Numerov algorithm'''
    ZZ, spin, eneigval, alpha_c, a_1, a_2, a_3, a_4, r_c = GetAtomParams(atom,nn,ll,jj)
    
    hh = 0.01 #Choose this to be small so that O(h**6) may be ignored
    h2 = hh**2 #For efficiecy
    r_i = hh/2 #Define starting point
    r_o = 2*nn*(nn+15) #Define end point
    x_i = np.log(r_i) #Change of variables
    x_o = np.log(r_o)
    xx = np.arange(x_i,x_o+hh,hh)
    rr = np.exp(xx) #For efficiency
    
    #Set up Schrodinger equation from Gallagher (2005), note Pritchard (2012)
    
    LdotS = (jj*(jj+1)-ll*(ll+1)-spin*(spin+1))/2 #Spin-orbit coupling
    
    spinorbitpotential = np.divide((finestrucconst**2)*LdotS,np.multiply(np.power(rr,3),2)) #Fine structure splitting
    
    radialcharge = np.add(1,np.add(np.multiply(np.exp(np.multiply(rr,-a_1)),(ZZ-1)),np.multiply(rr,np.multiply(np.exp(np.multiply(rr,-a_2)),np.add(np.multiply(rr,a_4),a_3))))) #Effective nuclear charge
    
    coulombpotential = np.multiply(np.add(np.divide(radialcharge,rr),np.multiply(np.divide(alpha_c,np.multiply(np.power(rr,4),2)),np.subtract(1,np.exp(np.multiply(np.power(np.divide(rr,r_c),6),-1))))),-1) #Coulomb potential
    
    totalpotential = np.add(spinorbitpotential,coulombpotential) #Total potential
    
    cenfugterm = (ll + 1/2)**2 #Centifugal term for Schrodinger equation
    
    #Apply Numerov method
    
    G_x = np.add(np.multiply(np.multiply(np.exp(np.multiply(xx,2)),2),np.subtract(totalpotential,eneigval)),cenfugterm) #Coefficient in differential equation
    
    T_x = np.multiply(G_x,h2/12) #For effiency
    
    Ysoln_vec = np.zeros(np.shape(xx)) #Generate a place holder array for solutions
    Ysoln_vec[len(Ysoln_vec.tolist())-1] = -1e-10 #The function must go to zero at infinity
    Ysoln_vec[len(Ysoln_vec.tolist())-2] = -2e-10
    
    #To perform the iteration with convenient indexing the vectors xx, T_x and Ysoln_vec are flipped
    
    fYsoln = np.flip(Ysoln_vec,axis=0)
    fT_x = np.flip(T_x,axis=0)
    
    itr = 2
    while (itr<len(xx.tolist())):
        fYsoln[itr] = ((2 + 10*fT_x[itr-1])*fYsoln[itr-1] - (1 - fT_x[itr-2])*fYsoln[itr-2])/(1-fT_x[itr])
        itr+=1
        
    Ysoln_vec = np.flip(fYsoln,axis=0) #Return solutions to proper ordering
    
    #Normalise (method adapted from Gallagher)
    normconstvec = np.zeros(len(rr.tolist())-1)
    itr = 1
    while (itr<len(rr.tolist())):
        deltar = rr[itr]-rr[itr-1]
        normconstvec[itr-1] = (Ysoln_vec[itr]**2)*rr[itr]*deltar
        itr+=1
        
    normconst = np.sqrt(normconstvec.sum())
    normY_sol = np.divide(Ysoln_vec,normconst)
    
    filename = filenamemake(atom,nn,ll,jj)
    
    np.savetxt(filename, np.array([rr, normY_sol]).T,header='rr normY_sol')
    
    return normY_sol, rr

def GetAtomParams(atom,nn,ll,jj):
    if (atom=='87Rb'):
        from Rb87Numbers import ZZ, spin, eneigval, alpha_c, a_1, a_2, a_3, a_4, r_c
    return ZZ, spin, eneigval, alpha_c, a_1, a_2, a_3, a_4, r_c

def filenamemake(atom,nn,ll,jj):
    '''Function for constructing filenames'''
    jstring = str(jj).replace('.','_')
    filename = 'Wavefunctions/' + atom + str(nn) + 'n' + str(ll) + 'l' + jstring + 'j.txt'
    return filename