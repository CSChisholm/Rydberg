#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 17:19:48 2019

@author: Craig Chisholm
"""

#Functions for Rydberg library

import numpy as np
from SIunits import *
import AtomData as atoms
import MatlabHacks as mfuncs

def numerovfunc(atom,nn,ll,jj):
    '''Function for getting solution to the radial Schrodinger equation using Numerov algorithm'''
    ZZ, spin, eneigval, alpha_c, a_1, a_2, a_3, a_4, r_c, ionlim, delta_nlj = GetAtomParams(atom,nn,ll,jj)
    
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
    
    radialcharge = np.add(1,np.subtract(np.multiply(np.exp(np.multiply(rr,-a_1)),(ZZ-1)),np.multiply(rr,np.multiply(np.exp(np.multiply(rr,-a_2)),np.add(np.multiply(rr,a_4),a_3))))) #Effective nuclear charge
    
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

def Radiative_Lifetimes(atom,nn,ll,jj):
    '''Function for calculating Rydberg state lifetimes'''
    temp = 300 #Approximation of room temperature in K
    
    ZZ, spin, eneigval, alpha_c, a_1, a_2, a_3, a_4, r_c, ionlim, delta_nlj = GetAtomParams(atom,nn,ll,jj)
    if (atom=='87Rb'):
        loadparam, stateenvec, Rdecayratevec, qvec1, qvec4 = atoms.Rb87decaypaths(nn,ll,jj)
    
    for qq in range(qvec1,qvec4+1):
        #Get matrix element
        matrixelementSI = matrix_elements(atom,nn,ll,jj,loadparam[qq,0],loadparam[qq,1],loadparam[qq,2])[1]
        #Compute energies
        enlow = stateenvec[qq]-ionlim #Energy of lower state in eV
        eneigvaleV = eneigval*atomenergy #Energy of state |n,l,j> in eV
        
        #Perform calculation
        omega = np.absolute((eneigvaleV-enlow)/rpceV) #Angular frequency of transimission
        term1 = (omega**3)/(3*np.pi*vacpmtvty*rpcJ*(lightc**3)) #splitting terms to make reading easier
        term2 = (2*jj+1)/(2*loadparam[qq,2]+1)
        Rdecayratevec[qq] = term1*term2*(matrixelementSI**2) #Radiative decay rate
        
    Rdecayrate = Rdecayratevec.sum()
    
    #Account for black body radiation Beterov et al. (2009)
    
    neff = nn - delta_nlj #Effective principal quantum number
    
    if (atom=='87Rb'):
        BBdecayrate = atoms.Rb87blackbody(neff,temp,ll,jj)
    
    decayrate = Rdecayrate + BBdecayrate #Total decay rate
    lifetime = 1/decayrate
    
    return lifetime

def BlockadeShift(atom,nn,ll,jj,mj):
    '''Function for calculating Rydberg blockade shifts.'''
    
    pceV = 2*np.pi*rpceV #Planck's constant in eV
    
    if (atom=='87Rb'):
        ZZ, spin, eneigval, alpha_c, a_1, a_2, a_3, a_4, r_c, ionlim, delta_nlj = GetAtomParams(atom,nn,ll,jj)
    
    #Create ranges for theta and R (R in atomic units)
    theta = np.arange(0,np.pi/2+0.005,0.01)
    RRSI = np.multiply(np.arange(4,10.0005,0.001),1e-6)
    RR + np.divide(RRSI,bohrrad)
    
    #Set defect threshold
    defectmax = 100e9 #Maximum energy defect in Hz
    entol = defectmax*pceV/atomeenergy #Converted to atomic units
    ntol = 4 #Maximum change in n
    
    #Determine the single particle state energy
    singen = eneigval
    
    #Create set of allowable interacting state quantum numbers
    
    
    return RRSI, theta, blockadeshiftGHzmesh, C_6val

def matrix_elements(atom,nn1,ll1,jj1,nn2,ll2,jj2):
    '''Function for calculating dipole matrix elements based on results from numerovfunc()'''
    matrixelement = radiel(atom,nn1,ll1,jj1,nn2,ll2,jj2)
    matrixelementSI = matrixelement*bohrrad*eleccharge
    
    return matrixelement, matrixelementSI

def radiel(atom,nn1,ll1,jj1,nn2,ll2,jj2):
    '''function to evaluate radial matrix elements'''
    #first load radial wavefunctions for calculating the radial part of the matrix element
    try:
        rscale1, radial1 = np.loadtxt(filenamemake(atom,nn1,ll1,jj1), unpack=True)
    except FileNotFoundError:
        radial1, rscale1 = numerovfunc(atom,nn1,ll1,jj1)
    
    try:
        rscale2, radial2 = np.loadtxt(filenamemake(atom,nn2,ll2,jj2), unpack=True)
    except FileNotFoundError:
        radial2, rscale2 = numerovfunc(atom,nn2,ll2,jj2)
    
    #Calculate the radial matrix element
    if (nn1 >= nn2):
        Yk1 = radial1
        Yk2 = radial2
        rscale = rscale1
    else:
        Yk1 = radial2
        Yk2 = radial1
        rscale = rscale2
    
    #Resize smaller solution vector by attaching zeros to the end such that the two solution vectors are the same length
    if not(len(Yk1.tolist())==len(Yk2.tolist())):
        szero = np.zeros(len(Yk1.tolist())-len(Yk2.tolist())).tolist()
        Yk2conc = np.array(Yk2.tolist()+szero)
    else:
        Yk2conc = np.copy(Yk2)
    
    #Solve the matrix elements using method adapted from Zimmerman et al. (1979)
    
    deltar = np.subtract(rscale[1:(len(rscale.tolist())-1)],rscale[0:(len(rscale.tolist())-2)])
    numervec = np.multiply(Yk1[1:(len(Yk1.tolist())-1)],np.multiply(Yk2conc[1:(len(Yk2conc.tolist())-1)],np.multiply(np.power(rscale[1:(len(rscale.tolist())-1)],2),deltar)))
    matrixelement = numervec.sum()
    
    return matrixelement

def GetAtomParams(atom,nn,ll,jj):
    if (atom=='87Rb'):
        ZZ, spin, eneigval, alpha_c, a_1, a_2, a_3, a_4, r_c, ionlim, delta_nlj = atoms.Rb87Numbers(nn,ll,jj)
    return ZZ, spin, eneigval, alpha_c, a_1, a_2, a_3, a_4, r_c, ionlim, delta_nlj

def filenamemake(atom,nn,ll,jj):
    '''Function for constructing filenames'''
    jstring = str(jj).replace('.','_')
    filename = 'Wavefunctions/' + atom + str(nn) + 'n' + str(ll) + 'l' + jstring + 'j.txt'
    return filename