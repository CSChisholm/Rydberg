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
import WignerFuncs as Wigner

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
    RR = np.divide(RRSI,bohrrad)
    
    #Set defect threshold
    defectmax = 100e9 #Maximum energy defect in Hz
    entol = defectmax*pceV/atomenergy #Converted to atomic units
    ntol = 4 #Maximum change in n
    
    #Determine the single particle state energy
    singen = eneigval
    
    #Create set of allowable interacting state quantum numbers
    nvec = mfuncs.rectpulse(np.arange(nn-ntol,nn+ntol+1),ll+2).T.ravel()
    lcands = np.arange(ll-1,ll+2,2)
    lcands = lcands[lcands>=0]
    if (lcands[0]==0):
        smalllvec = np.array([0] + mfuncs.rectpulse(lcands[1:len(lcands.tolist())],2).T.ravel().tolist())
        jmaker = np.array([0.5] + mfuncs.repmat(np.array([-0.5,0.5]),1,len(lcands[1:len(lcands.tolist())])).ravel().tolist())
    else:
        smalllvec = mfuncs.rectpulse(lcands,2).ravel()
        jmaker = mfuncs.repmat(np.array([-0.5,0.5]),1,len(lcands.tolist())).ravel()
    smalljvec = np.add(smalllvec,jmaker)
    lvec = mfuncs.repmat(smalllvec,1,int(len(nvec.tolist())/len(smalllvec.tolist()))).ravel()
    jvec = mfuncs.repmat(smalljvec,1,int(len(nvec.tolist())/len(smalljvec.tolist()))).ravel()
    truncspace1 = np.array([nvec.tolist(),lvec.tolist(),jvec.tolist()])
    pairs1 = mfuncs.combvec(truncspace1,truncspace1)
    Sp1 = np.shape(pairs1)
    
    #Check which of these states satisfy the infinite separation energy defect condition and selection rules for j
    pindex = np.zeros(Sp1[1])
    defects = np.zeros(Sp1[1])
    
    for kk in range(0,Sp1[1]):
        energy1 = envalfunc(atom,pairs1[0,kk],pairs1[1,kk],pairs1[2,kk])
        energy2 = envalfunc(atom,pairs1[3,kk],pairs1[4,kk],pairs1[5,kk])
        defects[kk] = energy1 + energy2 - 2*singen #Energy defect
        j1check = (pairs1[2,kk]==np.arange(jj-1,jj+1.5)).sum()
        j2check = (pairs1[5,kk]==np.arange(jj-1,jj+1.5)).sum()
        if ((np.absolute(defects[kk])<=entol) and (j1check==1) and (j2check==1)):
            pindex[kk] = 1
    pindex = np.where(pindex==1)[0]
    pairs2L = []
    for row in pairs1:
        pairs2L.append(row[pindex].tolist())
    pairs2L.append(defects[pindex].tolist())
    pairs2 = np.array(pairs2L)
    
    matel2part = np.zeros((len(theta.tolist()),np.shape(pairs2)[1])) #Vector to store matrix elements in
    
    #Call function to calculate matrix elements
    for kk in range(0,np.shape(pairs2)[1]):
        matel2part[:,kk] = matrixel2p(atom,nn,ll,jj,mj,pairs2[0,kk],pairs2[1,kk],pairs2[2,kk],pairs2[3,kk],pairs2[4,kk],pairs2[5,kk],pairs2[6,kk],theta)
    
    #Compute the blockade shift in atomic units
    summation = np.multiply(np.sum(matel2part,axis=1),-1)
    kindex = np.where(np.absolute(summation)==np.absolute(summation).max())[0][0] #Find the index of maximum blockade shift
    
    blockadeshiftau = np.divide(summation[kindex],np.power(RR,6))
    
    RRmesh, summationmesh = np.meshgrid(RR,summation)
    blockadeshiftaumesh = np.divide(summationmesh,np.power(RRmesh,6))
    
    #Convert units to GHz
    encon = (atomenergy/pceV)*1e-9 #Factor from atomic energy units to GHz
    blockadeshiftGHz = np.multiply(blockadeshiftau,encon)
    blockadeshiftGHzmesh = np.multiply(blockadeshiftaumesh,encon)
    
    #C6
    C_6val = float(blockadeshiftGHz[len(blockadeshiftGHz.tolist())-1]*(RRSI[len(RRSI.tolist())-1]**6)) #GHz/m^6
    
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

def matrixel2p(atom,nni,lli,jji,mji,nn1,ll1,jj1,nn2,ll2,jj2,defect,theta):
    '''Function for calculating the two particle matrix elements required for second order perturbation theory evaluation of the blockade shift'''
    
    sint = np.sin(theta) #For efficiency
    cost = np.cos(theta)
    sin2t = np.power(sint,2)
    cos2t = np.power(cost,2)
    
    #first calculate angular factors as vectors [minus, zero, plus] for each single particle interaction state
    
    angvec1 = angcalc(lli,jji,mji,ll1,jj1)
    angvec2 = angcalc(lli,jji,mji,ll2,jj2)
    
    #Radial matrix elements
    matelem1i = radiel(atom,nni,lli,jji,nn1,ll1,jj1)
    matelem2i = radiel(atom,nni,lli,jji,nn2,ll2,jj2)
    
    #Calculate terms in two particle matrix element using Wigner-Eckhart theorem (Pritchard, 2012), (Reinhard et al., 2007).
    line1 = np.add(np.add(np.multiply(angvec1[2],angvec2[2]),np.multiply(angvec1[0],angvec2[2])),np.multiply(np.multiply(angvec1[1],angvec2[1]),np.subtract(1,np.multiply(cos2t,3))))
    line2 = np.multiply(np.add(np.add(np.multiply(angvec1[2],angvec2[2]),np.multiply(angvec1[2],angvec2[0])),np.add(np.multiply(angvec1[0],angvec2[2]),np.multiply(angvec1[0],angvec2[0]))),np.multiply(sin2t,-3/2))
    line3 = np.multiply(np.add(np.add(np.multiply(angvec1[2],angvec2[1]),np.multiply(angvec1[0],angvec2[1])),np.add(np.multiply(angvec1[1],angvec2[2]),np.multiply(angvec1[1],angvec2[0]))),np.multiply(np.multiply(sint,cost),-3/np.sqrt(2)))
    
    matrixelement = np.multiply(np.multiply(matelem1i,matelem2i),np.add(line1,np.add(line2,line3)))
    
    mel2p = np.divide(np.power(matrixelement,2),defect)
    
    return mel2p

def angcalc(lli,jji,mji,ll1,jj1):
    angvec1 = np.zeros(3)
    for qq in range(-1,2):
        mj1 = mji - qq
        tf1 = (mj1==np.arange(-jj1,jj1+0.5)).sum()
        if (tf1==1):
            ang11 = (-1)**(jji-mji+espin+jj1+1)
            ang12 = np.sqrt((2*jji+1)*(2*jj1+1)*(2*lli+1)*(2*ll1+1))
            ang13 = Wigner.Wigner6j(jji,1,jj1,ll1,espin,lli)
            ang14 = Wigner.Wigner3j(jji,1,jj1,-mji,qq,mj1)
            ang15 = Wigner.Wigner3j(lli,1,ll1,0,0,0)
            angvec1[qq+1] = ang11*ang12*ang13*ang14*ang15 #Shift for Python indexing
    return angvec1

def envalfunc(atom,nn,ll,jj):
    return GetAtomParams(atom,nn,ll,jj)[2]

def GetAtomParams(atom,nn,ll,jj):
    if (atom=='87Rb'):
        ZZ, spin, eneigval, alpha_c, a_1, a_2, a_3, a_4, r_c, ionlim, delta_nlj = atoms.Rb87Numbers(nn,ll,jj)
    return ZZ, spin, eneigval, alpha_c, a_1, a_2, a_3, a_4, r_c, ionlim, delta_nlj

def filenamemake(atom,nn,ll,jj):
    '''Function for constructing filenames'''
    jstring = str(jj).replace('.','_')
    filename = 'Wavefunctions/' + atom + str(nn) + 'n' + str(ll) + 'l' + jstring + 'j.txt'
    return filename