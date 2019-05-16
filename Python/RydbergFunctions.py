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
        ZZ, spin, eneigval, alpha_c, a_1, a_2, a_3, a_4, r_c = Rb87Numbers(nn,ll,jj)
    return ZZ, spin, eneigval, alpha_c, a_1, a_2, a_3, a_4, r_c

def filenamemake(atom,nn,ll,jj):
    '''Function for constructing filenames'''
    jstring = str(jj).replace('.','_')
    filename = 'Wavefunctions/' + atom + str(nn) + 'n' + str(ll) + 'l' + jstring + 'j.txt'
    return filename

def Rb87Numbers(nn,ll,jj):
    #Numbers for 87Rb

    #Low level parameters
    ionlim = 4.1771270 #Rb groundstate energy in eV
    D1trans = 1.56051536 #D1 line transition energy in eV
    D2trans = 1.58999085 #D2 line transition energy in eV
    n4l2j5_2en = 2.40116159
    n4l2j3_2en = 2.40121692
    n6l0j1_2en = 2.49759249
    n6l1j1_2en = 2.94203794
    n6l1j3_2en = 2.95165365
    n5l0j1_2en = 0 #5S_{1/2} energy in eV
    #Other parameters
    ZZ = 37 #Nucler charge of Rubidium
    spin = 1/2 #Electron spin
    #Quantum defects from Li et al. (2003) and Han et al. (2006)
    if (ll==0):
        delta_0 = 3.1311804
        delta_2 = 0.1784
    elif ((ll==1) and (jj==1/2)):
        delta_0 = 2.65448849
        delta_2 = 0.2900
    elif ((ll==1) and (jj==3/2)):
        delta_0 = 2.6416737
        delta_2 = 0.2950
    elif ((ll==2) and (jj==3/2)):
        delta_0 = 1.34809171
        delta_2 = -0.60286
    elif ((ll==2) and (jj==5/2)):
        delta_0 = 1.34646572
        delta_2 = -0.59600
    elif ((ll==3) and (jj==5/2)):
        delta_0 = 0.0165192
        delta_2 = -0.085
    elif ((ll==3) and (jj==7/2)):
        delta_0 = 0.0165437
        delta_2 = -0.086
    else:
        delta_0 = 0
        delta_2 = 0
        
    #Approximation for quantum efect assuming n > 20
    delta_nlj = delta_0 + delta_2/((nn-delta_0)**2)
    #Model parameters from Marinescu et al. (1994)
    Ry = 1
    alpha_c = 9.0760
    if (ll==0):
        a_1 = 3.69628474
        a_2 = 1.64915255
        a_3 = -9.86069196
        a_4 = 0.19579987
        r_c = 1.66242117
    elif (ll==1):
        a_1 = 4.44088978
        a_2 = 1.92828831
        a_3 = -16.79597770
        a_4 = -0.81633314
        r_c = 1.50195124
    elif (11==2):
        a_1 = 3.78717363
        a_2 = 1.57027864
        a_3 = -11.65588970
        a_4 = 0.52942835
        r_c = 4.86851938
    else:
        a_1 = 2.39848933
        a_2 = 1.76810544
        a_3 = -12.07106780
        a_4 = 0.77256589
        r_c = 4.79831327
    
    #Get energy eigen value                         
    if ((nn==5) and (ll==0)):
        eneigval = (n5l0j1_2en-ionlim)/atomenergy
    elif ((nn==5) and (ll==1) and (jj==1/2)):
        eneigval = (D1trans-ionlim)/atomenergy
    elif ((nn==5) and (ll==1) and (jj==3/2)):
        eneigval = (D2trans-ionlim)/atomenergy
    elif ((nn==4) and (ll==2) and(jj==5/2)):
        eneigval = (n4l2j5_2en-ionlim)/atomenergy
    elif ((nn==4) and (ll==2) and (jj==3/2)):
        eneigval = (n4l2j3_2en-ionlim)/atomenergy
    elif ((nn==6) and (ll==0)):
        eneigval = (n6l0j1_2en-ionlim)/atomenergy
    elif ((nn==6) and (ll==1) and (jj==1/2)):
        eneigval = (n6l1j1_2en-ionlim)/atomenergy
    elif ((nn==6) and (ll==1) and (jj==3/2)):
        eneigval = (n6l1j3_2en-ionlim)/atomenergy
    else:
       eneigval = -Ry/(2*((nn-delta_nlj)**2))
       
    return ZZ, spin, eneigval, alpha_c, a_1, a_2, a_3, a_4, r_c