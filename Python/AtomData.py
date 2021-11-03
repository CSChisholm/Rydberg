#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 10:52:03 2019

@author: Craig Chisholm
"""

#Atom specific parameters

import numpy as np
import SIunits as SIunits

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
        eneigval = (n5l0j1_2en-ionlim)/SIunits.atomenergy
    elif ((nn==5) and (ll==1) and (jj==1/2)):
        eneigval = (D1trans-ionlim)/SIunits.atomenergy
    elif ((nn==5) and (ll==1) and (jj==3/2)):
        eneigval = (D2trans-ionlim)/SIunits.atomenergy
    elif ((nn==4) and (ll==2) and(jj==5/2)):
        eneigval = (n4l2j5_2en-ionlim)/SIunits.atomenergy
    elif ((nn==4) and (ll==2) and (jj==3/2)):
        eneigval = (n4l2j3_2en-ionlim)/SIunits.atomenergy
    elif ((nn==6) and (ll==0)):
        eneigval = (n6l0j1_2en-ionlim)/SIunits.atomenergy
    elif ((nn==6) and (ll==1) and (jj==1/2)):
        eneigval = (n6l1j1_2en-ionlim)/SIunits.atomenergy
    elif ((nn==6) and (ll==1) and (jj==3/2)):
        eneigval = (n6l1j3_2en-ionlim)/SIunits.atomenergy
    else:
       eneigval = -Ry/(2*((nn-delta_nlj)**2))
       
    return ZZ, spin, eneigval, alpha_c, a_1, a_2, a_3, a_4, r_c, ionlim, delta_nlj

def Rb87decaypaths(nn,ll,jj):
    '''Radiative decay pathways for 87Rb'''
    loadparam = np.array([[5,1,0.5],[5,1,1.5],[6,1,1.5],[6,1,0.5],[5,0,0.5],[4,2,1.5],[6,0,0.5],[4,2,2.5]]) #Values to load based on allowed decay pathways
    
    D1trans = 1.56051536 #D1 line transition energy in eV
    D2trans = 1.58999085 #D2 line transition energy in eV
    n4l2j5_2en = 2.40116159
    n4l2j3_2en = 2.40121692
    n6l0j1_2en = 2.49759249
    n6l1j1_2en = 2.94203794
    n6l1j3_2en = 2.95165365
    n5l0j1_2en = 0 #5S_{1/2} energy in eV
    
    stateenvec = np.array([D1trans, D2trans, n6l1j3_2en, n6l1j1_2en, n5l0j1_2en,    n4l2j3_2en, n6l0j1_2en, n4l2j5_2en]) #Energy levels
    Rdecayratevec = np.zeros(8) #Create a vector to store decay rates in
    
    if ((ll==0) or ((ll==2) and (jj==1.5))):
        qvec1 = 0
        qvec4 = 3
    elif ((ll==1) and (jj==1.5)):
        qvec1 = 4
        qvec4 = 7
    elif ((ll==1) and (jj==0.5)):
        qvec1 = 4
        qvec4 = 6
    elif ((ll==2) and (jj==2.5)):
        qvec1 = 1
        qvec4 = 2
    
    return loadparam, stateenvec, Rdecayratevec, qvec1, qvec4

def Rb87blackbody(neff,temp,ll,jj):
    '''Function for calculating black body decay rate for 87Rb Rydberg states'''
    
    #Numbers from Beterov et a;. (2009)
    if (ll==0):
        AA = 0.134
        BB = 0.251
        CC = 2.567
        DD = 4.426
    elif ((ll==1) and (jj==0.5)):
        AA = 0.053
        BB = 0.128
        CC = 2.183
        DD = 3.989
    elif ((ll==1) and (jj==1.5)):
        AA = 0.046
        BB = 0.109
        CC = 2.085
        DD = 3.901
    elif ((ll==2) and (jj==1.5)):
        AA = 0.033
        BB = 0.084
        CC = 1.912
        DD = 3.716
    elif ((ll==2) and (jj==2.5)):
        AA = 0.032
        BB = 0.082
        CC = 1.898
        DD = 3.703
    
    term3 = AA/(neff**DD)
    term4 = 2.14e10/(np.exp(315780*BB/((neff**CC)*temp))-1)
    BBdecayrate = term3*term4
    
    return BBdecayrate