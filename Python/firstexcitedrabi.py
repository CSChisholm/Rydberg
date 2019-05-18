#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 16:36:53 2019

@author: Craig Chisholm
"""

#Matrix elements with 5P_{3/2} and laser intensities with a user defined Rabi frequency

import numpy as np
import RydbergFunctions as Rydberg
import WignerFuncs as Wigner
from SIunits import *

def firstexcitedrabiRb(nn,ll,jj):
    
    #Define key constants
    atom = '87Rb'
    spin = 0.5 #electron spin
    ng = 5 #Principle qunatum number
    lp = 1 #Azimuthal quantum number for p state
    J3_2 = 1.5 #Total angular momentum of nP_{3/2} state
    
    #Define parameters relating to the decomposition of the hyperfine basis
    mjvec = np.array([0.5, 1.5])
    ClebschGordg = np.multiply(1/np.sqrt(2),np.array([1, -1]))
    matpart = np.array([0, 0])
    
    #Choose Rabi frequency
    RabifreqM = 1 #MHz
    Rabifreq = 2*np.pi*RabifreqM*1e6
    
    #Choose polarisation of light
    qq = 1
    
    beamwaist = 30 #um
    
    matrixelement = Rydberg.radiel(atom,ng,lp,J3_2,nn,ll,jj)
    
    for kk in range(0,2):   
        #Calculate the <Jm_J|er|J'm_J'> matrix element
        m_jprime = mjvec[kk]-qq
        angfac1 = (-1)**(J3_2-mjvec[kk]+spin+jj+1)
        angfac2 = np.sqrt((2*J3_2+1)*(2*jj+1)*(2*lp+1)*(2*ll+1))
        angfac3 = Wigner.Wigner6j(J3_2,1,jj,ll,spin,lp)
        angfac4 = Wigner.Wigner3j(J3_2,1,jj,-mjvec[kk],qq,m_jprime)
        angfac5 = Wigner.Wigner3j(lp,1,ll,0,0,0)
        angular = angfac1*angfac2*angfac3*angfac4*angfac5
        matpart[kk] = matrixelement*angular*ClebschGordg[kk]
    
    matrixfac = matpart.sum()
    
    #Calculate required intensity
    term01 = 1/(matrixfac**2)
    term02 = vacpmtvty*(rpcJ**2)*lightc/(2*(eleccharge**2)*(bohrrad**2))
    intensity = term01*term02*(Rabifreq**2) #Intensity in W/m^2
    
    power = (np.pi/2)*intensity*((beamwaist*1e-6)**2)
    return power