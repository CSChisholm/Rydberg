#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 13:35:07 2019

@author: Craig Chisholm
"""

import numpy as np

#Reverse engineered versions of inbuilt Matlab functions

def repmat(xx,m,n):
    outArray = np.tile(xx,(m,n))
    return outArray

def rectpulse(XX,Nsamp):
    '''Function to emulate Matlab's rectpulse'''
    UU = repmat(XX,Nsamp,1)
    if not(ndims(XX)==1):
        sX = np.shape(XX)
    else:
        sX = (1,len(XX.tolist()))
    YY = np.zeros(np.shape(UU))
    for kk in range(0,sX[0]):
        YY[(kk*Nsamp):((kk+1)*Nsamp),:] = UU[kk:np.shape(UU)[0]:sX[0],:]
    return YY

def ndims(xx):
    return len(np.shape(xx))

def combvec(*kargs):
    '''Function for mimicking Matlab's combvec.'''
    numargs = len(kargs)
    Mvec = np.zeros(numargs)
    Nvec = np.zeros(numargs)
    kk = 0
    for arg in kargs:
        if (ndims(arg)==1):
            Mvec[kk] = 1
            Mvec[kk] = len(arg.tolist())
        else:
            sA = np.shape(arg)
            Mvec[kk] = sA[0]
            Nvec[kk] = sA[1]
        kk+=1
    Nprod = int(Nvec.prod())
    outmatrix = np.zeros((int(Mvec.sum()),Nprod))
    Mstart = 0
    kk = 1
    for arg in kargs:
        Mend = int(Mvec[0:kk].sum())
        val1 = rectpulse(arg.T,int(Nprod/Nvec[(kk-1):numargs].prod())).T
        outmatrix[Mstart:Mend,:] = repmat(val1,1,int(Nprod/Nvec[0:kk].prod()))
        Mstart = Mend
        kk+=1
    return outmatrix