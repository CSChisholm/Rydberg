#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 13:35:07 2019

@author: Craig Chisholm
"""

import numpy as np

#Reverse engineered versions of inbuilt Matlab functions

def repmat(xx,m,n):
    tileArray = np.tile(xx,(m,n))
    if (m==1):
        outArray = tileArray[0]
    else:
        outArray = tileArray
    return outArray

def rectpulse(XX,Nsamp):
    '''Function to emulate Matlab's rectpulse'''
    UU = repmat(XX,1,Nsamp)
    sX = np.shape(XX)
    YY = np.zeros(np.shape(UU))
    for kk in range(0,sX[0]):
        YY[(kk*Nsamp):((kk+1)*Nsamp)] = UU[kk:(np.shape(UU)[0]-1):sX[0]]
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
            Mvec[kk] = np.shape(arg)[0]
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
        Mend = int(Mvec[0:kk].sum()) - 1
        val1 = rectpulse(arg.T,int(Nprod/Nvec[kk:(numargs-1)].prod())).T
        outmatrix[Mstart:Mend,:] = np.tile(val1,(int(Nprod/Nvec[0:kk].prod()),1)).T
        Mstart = Mend + 1
        kk+=1
    return outmatrix