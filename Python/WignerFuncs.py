#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 13:35:04 2019

@author: Craig Chisholm
"""

import numpy as np

def factorial(a):
    if isinstance(a,list):
        return [fac(int(x)) for x in a]
    else:
        return fac(int(a))

def fac(nn):
    if (nn in [0,1]):
        return 1
    else:
        return nn*fac(nn-1)

def triangle_coeff(a,b,c):
    '''Calculateing triange coefficients for Angular momenta.'''
    tri = factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/(factorial(a+b+c+1))
    return tri

def fung(t, j1,j2,j3,J1,J2,J3):
    '''Calculating the denominator in Racah Formula'''
    r = np.multiply(factorial(np.add(t,-j1-j2-j3)),np.multiply(factorial(np.add(t,-j1-J2-J3)),np.multiply(factorial(np.add(t,-J1-j2-J3)),np.multiply(factorial(np.add(t,-J1-J2-j3)),np.multiply(factorial(np.subtract(j1+j2+J1+J2,t)),np.multiply(factorial(np.subtract(j2+j3+J2+J3,t)),factorial(np.subtract(j3+j1+J3+J1,t))))))))
    return r

def Wigner6jcheck(j1,j2,j3,J1,J2,J3):
    '''Function to check Wigner-6j conditions'''
    Triad = np.array([[j1,j2,j3],[j1,J2,J3],[J1,j2,J3],[J1,J2,j3]])
    tfvec = np.zeros(8)
    
    #Check integer condition
    for kk in range(0,4):
        tfvec[kk] = np.floor(Triad[kk,:].sum())==(Triad[kk,:].sum())
    #Check triangle inequalities
    for kk in range(0,4):
        if ((np.absolute(Triad[kk,0]-Triad[kk,1])<=Triad[kk,2]) and (Triad[kk,2]<=(Triad[kk,0]+Triad[kk,1]))):
            tfvec[kk+4] = 1
    tfsum = tfvec.sum()
    check = tfsum==8
    return check

def Wigner6j(j1,j2,j3,J1,J2,J3):
    '''Wigner 6j-symbol calculator.
Originally written by Amita B Deb, Clarendon Lab. 2007 (Matlab).

Calculates { j1, j2 ,j3}  using Racah formula. See: Sobelman: Atomic Spectra and Radiative Transitions.
              J1  J2  J3'''
    check = Wigner6jcheck(j1,j2,j3,J1,J2,J3)
    if (check):
        #Finding Triangular coefficients
        tri1 = triangle_coeff(j1,j2,j3)
        tri2 = triangle_coeff(j1,J2,J3)
        tri3 = triangle_coeff(J1,j2,J3)
        tri4 = triangle_coeff(J1,J2,j3)
        #Finding the range of summation in Racah formula
        a = [j1+j2+j3,j1+J2+J3,J1+j2+J3,J1+J2+j3]
        rangei = max(a)
        k = [j1+j2+J1+J2,j2+j3+J2+J3,j3+j1+J3+J1]
        rangef = min(k)
        
        tt = np.arange(rangei,rangef+0.5)
        Wigner = np.multiply(np.sqrt(tri1*tri2*tri3*tri4),np.multiply(np.power(-1,tt),np.divide(factorial(np.add(tt,1)),fung(tt,j1,j2,j3,J1,J2,J3))))
    else:
        Wigner = 0
    return Wigner

def tfunction3j(tt,aa,bb,cc,alpha,beta):
    xoft3j = np.multiply(np.multiply(np.multiply(np.multiply(np.multiply(factorial(tt),factorial(np.add(tt,cc-bb+alpha))),factorial(np.add(tt,cc-aa-beta))),factorial(np.subtract(aa+bb-cc,tt))),factorial(np.subtract(aa-alpha,tt))),factorial(np.subtract(bb+beta,tt)))
    return xoft3j

def Wigner3j(aa,bb,cc,alpha,beta,gamma):
    '''Wigner 3j-symbol calculator'''
    #Find triangle coefficient using A. B. Deb (2007) script
    triang = triangle_coeff(aa, bb, cc)
    
    #Second term
    sqrtterm = np.sqrt(factorial(aa+alpha)*factorial(aa-alpha)*factorial(bb+beta)*factorial(bb-beta)*factorial(cc+gamma)*factorial(cc-gamma))
    
    #Finding the range of summation in the Racah formula
    tmin = [bb-cc-alpha,aa+beta-cc,0]
    tmax = [aa+bb-cc,aa-alpha,bb+beta]
    rangei = max(tmin)
    rangef = min(tmax)
    
    #Sum over the function of t
    tt = np.arange(rangei,rangef+0.5)
    
    sumscalar = np.divide(np.power(-1,tt),tfunction3j(tt,aa,bb,cc,alpha,beta))
    
    Wigner = ((-1)**(aa-bb-gamma))*np.sqrt(triang)*sqrtterm*(sumscalar.sum())
    
    return Wigner

def ClebschGord(jj1,jj2,JJ,mm1,mm2,MM):
    '''Function to calculate Clebsch-Gordan coeficients'''
    ClebschGord = ((-1)**(jj1-jj2+MM))*np.sqrt(2*JJ+1)*Wigner3j(jj1,jj2,JJ,mm1,mm2,-MM)
    return ClebschGord

def krondelt(ii,jj):
    '''Kronecker delta function'''
    if (ii==jj):
        krondelt = 1;
    else:
        krondelt = 0;
    return krondelt