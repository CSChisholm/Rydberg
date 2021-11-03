#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 17:16:25 2019

@author: craig
"""

#Script to load various SI values which are needed to convert back from
#atomic rydberg units

import scipy.constants as constants

finestrucconst = constants.alpha #Fine structure constant
atomenergy = 27.2 #Atomic energy units in eV
Ry = 1
rpcJ = constants.hbar #Reduced Planck constant in J s
lightc = constants.c #Speed of light in m/s
vacpmtvty = constants.epsilon_0 #Vacuum permitivity in F/m
espin = 0.5 #electron spin
boltzconst = constants.k #Boltzmann's constant in J/K
bohrrad = constants.physical_constants['Bohr radius'][0] #Bohr radius
eleccharge = constants.e #electron charge
rpceV = constants.hbar/constants.e #Reduced Planck constant in eV s