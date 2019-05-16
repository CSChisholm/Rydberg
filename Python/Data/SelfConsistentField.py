#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 20:03:21 2019

@author: Craig Chisholm
"""

import numpy as np

#Data from Callaway and Morgan

rr = np.array(np.arange(0,0.021,0.005).tolist() + np.arange(0.03,0.101,0.01).tolist() + np.arange(0.12,0.301,0.02).tolist() + np.arange(0.35,0.601,0.05).tolist() + np.arange(0.7,1.201,0.1).tolist() + np.arange(1.4,4.401,0.2).tolist() + [4.8, 5.2] + np.arange(5.6,8.401,0.4).tolist() + np.arange(9.2,13.201,0.8).tolist() + np.arange(14.8,26.01,1.6).tolist())

R_a = np.array([0, 0.02055, 0.03362, 0.04069, 0.043, 0.03745, 0.02353, 0.00576, -0.01281, -0.03018, -0.04517, -0.05711, -0.06573, -0.07309, -0.06882, -0.05565, -0.03648, -0.01396, 0.00967, 0.03269, 0.05386, 0.07231, 0.08754, 0.11040, 0.11257, 0.09803, 0.07170, 0.03828, 0.00167, -0.06995, -0.12867, -0.1691, -0.1907, -0.1950, -0.1848, -0.1333, -0.0580, 0.0252, 0.1064, 0.1807, 0.2464, 0.3032, 0.3510, 0.3902, 0.4215, 0.4456, 0.4632, 0.4751, 0.4819, 0.4844, 0.4831, 0.4713, 0.4503, 0.4233, 0.3927, 0.3604, 0.3276, 0.2953, 0.2643, 0.2350, 0.2078, 0.1601, 0.1213, 0.0907, 0.0669, 0.0489, 0.0354, 0.0180, 0.0089, 0.0044, 0.0021, 0.001, 0.0005, 0.0002, 0.0001])

#Tidy up for plotting

truncrr = rr[1:(len(R_a.tolist())-1)]
plotscale2 = np.sqrt(truncrr)
truncR_a = R_a[1:(len(R_a.tolist())-1)]
sqrtrR = np.divide(truncR_a,plotscale2)

np.savetxt('SCFcalc.txt', np.array([plotscale2, sqrtrR]).T,header='plotscale2 sqrtrR')