# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 18:50:24 2021

@author: Anatole Storck
"""

"""
This file return isotropic kick velocities for the DNSBs. While this file is not
used in the thesis, it could be helpful nonetheless.
"""

import math
import random
import numpy as np
import matplotlib.pyplot as plt
from pynverse import inversefunc

"""The method"""
#Sample size
X = np.arange(0, 3000, 1)
#Spherical
R = 100
Phi = []
Theta = []
#Creating distribution function in order to fix over densities at the poles
def cos(x):
    return math.cos(x)
for i in X:
    Phi.append(2*math.pi*random.random())
    Theta.append(inversefunc(cos, y_values=(2*random.random() - 1)))
#Converting to Cartesian
X = R*np.sin(Theta)*np.cos(Phi)
Y = R*np.sin(Theta)*np.sin(Phi)
Z = R*np.cos(Theta)
#Converting to Cylindrical
S = R*np.sin(Theta)
Z = R*np.cos(Theta)
Phi = Phi
#Plotting
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(X, Y, Z, c='k')
ax.set_box_aspect([1,1,1])
plt.show()
