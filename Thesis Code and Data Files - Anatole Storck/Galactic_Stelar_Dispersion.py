# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 14:25:11 2021

@author: Anatole Storck
"""

"""
This file creates a random distribution of distances in the galactic disk
that are consistent with the stellar distribution.
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from pynverse import inversefunc
import random

#Formating
plt.rcParams['axes.linewidth'] = 2.15
plt.rcParams["font.family"] = 'Courier New'
plt.rcParams.update({'font.size': 15})

"""Constants"""
#From McMillan
sigma_0 = 887
#Normalizing factor is taken after doing C(R) and taking the limit as R --> inf
norm_Factor = 40628.52132884386

"""The method"""
#For the surface density
Y_1 = []
#For N(R<r)/N
Y_2 = []
#Stellar surface density and its cumulative distribution function
def Sigma(R):
    return sigma_0*math.exp(-1*(R/2.7))
def C(R):
    if R > 0:
        return (2*math.pi*sigma_0*(math.exp(-1*(R/2.7))*(-2.7*R - 7.29) + 7.29)
                /norm_Factor)
    #Correcting for inversefunc's antics
    else:
        return R
#R from 0 to 15 pc
X = np.arange(0, 20, 0.1)
for i in X:
    Y_2.append(C(i))
    Y_1.append(Sigma(i))
#Creating C_inv(y) = x
C_inv = inversefunc(C, y_values=random.random())
#creating function to return random positions (consistent with given density)
def Stellar_distribution(k):
    x = []
    y = []
    for h in range(0, k):
        rand = random.random()
        x.append(inversefunc(C, y_values=rand) * 1e3)
        y.append(rand)
        if (x[h] < 0):
            print(x[h], y[h], h)
    return (x, y)
        
"""Testing the method"""
# =============================================================================
# (x, y) = Stellar_distribution(300)
# 
# x[:] = [b/1000 for b in x]
# 
# fig, (ax, ax2) = plt.subplots(ncols=2)
# fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2)
# ax.plot(X, Y_1, 'k')
# ax.minorticks_on()
# ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# ax.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
#                 which = 'minor', top=True, right=True)
# ax.set_xlim(0, 15)
# ax.set_ylim(0,900)
# ax.set_box_aspect(1)
# ax.set(xlabel=u'$\mathregular{R/kpc}$', ylabel=u'$\mathregular{\Sigma (R)}$')
# plt.gcf().set_size_inches(15, 15)
# ax2.plot(X, Y_2, 'k')
# ax2.minorticks_on()
# ax2.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# ax2.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
#                 which = 'minor', top=True, right=True)
# ax2.set_xlim(0, 15)
# ax2.set_ylim(0,1)
# ax2.set_box_aspect(1)
# ax2.set(xlabel=u'$\mathregular{R/kpc}$', ylabel=u'$\mathregular{C(R)}$')
# plt.gcf().set_size_inches(15, 15)
# plt.show()
#     
# fig, ax3 = plt.subplots(ncols=1)
# ax3.scatter(x, y, color='k', lw = 0.5, s = 1, alpha = 1)
# ax3.minorticks_on()
# ax3.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# ax3.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
#                 which = 'minor', top=True, right=True)
# ax3.set_xlim(0, 15)
# ax3.set_ylim(0, 1)
# ax3.set_box_aspect(1)
# ax3.set(xlabel=u'$\mathregular{R/kpc}$', ylabel=u'$\mathregular{N(R<r)/N}$')
# plt.gcf().set_size_inches(8, 8)
# plt.show()
# 
# =============================================================================
