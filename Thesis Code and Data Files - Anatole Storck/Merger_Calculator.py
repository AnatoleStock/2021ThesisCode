# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 15:50:14 2021

@author: Anatole Storck
"""

"""
This file uses the Runge-Kutta method to integrate orbits with initial values
in a binary system. Specifically: semi-major axis, eccentricity, masses.
This file is similar compared to the Binary_Orbit_Integrator except it asses
when the binary mergers. We use this to calculate merger times.
"""

import math
import matplotlib.pyplot as plt
#Plot formating stuff
plt.rcParams['axes.linewidth'] = 2.15
plt.rcParams["font.family"] = 'Courier New'
plt.rcParams.update({'font.size': 16})
#R_solar --> 696340 km
#with units of R_solar, years, and solar masses
G = 390983166.6
#In R_solar per year
c = 13577063.75
#From Peter's 1964 thesis (Eq. 5.41 and 5.42)
def a_dot(a, e, m1, m2):
    term_1 = -1*(64/5)*G**3
    term_2 = (m1*m2*(m1+m2))/(c**5*a**3)
    term_3 = (1-e**2)**(-7/2)
    term_4 = 1 + (73/24)*e**2 + (37/96)*e**4
    return term_1*term_2*term_3*term_4
def e_dot(a, e, m1, m2):
    term_1 = -1*(304/15)*G**3*e
    term_2 = (m1*m2*(m1+m2))/(c**5*a**4)
    term_3 = (1-e**2)**(-5/2)
    term_4 = 1 + (121/304)*e**2
    return term_1*term_2*term_3*term_4
#Creating our integrator using the Runge-Kutta method
def RK4(a, e, m1, m2, t_birth):
    A = [] #Semi-major axis
    E = [] #Eccentricity
    T = [t_birth*1e6] #from Myr to yr since that's the units of RK4
    while a > (50/696340):
        #choosing a constant of 1e-2.4
        delta_t = abs(a/a_dot(a, e, m1, m2))*(10**(-3.3))
        A.append(a)
        E.append(e)
        T.append(T[len(T)-1] + delta_t)
        k1 = (a_dot(a, e, m1, m2), e_dot(a, e, m1, m2))
        k2 = (a_dot(a + delta_t*(k1[0]/2), e + delta_t*(k1[1]/2), m1, m2),
              e_dot(a + delta_t*(k1[0]/2), e + delta_t*(k1[1]/2), m1, m2))
        k3 = (a_dot(a + delta_t*(k2[0]/2), e + delta_t*(k2[1]/2), m1, m2),
              e_dot(a + delta_t*(k2[0]/2), e + delta_t*(k2[1]/2), m1, m2))
        k4 = (a_dot(a + delta_t*k3[0], e + delta_t*k3[1], m1, m2),
              e_dot(a + delta_t*k3[0], e + delta_t*k3[1], m1, m2))
        (a, e) = (a + (delta_t/6)*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]),
                  e + (delta_t/6)*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]))
    #Changing back to Myr and defining t1 and t2 (entering/leave LISA Band)
    for i in range(0, len(T)):
        T[i] = T[i]/1e6
    T = T[:-1]
    return (A, E, T)

#TESTING THE METHOD
# =============================================================================
# (x, y, T) = RK4(2, 0.5, 1.4, 1.4, -3e3)
# 
# fig, (ax, ax2) = plt.subplots(ncols=2)
# fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2)
# #Plotting semimajor-axis
# ax.plot(T, x, 'k')
# ax.minorticks_on()
# ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# ax.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
#                 which = 'minor', top=True, right=True)
# ax.set_xlim(T[0])
# ax.set_ylim(0, x[0])
# ax.set_box_aspect(1)
# ax.set(xlabel=u'$\mathregular{t/Myr}$', ylabel=u'$\mathregular{a/R_\odot}$')
# plt.gcf().set_size_inches(15, 15)
# #Plotting eccentricity
# ax2.plot(T, y, 'k')
# ax2.minorticks_on()
# ax2.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# ax2.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
#                 which = 'minor', top=True, right=True)
# ax2.set_xlim(T[0])
# ax2.set_ylim(0, y[0])
# ax2.set_box_aspect(1)
# ax2.set(xlabel=u'$\mathregular{t/Myr}$', ylabel=u'$\mathregular{e}$')
# plt.gcf().set_size_inches(15, 15)
# plt.show()
# =============================================================================
