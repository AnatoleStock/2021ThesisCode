# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 15:50:14 2021

@author: Anatole Storck
"""

"""
This file uses the Runge-Kutta method to integrate orbits with initial values
in a binary system. Specifically: semi-major axis, eccentricity, masses.

INPUTS of RK4: a (semi-major axis), e (eccentricity), m1 (mass 1), m2 (mass 2),
t_birth (birth time of the DNSB).

OUTPUTS of RK4: A (full semi-major axis evolution), E (full eccentricity evolution),
T (full time evolution), t1 (time in T where the DNSB enters lisa band),
t2 (time in T where the DNSB leaves the LISA band).
"""

import math
import matplotlib.pyplot as plt
#Plot formating
plt.rcParams['axes.linewidth'] = 2.15
plt.rcParams["font.family"] = 'Courier New'
plt.rcParams.update({'font.size': 16})

"""Constants"""
#(1 R_solar --> 696340 km)
#With units of R_solar, years, and solar masses
G = 390983166.6
#In R_solar per year
c = 13577063.75

"""Post-Keplarian equations"""
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
#Making our orbital period in seconds using 3.154e7
def P(a, m1, m2):
    return (4*math.pi**2*a**3/(G*(m1 + m2)))**(0.5) * (3.154e7)
def f_gw(a, m1, m2):
    return 2/P(a, m1, m2)
#%%
"""RK4 Integrator for the LISA band"""
def RK4_LISA(a, e, m1, m2, t_birth):
    A = [] #Semi-major axis
    E = [] #Eccentricity
    T = [t_birth*1e6] #from Myr to yr since that's the units of RK4
    t = [] #t is the from the start to end of LISA band
    #Making our LISA Band endpoint requirement
    while f_gw(a, m1, m2) < 1:
        #Making our LISA Band entry-point requirement
        if f_gw(a, m1, m2) > 1e-5:
            t.append(T[len(T)-1])
        #choosing a constant of 1e-2.4
        delta_t = abs(a/a_dot(a, e, m1, m2))*(10**(-3.3))
        A.append(a)
        E.append(e)
        T.append(T[len(T)-1] + delta_t)
        #The Runge-Kutta method
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
        T[i] = T[i]/1e6#yr --> Myr
    T = T[:-1]
    if len(t) == 0:
        t1 = 0
        t2 = 0
    else:
        t1 = t[0]/1e6#yr --> Myr
        t2 = t[-1]/1e6#yr --> Myr
    return (A, E, T, t1, t2)
#%%
"""RK4 Integrator for the merger time"""
def RK4_Merging(a, e, m1, m2, t_birth):
    A = [] #Semi-major axis
    E = [] #Eccentricity
    T = [t_birth*1e6] #converting from Myr to yr since that's the units of RK4
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
    #Changing back to Myr
    for i in range(0, len(T)):
        T[i] = T[i]/1e6
    T = T[:-1]
    return (A, E, T)
#%%
"""TESTING THE METHOD"""
(x1, y1, T1, t1, t2) = RK4_LISA(2, 0.5, 1.4, 1.4, 0)
(x2, y2, T2, t1, t2) = RK4_LISA(1.95, 0.5, 1.4, 1.4, 0)
(x3, y3, T3, t1, t2) = RK4_LISA(1.9, 0.5, 1.4, 1.4, 0)
(x4, y4, T4, t1, t2) = RK4_LISA(1.85, 0.5, 1.4, 1.4, 0)
(x5, y5, T5, t1, t2) = RK4_LISA(1.80, 0.5, 1.4, 1.4, 0)
(x6, y6, T6, t1, t2) = RK4_LISA(1.75, 0.5, 1.4, 1.4, 0)
(x7, y7, T7, t1, t2) = RK4_LISA(1.70, 0.5, 1.4, 1.4, 0)
(x8, y8, T8, t1, t2) = RK4_LISA(1.65, 0.5, 1.4, 1.4, 0)
(x9, y9, T9, t1, t2) = RK4_LISA(1.60, 0.5, 1.4, 1.4, 0)
(x10, y10, T10, t1, t2) = RK4_LISA(1.55, 0.5, 1.4, 1.4, 0)

fig, (ax, ax2) = plt.subplots(2, sharex=True)
fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, hspace=0.1)
#Plotting semimajor-axis
ax.plot(T1, x1)
ax.plot(T2, x2)
ax.plot(T3, x3)
ax.plot(T4, x4)
ax.plot(T5, x5)
ax.plot(T6, x6)
ax.plot(T7, x7)
ax.plot(T8, x8)
ax.plot(T9, x9)
ax.plot(T10, x10)
ax.minorticks_on()
ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax.set_xlim(T1[0])
ax.set_ylim(0, x1[0])
ax.set_box_aspect(1)
ax.set(ylabel=u'$\mathregular{a/R_\odot}$')
ax.set_title('Semi-major axis evolution')
plt.gcf().set_size_inches(15, 15)
#Plotting eccentricity
ax2.plot(T1, y1, 'k')
ax2.plot(T2, y2, 'k')
ax2.plot(T3, y3, 'k')
ax2.minorticks_on()
ax2.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax2.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax2.set_xlim(T1[0])
ax2.set_ylim(0, y1[0])
ax2.set_box_aspect(1)
ax2.set(xlabel=u'$\mathregular{t/Myr}$', ylabel=u'$\mathregular{e}$')
ax2.set_title('Eccentricity evolution')
plt.gcf().set_size_inches(15, 15)
#fig.savefig('BOI_plot.png', dpi=300, bbox_inches='tight')
plt.show()
