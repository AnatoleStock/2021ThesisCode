# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 11:18:30 2021

@author: Anatole Storck
"""

"""
This file uses the Runge-Kutta method to integrate orbits with initial values
in the Milky Way. Specifically: orbital velocity, position, mass (cylindrical coord.)

INPUTS of RK4: R (distance, in galactic plane, from galactic center),
z (distance out of the galactic plane), phi (angle around the galactic plane),
VR (initial velocity in R direction), Vz (initial velocity in z direction),
V_phi0 (initial velocity in phi direction).

OUTPUTS of RK4: P (tuple of evoluton of R, z, and phi),
V (tuple of evolution of VR and Vz), E (evolution of specific energy of system),
T (evolution of time).
"""

"""Preamble"""
import numpy as np
import math
import matplotlib.pyplot as plt
#Plot formating
plt.rcParams['axes.linewidth'] = 2.15
plt.rcParams["font.family"] = 'Courier New'
plt.rcParams.update({'font.size': 16})

"""Constants for the GalPot"""
#conversion for (km/pc) --> 3.086e13
#pc solarmasses**-1 (km/s)**2
G = 4.301e-3
#solarmasses
M_s = 1.12e10
#pc
a_s = 0.0
b_s = 0.277e3
#solarmasses
M_d = 8.07e10
#pc
a_d = 3.7e3
b_d = 0.20e3
#solarmasses
M_c = 5.0e10
#pc
r_c = 6.0e3

""""GalPot and related equations"""
#Galactic Potential Components (Repetto 2018)
def GalPot_s(R, z):
    return -1*(G*M_s)/((R**2+(a_s+(z**2+b_s**2)**0.5)**2)**(0.5))
def GalPot_d(R, z):
    return -1*(G*M_d)/((R**2+(a_d+(z**2+b_d**2)**0.5)**2)**(0.5))
def GalPot_h(R, z):
    return ((G*M_c/r_c)*
            (0.5*math.log(1 + ((R**2 + z**2)/r_c**2), math.e)
             + (r_c/(R**2 + z**2)**0.5)*math.atan((R**2 + z**2)**0.5/r_c)))
#Full Potential
def GalPot(R, z):
    return GalPot_s(R, z) + GalPot_d(R, z) + GalPot_h(R, z)
#Partial Derivaties of the Galactic Potential
def GalPot_s_dR(R, z):
    return GalPot_s(R, z)*(-1)*R*(R**2+(a_s+(z**2+b_s**2)**0.5)**2)**(-1)
def GalPot_d_dR(R, z):
    return GalPot_d(R, z)*(-1)*R*(R**2+(a_d+(z**2+b_d**2)**0.5)**2)**(-1)
def GalPot_h_dR(R, z):
    return (G*M_c*R*(1/(r_c*(R**2 + z**2))
                     - (math.atan((R**2+z**2)**0.5/r_c)/(R**2+z**2)**1.5)))
def GalPot_s_dz(R, z):
    return (GalPot_s(R, z)*(-1)*z*(a_s+(z**2+b_s**2)**(0.5))*
            ((z**2+b_s**2)**(0.5)*(R**2+(a_s+(z**2+b_s**2)**(0.5))**2))**(-1))
def GalPot_d_dz(R, z):
    return (GalPot_d(R, z)*(-1)*z*(a_d+(z**2+b_d**2)**(0.5))*
            ((z**2+b_d**2)**(0.5)*(R**2+(a_d+(z**2+b_d**2)**(0.5))**2))**(-1))
def GalPot_h_dz(R, z):
    return (G*M_c*z*(1/(r_c*(R**2 + z**2))
                     - (math.atan((R**2+z**2)**0.5/r_c)/(R**2+z**2)**1.5)))
#Full Partial Derivatives
def GalPot_dR(R, z):
    return GalPot_s_dR(R, z) + GalPot_d_dR(R, z) + GalPot_h_dR(R, z)
def GalPot_dz(R, z):
    return GalPot_s_dz(R, z) + GalPot_d_dz(R, z) + GalPot_h_dz(R, z)
#Defining our Velocity Derivatives and V_phi
def VR_dot(R, z, J_z):
    return ((J_z**2/R**3) - GalPot_dR(R, z))
def Vz_dot(R, z):
    return -1*GalPot_dz(R, z)
def V_phi(R, J_z):
    return J_z/R
#Circular Velocity
def Vc(R, z):
    return (R*GalPot_dR(R, z))**(0.5)
#Specific Energy Calculator
def Epsilon(Vr, Vz, R, z, J_z):
    return 0.5*(Vr**2 + Vz**2 + V_phi(R, J_z)**2) + GalPot(R, z)

"""RK4 Integrator"""
def Integrator(R, z, phi, VR, Vz, Vphi_0, T_birth):    
    T = [T_birth]
    #P is in R, z, phi
    P = ([], [], [])
    #V is in VR and Vz
    V = ([], [])
    #Specific energy
    E = []
    #Calculating initial J_z since it is constant
    J_z = Vphi_0*R
    #Numerator: abs(V_phi(R, J_z)
    #for i in range(0, 3000):
    while (T[-1] < 0):
        #Time-step condition
        if (abs(V_phi(R, J_z)/R)**0.85 > 8e-3):
            delta_t = 8e-3/abs(V_phi(R, J_z)/R)**0.85
        else:
            delta_t = 1
        #Alternative time-step condition
# =============================================================================
#         if ((1/R**(1/3))*delta_t > 1e-6):
#             delta_t = 1e-3/(1/R**(1/3))
# =============================================================================
        T.append(T[len(T)-1] + delta_t)
        P[0].append(R)
        P[1].append(z)
        P[2].append(phi)
        V[0].append(VR)
        V[1].append(Vz)
        E.append(Epsilon(VR, Vz, R, z, J_z))
        #Variable order for RK4 is (R, z, phi, Vr, Vz)
        k1 = (VR, Vz, V_phi(R, J_z)/R, VR_dot(R, z, J_z), Vz_dot(R, z))
        k2 = (VR + delta_t*(k1[3]/2),
              Vz + delta_t*(k1[4]/2),
              V_phi(R + delta_t*(k1[0]/2), J_z)/(R + delta_t*(k1[0]/2)),
              VR_dot(R + delta_t*(k1[0]/2), z + delta_t*(k1[1]/2), J_z),
              Vz_dot(R + delta_t*(k1[0]/2), z + delta_t*(k1[1]/2)))
        k3 = (VR + delta_t*(k2[3]/2),
              Vz + delta_t*(k2[4]/2),
              V_phi(R + delta_t*(k2[0]/2), J_z)/(R + delta_t*(k2[1]/2)),
              VR_dot(R + delta_t*(k2[0]/2), z + delta_t*(k2[1]/2), J_z),
              Vz_dot(R + delta_t*(k2[0]/2), z + delta_t*(k2[1]/2)))
        k4 = (VR + delta_t*(k3[3]),
              Vz + delta_t*(k3[4]),
              V_phi(R + delta_t*(k3[0]), J_z)/(R + delta_t*(k3[0])),
              VR_dot(R + delta_t*(k3[0]), z + delta_t*(k3[1]), J_z),
              Vz_dot(R + delta_t*(k3[0]), z + delta_t*(k3[1])))
        (R, z, phi, VR, Vz) = (R+(delta_t/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),
                               z+(delta_t/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1]),
                               phi+(delta_t/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2]),
                               VR+(delta_t/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3]),
                               Vz+(delta_t/6)*(k1[4]+2*k2[4]+2*k3[4]+k4[4]),)
    T = T[:-1]
    return (P, V, E, T)
#%%
"""TESTING THE METHOD"""
(P, V, E, T) = Integrator(3787, 0, 0.62825,
                          41.4, -11.1, 40.4, -932)

#Plotting
fig, (ax1, ax2) = plt.subplots(ncols=2)
#Plotting R
ax1.plot(T, P[0], 'k')
ax1.minorticks_on()
ax1.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax1.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax1.set_box_aspect(1)
ax1.set(xlabel=u'$\mathregular{t/Myr}$', ylabel=u'$\mathregular{R/pc}$')
plt.gcf().set_size_inches(15, 15)
#Plotting Z
ax2.plot(T, P[1], 'k')
ax2.minorticks_on()
ax2.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax2.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax2.set_box_aspect(1)
ax2.set(xlabel=u'$\mathregular{t/Myr}$', ylabel=u'$\mathregular{z/pc}$')
plt.gcf().set_size_inches(15, 15)
plt.show()
fig, (ax3, ax4) = plt.subplots(ncols=2)
#Plotting z against R
ax3.scatter(P[0], P[1], color='k', lw = 0.5, s = 0.5, alpha = 0.5)
ax3.minorticks_on()
ax3.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax3.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax3.set_box_aspect(1)
ax3.set(xlabel=u'$\mathregular{R/pc}$', ylabel=u'$\mathregular{z/pc}$')
plt.gcf().set_size_inches(15, 15)
#Plotting phi
ax4.scatter(T, E, color='k', lw = 0.5, s = 0.5, alpha = 0.5)
ax4.minorticks_on()
ax4.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax4.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax4.set_box_aspect(1)
ax4.set(xlabel=u'$\mathregular{t/Myr}$', ylabel=u'$\mathregular{\epsilon}$')
plt.gcf().set_size_inches(15, 15)
plt.show()
fig, (ax5, ax6) = plt.subplots(ncols=2)
#Plotting Vr
ax5.scatter(T, V[0], color='k', lw = 0.5, s = 0.5, alpha = 0.5)
ax5.minorticks_on()
ax5.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax5.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax5.set_box_aspect(1)
ax5.set(xlabel=u'$\mathregular{t/Myr}$', ylabel=u'$\mathregular{Vr/(km/s)}$')
plt.gcf().set_size_inches(15, 15)
#Plotting Vz
ax6.scatter(T, V[1], color='k', lw = 0.5, s = 0.5, alpha = 0.5)
ax6.minorticks_on()
ax6.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax6.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax6.set_box_aspect(1)
ax6.set(xlabel=u'$\mathregular{t/Myr}$', ylabel=u'$\mathregular{Vz/(km/s)}$')
plt.gcf().set_size_inches(15, 15)
plt.show()

fig, ax7 = plt.subplots(ncols=1)
#Plotting Energy
ax7.scatter(T, E, color='k', lw = 1.5, s = 0.5, alpha = 0.5)
ax7.minorticks_on()
ax7.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax7.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
#ax7.set_box_aspect(1)
ax7.set(xlabel=u'$\mathregular{t/Myr}$', ylabel=u'$\mathregular{\epsilon}$')
plt.gcf().set_size_inches(15, 8)
plt.show()
#3D Plotting
X, Y, Z = P[0]*np.cos(P[2]), P[0]*np.sin(P[2]), P[1]

fig = plt.figure()
ax8 = plt.axes(projection='3d')
ax8.scatter(X, Y, Z, c='k', s=0.05)
ax8.set_box_aspect([1,1,1])
plt.show()
fig, (ax9, ax10) = plt.subplots(ncols=2)
#Plotting RvsVr
ax9.scatter(P[0], V[1], color='k', lw = 0.5, s = 0.5, alpha = 0.5)
ax9.minorticks_on()
ax9.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax9.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax9.set_box_aspect(1)
ax9.set(xlabel=u'$\mathregular{R/pc}$', ylabel=u'$\mathregular{Vz/(km/s)}$')
plt.gcf().set_size_inches(15, 15)
#Plotting zvsVz
ax10.scatter(P[1], V[0], color='k', lw = 0.5, s = 0.5, alpha = 0.5)
ax10.minorticks_on()
ax10.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax10.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax10.set_box_aspect(1)
ax10.set(xlabel=u'$\mathregular{z/pc}$', ylabel=u'$\mathregular{Vr/(km/s)}$')
plt.gcf().set_size_inches(15, 15)
plt.show()
