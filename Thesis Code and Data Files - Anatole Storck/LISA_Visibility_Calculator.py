# -*- coding: utf-8 -*-
"""
Created in April 2021

@author: Anatole Storck
"""

"""
This file looks at the visibility of DNSBs to LISA. Calculating their SNR and
plotting various relationships.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import copy
import pandas as pd
import cmocean
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

"""Preamble"""
plt.rcParams['axes.linewidth'] = 2.15
plt.rcParams["font.family"] = 'Courier New'
plt.rcParams.update({'font.size': 16})
cmap = cmocean.cm.thermal
newcmap = cmocean.tools.crop_by_percent(cmap, 4, which='max', N=None)
color = newcmap
cmap_r = cmocean.cm.haline_r
newcmap_r = cmocean.tools.crop_by_percent(cmap_r, 15, which='min', N=None)
color_r = newcmap_r
color1 = cmocean.cm.solar
color2 = cmocean.cm.haline
#Constants
G = 4.30091e-3*3.086e13
c = 299792.458
#Importing Files (Dat comes from Ross' data files while Val are calculated and given values)
Vis_DNSB_Val = pd.read_csv(r'C:\Users\storc\OneDrive\Thesis Work\Data\Visible_DNSB_Values(300_2_upper200).csv')
Vis_DNSB_Dat = pd.read_csv(r'C:\Users\storc\OneDrive\Thesis Work\Data\Visible_DNSB_Data(300_2_upper200).csv')
#Fin are final positions R and z
Vis_DNSB_Fin = pd.read_csv(r'C:\Users\storc\OneDrive\Thesis Work\Data\Final_Position(300_2_upper200).csv')
#length of the list
List_len = len(Vis_DNSB_Val['Initial Position (R)'])
#converting from Cyl. Gal. to Helio Gal.
R_dot = 8.330 #kpc
R = copy.deepcopy(Vis_DNSB_Fin['Final Position (R)'])
z = copy.deepcopy(Vis_DNSB_Fin['Final Position (z)'])
Phi = copy.deepcopy(Vis_DNSB_Fin['Final Angle (Phi)'])
#To Cart. Gal
X = R*np.cos(Phi)
Y = R*np.sin(Phi)
Z = z
#Translating the system to the solar position
X_sol = X - 8.330
#Helio Gal.
r = np.sqrt(np.square(X_sol) + np.square(Y) + np.square(Z))
l = np.arctan2(-Y,-1*X_sol)
b = np.arcsin(Z/(np.sqrt(np.square(X_sol) + np.square(Y) + np.square(Z))))
#To degrees
l_deg = l*(180/np.pi)
b_deg = b*(180/np.pi)

"""
9 DNS systems Properties
"""
RP1 = 'J0453+1559'
RP1_m1 = 1.559
RP1_m2 = 1.174
RP1_P = 4.072*(86400) #s
RP1_d = 1.07*1e3
RP2 = 'J0737-3039'
RP2_m1 = 1.338
RP2_m2 = 1.249
RP2_P = 0.102*(86400) #s
RP2_d = 1.15*1e3 
RP3 = 'B1534+12'
RP3_m1 = 1.333
RP3_m2 = 1.346
RP3_P = 0.421*(86400) #s
RP3_d = 1.05*1e3 
RP4 = 'J1756-2251'
RP4_m1 = 1.341
RP4_m2 = 1.230
RP4_P = 0.320*(86400) #s
RP4_d = 0.73*1e3
RP5 = 'J1811-1736'
RP5_m1 = 1.64 #>
RP5_m2 = 0.93 #<
RP5_P = 18.779*(86400) #s
RP5_d = 5.93*1e3 
RP6 = 'J1829+2456'
RP6_m1 = 1.38 #>
RP6_m2 = 1.22 #<
RP6_P = 1.176*(86400) #s
RP6_d = 0.74*1e3 
RP7 = 'J1906+0746'
RP7_m1 = 1.291
RP7_m2 = 1.322
RP7_P = 0.166*(86400) #s
RP7_d = 7.40*1e3 
RP8 = 'B1913+16'
RP8_m1 = 1.440
RP8_m2 = 1.389
RP8_P = 0.323*(86400) #s
RP8_d = 9.80*1e3 
RP9 = 'J1930-1852'
RP9_m1 = 1.32 #>
RP9_m2 = 1.30 #<
RP9_e = 0.399
RP9_P = 45.060*(86400) #s
RP9_d = 1.5*1e3 

"""
10 DWD system properties
"""
WD1 = 'HM Cnc'
WD1_m1 = 0.55
WD1_m2 = 0.27
WD1_P = 321.529 #s
WD1_d = 5*1e3 #kpc
WD2 = 'V407 Vul'
WD2_m1 = 0.8
WD2_m2 = 0.177
WD2_P = 569.395 #s
WD2_d = 1.786*1e3 #kpc
WD3 = 'ES Cet'
WD3_m1 = 0.8
WD3_m2 = 0.161
WD3_P = 620.21 #s
WD3_d = 1.584*1e3 #kpc
WD4 = 'SDSS J135154.46-064309.0'
WD4_m1 = 0.8
WD4_m2 = 0.1
WD4_P = 943.84 #s
WD4_d = 1.317*1e3 #kpc
WD5 = 'AM CVn'
WD5_m1 = 0.68
WD5_m2 = 0.125
WD5_P = 1028.73 #s
WD5_d = 0.299*1e3 #kpc
WD6 = 'SDSS J190817.07+394036.4'
WD6_m1 = 0.8
WD6_m2 = 0.085
WD6_P = 1085.7 #s
WD6_d = 1.044*1e3 #kpc
WD7 = 'HP Lib'
WD7_m1 = 0.49
WD7_m2 = 0.048
WD7_P = 1102.70 #s
WD7_d = 0.276*1e3 #kpc
WD8 = 'SDSS J065133.34+284423.4'
WD8_m1 = 0.247
WD8_m2 = 0.49
WD8_P = 765.5 #s
WD8_d = 0.933*1e3 #kpc
WD9 = 'SDSS J093506.92+441107.0'
WD9_m1 = 0.312
WD9_m2 = 0.75
WD9_P = 1188.0 #s
WD9_d = 0.645*1e3 #kpc
WD10 = 'CR Boo'
WD10_m1 = 0.67
WD10_m2 = 0.044
WD10_P = 1471.3 #s
WD10_d = 0.337*1e3 #kpc
WD11 = 'V803 Cen'
WD11_m1 = 0.78
WD11_m2 = 0.059
WD11_P = 1596.4 #s
WD11_d = 0.347*1e3 #kpc

"""Equations used for LISA curve, background, PSD, h_c, SNR"""
#From LISA SciDoc
def R(f):
    return 1 + (f/(25e-3))**2
def S_II(f):
    return 3.6e-41
def S_I(f):
    return 5.76e-48*(1 + (f/(0.4e-3))**-2)
def S_h(f):
    return (10/3)*((S_I(f)/(2*np.pi*f)**4) + S_II(f))*R(f)
def srqt_S_h(f):
    return np.sqrt(S_h(f))
def h_s(f):
    return np.sqrt(f*S_h(f))
#Creating a frequency range
F_background = np.logspace(-5, -2.2, 1000)
def F1(f):
    return 7e-20 *(f/7.5e-5)**(-0.734862392)
def F2(f):
    return 2e-22 *(f/3e-3)**(-8.82746912)
def F3(f):
    return ((1/F1(f)) + (1/F2(f)))**(-1)
def F3PSD(f):
    return F3(f)*f**(-1/2)
#Creating another frequency range
F = np.logspace(-5, 1, num=1000)
def M_c(m1, m2):
    return (m1*m2)**(3/5)/(m1+m2)**(1/5)
def N_cycle(f_gw, T_obs):
    return f_gw * (T_obs*3.154e7)
def h(m1, m2, d, f_gw):
    return (2*(G*M_c(m1, m2))**(5/3)/(c**4*d*3.086e13))*(np.pi*f_gw)**(2/3)
def h_c(m1, m2, d, f_gw, T_obs):
    return h(m1, m2, d, f_gw)*np.sqrt(N_cycle(f_gw, T_obs))
def sqrt_S_h(m1, m2, d, f_gw, T_obs):
    return h_c(m1, m2, d, f_gw, T_obs)*f_gw**(-1/2)

def SNR(m1, m2, d, f_gw, T_obs):
    return (h_c(m1, m2, d, f_gw, T_obs)/max(h_s(f_gw), F3(f_gw)))
#Creating several lists
DNSB_f_gw = []
DNSB_h_c = []
DNSB_sqrt_S_h = []
DNSB_SNR = []
DNSB_SNR_noBackground = []
#Appending several lists
for i in range(0, List_len):
    f = 2/Vis_DNSB_Val.at[i, 'Period at present (P)']
    DNSB_f_gw.append(f)
    DNSB_h_c.append(h_c(Vis_DNSB_Dat.at[i, 'Mass1 (m1)'],
                        Vis_DNSB_Dat.at[i, 'Mass2 (m2)'], r[i]*1e3, f, 4))
    DNSB_sqrt_S_h.append(sqrt_S_h(Vis_DNSB_Dat.at[i, 'Mass1 (m1)'],
                        Vis_DNSB_Dat.at[i, 'Mass2 (m2)'], r[i]*1e3, f, 4))
    DNSB_SNR.append(SNR(Vis_DNSB_Dat.at[i, 'Mass1 (m1)'],
                        Vis_DNSB_Dat.at[i, 'Mass2 (m2)'], r[i]*1e3, f, 4))
    DNSB_SNR_noBackground.append(h_c(Vis_DNSB_Dat.at[i, 'Mass1 (m1)'],
                        Vis_DNSB_Dat.at[i, 'Mass2 (m2)'], r[i]*1e3, f, 4)/
                                 h_s(f))

"""
PLOTTING THE LISA CURVE
"""
#Plotting
# =============================================================================
# fig, (ax, ax1) = plt.subplots(ncols=2)
# fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2)
# # SIGNAL AND SENSITIVITY CURVES FOR THE ENTIRE OBSERVATION TIME
# ax.loglog(F, srqt_S_h(F), 'k', linestyle='solid', linewidth=1, label='Lisa curve')
# ax.loglog(F_background, F3PSD(F_background), 'k', linestyle=(0, (5, 5)), linewidth=1, label='Background')
# ax.minorticks_on()
# ax.tick_params(direction='in', length=10, width=1.5,
#                 which = 'major', top=True, right=True, pad=7)
# ax.tick_params(direction='in', length=5, width=1.5,
#                 which = 'minor', top=True, right=True)
# ax.set_xlim(1e-5, 1)
# ax.set_ylim(1e-21, 1e-16)
# ax.set_xlabel(u'$\mathregular{Frequency/Hz}$')
# ax.set_ylabel(r'$\mathregular{Strain\:\:\:\:Sensitivity\:\:\:\:\left[\:\sqrt{Hz^{-1}}\:\:\:\right]}$')
# ax.set_box_aspect(1)
# ax.legend(loc="upper right")
# plt.gcf().set_size_inches(15, 15)
# plt.tight_layout()
# 
# # SIGNAL AND SENSITIVITY CURVES FOR THE ENTIRE OBSERVATION TIME
# ax1.loglog(F, h_s(F), 'k', linestyle='solid', linewidth=1, label='Lisa curve')
# ax1.loglog(F_background, F3(F_background), 'k', linestyle=(0, (5, 5)), linewidth=1, label='Background')
# ax1.minorticks_on()
# ax1.tick_params(direction='in', length=10, width=1.5,
#                 which = 'major', top=True, right=True, pad=7)
# ax1.tick_params(direction='in', length=5, width=1.5,
#                 which = 'minor', top=True, right=True)
# ax1.set_xlim(1.01e-5, 1)
# ax1.set_ylim(1e-22, 1e-18)
# ax1.set_xlabel(u'$\mathregular{Frequency/Hz}$')
# ax1.set_ylabel(u'$\mathregular{Characteristic\:\:\:\:Strain}$')
# ax1.set_box_aspect(1)
# ax1.legend(loc="upper right")
# plt.gcf().set_size_inches(15, 15)
# plt.tight_layout()
# fig.savefig('LISA_curve.png', dpi=250, bbox_inches='tight')
# plt.show()
# =============================================================================

"""
PLOTTING THE LISA CURVE WITH ALL BINARIES
"""
# =============================================================================
# fig, (ax2, ax3) = plt.subplots(ncols=2)
# fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2)
# ax2.grid(which='both')
# ax2.set_axisbelow(True)
# # SIGNAL AND SENSITIVITY CURVES FOR THE ENTIRE OBSERVATION TIME
# ax2.loglog(F, srqt_S_h(F), 'b', linestyle='solid', linewidth=1, label='Lisa curve')
# ax2.plot(F_background, F3PSD(F_background), 'r', linestyle='solid', linewidth=1, label='Background')
# ax2.scatter(DNSB_f_gw, DNSB_sqrt_S_h, color='k', lw = 0.35, s = 7, alpha = 0.75, marker='^')
# ax2.scatter(1, 1, color='k', s=65, marker='d', label='Observered DNSBs')
# ax2.scatter(1, 1, color='k', s=65, marker='*', label='Verification DWDs')
# ax2.scatter(1, 1, color='k', s=65, marker='^', label='Resolvable DNSBs')
# #DNSB
# ax2.scatter(2/RP1_P, sqrt_S_h(RP1_m1, RP1_m2, RP1_d, 2/RP1_P, 4), s=65, marker='d')
# ax2.scatter(2/RP2_P, sqrt_S_h(RP2_m1, RP2_m2, RP2_d, 2/RP2_P, 4), s=65, marker='d')
# ax2.scatter(2/RP3_P, sqrt_S_h(RP3_m1, RP3_m2, RP3_d, 2/RP3_P, 4), s=65, marker='d')
# ax2.scatter(2/RP4_P, sqrt_S_h(RP4_m1, RP4_m2, RP4_d, 2/RP4_P, 4), s=65, marker='d')
# ax2.scatter(2/RP5_P, sqrt_S_h(RP5_m1, RP5_m2, RP5_d, 2/RP5_P, 4), s=65, marker='d')
# ax2.scatter(2/RP6_P, sqrt_S_h(RP6_m1, RP6_m2, RP6_d, 2/RP6_P, 4), s=65, marker='d')
# ax2.scatter(2/RP7_P, sqrt_S_h(RP7_m1, RP7_m2, RP7_d, 2/RP7_P, 4), s=65, marker='d')
# ax2.scatter(2/RP8_P, sqrt_S_h(RP8_m1, RP8_m2, RP8_d, 2/RP8_P, 4), s=65, marker='d')
# ax2.scatter(2/RP9_P, sqrt_S_h(RP9_m1, RP9_m2, RP9_d, 2/RP9_P, 4), s=65, marker='d')
# #DWDB
# ax2.scatter(2/WD1_P, sqrt_S_h(WD1_m1, WD1_m2, WD1_d, 2/WD1_P, 4), s=70, marker='*')
# ax2.scatter(2/WD2_P, sqrt_S_h(WD2_m1, WD2_m2, WD2_d, 2/WD2_P, 4), s=70, marker='*')
# ax2.scatter(2/WD3_P, sqrt_S_h(WD3_m1, WD3_m2, WD3_d, 2/WD3_P, 4), s=70, marker='*')
# ax2.scatter(2/WD4_P, sqrt_S_h(WD4_m1, WD4_m2, WD4_d, 2/WD4_P, 4), s=70, marker='*')
# ax2.scatter(2/WD5_P, sqrt_S_h(WD5_m1, WD5_m2, WD5_d, 2/WD5_P, 4), s=70, marker='*')
# ax2.scatter(2/WD6_P, sqrt_S_h(WD6_m1, WD6_m2, WD6_d, 2/WD6_P, 4), s=70, marker='*')
# ax2.scatter(2/WD7_P, sqrt_S_h(WD7_m1, WD7_m2, WD7_d, 2/WD7_P, 4), s=70, marker='*')
# ax2.scatter(2/WD8_P, sqrt_S_h(WD8_m1, WD8_m2, WD8_d, 2/WD8_P, 4), s=70, marker='*')
# ax2.scatter(2/WD9_P, sqrt_S_h(WD9_m1, WD9_m2, WD9_d, 2/WD9_P, 4), s=70, marker='*')
# ax2.scatter(2/WD10_P, sqrt_S_h(WD10_m1, WD10_m2, WD10_d, 2/WD10_P, 4), s=70, marker='*')
# ax2.scatter(2/WD11_P, sqrt_S_h(WD11_m1, WD11_m2, WD11_d, 2/WD11_P, 4), s=70, marker='*')
# ax2.minorticks_on()
# ax2.tick_params(direction='in', length=10, width=1.5,
#                 which = 'major', top=True, right=True, pad=7)
# ax2.tick_params(direction='in', length=5, width=1.5,
#                 which = 'minor', top=True, right=True)
# ax2.set_xlim(1.01e-5, 2e-1)
# ax2.set_ylim(1e-21, 1e-16)
# ax2.set_xlabel(u'$\mathregular{Frequency/Hz}$')
# ax2.set_ylabel(u'$\mathregular{Strain\:\:\:\:Sensitivity\:\:\:\:[\sqrt{Hz}^{-1}]}$')
# ax2.set_box_aspect(1)
# ax2.legend(loc="lower right", fontsize="small")
# plt.gcf().set_size_inches(15, 15)
# plt.tight_layout()
# 
# ax3.grid(which='both')
# ax3.set_axisbelow(True)
# # SIGNAL AND SENSITIVITY CURVES FOR THE ENTIRE OBSERVATION TIME
# ax3.loglog(F, h_s(F), 'b', linestyle='solid', linewidth=1, label='Lisa curve')
# ax3.plot(F_background, F3(F_background), 'r', linestyle='solid', linewidth=1, label='Background')
# ax3.scatter(DNSB_f_gw, DNSB_h_c, color='k', lw = 0.35, s = 7, alpha = 0.75, marker='^')
# ax3.scatter(1, 1, color='k', s=65, marker='d', label='Observered DNSBs')
# ax3.scatter(1, 1, color='k', s=65, marker='*', label='Verification DWDs')
# ax3.scatter(1, 1, color='k', s=65, marker='^', label='Resolvable DSNBs')
# #DNSB
# ax3.scatter(2/RP1_P, h_c(RP1_m1, RP1_m2, RP1_d, 2/RP1_P, 4), s=65, marker='d')
# ax3.scatter(2/RP2_P, h_c(RP2_m1, RP2_m2, RP2_d, 2/RP2_P, 4), s=65, marker='d')
# ax3.scatter(2/RP3_P, h_c(RP3_m1, RP3_m2, RP3_d, 2/RP3_P, 4), s=65, marker='d')
# ax3.scatter(2/RP4_P, h_c(RP4_m1, RP4_m2, RP4_d, 2/RP4_P, 4), s=65, marker='d')
# ax3.scatter(2/RP5_P, h_c(RP5_m1, RP5_m2, RP5_d, 2/RP5_P, 4), s=65, marker='d')
# ax3.scatter(2/RP6_P, h_c(RP6_m1, RP6_m2, RP6_d, 2/RP6_P, 4), s=65, marker='d')
# ax3.scatter(2/RP7_P, h_c(RP7_m1, RP7_m2, RP7_d, 2/RP7_P, 4), s=65, marker='d')
# ax3.scatter(2/RP8_P, h_c(RP8_m1, RP8_m2, RP8_d, 2/RP8_P, 4), s=65, marker='d')
# ax3.scatter(2/RP9_P, h_c(RP9_m1, RP9_m2, RP9_d, 2/RP9_P, 4), s=65, marker='d')
# #DWDB
# ax3.scatter(2/WD1_P, h_c(WD1_m1, WD1_m2, WD1_d, 2/WD1_P, 4), s=70, marker='*')
# ax3.scatter(2/WD2_P, h_c(WD2_m1, WD2_m2, WD2_d, 2/WD2_P, 4), s=70, marker='*')
# ax3.scatter(2/WD3_P, h_c(WD3_m1, WD3_m2, WD3_d, 2/WD3_P, 4), s=70, marker='*')
# ax3.scatter(2/WD4_P, h_c(WD4_m1, WD4_m2, WD4_d, 2/WD4_P, 4), s=70, marker='*')
# ax3.scatter(2/WD5_P, h_c(WD5_m1, WD5_m2, WD5_d, 2/WD5_P, 4), s=70, marker='*')
# ax3.scatter(2/WD6_P, h_c(WD6_m1, WD6_m2, WD6_d, 2/WD6_P, 4), s=70, marker='*')
# ax3.scatter(2/WD7_P, h_c(WD7_m1, WD7_m2, WD7_d, 2/WD7_P, 4), s=70, marker='*')
# ax3.scatter(2/WD8_P, h_c(WD8_m1, WD8_m2, WD8_d, 2/WD8_P, 4), s=70, marker='*')
# ax3.scatter(2/WD9_P, h_c(WD9_m1, WD9_m2, WD9_d, 2/WD9_P, 4), s=70, marker='*')
# ax3.scatter(2/WD10_P, h_c(WD10_m1, WD10_m2, WD10_d, 2/WD10_P, 4), s=70, marker='*')
# ax3.scatter(2/WD11_P, h_c(WD11_m1, WD11_m2, WD11_d, 2/WD11_P, 4), s=70, marker='*')
# ax3.minorticks_on()
# ax3.tick_params(direction='in', length=10, width=1.5,
#                 which = 'major', top=True, right=True, pad=7)
# ax3.tick_params(direction='in', length=5, width=1.5,
#                 which = 'minor', top=True, right=True)
# ax3.set_xlim(1.01e-5, 2e-1)
# ax3.set_ylim(1e-22, 1e-18)
# ax3.set_xlabel(u'$\mathregular{Frequency/Hz}$')
# ax3.set_ylabel(u'$\mathregular{Characteristic\:\:\:\:Strain}$')
# ax3.set_box_aspect(1)
# ax3.legend(loc="lower right", fontsize="small")
# plt.gcf().set_size_inches(15, 15)
# plt.tight_layout()
# fig.savefig('LISA_curve_DNSBs.png', dpi=250, bbox_inches='tight')
# plt.show()
# =============================================================================

"""
PLOT OF THE MASS DEPENDENCY ON THE PSD
"""
# =============================================================================
# #CMAP
# DNSB_Mtot = []
# for i in range(0, List_len):
#     DNSB_Mtot.append(Vis_DNSB_Dat.at[i, 'Mass1 (m1)'] + Vis_DNSB_Dat.at[i, 'Mass2 (m2)']) 
# pos_color_log = []
# pos_color = []
# mass_color = []
# for i in range(0, List_len):
#     pos_color_log.append(math.log(r[i], 10))
#     pos_color.append(r[i])
#     mass_color.append(1)
# #Plotting
# fig, (ax9, ax10) = plt.subplots(ncols=2)
# fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2)
# im9 = ax9.scatter(DNSB_Mtot, DNSB_sqrt_S_h,
#                   cmap=cmocean.cm.thermal, c=pos_color, vmin = min(pos_color), vmax = 50,
#                   lw = 1, s = 2, alpha= 1)
# ax9_divider = make_axes_locatable(ax9)
# ax9.minorticks_on()
# ax9.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# ax9.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
#                 which = 'minor', top=True, right=True)
# ax9.set_box_aspect(1)
# ax9.set(xlabel=u'$\mathregular{M_{total}/M_\odot}$', ylabel=u'$\mathregular{Strain\:\:\:\:Sensitivity\:\:\:\:[\sqrt{Hz}^{-1}]}$')
# ax9.set_xlim(2.5, 5.05)
# ax9.set_ylim(2e-21, 4e-17)
# ax9.set_yscale('log')
# plt.gcf().set_size_inches(15, 15)
# cax9 = ax9_divider.append_axes("top", size="1.3%", pad=-6.59)
# cb9 = plt.colorbar(im9, cax=cax9, orientation='horizontal', fraction=0.0476, pad=0)
# cb9.ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# cb9.set_alpha(1)
# cb9.draw_all()
# cax9.xaxis.set_ticks_position("top")
# ax9.set_title(u'$\mathregular{r/kpc}$', y=1.105)
# 
# im10 = ax10.scatter(r, DNSB_sqrt_S_h,
#                   cmap='viridis_r', c=DNSB_Mtot, vmin = min(DNSB_Mtot), vmax = max(DNSB_Mtot),
#                   lw = 1, s = 2, alpha= 1)
# ax10_divider = make_axes_locatable(ax10)
# ax10.minorticks_on()
# ax10.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# ax10.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
#                 which = 'minor', top=True, right=True)
# ax10.set_box_aspect(1)
# ax10.set(xlabel=u'$\mathregular{r/kpc}$', ylabel=u'$\mathregular{Strain\:\:\:\:Sensitivity\:\:\:\:[\sqrt{Hz}^{-1}]}$')
# ax10.set_xlim(5e-2, 4e2)
# ax10.set_ylim(2e-21, 4e-17)
# ax10.loglog()
# ax10.invert_xaxis()
# plt.gcf().set_size_inches(15, 15)
# cax10 = ax10_divider.append_axes("top", size="1.3%", pad=-6.59)
# cb10 = plt.colorbar(im10, cax=cax10, orientation='horizontal', fraction=0.0476, pad=0)
# cb10.ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# cb10.set_alpha(1)
# cb10.draw_all()
# cax10.xaxis.set_ticks_position("top")
# ax10.set_title(u'$\mathregular{M_{total}/M_\odot}$', y=1.105)
# fig.savefig('Mass_distance.png', dpi=250, bbox_inches='tight')
# plt.show()
# =============================================================================

"""
PLOTTING THE RB IN THE LISA CURVE 
"""
# =============================================================================
# #Plotting only the resolvable binaries
# fig, (ax4, ax5) = plt.subplots(ncols=2)
# fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2)
# ax4.grid(which='both')
# ax4.set_axisbelow(True)
# # SIGNAL AND SENSITIVITY CURVES FOR THE ENTIRE OBSERVATION TIME
# ax4.loglog(F, h_s(F), 'b', linestyle='solid', linewidth=1, label='Lisa curve')
# ax4.scatter(1, 1, color='k', s=65, marker='*', label='Verification DWDs')
# ax4.scatter(1, 1, color='k', s=65, marker='^', label='Resolvable DNSBs')
# #Counting the resolvable binaries
# SNR_5_noBack = 0
# SNR_5_Back = 0
# for i in range(0, List_len):
#     if(DNSB_SNR_noBackground[i] > 5):
#         ax4.scatter(DNSB_f_gw[i], DNSB_h_c[i], color='k', lw = 0.6, s = 40, alpha = 0.75, marker='^')
#         SNR_5_noBack += 1
# #DWDB
# ax4.scatter(2/WD1_P, h_c(WD1_m1, WD1_m2, WD1_d, 2/WD1_P, 4), s=55, marker='*')
# ax4.scatter(2/WD2_P, h_c(WD2_m1, WD2_m2, WD2_d, 2/WD2_P, 4), s=55, marker='*')
# ax4.scatter(2/WD3_P, h_c(WD3_m1, WD3_m2, WD3_d, 2/WD3_P, 4), s=55, marker='*')
# ax4.scatter(2/WD4_P, h_c(WD4_m1, WD4_m2, WD4_d, 2/WD4_P, 4), s=55, marker='*')
# ax4.scatter(2/WD5_P, h_c(WD5_m1, WD5_m2, WD5_d, 2/WD5_P, 4), s=55, marker='*')
# ax4.scatter(2/WD6_P, h_c(WD6_m1, WD6_m2, WD6_d, 2/WD6_P, 4), s=55, marker='*')
# ax4.scatter(2/WD7_P, h_c(WD7_m1, WD7_m2, WD7_d, 2/WD7_P, 4), s=55, marker='*')
# ax4.scatter(2/WD8_P, h_c(WD8_m1, WD8_m2, WD8_d, 2/WD8_P, 4), s=55, marker='*')
# ax4.scatter(2/WD9_P, h_c(WD9_m1, WD9_m2, WD9_d, 2/WD9_P, 4), s=55, marker='*')
# ax4.scatter(2/WD10_P, h_c(WD10_m1, WD10_m2, WD10_d, 2/WD10_P, 4), s=55, marker='*')
# ax4.scatter(2/WD11_P, h_c(WD11_m1, WD11_m2, WD11_d, 2/WD11_P, 4), s=55, marker='*')
# ax4.minorticks_on()
# ax4.tick_params(direction='in', length=10, width=1.5,
#                 which = 'major', top=True, right=True, pad=7)
# ax4.tick_params(direction='in', length=5, width=1.5,
#                 which = 'minor', top=True, right=True)
# ax4.set_xlim(1.01e-4, 5e-2)
# ax4.set_ylim(6e-22, 7e-18)
# ax4.set_xlabel(u'$\mathregular{Frequency/Hz}$')
# ax4.set_ylabel(u'$\mathregular{Characteristic\:\:\:\:Strain}$')
# ax4.set_box_aspect(1)
# ax4.legend(loc="lower right", fontsize="small")
# plt.gcf().set_size_inches(15, 15)
# plt.tight_layout()
# 
# ax5.grid(which='both')
# ax5.set_axisbelow(True)
# # SIGNAL AND SENSITIVITY CURVES FOR THE ENTIRE OBSERVATION TIME
# ax5.loglog(F, h_s(F), 'b', linestyle='solid', linewidth=1, label='Lisa curve')
# ax5.plot(F_background, F3(F_background), 'r', linestyle='solid', linewidth=1, label='Background')
# ax5.scatter(1, 1, color='k', s=65, marker='*', label='Verification DWDs')
# ax5.scatter(1, 1, color='k', s=65, marker='^', label='Resolvable DNSBs')
# for i in range(0, List_len):
#     if(DNSB_SNR[i] > 5):
#         ax5.scatter(DNSB_f_gw[i], DNSB_h_c[i], color='k', lw = 0.6, s = 40, alpha = 0.75, marker='^')
#         SNR_5_Back += 1
# #DWDB
# ax5.scatter(2/WD1_P, h_c(WD1_m1, WD1_m2, WD1_d, 2/WD1_P, 4), s=55, marker='*')
# ax5.scatter(2/WD2_P, h_c(WD2_m1, WD2_m2, WD2_d, 2/WD2_P, 4), s=55, marker='*')
# ax5.scatter(2/WD3_P, h_c(WD3_m1, WD3_m2, WD3_d, 2/WD3_P, 4), s=55, marker='*')
# ax5.scatter(2/WD4_P, h_c(WD4_m1, WD4_m2, WD4_d, 2/WD4_P, 4), s=55, marker='*')
# ax5.scatter(2/WD5_P, h_c(WD5_m1, WD5_m2, WD5_d, 2/WD5_P, 4), s=55, marker='*')
# ax5.scatter(2/WD6_P, h_c(WD6_m1, WD6_m2, WD6_d, 2/WD6_P, 4), s=55, marker='*')
# ax5.scatter(2/WD7_P, h_c(WD7_m1, WD7_m2, WD7_d, 2/WD7_P, 4), s=55, marker='*')
# ax5.scatter(2/WD8_P, h_c(WD8_m1, WD8_m2, WD8_d, 2/WD8_P, 4), s=55, marker='*')
# ax5.scatter(2/WD9_P, h_c(WD9_m1, WD9_m2, WD9_d, 2/WD9_P, 4), s=55, marker='*')
# ax5.scatter(2/WD10_P, h_c(WD10_m1, WD10_m2, WD10_d, 2/WD10_P, 4), s=55, marker='*')
# ax5.scatter(2/WD11_P, h_c(WD11_m1, WD11_m2, WD11_d, 2/WD11_P, 4), s=55, marker='*')
# ax5.minorticks_on()
# ax5.tick_params(direction='in', length=10, width=1.5,
#                 which = 'major', top=True, right=True, pad=7)
# ax5.tick_params(direction='in', length=5, width=1.5,
#                 which = 'minor', top=True, right=True)
# ax5.set_xlim(1.01e-4, 5e-2)
# ax5.set_ylim(6e-22, 7e-18)
# ax5.set_xlabel(u'$\mathregular{Frequency/Hz}$')
# ax5.set_ylabel(u'$\mathregular{Characteristic\:\:\:\:Strain}$')
# ax5.set_box_aspect(1)
# ax5.legend(loc="lower right", fontsize="small")
# plt.gcf().set_size_inches(15, 15)
# plt.tight_layout()
# fig.savefig('LISA_curve_RB.png', dpi=250, bbox_inches='tight')
# plt.show()
# print('No Back: ' + str(SNR_5_noBack))
# print('')
# print('Back: ' + str(SNR_5_Back))
# =============================================================================

"""
PLOTTING THE RB IN THE TOP DOWN VIEW AND ELLIPTICAL PROJECTION
"""

#Plotting the positions of these resolvable binaries
# =============================================================================
# pos_color_log = []
# pos_color = []
# for i in range(0, List_len):
#     pos_color_log.append(math.log(r[i], 10))
#     pos_color.append(r[i])
# 
# fig, (ax9, ax10) = plt.subplots(ncols=2)
# fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2)
# im9 = ax9.scatter(l, b,
#                   cmap='gist_ncar_r', c=pos_color_log, vmin = min(pos_color_log), vmax = max(pos_color_log),
#                   lw = 0.5, s = 4, alpha= 0.4)
# for i in range(0, List_len):
#     if(DNSB_SNR[i] > 5):
#         ax9.scatter(l[i], b[i], s=80, marker='^')
# ax9.scatter(0, 0, s=120, marker='x', color='k', label = 'Galactic center')
# ax9_divider = make_axes_locatable(ax9)
# ax9.minorticks_on()
# ax9.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# ax9.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
#                 which = 'minor', top=True, right=True)
# ax9.set_box_aspect(1)
# ax9.set(xlabel=u'$\mathregular{l/rad}$', ylabel=u'$\mathregular{b/rad}$')
# ax9.set_xlim(-1*math.pi, math.pi)
# ax9.set_ylim(-math.pi/2, math.pi/2)
# plt.gcf().set_size_inches(15, 15)
# cax9 = ax9_divider.append_axes("top", size="1.3%", pad=-6.59)
# cb9 = plt.colorbar(im9, cax=cax9, orientation='horizontal', fraction=0.0476, pad=0)
# cb9.ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# cb9.set_alpha(1)
# cb9.draw_all()
# ax9.legend(loc="lower left")
# cax9.xaxis.set_ticks_position("top")
# ax9.set_title('log(r)/kpc', y=1.105)
# 
# im10 = ax10.scatter(X_sol, Y,
#                   cmap='tab20c', c=pos_color, vmin = min(pos_color), vmax = 50,
#                   lw = 1, s = 4, alpha= 1)
# for i in range(0, List_len):
#     if(DNSB_SNR[i] > 5):
#         ax10.scatter(X_sol[i], Y[i], s=80, marker='^')
# ax10.scatter(-8.33, 0, s=120, marker='x', color='k', label = 'Galactic center')
# ax10.scatter(0, 0, s=200, marker='*', color='k', label = 'Solar System')
# ax10_divider = make_axes_locatable(ax10)
# ax10.minorticks_on()
# ax10.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# ax10.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
#                 which = 'minor', top=True, right=True)
# ax10.set_box_aspect(1)
# ax10.set(xlabel=u'$\mathregular{x/kpc}$', ylabel=u'$\mathregular{y/kpc}$')
# ax10.set_xlim(-10-8.1, 20-8.1)
# ax10.set_ylim(-15, 15)
# plt.gcf().set_size_inches(15, 15)
# cax10 = ax10_divider.append_axes("top", size="1.3%", pad=-6.59)
# cb10 = plt.colorbar(im10, cax=cax10, orientation='horizontal', fraction=0.0476, pad=0)
# cb10.ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# cb10.set_alpha(1)
# cb10.draw_all()
# ax10.legend(loc="lower left")
# cax10.xaxis.set_ticks_position("top")
# ax10.set_title('r/kpc', y=1.105)
# plt.show()
# 
# =============================================================================

"""
PLOTTING THE ELLIPTICAL PROJECTION
"""
# =============================================================================
# #CMAP
# pos_color_log = []
# pos_color = []
# for i in range(0, List_len):
#     pos_color_log.append(math.log(r[i], 10))
#     pos_color.append(r[i])
# #Projection b and l range
# circ_long = np.linspace(-np.pi, np.pi, 13)[1:-1]
# circ_lat = np.linspace(-np.pi / 2, np.pi / 2, 7)[1:-1]
# radius = 10 * np.pi / 180.
# #Plotting
# fig = plt.figure(figsize=(5, 4))
# plt.subplots_adjust(hspace=0, wspace=0.12,
#                     left=0.08, right=0.95,
#                     bottom=0.05, top=1.0)
# ax = plt.subplot(221, projection='aitoff')
# ax.grid(True, which='both', alpha = 0.7)
# ax.set_axisbelow(True)
# im = ax.scatter(l, b,
#                   cmap='gnuplot2', c=pos_color_log, vmin = -0.5,
#                   vmax = 2.5,
#                   lw = 0.5, s = 4.5, alpha= 0.45)
# #Plotting the resolvable DNSBs
# # =============================================================================
# # for i in range(0, List_len):
# #     if(DNSB_SNR[i] > 5):
# #         ax.scatter(l[i], b[i], s=90, marker='^')
# # ax.scatter(220, 220, s=120, marker='^', color='k', label = 'Resolvable DNSBs')
# # =============================================================================
# ax.scatter(0, 0, s=120, marker='x', color='k', label = 'Galactic center')
# ax_divider = make_axes_locatable(ax)
# ax.set(xlabel=u'$\mathregular{l/deg}$', ylabel=u'$\mathregular{b/deg}$')
# plt.gcf().set_size_inches(45, 45)
# cb = plt.colorbar(im, orientation='vertical', fraction=0.022, pad=0)
# cb.ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# cb.set_alpha(1)
# cb.draw_all()
# cb.set_label(r'log(r)/kpc')
# ax.legend(loc=(0.825, 0.02))
# fig.savefig('Elliptical_projection_DNSB.png', dpi=250, bbox_inches='tight')
# plt.show()
# =============================================================================

"""
PLOTTING THE TOP-DOWN VIEW
"""
# =============================================================================
# #CMAP
# pos_color_log = []
# pos_color = []
# for i in range(0, List_len):
#     pos_color_log.append(math.log(r[i], 10))
#     pos_color.append(r[i])
# #Plotting
# fig, ((ax11, ax12), (ax9, ax10)) = plt.subplots(2, 2)
# fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2, hspace=-0.5)
# 
# im11 = ax11.scatter(X_sol, Z,
#                   cmap='tab20c', c=pos_color, vmin = min(pos_color), vmax = 50,
#                   lw = 1, s = 2, alpha= 1)
# for i in range(0, List_len):
#     if(DNSB_SNR_noBackground[i] > 5):
#         ax11.scatter(X_sol[i], Z[i], s=80, marker='^')
# ax11.scatter(-8.33, 0, s=120, marker='x', color='k', label = 'Galactic center')
# ax11.scatter(0, 0, s=200, marker='*', color='k', label = 'Solar System')
# ax11.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# ax11.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
#                 which = 'minor', top=True, right=True)
# ax11_divider = make_axes_locatable(ax11)
# cax11 = ax11_divider.append_axes("top", size="9.45%", pad=0)
# cb11 = plt.colorbar(im11, cax=cax11, orientation='horizontal', fraction=0.0476, pad=0)
# cb11.ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# cb11.set_alpha(1)
# cb11.draw_all()
# ax11.axes.get_xaxis().set_ticklabels([])
# ax11.minorticks_on()
# ax11.set_xlim(-50.1-8.1, 50.1-8.1)
# ax11.set_ylim(-15, 15)
# ax11.set_aspect(1)
# cax11.xaxis.set_ticks_position("top")
# ax11.set(ylabel=u'$\mathregular{z/kpc}$')
# ax11.set_title(u'$\mathregular{r/kpc}$', y=1.25)
# 
# im12 = ax12.scatter(X_sol, Z,
#                   cmap='tab20c', c=pos_color, vmin = min(pos_color), vmax = 50,
#                   lw = 1, s = 2, alpha= 1)
# for i in range(0, List_len):
#     if(DNSB_SNR_noBackground[i] > 5):
#         ax12.scatter(X_sol[i], Z[i], s=80, marker='^')
# ax12.scatter(-8.33, 0, s=120, marker='x', color='k', label = 'Galactic center')
# ax12.scatter(0, 0, s=200, marker='*', color='k', label = 'Solar System')
# ax12.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# ax12.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
#                 which = 'minor', top=True, right=True)
# ax12_divider = make_axes_locatable(ax12)
# cax12 = ax12_divider.append_axes("top", size="2.63%", pad=-2.41)
# cb12 = plt.colorbar(im12, cax=cax12, orientation='horizontal', fraction=0.0476, pad=0)
# cb12.ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# cb12.set_alpha(1)
# cb12.draw_all()
# ax12.axes.get_xaxis().set_ticklabels([])
# ax12.minorticks_on()
# ax12.set_xlim(-10.1-8.1, 18.1-8.1)
# ax12.set_ylim(-15, 15)
# ax12.set_aspect(0.285)
# cax12.xaxis.set_ticks_position("top")
# ax12.set(ylabel=u'$\mathregular{z/kpc}$')
# ax12.set_title(u'$\mathregular{r/kpc}$', y=1.25)
# ax12_divider.append_axes("top", size="0%", pad=2.4)
# 
# ax9.scatter(X_sol, Y,
#                   cmap='tab20c', c=pos_color, vmin = min(pos_color), vmax = 50,
#                   lw = 1, s = 2, alpha= 1)
# for i in range(0, List_len):
#     if(DNSB_SNR_noBackground[i] > 5):
#         ax9.scatter(X_sol[i], Y[i], s=80, marker='^')
# ax9.scatter(-8.33, 0, s=120, marker='x', color='k', label = 'Galactic center')
# ax9.scatter(0, 0, s=200, marker='*', color='k', label = 'Solar System')
# ax9.scatter(220, 220, s=120, marker='^', color='k', label = 'Resolvable DNSBs')
# ax9.minorticks_on()
# ax9.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# ax9.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
#                 which = 'minor', top=True, right=True)
# ax9.set_box_aspect(1)
# ax9.set(xlabel=u'$\mathregular{x/kpc}$', ylabel=u'$\mathregular{y/kpc}$')
# ax9.set_xlim(-50-8.1, 50-8.1)
# ax9.set_ylim(-50, 50)
# plt.gcf().set_size_inches(15, 15)
# ax9.legend(loc="lower left", fontsize="small")
# 
# ax10.scatter(X_sol, Y,
#                   cmap='tab20c', c=pos_color, vmin = min(pos_color), vmax = 50,
#                   lw = 1, s = 2, alpha= 1)
# for i in range(0, List_len):
#     if(DNSB_SNR_noBackground[i] > 5):
#         ax10.scatter(X_sol[i], Y[i], s=80, marker='^')
# ax10.scatter(-8.33, 0, s=120, marker='x', color='k', label = 'Galactic center')
# ax10.scatter(0, 0, s=200, marker='*', color='k', label = 'Solar System')
# ax10.scatter(220, 220, s=120, marker='^', color='k', label = 'Resolvable DNSBs')
# ax10.minorticks_on()
# ax10.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
#                 which = 'major', top=True, right=True, pad=7)
# ax10.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
#                 which = 'minor', top=True, right=True)
# ax10.set_box_aspect(1)
# ax10.set(xlabel=u'$\mathregular{x/kpc}$', ylabel=u'$\mathregular{y/kpc}$')
# ax10.set_xlim(-10.1-8.1, 18.1-8.1)
# ax10.set_ylim(-15, 15)
# plt.gcf().set_size_inches(15, 15)
# ax10.legend(loc="lower left", fontsize="small")
# fig.savefig('RB_TopDown.png', dpi=300, bbox_inches='tight')
# plt.show()
# =============================================================================
