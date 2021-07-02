# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 14:02:18 2021

@author: Anatole Storck
"""

"""
This file takes the treated DNSB population data and creates several plots
"""

import math
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import cmocean

"""Preamble"""
#Plotting formating and coloring
plt.rcParams['axes.linewidth'] = 2.15
plt.rcParams["font.family"] = 'Courier New'
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'axes.titlesize': 16})
cmap = cmocean.cm.thermal
newcmap = cmocean.tools.crop_by_percent(cmap, 4, which='max', N=None)
color = newcmap
cmap_r = cmocean.cm.haline_r
newcmap_r = cmocean.tools.crop_by_percent(cmap_r, 15, which='min', N=None)
color_r = newcmap_r
#Importing Files (Dat comes from Ross' data files while Val are calculated and given values)
Vis_DNSB_Val = pd.read_csv(r'Data\Visible_DNSB_Values(MergerRate266).csv')
Vis_DNSB_Dat = pd.read_csv(r'Data\Visible_DNSB_Data(MergerRate266).csv')
#Fin are final positions R and z
Vis_DNSB_Fin = pd.read_csv(r'Data\Final_Position(MergerRate266).csv')
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
#%%
"""PLOTTING t1 AND t2 AS A FUNCTION OF a AND e"""
Tot_DNSB_Dat = pd.read_csv(r'Data\nsns.doubleAcc.Dat', sep=' ')
Unique_a = set(Vis_DNSB_Dat['Semi-major (a)'])
#Creating the CMAP range
t1_color = []
t2_color = []
for i in range(0, List_len):
    t1_color.append(math.log(max(abs(Vis_DNSB_Val.at[i, 'Birth Time (t_birth)']-
                                     Vis_DNSB_Val.at[i, 'First Visible (t1)']), 0.11), 10))
    t2_color.append(math.log(Vis_DNSB_Val.at[i, 'Last Visible (t2)'], 10))
#Plotting
fig, (ax1, ax2) = plt.subplots(ncols=2)
fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2)
ax1.scatter(Tot_DNSB_Dat['#a/Rsun'], Tot_DNSB_Dat['e'], c='gray', lw = 1, s = 5, alpha= 0.4)
im1 = ax1.scatter(Vis_DNSB_Dat['Semi-major (a)'], Vis_DNSB_Dat['Ecc (e)'],
                  cmap=color_r, c=t1_color,vmin = min(t1_color), vmax = max(t1_color),
                  lw = 1, s = 8, alpha= 1)
ax1_divider = make_axes_locatable(ax1)
ax1.minorticks_on()
ax1.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax1.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax1.set_box_aspect(1)
ax1.set(xlabel=u'$\mathregular{a/R_\odot}$', ylabel=u'$\mathregular{e}$')
ax1.set_xscale('log')
ax1.set_xlim(2e-2, 1e4)
ax1.set_ylim(0, 1)
plt.gcf().set_size_inches(15, 15)
cax1 = ax1_divider.append_axes("top", size="1.3%", pad=-6.59)
cb1 = plt.colorbar(im1, cax=cax1, orientation='horizontal', fraction=0.0476, pad=0)
cb1.ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
cb1.set_alpha(1)
cb1.draw_all()
cax1.xaxis.set_ticks_position("top")
ax1.set_title(u'$\mathregular{log(t_1 - t_{birth})/Myr}$', y=1.105)

ax2.scatter(Tot_DNSB_Dat['#a/Rsun'], Tot_DNSB_Dat['e'], c='gray', lw = 1, s = 5, alpha= 0.4)
im2 = ax2.scatter(Vis_DNSB_Dat['Semi-major (a)'], Vis_DNSB_Dat['Ecc (e)'],
                  cmap=color, c=t2_color,vmin = min(t2_color), vmax = max(t2_color),
                  lw = 1, s = 8, alpha= 1)
ax2_divider = make_axes_locatable(ax2)
ax2.minorticks_on()
ax2.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax2.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax2.set_box_aspect(1)
ax2.set(xlabel=u'$\mathregular{a/R_\odot}$', ylabel=u'$\mathregular{e}$')
ax2.set_xscale('log')
ax2.set_xlim(2e-2, 1e4)
ax2.set_ylim(0, 1)
plt.gcf().set_size_inches(15, 15)
cax2 = ax2_divider.append_axes("top", size="1.3%", pad=-6.59)
cb2 = plt.colorbar(im2, cax=cax2, orientation='horizontal', fraction=0.0476, pad=0)
cb2.ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
cb2.set_alpha(1)
cb2.draw_all()
cax2.xaxis.set_ticks_position("top")
ax2.set_title(u'$\mathregular{log(t_2)/Myr}$', y=1.105)
#fig.savefig('DNSB_t1t2.png', dpi=250, bbox_inches='tight')
plt.show()
#%%
"""PLOTING FINAL POSITIONS"""
z_initial = []
R_tbirth_color = []
for i in range (0, List_len):
    z_initial.append(0)#assumption is all binaries start in galactic plane
    R_tbirth_color.append(Vis_DNSB_Val.at[i, 'Initial Position (R)']/1000)
#Changing from pc to kpc
R_0 = copy.deepcopy(Vis_DNSB_Val['Initial Position (R)'])
R_f = copy.deepcopy(Vis_DNSB_Fin['Final Position (R)'])
z = copy.deepcopy(Vis_DNSB_Fin['Final Position (z)'])
number = 1 #since final R and z are given in kpc already
R_0[:] = [x / 1000 for x in R_0]
R_f[:] = [x / number for x in R_f]
z[:] = [x / number for x in z]
#Plotting
fig, (ax7, ax8) = plt.subplots(ncols=2)
fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2)
im7 = ax7.scatter(R_0, z_initial, 
                  cmap='plasma_r', c=R_tbirth_color, vmin = min(R_tbirth_color),
                  vmax = max(R_tbirth_color), lw = 1.2, s = 7)
ax7_divider = make_axes_locatable(ax7)
ax7.minorticks_on()
ax7.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax7.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax7.set_xlim(0, 75)
ax7.set_ylim(-20, 20)
ax7.set_box_aspect(1)
ax7.set(xlabel=u'$\mathregular{R/kpc}$', ylabel=u'$\mathregular{z/kpc}$')
plt.gcf().set_size_inches(15, 15)
cax7 = ax7_divider.append_axes("top", size="1.3%", pad=-6.59)
cb7 = plt.colorbar(im7, cax=cax7, orientation='horizontal', fraction=0.0476, pad=0)
cb7.ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
cax7.xaxis.set_ticks_position("top")
ax7.set_title(u'$\mathregular{R/kpc}$', y=1.105)

im8 = ax8.scatter(R_f, z,
                  cmap='plasma_r', c=R_tbirth_color, vmin = min(R_tbirth_color),
                  vmax = max(R_tbirth_color), lw = 1, s = 5, alpha=0.75)
ax8_divider = make_axes_locatable(ax8)
ax8.minorticks_on()
ax8.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax8.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax8.set_xlim(0, 75)
ax8.set_ylim(-20, 20)
ax8.set_box_aspect(1)
ax8.set(xlabel=u'$\mathregular{R/kpc}$', ylabel=u'$\mathregular{z/kpc}$')
plt.gcf().set_size_inches(15, 15)
cax8 = ax8_divider.append_axes("top", size="1.3%", pad=-6.59)
cb8 = plt.colorbar(im8, cax=cax8, orientation='horizontal', fraction=0.0476, pad=0)
cb8.ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
cb8.set_alpha(1)
cb8.draw_all()
cax8.xaxis.set_ticks_position("top")
ax8.set_title(u'$\mathregular{R/kpc}$', y=1.105)
#fig.savefig('DNSB_dispersion.png', dpi=250, bbox_inches='tight')
plt.show()
#%%
"""PLOTTING KICK VS MERGER TIME"""
V_kick = []
t_merge = Vis_DNSB_Val['Last Visible (t2)'] - Vis_DNSB_Val['First Visible (t1)']
for i in range(0, List_len):
    V_kick.append(math.sqrt(Vis_DNSB_Dat.at[i, 'Kick VelocityX (Vx)']**2 +
                            Vis_DNSB_Dat.at[i, 'Kick VelocityY (Vy)']**2 +
                            Vis_DNSB_Dat.at[i, 'Kick VelocityZ (Vz)']**2))
    
fig, ax3 = plt.subplots(ncols=1)
fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2)
ax3.scatter(V_kick, t_merge, color='k', lw = 0.5, s = 1, alpha = 0.5)
ax3.minorticks_on()
ax3.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax3.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
#ax3.set_xlim(0, 300)
ax3.set_yscale('log')
ax3.set_box_aspect(1)
ax3.set(xlabel=u'$\mathregular{v_{kick}/(kms^{-1})}$',
        ylabel=u'$\mathregular{t_{merge}/Myr}$')
plt.gcf().set_size_inches(7.5, 7.5)
plt.show()
