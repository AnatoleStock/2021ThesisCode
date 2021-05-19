# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 07:42:06 2021

@author: Anatole Storck
"""

"""
This file aims to to scale the population by a certain factor, then calculating
the final positions using the Galactic_Orbit_Integrator.
"""

import random
import numpy as np
import math
import pandas as pd
import Binary_Orbit_Integrator as bi
import Galactic_Orbit_Integrator as gal
import Galactic_Stelar_Dispersion as de
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import cmocean

#Plot formating stuff
plt.rcParams['axes.linewidth'] = 2.15
plt.rcParams["font.family"] = 'Courier New'
plt.rcParams.update({'font.size': 16})
cmap = cmocean.cm.solar
newcmap = cmocean.tools.crop_by_percent(cmap, 7, which='max', N=None)
color = newcmap
#Importing the data
NSNS_data = pd.read_csv(r'C:\Users\storc\OneDrive\Thesis Work\Data\nsns.doubleAcc.DAT', sep=' ')
Master_t1t2 = pd.read_csv(r'C:\Users\storc\OneDrive\Thesis Work\Data\Master_All_DNSB_t1t2.csv')
scale_factor = 266
#Elements in order: a, e, m1, m2, Vx, Vy, Vz
NSNS_full = ([], [], [], [], [], [], [])
for k in range(0, scale_factor):
    for i in range(0, 2857):
        #Exporting the data to lists
        NSNS_full[0].append(NSNS_data.at[i, '#a/Rsun'])
        NSNS_full[1].append(NSNS_data.at[i, 'e'])
        NSNS_full[2].append(NSNS_data.at[i, 'm1/Msun'])
        NSNS_full[3].append(NSNS_data.at[i, 'm2/Msun'])
        NSNS_full[4].append(NSNS_data.at[i, 'vx/kms'])
        NSNS_full[5].append(NSNS_data.at[i, 'vy/kms'])
        NSNS_full[6].append(NSNS_data.at[i, 'vz/kms'])
#Elements in order: a, e, m1, m2, t_birth, t1, t2
NSNS_visible_binary = ([], [], [], [], [], [], [])
#Elements in order: R, Phi, Vx, Vy, Vz
NSNS_visible_galactic = ([], [], [], [], [])
#Elements at present day: a, e, P
NSNS_visible_present = ([], [], [])
#Since merger time does not change for different birth times
#we can calculate the binary integrator once and apply it to
#the rest of the iterations
t1_local = Master_t1t2['T1']
t2_local = Master_t1t2['T2']
#Creating basis for the 2857 unique binary configurations from Ross' data
# =============================================================================
# #NO LONGER NEEDED WITH THE MASTER FILE FOR ALL T1s and T2s
# for i in range(0, 2857):
#         (a, e, T, t1, t2) = bi.RK4(NSNS_full[0][i], NSNS_full[1][i],
#                                NSNS_full[2][i], NSNS_full[3][i], 0)
#         t1_local.append(t1)
#         t2_local.append(t2)
#         print ('binary: ' + str(i))
#t1 and t2 are local since the binary integrator starts the system at 0 yrs
# =============================================================================
t1_abs = []
t2_abs = []
#counting the amount of visible binaries
count = 0
for k in range(0, scale_factor):
    for i in range(0, 2857):
        #Assigning a birth time (-10e3 to 0e3 Myrs) for each binary
        t_birth = -1*random.random()*(10e3)
        #Changing t1 and t2 to our reference time of 0 yrs at present day
        t1_abs = t_birth + t1_local[i]
        t2_abs = t_birth + t2_local[i]
        if(t1_abs < 0 < t2_abs):
            count += 1
            #Doing the integrator to find t_present (since t_birth shifts it)
            (a, e, T, t1, t2) = bi.RK4(NSNS_full[0][i], NSNS_full[1][i],
                                NSNS_full[2][i], NSNS_full[3][i], t_birth)
            present_index = T.index(min(T, key=abs)) #We want the index at t_present
            NSNS_visible_present[0].append(a[present_index])
            NSNS_visible_present[1].append(e[present_index])
            NSNS_visible_present[2].append(bi.P(a[present_index],
                                                NSNS_data.at[i, 'm1/Msun'],
                                                NSNS_data.at[i, 'm2/Msun']))
            NSNS_visible_binary[0].append(NSNS_data.at[i, '#a/Rsun'])
            NSNS_visible_binary[1].append(NSNS_data.at[i, 'e'])
            NSNS_visible_binary[2].append(NSNS_data.at[i, 'm1/Msun'])
            NSNS_visible_binary[3].append(NSNS_data.at[i, 'm2/Msun'])
            NSNS_visible_binary[4].append(t_birth)
            NSNS_visible_binary[5].append(t1_abs)
            NSNS_visible_binary[6].append(t2_abs)
            NSNS_visible_galactic[2].append(NSNS_data.at[i, 'vx/kms'])
            NSNS_visible_galactic[3].append(NSNS_data.at[i, 'vy/kms'])
            NSNS_visible_galactic[4].append(NSNS_data.at[i, 'vz/kms'])
            print(count)
            print('')
#Creating random positions and adding angle and R to list
R_random = de.Stellar_distribution(count)[0]
for i in range(0, count):
    NSNS_visible_galactic[0].append(R_random[i])
    NSNS_visible_galactic[1].append(2*math.pi*random.random())

"""Exporting the calculated values"""
#Creating a permanent data file for randomly generated values
Permanent_DNSB_Values = {'Initial Position (R)': NSNS_visible_galactic[0],
                         'Initial Angle (Phi)': NSNS_visible_galactic[1],
                         'Birth Time (t_birth)': NSNS_visible_binary[4],
                         'First Visible (t1)': NSNS_visible_binary[5],
                         'Last Visible (t2)': NSNS_visible_binary[6],
                         'Semi-major at present (a)': NSNS_visible_present[0],
                         'Ecc at present (e)': NSNS_visible_present[1],
                         'Period at present (P)': NSNS_visible_present[2]}
df_Values = pd.DataFrame(Permanent_DNSB_Values, columns= ['Initial Position (R)',
                                                   'Initial Angle (Phi)',
                                                   'Birth Time (t_birth)',
                                                   'First Visible (t1)',
                                                   'Last Visible (t2)',
                                                   'Semi-major at present (a)',
                                                   'Ecc at present (e)',
                                                   'Period at present (P)'])
df_Values.to_csv (r'C:\Users\storc\Desktop\Visible_DNSB_Values(MergerRate266).csv',
                  index = False, header=True)

#Creating a permanent data file for elements of binaries that are visible
Visible_DNSB_Data = {'Semi-major (a)': NSNS_visible_binary[0],
                     'Ecc (e)': NSNS_visible_binary[1],
                     'Mass1 (m1)': NSNS_visible_binary[2],
                     'Mass2 (m2)': NSNS_visible_binary[3],
                     'Kick VelocityX (Vx)': NSNS_visible_galactic[2],
                     'Kick VelocityY (Vy)': NSNS_visible_galactic[3],
                     'Kick VelocityZ (Vz)': NSNS_visible_galactic[4]}
df_Data = pd.DataFrame(Visible_DNSB_Data, columns= ['Semi-major (a)',
                                                   'Ecc (e)',
                                                   'Mass1 (m1)',
                                                   'Mass2 (m2)',
                                                   'Kick VelocityX (Vx)',
                                                   'Kick VelocityY (Vy)',
                                                   'Kick VelocityZ (Vz)'])
df_Data.to_csv (r'C:\Users\storc\Desktop\Visible_DNSB_Data(MergerRate266).csv', 
                  index = False, header=True)

"""Plotting data related to entering and leaving the LISA band"""
fig, (ax3, ax4) = plt.subplots(ncols=2)
fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2)
cmap_r = cmocean.cm.solar_r
newcmap_r = cmocean.tools.crop_by_percent(cmap_r, 7, which='min', N=None)
color_r = newcmap_r
t1_color = []
t2_color = []
for i in range(0, len(NSNS_visible_binary[5])):
    if (NSNS_visible_binary[5][i]-NSNS_visible_binary[4][i]) == 0:
        t1_color.append(-1)
    elif (math.log(abs(NSNS_visible_binary[5][i]-NSNS_visible_binary[4][i]), 10)) < -1:
        t1_color.append(-1)
    else:
        t1_color.append(math.log(abs(NSNS_visible_binary[5][i]-NSNS_visible_binary[4][i]), 10))
    t2_color.append(math.log(NSNS_visible_binary[6][i]-NSNS_visible_binary[4][i], 10))

im3 = ax3.scatter(NSNS_visible_binary[0], NSNS_visible_binary[1], cmap=color_r,
            c=t1_color,vmin = min(t1_color), vmax = max(t1_color),
            lw = 1.2, s = 45)
ax3_divider = make_axes_locatable(ax3)
ax3.minorticks_on()
ax3.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax3.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax3.set_box_aspect(1)
ax3.set(xlabel=u'$\mathregular{a/R_\odot}$', ylabel=u'$\mathregular{e}$')
ax3.set_xscale('log')
ax3.set_ylim(0, 1)
plt.gcf().set_size_inches(15, 15)
cax3 = ax3_divider.append_axes("top", size="1.3%", pad=-6.59)
cb3 = plt.colorbar(im3, cax=cax3, orientation='horizontal', fraction=0.0476, pad=0)
cb3.ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
cax3.xaxis.set_ticks_position("top")

im4 = ax4.scatter(NSNS_visible_binary[0], NSNS_visible_binary[1], cmap=color,
            c=t2_color,vmin = min(t2_color), vmax = max(t2_color),
            lw = 1.2, s = 45)
ax4_divider = make_axes_locatable(ax4)
ax4.minorticks_on()
ax4.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
ax4.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                which = 'minor', top=True, right=True)
ax4.set_box_aspect(1)
ax4.set(xlabel=u'$\mathregular{a/R_\odot}$', ylabel=u'$\mathregular{e}$')
ax4.set_xscale('log')
ax4.set_ylim(0, 1)
plt.gcf().set_size_inches(15, 15)
cax4 = ax4_divider.append_axes("top", size="1.3%", pad=-6.59)
cb4 = plt.colorbar(im4, cax=cax4, orientation='horizontal', fraction=0.0476, pad=0)
cb4.ax.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                which = 'major', top=True, right=True, pad=7)
cax4.xaxis.set_ticks_position("top")


#Importing the just-exported data
Vis_DNSB_Val = pd.read_csv(r'C:\Users\storc\Desktop\Visible_DNSB_Values(MergerRate266).csv')
Vis_DNSB_Dat = pd.read_csv(r'C:\Users\storc\Desktop\Visible_DNSB_Data(MergerRate266).csv')
List_len = len(Vis_DNSB_Val['Initial Position (R)'])

"""Running the Galactic integrator to get final positions"""
#counting the orbits
COUNT = 0
Final_Position = ([], [], [])
for i in range(0, List_len):
    (P, V, E, T) = gal.Integrator(Vis_DNSB_Val.at[i, 'Initial Position (R)'],
                                  0, Vis_DNSB_Val.at[i, 'Initial Angle (Phi)'],
                                  Vis_DNSB_Dat.at[i, 'Kick VelocityX (Vx)'],
                                  Vis_DNSB_Dat.at[i, 'Kick VelocityZ (Vz)'],
                                  (Vis_DNSB_Dat.at[i, 'Kick VelocityY (Vy)']
                                  + gal.Vc(Vis_DNSB_Val.at[i, 'Initial Position (R)'], 0)),
                                  Vis_DNSB_Val.at[i, 'Birth Time (t_birth)'])
    #To cartesian for top down view and going to kpc
    P[0][:] = [x / 1000 for x in P[0]]
    P[1][:] = [x / 1000 for x in P[1]]
    if (COUNT % 1000 == 0):
        X = P[0]*np.cos(P[2])
        Y = P[0]*np.sin(P[2])
        fig, (ax5, ax6) = plt.subplots(ncols=2)
        fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.2)
        #Plotting z against R
        ax5.scatter(P[0], P[1], color='k', lw = 0.5, s = 0.5, alpha = 0.5)
        ax5.minorticks_on()
        ax5.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                        which = 'major', top=True, right=True, pad=7)
        ax5.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                        which = 'minor', top=True, right=True)
        ax5.set_box_aspect(1)
        ax5.set(xlabel=u'$\mathregular{R/kpc}$', ylabel=u'$\mathregular{z/kpc}$')
        plt.gcf().set_size_inches(15, 15)
        #ax5.set_title('Binary #' + str(i+1))
        #Plotting phi
        ax6.scatter(X, Y, color='k', lw = 0.5, s = 0.5, alpha = 0.5)
        ax6.minorticks_on()
        ax6.tick_params(direction='in', length=10, width=1.5, grid_alpha=0.55,
                        which = 'major', top=True, right=True, pad=7)
        ax6.tick_params(direction='in', length=5, width=1.5, grid_alpha=0.55,
                        which = 'minor', top=True, right=True)
        ax6.set_box_aspect(1)
        ax6.set(xlabel=u'$\mathregular{x/kpc}$', ylabel=u'$\mathregular{y/kpc}$')
        plt.gcf().set_size_inches(15, 15)
        plt.show()
    Final_Position[0].append(P[0][-1])
    Final_Position[1].append(P[1][-1])
    Final_Position[2].append(P[2][-1])
    COUNT += 1
    print(COUNT)

#CREATING FINAL POSITIONS FILE
Final_P = {'Final Position (R)': Final_Position[0],
           'Final Position (z)': Final_Position[1],
           'Final Angle (Phi)': Final_Position[2]}
df_Data = pd.DataFrame(Final_P, columns= ['Final Position (R)',
                                          'Final Position (z)',
                                          'Final Angle (Phi)'])
df_Data.to_csv (r'C:\Users\storc\Desktop\Final_Position(MergerRate266).csv', 
                  index = False, header=True)