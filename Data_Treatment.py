# -*- coding: utf-8 -*-
"""
Created on Wed May  5 12:46:36 2021

@author: Anatole Storck
"""

"""
This file is to treat the data from any failings that the galactic integrator
might have procured, as well as removing DNSBs too far away from LISA.
"""

import numpy as np
import pandas as pd

"""Removing anomaly (negative R) from the list and putting upper bound (200 kpc) for DNSB population"""
#Importing Files (Dat comes from Ross' data files while Val are calculated and given values)
Vis_DNSB_Val = pd.read_csv(r'Data\Visible_DNSB_Values(300).csv')
Vis_DNSB_Dat = pd.read_csv(r'Data\Visible_DNSB_Data(300).csv')
#Fin are final positions R and z
Vis_DNSB_Fin = pd.read_csv(r'Data\Final_Position(300).csv')
#length of the list
List_len = len(Vis_DNSB_Val['Initial Position (R)'])
print(List_len)

drop_list1 = []

for i in range(0, List_len):
    if(Vis_DNSB_Fin.at[i, 'Final Position (R)'] < 0):
        drop_list1.append(i)
    elif(np.sqrt(np.square(Vis_DNSB_Fin.at[i, 'Final Position (R)'])
                 + np.square(Vis_DNSB_Fin.at[i, 'Final Position (z)'])) > 778):
        drop_list1.append(i)
        
Vis_DNSB_Fin.drop(labels=drop_list1, axis=0, inplace=True)
Vis_DNSB_Dat.drop(labels=drop_list1, axis=0, inplace=True)
Vis_DNSB_Val.drop(labels=drop_list1, axis=0, inplace=True)
    
List_len1 = len(Vis_DNSB_Val['Initial Position (R)'])
print(List_len1)
#%%
"""Exporting the treated data file"""
Vis_DNSB_Fin.to_csv (r'Desktop\Final_Position(300_2_Andromeda).csv', 
                  index = False, header=True)
Vis_DNSB_Dat.to_csv (r'Desktop\Visible_DNSB_Data(300_2_Andromeda).csv', 
                  index = False, header=True)
Vis_DNSB_Val.to_csv (r'Desktop\Visible_DNSB_Values(300_2_Andromeda).csv', 
                  index = False, header=True)
    