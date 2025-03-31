#!/usr/bin/env python
# coding: utf-8

# In[1]:


#python task 1 convert UBVRI to ugriz
#TJ define all the properties found from the standard star link
name = "PG1633+099A"
ra = ((16)+(35/60)+(26/3600))*15
dec = (9)+(47/60)+(53/3600)
V = 15.256
BmV = 0.873 #TJ (B - V)
UmB = 0.320 #TJ (U - B)
VmR = 0.505 #TJ (V - R)
RmI = 0.511 #TJ (R - I)
#TJ define properties from SDSS
SDSS_g = 15.7
SDSS_u = 17.30
SDSS_z = 14.55
#TJ start using equations for transformation.
umg = 1.28*(UmB) + 1.13
gmr = 1.02*(BmV) - 0.22
rmz = 1.72*(RmI) - 0.41
g = V + 0.60*(BmV) - 0.12
z = g - gmr - rmz #TJ g - (g - r) - (r - z) = z
print(f'computed g-mag : {g}, compared to SDSS g-mag : {SDSS_g}')
print(f'computed z-mag : {z}, compared to SDSS z-mag : {SDSS_z}')
print('well.... close enough?')


# In[2]:


#python task 2 
from astropy.table import Table
import numpy as np
path = '/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-240p005-250p010.fits'
table = Table.read(path)
distance = np.sqrt((table["RA"] - ra)**2 + (table["DEC"] - dec)**2)
nearest_index = np.argmin(distance)
star = table[nearest_index]
star_g_mag = round(22.5 - 2.5*np.log10(star["FLUX_G"]), 2)
star_r_mag = round(22.5 - 2.5*np.log10(star["FLUX_R"]), 2)
star_z_mag = round(22.5 - 2.5*np.log10(star["FLUX_Z"]), 2)
SDSS_grz_values = [15.7, 15.19, 14.55]
print(f'compare calculated g, r, z = {[star_g_mag, star_r_mag, star_z_mag]} \nto SDSS values g, r, z = {SDSS_grz_values}')
print('these are.... kinda close? I am not too happy about how far apart the r_mag values are.')
print(f'Star WISE magnitudes : \nW1: {star["FLUX_W1"]} \nW2: {star["FLUX_W2"]} \nW3: {star["FLUX_W3"]} \nW4: {star["FLUX_W4"]}')
print('based on the magnitude for W4, this star was NOT detected in this band')


# In[ ]:





# In[ ]:




