#!/usr/bin/env python
# coding: utf-8

# In[1]:


# TJ needed imports
import sys
import os
from astropy.coordinates import SkyCoord
from astropy.table import Table
import numpy as np
from astropy import units as u
import glob

#TJ change directory to include the entire ASTRO5160 directory
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '..'))) 
from Homeworks.useful_functions import *


# In[2]:


#python task #1 Find object closest to (188.53667, 21.04572)
sweep_directory = '/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0' #TJ assign necessary sweep file directory
loc = SkyCoord(188.53667, 21.04572, unit = u.deg)

sweep_file = find_sweep_file(sweep_directory, loc) #TJ grab the only sweep file we need
sweep_data = read_partial_fits(sweep_file, ["RA", "DEC","ALLMASK_G", "ALLMASK_R", "ALLMASK_Z", "TYPE"]) #TJ read in the columns that we need
sweep_obj = find_matching_object(loc, sweep_data) #TJ find the closest obejct to our given location
print(f'object located at ({sweep_obj["RA"]} {sweep_obj["DEC"]}) is classified as an {sweep_obj["TYPE"]} which is an expontial galaxy')
print('This means the light profile is extended over several pixels (not a point source)')


# In[3]:


#Python task #2 Is this object saturated in any bands?
print(f'This object is saturated in all the G, R, and Z bands')
print('Using the viewer, this object is saturated, evidenced by the bright white/red discolored center')
print('This object is likely a galaxy because it is larger than a standard PSF would indicate, but I suppose it could be a quasar that just \nsaturates the nearby CCDs and simply appears to be larger than a point source.')


# In[4]:


#TJ Python task #3 find all PSF objects within 3 degrees of 180,30 with r_mag < 20, match to qso dataset
center = SkyCoord(180, 30, unit=u.deg) #TJ define center of field
needed_ra = 179, 179, 181, 181 #TJ assign random objects that will force the sweep files to grab all 4 of the sweep files that share this corner
needed_decs = 29, 31, 31, 29
dummy_locations = SkyCoord(needed_ra*u.deg, needed_decs*u.deg, unit = u.deg) #TJ create dummy locations to get all the neighboring sweep files
sweep_files = find_sweep_files_with(sweep_directory, dummy_locations)
sweep_columns_to_keep = ["RA", "DEC", "FLUX_G", "FLUX_R", "FLUX_Z", "FLUX_W1", "FLUX_W2", "TYPE"] #TJ only read in columns we need
sweep_data = vstack([read_partial_fits(file, sweep_columns_to_keep) for file in sweep_files]) #TJ subsequent stack sweep files as new rows on bottom
sweep_objs = SkyCoord(sweep_data["RA"], sweep_data["DEC"], unit=u.deg) #TJ create skycoord objects
near_center = sweep_data[center.separation(sweep_objs) < 3*u.deg] #TJ filter by distance to center less than 3 degrees
for flux_col in ['FLUX_W1', 'FLUX_W2', "FLUX_G", 'FLUX_R', "FLUX_Z"]: #TJ prevent attempting to take log10 of negative number by replacing with a small positive number
    near_center[flux_col] = np.where(near_center[flux_col] <= 0, 1e-10, near_center[flux_col]) 
near_center = add_mag_column(near_center, ['FLUX_W1', 'FLUX_W2', 'FLUX_G', 'FLUX_R', "FLUX_Z"]) #TJ add magnitude columns based on fluxes
valid_matches = ((near_center['R_mag'] < 20) & (near_center["TYPE"] == "PSF")) #TJ filter by r_mag and type
psfobjs = near_center[valid_matches]
print(f'Keeping {len(psfobjs)} objects. This should be 39655.')

qso_file_path = '/d/scratch/ASTR5160/week10/qsos-ra180-dec30-rad3.fits'
qso_table = Table.read(qso_file_path)
#TJ the first output is the psfobjs that have matches, second is qso_table objects, third is the unmatched rows from the psfobjs
qsos, _, not_qsos = cross_match_objects(psfobjs, qso_table, max_sep=1/3600)
print(f'Found {len(qsos)} matching qso objects. This should be 275.')

psfobjs.write('../Data_files/Class21_psfobjs.fits', format='fits', overwrite=True)
qsos.write('../Data_files/Class21_qsos.fits', format='fits', overwrite=True)
not_qsos.write('../Data_files/Class21_not_qsos.fits', format='fits', overwrite=True)


# In[ ]:




