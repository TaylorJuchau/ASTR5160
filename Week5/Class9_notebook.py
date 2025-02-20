#!/usr/bin/env python
# coding: utf-8

# In[17]:


#Import needed packages
import healpy as hp
from numpy.random import random
import numpy as np
import matplotlib.pyplot as plt


# In[19]:


#Python task #1 populate sphere with 1mil points
ra = 360.*(random(1000000)) #TJ generate 1 million right ascensions
dec = (180/np.pi)*np.arcsin(1.-random(1000000)*2.) #TJ generate 1 million declinations
#TJ plot these 1 million points to verify they uniformly paint a sphere
fig = plt.figure()
ax = fig.add_subplot(111, projection="aitoff")
ax.scatter(ra, dec, s=0.01)
xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
ax.set_xticklabels(xlab, weight=600)
ax.grid(color='blue', linestyle='--', linewidth=2)
ax.set_xlabel("Right Ascension", fontsize=12)
ax.set_ylabel("Declination", fontsize=12)
ax.set_title("Aitoff Projection of random dots", fontsize=12)


# In[20]:


#Python task #2 check which pixel each point is at
nside = 1 #TJ Assign nside number for pixel size
in_pixel = hp.ang2pix(nside, ra, dec, lonlat=True)


# In[37]:


#Python task #3 Check that each pixel has the same number of points in it
pixel_index, counts = np.unique(in_pixel, return_counts=True) #TJ extract pixel indices, and how many points are in each of those pixels
print('Coefficient of variation: ', np.std(counts)/np.mean(counts)) #TJ take the ratio of the standard deviation to the mean value.
print('This is very small, indicating that each pixel has about the same number of points in it')


# In[63]:


#Python task #4 plot pixels 2, 5 and 8
ii = in_pixel == 2 #TJ find indices that are in pixel #2
v = in_pixel == 5 #TJ find indices that are in pixel #5
viii = in_pixel == 8 #TJ find indices that are in pixel #8
plt.scatter(ra, dec, s=0.001, color = 'black')
plt.scatter(ra[ii], dec[ii], s = 0.1, color = 'red', label = 'pixel 2')
plt.scatter(ra[v], dec[v], s = 0.1, color = 'blue', label = 'pixel 5')
plt.scatter(ra[viii], dec[viii], s = 0.1, color = 'green', label = 'pixel 8')
plt.xlabel('RA')
plt.ylabel('Declination')
plt.title('visualize pixels 2, 5, and 8')
plt.legend(loc = 'upper left')
plt.show()


fig = plt.figure()
ax = fig.add_subplot(111, projection="aitoff")
ax.scatter((ra-180)*np.pi/180, dec*np.pi/180, s=0.001, color = 'black')
ax.scatter((ra[ii]-180)*np.pi/180, dec[ii]*np.pi/180, s = 0.1, color = 'red', label = 'pixel 2')
ax.scatter((ra[v]-180)*np.pi/180, dec[v]*np.pi/180, s = 0.1, color = 'blue', label = 'pixel 5')
ax.scatter((ra[viii]-180)*np.pi/180, dec[viii]*np.pi/180, s = 0.1, color = 'green', label = 'pixel 8')
xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
ax.set_xticklabels(xlab, weight=600)
ax.grid(color='blue', linestyle='--', linewidth=2)
ax.set_xlabel("Right Ascension", fontsize=12)
ax.set_ylabel("Declination", fontsize=12)
ax.set_title("Aitoff Projection of random dots", fontsize=12)


# In[60]:


#Python task #5 go to nside tier 2
smaller_nside = 2 #TJ Assign nside number for pixel size
in_smaller_pixel = hp.ang2pix(smaller_nside, ra, dec, lonlat=True)
nside1_pixel5 = in_pixel == 5
smaller_pixels, counts = np.unique(in_smaller_pixel[nside1_pixel5], return_counts=True)
print(f'nside==2 pixels that are inside nside==1 pixel #5 : {smaller_pixels}')


# In[ ]:




