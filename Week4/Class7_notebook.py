#!/usr/bin/env python
# coding: utf-8

# In[66]:


import numpy as np
from numpy.random import random
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy import wcs
import dustmaps
from dustmaps.config import config
from dustmaps.sfd import SFDQuery


# In[32]:


#Python task #1 
ra = 2*np.pi*(random(10000)-0.5) #TJ generates 10000 random RA's measured in radians between -pi and pi
dec = np.arcsin(1.-random(10000)*2.) #TJ creates 10000 random declinations in radians between +pi/2 and -pi/2
plt.scatter(ra,dec, s=0.1) #TJ plot in standard x,y axes
plt.xlabel('Right Ascension')
plt.ylabel('Declination')
plt.title('Random dots on a sphere')
print('There are more dots near the equator, since there is more surface area in that region')


# In[56]:


#TJ Python task #2 project onto a sphere
fig = plt.figure()
ax = fig.add_subplot(111, projection="aitoff")
ax.scatter(ra, dec, s=0.1)
xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
ax.set_xticklabels(xlab, weight=600)
ax.grid(color='blue', linestyle='--', linewidth=2)
ax.set_xlabel("Right Ascension", fontsize=12)
ax.set_ylabel("Declination", fontsize=12)
ax.set_title("Aitoff Projection of random dots", fontsize=12)


fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection="lambert")
ax2.scatter(ra, dec, s=0.1)
xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
ax2.set_xticklabels(xlab, weight=600)

ax2.grid(color='blue', linestyle='--', linewidth=2)
ax2.set_xlabel("Right Ascension", fontsize=12)
ax2.set_ylabel("Declination", fontsize=12)
ax2.set_title("Lambert Projection of random dots", fontsize=12)
fig2.show()



# In[61]:


#Python task #3 Create dust grid map #I cant do this from my computer, I can do it on a separate notebook on campus if you'd like
mesh_array = np.meshgrid(np.linspace(0.5,359.5,360),np.linspace(-89.5,89.5,180)) #TJ make the 1by1 grid of locations
dust_map_1 = sfd(SkyCoord(mesh_array[0]*u.deg,mesh_array[1]*u.deg))
w = wcs.WCS(naxis=2)
w.wcs.ctype = ["RA---AIT", "DEC--AIT"]
– x, y = w.wcs_world2pix(ra, dec, 1)
Python tasks


# In[62]:





# In[ ]:




