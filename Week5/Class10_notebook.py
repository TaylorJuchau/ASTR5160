#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u


# In[24]:


#Python task #1 cap bounded by 5h right ascension
def define_cap_bounded_by_ra(right_ascension):
    '''Define a cap bounded by a given Right Ascension
    
    Parameters
    -------------
    RA : right ascension in hours

    Returns
    -------------
    4-d vector representing the x,y,z pointing to the center of the cap, and the size of the cap in 1-cos(theta), 
    here assumed to be half the hemisphere
    
    '''
    RA = (right_ascension+6)*u.hourangle #TJ set cap center to be offset by 90 degrees
    Dec = 0*u.deg #TJ set declination to 0
    center = SkyCoord(ra = RA, dec = Dec) #TJ get skycoord of center of cap
    cap_xyz = center.cartesian #TJ convert to x, y, z values
    h = 1 #TJ set height of cap to 2
    return (float(cap_xyz.x), float(cap_xyz.y), float(cap_xyz.z), 1)
if __name__ == "__main__":
    answer = define_cap_bounded_by_ra(5)
    print(answer)


# In[41]:


def define_cap_bounded_by_dec(declination):
    '''Define a cap bounded by a given declination in degrees
    
    Parameters
    -------------
    Declination : declination of the bottom of the circumpolar cap

    Returns
    -------------
    4-d vector representing the x,y,z pointing to the center of the cap, and the size of the cap in 1-cos(theta) set by
    the provided declination
    
    '''
    RA = 0*u.hourangle
    Dec = 90*u.deg
    center = SkyCoord(ra = RA, dec = Dec)
    cap_xyz = center.cartesian
    h = np.sin(declination*u.deg)
    x, y, z = float(cap_xyz.x), float(cap_xyz.y), float(cap_xyz.z)
    if abs(cap_xyz.x) < 0.00000001:
        x = 0.
    if abs(cap_xyz.y) < 0.00000001:
        y = 0.
    if abs(cap_xyz.z) < 0.00000001:
        z = 0.
    return (x, y, z, float(1-h))
if __name__ == "__main__":
    answer = define_cap_bounded_by_dec(36)
    print(answer)


# In[51]:


def define_cap(Right_Ascension, Declination, Radius):
    '''Define a cap centered at provided Right Ascension and Declination, with radius in degrees
    
    Parameters
    -------------
    Right_Ascension : RA in hours
    Declination : declination of the center of cap in degrees
    Radius : sets how many degrees the cap's radius is

    Returns
    -------------
    4-d vector representing the x,y,z pointing to the center of the cap, and the size of the cap in 1-cos(theta)
    
    '''
    RA = Right_Ascension*u.hourangle
    Dec = Declination*u.deg
    center = SkyCoord(ra = RA, dec = Dec)
    cap_xyz = center.cartesian
    x, y, z = float(cap_xyz.x), float(cap_xyz.y), float(cap_xyz.z)
    if abs(cap_xyz.x) < 0.00000001:
        x = 0.
    if abs(cap_xyz.y) < 0.00000001:
        y = 0.
    if abs(cap_xyz.z) < 0.00000001:
        z = 0.
    h = np.cos(Radius*u.deg)
    return (x, y, z, float(1-h))
if __name__ == "__main__":
    answer = define_cap(5,36,1)
    print(answer)


# In[ ]:




