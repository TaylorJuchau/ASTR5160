#!/usr/bin/env python
# coding: utf-8

# In[20]:


#TJ import needed packages
import numpy as np
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
#Python task #2 convert RA and Dec to degrees
sample_dec_string = "67 36 4" #TJ just a random Declination in units of (degrees, minutes, seconds)
sample_RA_string = "13 25 15" #TJ just a random Right Ascension in units of (hours, minutes, seconds)
coordinates = SkyCoord(sample_RA_string, sample_dec_string, unit=(u.hourangle, u.deg)) #TJ convert strings into a SkyCoord coordinate
manually_converted = [15*(13+(25/60)+(15/3600)), (67+(36/60)+(4/3600))]
print(f'{coordinates}     manually converted: {manually_converted}') #TJ print both methods to check


# In[2]:


#Python task #3 Find current Julian date and modified Julian date
now = Time.now() #TJ create timestamp
now_jd = now.jd #TJ convert to Julian date
now_mjd = now.mjd #TJ convert to modified Julian date
print(f'jd : {now_jd},  mjd : {now_mjd}, difference : {now_jd-now_mjd}') #print to check relationship


# In[15]:


#Python task #4 Create array of Julian dates for days around today
julian_dates_around_today = [] #TJ initialize array
for i in np.arange(-10,11):
    time_for_new_day = Time.now()+i #TJ create timestamo for this new day
    julian_dates_around_today.append([f"today {i} days : ", time_for_new_day.mjd]) #TJ convert new day into MJD and label row
julian_dates_around_today #TJ display the array


# In[19]:


#Python task #5 set WIRO location
latitude = (41 + (5 / 60) + (49 / 3600)) *u.deg   #TJ latitude of WIRO in degrees
longitude = (-(105 + (58 / 60) + (33 / 3600)))*u.deg  #TJ longitude of WIRO in degrees
altitude = 2943 * u.m  #TJ altitude of WIRO in meters

WIRO_location = EarthLocation.from_geodetic(lat=latitude, lon=longitude, height=altitude) #TJ set WIRO location


# In[29]:


#python task #6 calculate air mass
utcoffset = -7 * u.hour  #TJ define offset from UTC to Mountain time
time_tonight = Time("2025-1-31 23:00:00") - utcoffset #TJ set time we are going to observe at tonight
object_RA = '12 0 0' #TJ set object RA in hours, min, sec
object_dec = '30 0 0' #TJ set object declination in degrees
object_coordinates = SkyCoord(object_RA, object_dec, unit=(u.hourangle, u.deg)) #create SkyCoord location
object_altaz_tonight = object_coordinates.transform_to(AltAz(obstime=time_tonight, location=WIRO_location)) #TJ get altitude and azimuth
print(f"object's Altitude tonight= {object_altaz_tonight.alt:.2}    object's Azimuth = {object_altaz_tonight.az:.2}")
object_airmass_tonight = object_altaz_tonight.secz #get airmass at this time
print(f'object airmass tonight at 11pm : {object_airmass_tonight}')
time_next_month = Time("2025-3-02 23:00:00") - utcoffset #TJ set time we are going to observe at in a month
object_altaz_next_month = object_coordinates.transform_to(AltAz(obstime=time_next_month, location=WIRO_location)) #TJ get altitude and azimuth
print(f"object's Altitude next month= {object_altaz_next_month.alt:.2}    object's Azimuth = {object_altaz_next_month.az:.2}")
object_airmass_next_month = object_altaz_next_month.secz #get airmass at this time
print(f'object airmass next month at 11pm : {object_airmass_next_month}')



# In[21]:





# In[ ]:




