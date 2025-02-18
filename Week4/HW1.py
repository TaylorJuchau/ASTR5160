#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import re
from tabulate import tabulate
import argparse


data_file_location = '/d/scratch/ASTR5160/week4/HW1quasarfile.txt' #Define data file location
locations = [] #TJ Initialize what will be the full list of 1111 locations
with open(data_file_location, "r") as file:
    lines = file.readlines()  #TJ Read out data from each line
    for line in lines:
        line_string = line.strip() #TJ strip off the /n from each line
        parts = re.split(r'(\+|-)', line_string, maxsplit=1)  #TJ Split at first + or -
        location_strings = [parts[0], parts[1] + parts[2]] #TJ Split ra and dec into separate strings
        #TJ parse each string to extract hour, min, second from the RA
        RA = f'{int(location_strings[0][:2])}'+'h'+f'{int(location_strings[0][2:4])}'+'m'+ f'{float(location_strings[0][4:])}'+'s'
        #TJ parse each string to extract degrees, min, sec from Dec, make sure to count the sign in front and return it to the 
        sign = -1 if location_strings[1][0] == '-' else 1 #record the sign of the declination
        Dec = f'{int(location_strings[1][1:3]) * sign}' + 'd' + f'{int(location_strings[1][3:5])}' + 'm' + f'{float(location_strings[1][5:])}'+'s'
        locations.append([RA, Dec])
locations = np.array(locations)
QSO_locations = SkyCoord(locations[:,0], locations[:,1], unit = (u.hourangle, u.deg), frame='icrs')
Observing_location = EarthLocation.of_site('kpno') #TJ set observation location to Kitt Peak, in cartesian coordinates with units of meters
utcoffset = -7 * u.hour  #TJ define offset from UTC to Mountain time
def find_lowest_airmass(argument):
    '''Find quasars that have the lowest air mass for every day of a given month, when observed from Kitt Peak at 11pm MST
    
    Parameters
    -------------
    month : class integer, representing a month of the year 2025 for which to run this analysis over.
    or
    month : class string, that is the three letter abreviation of a month, for examole "Jan" for January
    
    Returns
    -------------
    Table of values with columns for:
    day (should be one row for every day in the provided month)
    object location in hhmmss.ss ddmmss.ss
    object Right Ascension given in degrees
    object Declination in degrees
    airmass of object at 11pm MST when viewed from Kitt Peak Observatory
    '''
    year = 2025
    if type(argument) == int:
        if argument < 1 or argument > 13:
            print('integer month must be between 1-12')
            return None
        month = argument
    elif type(argument) == str and len(argument) == 3:
        arg_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        if argument not in arg_list:
            print('month not in recognizable format: acceptable arguments are integers between 1-12, or three letter abreviation of month, or full month. \nString arguments must start with a capital letter')
            return None
        month = arg_list.index(argument)+1
    elif type(argument) == str and len(argument) > 3:
        arg_list = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
        if argument not in arg_list:
            print('month not in recognizable format: acceptable arguments are integers between 1-12, or three letter abreviation of month, or full month. \nString arguments must start with a capital letter')
            return None
        month = arg_list.index(argument)+1
    else:
        print('month not in recognizable format: acceptable arguments are integers between 1-12, or three letter abreviation of month, or full month. \nString arguments must start with a capital letter')
        return None
    month_lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    month_length = month_lengths[month-1]
    start_time = Time(f"{year}-{month}-02 06:00:00") +utcoffset
    rows = [["Date", "Quasar Coordinates (hms.ss ◦ ′ ′′)", "RA (◦)", "Dec (◦)", "Airmass"]]
    for day in range(0,month_length):
        obs_time = start_time + day*24*u.hour
        QSO_altaz_tonight = QSO_locations.transform_to(AltAz(obstime=obs_time, location=Observing_location)) #TJ get altitude and azimuth at this night
        airmass_array = QSO_altaz_tonight.secz #get airmass at this time
        positive_airmass_indices = np.where(airmass_array > 0)[0]  #TJ Get indices of positive values for airmass
        prime_object_index = positive_airmass_indices[np.argmin(airmass_array[positive_airmass_indices])]
        time = obs_time.value
        prime_object_raw_location = lines[prime_object_index].strip()
        RA, Dec = QSO_locations[prime_object_index].ra, QSO_locations[prime_object_index].dec
        min_air_mass = airmass_array[prime_object_index]
        rows.append([time, prime_object_raw_location, RA, Dec, min_air_mass])
    print(tabulate(rows[1:], headers = rows[0]))
    
    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find lowest airmass quasars for each day of a given month.")
    parser.add_argument("month", type=str, help="Month as an integer (1-12) or name (e.g., 'Jan', 'January')")
    args = parser.parse_args()
    
    find_lowest_airmass(args.month)


# In[2]:


find_lowest_airmass(1)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




