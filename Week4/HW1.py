#!/usr/bin/env python
# coding: utf-8

# In[41]:


import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import re
from tabulate import tabulate
import argparse


data_file_location = '/d/scratch/ASTR5160/week4/HW1quasarfile.txt' #TJ Define default data file location

#TJ allow file to be run on my laptop by checking if I currently have access to the /d/ drive, if not, define data file location locally
try: 
    open(data_file_location, "r")
except FileNotFoundError:
    data_file_location = "C:/Users/tj360/ASTR5160/Data_files/HW1_data" #TJ local data files stored in this folder (not synced to github)
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
        locations.append([RA, Dec]) #TJ add to array of all locations in 'xxhxxmxx.xs', 'xxdxxmxx.xs' format
locations = np.array(locations) #TJ convert to numpy array for better manipulation later
QSO_locations = SkyCoord(locations[:,0], locations[:,1], unit = (u.hourangle, u.deg), frame='icrs') #TJ convert location array to SkyCoord array
Observing_location = EarthLocation.of_site('kpno') #TJ set observation location to Kitt Peak, in cartesian coordinates with units of meters

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
    
    year = 2025 #TJ set year to 2025, this can be easilly added as an argument, but assignment specifies that this is only for this year
    if type(argument) == int: #TJ check if argument is an integer, if it is, check that it is a valid month (between 1 and 12)
        if argument < 1 or argument > 12:
            print('integer month must be between 1-12')
            return None
        month = argument
    elif type(argument) == str and len(argument) == 3: #TJ check if argument is a string, and is in the accepted list of strings
        arg_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        if argument not in arg_list:
            print('month not in recognizable format: acceptable arguments are integers between 1-12, or three letter abreviation of month, or full month. \nString arguments must start with a capital letter')
            return None
        month = arg_list.index(argument)+1
    elif type(argument) == str and len(argument) > 3: #TJ check if month is spelled out entirely
        arg_list = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']
        if argument not in arg_list:
            print('month not in recognizable format: acceptable arguments are integers between 1-12, or three letter abreviation of month, or full month. \nString arguments must start with a capital letter')
            return None
        month = arg_list.index(argument)+1
    else:
        print('month not in recognizable format: acceptable arguments are integers between 1-12, or three letter abreviation of month, or full month. \nString arguments must start with a capital letter')
        return None
    month_lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] #TJ assign lengths of each month
    month_length = month_lengths[month-1] #TJ assign the appropriate number of days in the given month
    utcoffset = -7 * u.hour  #TJ define offset from UTC to Mountain Standard Time outside daylight savings
    start_time = Time(f"{year}-{month}-02 06:00:00") + utcoffset #TJ assign start date to be the 1st at 11pm MST
    rows = [["Date \nYYYY-MM-DD HH:MM:SS.SSS MST", "Quasar Coordinates \n(hhmmss.ss ◦ ′ ′′)", "RA \nDegrees (◦)", "Dec \nDegrees (◦)", "Airmass"]]
    for day in range(0,month_length):
        obs_time = start_time + day*24*u.hour #TJ observation time is an integer number of days after start time
        QSO_altaz_tonight = QSO_locations.transform_to(AltAz(obstime=obs_time, location=Observing_location)) #TJ get altitude and azimuth at this night
        airmass_array = QSO_altaz_tonight.secz #TJ get airmass at this time
        positive_airmass_indices = np.where(airmass_array > 0)[0]  #TJ Get indices of positive values for airmass
        prime_object_index = positive_airmass_indices[np.argmin(airmass_array[positive_airmass_indices])] #TJ extract lowest positive airmass
        time = Time(f"{year}-{month}-{day+1} 11:00:00") #TJ set obs_time manually because daylight savings time screws up just using obs_time.value
        prime_object_raw_location = lines[prime_object_index].strip() #TJ get raw location in given units for table
        RA, Dec = QSO_locations[prime_object_index].ra.deg, QSO_locations[prime_object_index].dec.deg #TJ get RA and Dec in normal person units
        min_air_mass = airmass_array[prime_object_index] #TJ get airmass for table
        rows.append([time, prime_object_raw_location, RA, Dec, min_air_mass]) #TJ add this day's data to table
    print(tabulate(rows[1:], headers = rows[0])) #TJ print table
    
    return None

#TJ make sure this part runs when called from the command line
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find lowest airmass quasars for each day of a given month.") #TJ add argparse object and description
    parser.add_argument("month", type=str, help="Month as an integer (1-12) or name (e.g., Jan, January)") #TJ assign argument type and help note
    args = parser.parse_args() #TJ assign arugments

    #TJ try to convert argument to integer, if that fails, assume it is a month name
    try:
        month_input = int(args.month)
    except ValueError:
        month_input = args.month  

    find_lowest_airmass(month_input) #TJ run the function


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




