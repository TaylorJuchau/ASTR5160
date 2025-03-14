#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data_file_location = "C:/Users/tj360/ASTR5160/Data_files/Class15_SQL_query.csv" 
#TJ local data files stored in this folder (not synced to github)


# In[7]:


#Python task 2 plot ra and dec using circles
data = pd.read_csv(data_file_location) #TJ get data from csv file
ra, dec = data.ra, data.dec #TJ extract ra and dec columns
plt.scatter(ra, dec, marker='o', s=20, facecolors='none', edgecolors='black')  #TJ assuming task means hollow circles
plt.xlabel('ra (degrees)')
plt.ylabel('dec (degrees)')
#TJ matplotlib kept displaying the x-axis as centred at 0, with ticks being +3e2, this fixes that
plt.ticklabel_format(useOffset=False, axis='x') 
plt.gca().invert_xaxis() #TJ invert x-axis so RA decreases to the right.
plt.title("circles representing objects within 2' of ra,dec = 300,-1")


plt.show()


# In[14]:


#Python task #3 make size represent g-magnitude
#TJ I know this says to bin the points, and I could do that, but I think its better to use a continuous size distribution
g = data.g #TJ extract 
size_array = 25 - 9*((g-np.mean(g))/(np.std(g))) #TJ create size array where the average size is 25, and
#TJ each point has its size changed by an amount proportional to its number of standard deviations from the mean
#TJ the minus sign indicates that if g is bigger than the mean g_magnitude, it appears smaller
plt.scatter(ra, dec, marker='o', s=size_array, facecolors='none', edgecolors='black')  #TJ assuming task means hollow circles
plt.xlabel('ra (degrees)')
plt.ylabel('dec (degrees)')
plt.ticklabel_format(useOffset=False, axis='x') #TJ fix offset problem
plt.gca().invert_xaxis() #TJ invert x-axis
plt.title("circles representing objects within 2' of ra,dec = 300,-1\nsize of circles represent g_magnitude")
plt.show()


# In[4]:


#Python task #4
print('yes! they do look pretty similar! it did make me realize I need to invert the x-axis though because there is\na bright star near the middle, but off to the left in the navigator, which corresponds to a big circle that\nwas previously off to the right on my scatter plot.')


# In[ ]:





# In[ ]:





# In[ ]:




