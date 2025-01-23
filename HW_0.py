#!/usr/bin/env python
# coding: utf-8

# In[1]:


#HW_0
#problem 1

#TJ import packages
import numpy as np
import random
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def make_line(m,b):
#TJ generates 10 datapoints along the line Y = mx + b, then generates random error with a standard deviation of
#TJ 0.5. Takes arguments for slope and y-intercept, returns 3 10x1 arrays for x, y, and error.
    #TJ Initialize arrays for x, y_exact, and the error on each y value
    x_points = []
    y_points = []
    y_error = []
    for i in range(0,10):
        x_points.append(random.uniform(0,10)) #TJ generate 10 random x values between 0 and 10
        correct_y = (x_points[-1]*m)+b #TJ generates correct y values
        y_points.append(correct_y + (random.gauss(0,0.5))) #TJ adds random error to correct y values (mean=0, std=0.5)
        y_error.append(y_points[-1] - correct_y) #TJ calculate error on y and record this value
    x_points = np.array(x_points) #TJ converts list of values into a manipulatable array
    y_points = np.array(y_points) #TJ converts list of values into a manipulatable array
    y_error = np.array(y_error) #TJ converts list of values into a manipulatable array
    
    return x_points, y_points, y_error    #TJ return desired arrays as outputs


# In[2]:


#problem 2

m = random.uniform(-10,10) #TJ define parameters for the original line, randomly generated to make sure plot looks good for all
b = random.uniform(-10,10)
data_array_unsorted = make_line(m,b) #TJ generate data to be analyzed

#Define basic linear function to use to fit to data
def linear_func(x,m,b):
#basic linear function, takes arguments for an x-value, slope and y-intercept and returns the y value
    return (m*x) + b

params, cov = curve_fit(linear_func, data_array_unsorted[0], data_array_unsorted[1]) #extract fitted values for m and b to variable params


# In[3]:


#problem 3
#TJ commented out so it doesnt run from the command line when this python file is called
'''
x_range = np.arange(0,11)
sorted_indices = np.argsort(data_array_unsorted[0]) #TJ original plot looks terrible, I need to reorder them.
x_data = data_array_unsorted[0][sorted_indices] #TJ create array of x_data that is sorted
y_data = data_array_unsorted[1][sorted_indices] #TJ create array of y_data that is sorted
error_bars = abs(data_array_unsorted[2][sorted_indices])
plt.plot(x_range, linear_func(x_range, m, b), label = f"orignal y = {np.round(m,1)}x+{np.round(b,1)} line", color = 'red', linestyle='dotted') # TJ plot after errors
plt.plot(x_data, y_data, label = 'original line w/ error added', color = 'blue') #TJ plot original y=mx+b line
plt.errorbar(x_data, y_data, yerr=0.5, xerr=0, label='')
plt.plot(x_range, linear_func(x_range, *params), linestyle='dashed', label = "reconsituted line from fitted paramters", color = 'purple')
plt.legend()
plt.xlim(0, 10)  #TJ Set x-axis to be 0 to 10, since thats the range the x-values were selected in
y_min = np.min(y_data)-2 #TJ set y limits to be slightly more extreme than the most extreme y-values
y_max = np.max(y_data)+2
plt.ylim(y_min,y_max)
plt.title('original line, line w/ error, and reconstituted line')
plt.xlabel('x')
plt.ylabel('y')
plt.show()'''


# In[4]:


#problem 4
#TJ commented out so it doesnt run from the command line when this python file is called

'''
plt.plot(x_range, linear_func(x_range, m, b), label = f"orignal y = {np.round(m,1)}x+{np.round(b,1)} line", color = 'red', linestyle='dotted') # TJ plot after errors
plt.plot(x_data, y_data, label = 'original line w/ error added', color = 'blue') #TJ plot original y=mx+b line
plt.errorbar(x_data, y_data, yerr=0.5, xerr=0, label='')
plt.plot(x_range, linear_func(x_range, *params), linestyle='dashed', label = "reconsituted line from fitted paramters", color = 'purple')
plt.legend()
plt.xlim(0, 10)  #TJ Set x-axis to be 0 to 10, since thats the range the x-values were selected in
y_min = np.min(y_data)-2 #TJ set y limits to be slightly more extreme than the most extreme y-values
y_max = np.max(y_data)+2
plt.ylim(y_min,y_max)
plt.title('original line, line w/ error, and reconstituted line')
plt.xlabel('x')
plt.ylabel('y')
plt.savefig("HW_0_figure.png", dpi=300, bbox_inches='tight')
plt.show()'''


# In[5]:


#problem 5

def do_HW_0(m,b):
    #TJ Takes arguements for slope and y intercept, generates 10 numbers between 0 and 10 to use as x points
    #   then generates y values according to y = mx+b, then tweaks each of those y-values by a number randomly 
    #   selected from a normal distribution centered at 0 with a standard deviation of 0.5, then fits parameters
    #   to the tweaked dataset to reconsitute a y=mx+b line. Plots three datasets (original, w/errors, and reconstituted)
    #   on the same plot. Saves plot as PNG to working directory with title "HW_0_figure.png". Returns nothing.
    data_array_unsorted = make_line(m,b)
    params, cov = curve_fit(linear_func, data_array_unsorted[0], data_array_unsorted[1]) #extract fitted values for m and b to variable params
    x_range = np.arange(0,11)
    sorted_indices = np.argsort(data_array_unsorted[0]) #TJ original plot looks terrible, I need to reorder them.
    x_data = data_array_unsorted[0][sorted_indices] #TJ create array of x_data that is sorted
    y_data = data_array_unsorted[1][sorted_indices] #TJ create array of y_data that is sorted
    error_bars = abs(data_array_unsorted[2][sorted_indices])
    plt.plot(x_range, linear_func(x_range, m, b), label = f"orignal y = {np.round(m,1)}x+{np.round(b,1)} line", color = 'red', linestyle='dotted') # TJ plot after errors
    plt.plot(x_data, y_data, label = 'original line w/ error added', color = 'blue') #TJ plot original y=mx+b line
    plt.errorbar(x_data, y_data, yerr=0.5, xerr=0, label='')
    plt.plot(x_range, linear_func(x_range, *params), linestyle='dashed', label = "reconsituted line from fitted paramters", color = 'purple')
    plt.legend()
    plt.xlim(0, 10)  #TJ Set x-axis to be 0 to 10, since thats the range the x-values were selected in
    y_min = np.min(y_data)-2 #TJ set y limits to be slightly more extreme than the most extreme y-values
    y_max = np.max(y_data)+2
    plt.ylim(y_min,y_max)
    plt.title('original line, line w/ error, and reconstituted line')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig("HW_0_figure.png", dpi=300, bbox_inches='tight')
    plt.show()
    return None


# In[6]:


if __name__ == "__main__":
    m = float(input("Slope: "))
    b= float(input("Y-intercept: "))
    do_HW_0(m,b)


# In[ ]:




