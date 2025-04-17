#!/usr/bin/env python
# coding: utf-8

# In[1]:


from scipy.stats import chi2
import numpy as np
import matplotlib.pyplot as plt


# In[23]:


#Python task #1
data_file_path = '/d/scratch/ASTR5160/week13/line.data'
show_plot = False #TJ change this to True to print plot to screen if you want
y_vals = np.loadtxt(data_file_path)
x_vals = np.linspace(0.5, 9.5, 10)
mean_y = [np.mean(y_vals[:,i]) for i in range(10)]
y_std = [np.std(y_vals[:,i]) for i in range(10)]
y_var = [np.var(y_vals[:,i], ddof=1) for i in range(10)] #TJ this will be the diagonal entries for the cov_matrix
#TJ this matrix is how well each of the 10 data bins correlate to each other, theres no reason for there to be 20 entries for any x_bin
cov_matrix = np.cov(y_vals, rowvar=False) #TJ this should be a 10by10 matrix
print('Verify that the first numbers (diagnonals of the cov matrix) are equal to the second numbers (variances)')
for i in range(10):
    print(cov_matrix[i][i], y_var[i])





#TJ if desired, show plot
if show_plot:
    for i in range(20):
        plt.scatter(x_vals, y_vals[i], s = 1)
    plt.errorbar(x_vals, mean_y, yerr=y_std, fmt='o', 
                 color='blue', ecolor='red', capsize=5,
                 label='Mean with std deviation')
    plt.title('some random data I found in a file')
    plt.xlabel('anonymous x axis')
    plt.ylabel('anonymous y axis')
    plt.show()


# In[47]:


#Python task #2 find the most anticorrelated bins
i, j = np.where(cov_matrix == np.min(cov_matrix)) #TJ search for most anticorrelated (lowest value in cov_matrix)
print(f'The values at [i,j] = [{i[0],j[0]}] and [{i[1],j[1]}] are both {cov_matrix[i[0]][j[0]]}, which is the lowest value in the covariance matrix')
print('This means that the datapoints in the 1st and 4th columns are the most ANTI-correlated')
np.fill_diagonal(cov_matrix, -np.inf) #TJ replace diagnonals with negative infinity to ignore them for consideration of most correlated
i, j = np.where(cov_matrix == np.max(cov_matrix)) #TJ search for most correlated (largest value in cov_matrix)
print(f'The values at [i,j] = [{i[0],j[0]}] and [{i[1],j[1]}] are both {cov_matrix[i[0]][j[0]]}, which is the largest [non-diagonal] value in the covariance matrix')
print('This means that the datapoints in the 5th and 8th columns are the most correlated')

