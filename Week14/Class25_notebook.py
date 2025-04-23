#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


# In[2]:


#Python task #1 read in data and get variance
data_file_path = '/d/scratch/ASTR5160/week13/line.data'
show_plot = False #TJ change this to True to print plot to screen if you want
y_vals = np.loadtxt(data_file_path)
x_vals = np.linspace(0.5, 9.5, 10)
y_var = [np.var(y_vals[:,i], ddof=1) for i in range(10)]
y_mean = [np.mean(y_vals[:,i]) for i in range(10)]


# In[3]:


#Python task #2 get ln of post-prob

def ln_prior(m, b):
    if (0 < b < 8) and (1 < m < 10):
        return 0.0  #TJ flat prior inside acceptable range
    else:
        return -np.inf #TJ infinitely deweight models outside acceptable range 



def ln_post_prob(m, b):
    '''find log of posterior probability of the m, b values for this particular dataset given at /d/scratch/ASTR5160/week13/line.data
    -------------

    Parameters
    -------------
    m :  type = float - represents the model's slope
    b :  type = float - represents the model's y-intercept
    
    Returns
    -------------
    log of posterior probability
    '''
    return ln_prior(m, b) - 0.5*sum((((y_mean- (x_vals*m + b))**2)/y_var)+np.log(2*np.pi*np.array(y_var)))


# In[4]:


#Python task #3 use MCMC chain of parameters

def met_hast(x_vals, y_vals, y_vars, m_0=3.5, b_0=5, steps=5000, step_size=0.3):
    '''create MCMC chain of parameters
    -------------

    Parameters
    -------------
    x_vals : type = list - list of independent variables
    y_vals : type = list - list of dependent variables (data)
    y_vars : type = list - list of variances in y_vals
    m_0 (optional, defaults to 3.5) : type = float - starting slope parameter
    b_0 (optional, defaults to 4) : type = float - starting y-intercept parameter
    steps (optional, defaults to 5000) : type = int - integer number of steps through parameter space we will take
    step_size (optional, defaults to 0.3) : type = float - standard deviation of gaussian used for generating new parameters
    **note** step size was varied until acceptance rate of ~30% was achieved. landed on 0.3 step size
    
    Returns
    -------------
    MCMC chain
    '''
    
    m_chain = [m_0] #TJ initialize chains
    b_chain = [b_0]
    accepted = [] #TJ initialize accepted array (the mean of this array is the acceptance rate)
    current_ln_post = ln_post_prob(m_0, b_0) #TJ this will be overwritten with a better value if one is found
    
    for i in range(steps):
        #TJ generate new m and b to compare against previously held value
        m_prop = np.random.normal(m_chain[-1], step_size)
        b_prop = np.random.normal(b_chain[-1], step_size)
        
        #TJ compute new log-posterior
        prop_ln_post = ln_post_prob(m_prop, b_prop)
        
        #TJ compute difference in post_prob values and accept or reject based on a randomly generated number
        ln_R = prop_ln_post - current_ln_post
        if np.random.rand() < np.exp(ln_R):
            m_chain.append(m_prop)
            b_chain.append(b_prop)
            current_ln_post = prop_ln_post
            accepted.append(1)
        else:
            m_chain.append(m_chain[-1])
            b_chain.append(b_chain[-1])
            accepted.append(0)
    
    return np.array(m_chain), np.array(b_chain), np.mean(accepted)
m_chain, b_chain, acc_rate = met_hast(x_vals, y_vals, y_var)
#TJ perform test, but ignore the first 500 values for the means, since these values will be skewed toward initial guesses
print(f'This resulted in an acceptance rate of {acc_rate} and a mean value of m = {np.mean(m_chain[500:]):.3f} +/- {np.std(m_chain[500:]):.3f}, b = {np.mean(b_chain[500:]):.3f} +/- {np.std(b_chain[500:]):.3f}')

