#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.table import Table
import corner
import emcee
import numpy as np
import matplotlib.pyplot as plt

#######TJ setting up the linear model######
def linear_model(theta, x): #TJ define linear function in general terms
    """Simple y = mx+b model that returns a y value for a given m, b and x value

    Parameters
    -------------
    theta :  type = list - list of length 2 that contains a slope "m" value as the first entry and a y-intercept "b" value as its second entry
    x : type = float - float representing an x value for which an associated y value will be returned
    
    Returns
    -------------
    a single y-value
    """
    m, b = theta
    return m * x + b

def log_prior_linear(theta): #TJ assign priors to eliminate impossible values that are definitely not correct
    """function to assign uniform prior weights to a bayesian analysis

    Parameters
    -------------
    theta :  type = list - list of length 2 that contains a slope "m" value as the first entry and a y-intercept "b" value as its second entry
    
    Returns
    -------------
    either zero or negative infinity depending on whether or not the passed m and b values are physically realistic to even consider.      
    """
    m, b = theta
    if -5 < m < 0 and -7 < b < 7:
        return 0.0  # log(1)
    return -np.inf

def log_likelihood_linear(theta, x, y, yerr): #TJ define log likelihood function for the linear fit
    """Uses the formula for log of the likelihood of the linear model

    Parameters
    -------------
    theta :  type = list - list of length 2 that contains a slope "m" value as the first entry and a y-intercept "b" value as its second entry
    x : type = list of floats - float representing an x value for an associated y value
    y : type = list of floats - float representing an y value for an associated x value
    yerr : type = list of floats - float representing the error on a given y value
    Returns
    -------------
    float to be compared to other outputs for different m,b values, closer to zero is better.   
    """
    model = linear_model(theta, x)
    return -0.5 * np.sum(((y - model) / yerr) ** 2 + np.log(2 * np.pi * yerr**2))

def log_posterior_linear(theta, x, y, yerr): #TJ define a function to add the priors to the log likelihood
    """scales the log of the likelihood by the priors assigned in the log_prior_linear function

    Parameters
    -------------
    theta :  type = list - list of length 2 that contains a slope "m" value as the first entry and a y-intercept "b" value as its second entry
    x : type = float - float representing an x value for which an associated y value will be returned
    y : type = list of floats - float representing an y value for an associated x value
    yerr : type = list of floats - float representing the error on a given y value
    Returns
    -------------
    float to be compared to other outputs for different m,b values, closer to zero is better.         
    """
    lp = log_prior_linear(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood_linear(theta, x, y, yerr)


#####TJ setting up the quadratic model######




def quadratic_model(theta, x): #TJ define linear function in general terms
    """Simple y = a2*x^2 + a1*x + a0 model that returns a y value for a given x value based on a2, a1, a0

    Parameters
    -------------
    theta :  type = list - list of length 3 that contains entries in order [a2, a1, a0]
    x : type = float - float representing an x value for which an associated y value will be returned
    
    Returns
    -------------
    a single y-value        
    """
    a2, a1, a0 = theta
    return a2 * x**2 + a1 * x + a0

def log_prior_quadratic(theta): #TJ assign priors to eliminate impossible values that are definitely not correct
    """function to assign uniform prior weights to a bayesian analysis

    Parameters
    -------------
    theta :  type = list - list of length 3 that contains entries in order [a2, a1, a0]
    
    Returns
    -------------
    either zero or negative infinity depending on whether or not the passed a2, a1, a0 values are physically realistic to even consider.      
    """
    a2, a1, a0 = theta
    if -2 < a2 < 2 and -5 < a1 < 5 and -7 < a0 < 7:
        return 0.0
    return -np.inf

def log_likelihood_quadratic(theta, x, y, yerr): #TJ define log likelihood function for the linear fit
    """Uses the formula for log of the likelihood of the quadratic model, I could have put both models under the same function
        and then just left which model it should use to be determined by an argument, but I have already copied and pasted it...
        ah... hindsight is 20/20 I guess

    Parameters
    -------------
    theta :  type = list - list of length 3 that contains entries in order [a2, a1, a0]
    x : type = list of floats - float representing an x value for an associated y value
    y : type = list of floats - float representing an y value for an associated x value
    yerr : type = list of floats - float representing the error on a given y value
    Returns
    -------------
    float to be compared to other outputs for different coefficient values, closer to zero is better.   
    """
    model = quadratic_model(theta, x)
    return -0.5 * np.sum(((y - model) / yerr) ** 2 + np.log(2 * np.pi * yerr**2))

def log_posterior_quadratic(theta, x, y, yerr): #TJ define a function to add the priors to the log likelihood
    """scales the log of the likelihood by the priors assigned in the log_prior_quadratic function

    Parameters
    -------------
    theta :  type = list - list of length 3 that contains entries in order [a2, a1, a0]
    x : type = float - float representing an x value for which an associated y value will be returned
    y : type = list of floats - float representing an y value for an associated x value
    yerr : type = list of floats - float representing the error on a given y value
    Returns
    -------------
    float to be compared to other outputs for different coefficient values, closer to zero is better.         
    """
    lp = log_prior_quadratic(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood_quadratic(theta, x, y, yerr)

##########TJ run the model fit procedure from emcee############


def run_emcee(log_posterior, ndim, nwalkers=32, nsteps=5000, initial_guess=None):
    """runs the emcee protocol for a set of data

    Parameters
    -------------
    log_posterior :  type = function - function determining the log_posterior for a given model (including priors)
    ndim : type = int - integer representing the number of parameters that the model needs to fit
    nwalkers (optional, defaults to 32) : type = int - integer number of walkers used when sampling the parameter space (should be at least twice ndim)
    nsteps (optional, defaults to 5000) : type = int - integer number of steps each walker will take on its path
    initial_guess (optional, defaults to none) : type = float - initial starting points for the walkers, will be randomly generated if left as none
    Returns
    -------------
    float to be compared to other outputs for different coefficient values, closer to zero is better.         
    """
    
    if initial_guess is None:
        initial_guess = np.random.randn(nwalkers, ndim) #TJ create starting points for walkers
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(x, y, yerr)) #TJ initialize the samplers
    sampler.run_mcmc(initial_guess, nsteps, progress=True) #TJ perform the emcee protocol
    return sampler

def summarize_samples(samples, labels):
    """summarizes the samples statistics to produce the 16, 50, 84 precentile values

    Parameters
    -------------
    samples :  type = np array - array containing 
    ndim : type = int - integer representing the number of parameters that the model needs to fit
    nwalkers (optional, defaults to 32) : type = int - integer number of walkers used when sampling the parameter space (should be at least twice ndim)
    nsteps (optional, defaults to 5000) : type = int - integer number of steps each walker will take on its path
    initial_guess (optional, defaults to none) : type = float - initial starting points for the walkers, will be randomly generated if left as none
    Returns
    -------------
    float to be compared to other outputs for different coefficient values, closer to zero is better.         
    """
    
    percentiles = np.percentile(samples, [16, 50, 84], axis=0) #TJ calculate the 16th, 50th and 84 percentile values
    summaries = {}
    for i, label in enumerate(labels):
        p16, p50, p84 = percentiles[:, i]
        summaries[label] = (p50, p50 - p16, p84 - p50)
    return summaries

if __name__ == "__main__":
    
    import argparse
    #TJ add argparse object and description
    parser = argparse.ArgumentParser(description="""Runs emcee package on data stored at /d/scratch/ASTR5160/final/dataxy.fits to fit both
    a linear model (m and b values) as well as a quadratic (a2, a1, a0 values). Then using 16, 50, and 84 percentile values, reports the best
    fit parameters as median plus or minus 1 sigma. Then reports on whether the linear fit is sufficient to model the data (it is not).
    
    """) 


    
    data_file_path = '/d/scratch/ASTR5160/final/dataxy.fits' #TJ assign datafile location
    data = Table.read(data_file_path)
    x = data['x']
    y = data['y']
    yerr = data['yerr']



    
    #TJ run the linear fit analysis
    sampler_linear = run_emcee(log_posterior_linear, ndim=2)
    samples_linear = sampler_linear.get_chain(discard=1000, thin=15, flat=True) #TJ discard the first 1000 steps to ensure convergence
    
    # Run quadratic fit
    sampler_quadratic = run_emcee(log_posterior_quadratic, ndim=3)
    samples_quadratic = sampler_quadratic.get_chain(discard=1000, thin=15, flat=True)
    lin_labels = ["m", "b"]
    summary_linear = summarize_samples(samples_linear, lin_labels )
    quad_labels = ["a2", "a1", "a0"]
    summary_quadratic = summarize_samples(samples_quadratic, quad_labels)
    
    
    
    xfit = np.linspace(min(x), max(x), 1000) #TJ make a higher resolution x-axis for fitted functions
    
    # Plot linear fit
    plt.errorbar(x, y, yerr=yerr, fmt='.k', label='Data')
    plt.plot(xfit, linear_model([summary_linear['m'][0], summary_linear['b'][0]], xfit), color = 'red', label = 'best fit linear')
    plt.plot(xfit, quadratic_model([summary_quadratic['a2'][0], summary_quadratic['a1'][0], summary_quadratic['a0'][0]], xfit), color = 'blue', label = 'best fit quadratic')
    
    plt.title('fitted functions compared to real data')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()
    
    fig = corner.corner(samples_linear, labels=lin_labels,
                        truths=[summary_linear[label][0] for label in lin_labels])
    fig.suptitle("Posterior Distributions for Linear Model")
    plt.show()
    
    fig = corner.corner(samples_quadratic, labels=quad_labels,
                        truths=[summary_quadratic[label][0] for label in quad_labels])
    fig.suptitle("Posterior Distributions for Quadratic Model")
    plt.show()
    #TJ this commented out section displays the values and their +/- 1sigma limits very neatly, but this does not work well from the command line
    
    print("fitted m and b values based on 16%, 50%, and 84% uncertainties :")
    labels = ["m", "b"]
    for i in range(samples_linear.shape[1]):
        mcmc = np.percentile(samples_linear[:, i], [16, 50, 84])
        lower = mcmc[1] - mcmc[0]
        upper = mcmc[2] - mcmc[1]
        print(f'{labels[i]} : {mcmc[1]} (+{upper} / - {lower})')
    print() #TJ blank line for clarity
    print("fitted a2, a1, and a0 values based on 16%, 50%, and 84% uncertainties :")
    labels = ["a2", "a1", "a0"]
    for i in range(samples_quadratic.shape[1]):
        mcmc = np.percentile(samples_quadratic[:, i], [16, 50, 84])
        lower = mcmc[1] - mcmc[0]
        upper = mcmc[2] - mcmc[1]
        print(f'{labels[i]} : {mcmc[1]} (+{upper} / - {lower})')
    print() #TJ leave blank space for clarity
    print('Based on the information from the probability distribution, it is unlikely that the linear model will be sufficient to fit the data.')
    print('If the a2 values had +/- 1sigma values that straddled zero, then I would feel ok just using a linear fit model, but as it stands,')
    print('the a2 values are NOT consistent with zero, and therefore we do actually need a quadratic fit.')



