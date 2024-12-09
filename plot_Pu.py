#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 13:01:00 2023

@author: ryota
"""


import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
import matplotlib.colors
import matplotlib as mpl
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import matplotlib.ticker as mticker
import scipy.stats as stats

L=84
N=64

#pairs = [[i, j] for i in range(L * N) for j in range(i + 2, L * N)]
npairs=14442625

def FmeanVar(H, lp):
    ensvec=[1,2]
    # Initialize an empty dictionary to hold all trajectory data
    all_trajectory_data = {}

    for ens in ensvec:
        #filepath =  f'/Volumes/My Passport/Data/LLPS_Glass_Simulation/H={H}/Pd_fraction_H={H}_lp={lp}_ens={ens}.dat'
        filepath =  f'/Volumes/My Passport/Data/LLPS_Glass_Simulation/H={H}/Pd_Nodivision__H={H}_lp={lp}_ens={ens}.dat'
        
        # Read the data from each file
        with open(filepath, 'r') as file:
            Fdata = file.readlines()
        
        # Temporary dictionary to hold data for the current file
        trajectory_data = {}
        
        for line in Fdata:
            if line.strip():  # skip empty lines
                parts = line.split()
                trajectory_index = int(parts[0].split('-')[0]) + (ens - 1) * 10  # Adjust the index based on ens value
                time_point = float(parts[1])
                value = float(parts[2])
                
                if trajectory_index not in trajectory_data:
                    trajectory_data[trajectory_index] = {'times': [], 'values': []}
                
                trajectory_data[trajectory_index]['times'].append(time_point)
                trajectory_data[trajectory_index]['values'].append(value)
        
        # Merge current file data into all trajectory data
        all_trajectory_data.update(trajectory_data)
    
    # Align trajectories by shifting time
    for traj in all_trajectory_data.values():
        start_time = traj['times'][0]
        traj['times'] = [t - start_time for t in traj['times']]
    
    # Determine the number of time points from the first trajectory
    num_time_points = len(all_trajectory_data[next(iter(all_trajectory_data))]['times'])
    
    # Compute the mean and variance across all trajectories for each time point
    aveFvec = []
    varFvec = []
    for i in range(num_time_points):
        values_at_time = [traj['values'][i] for traj in all_trajectory_data.values() if i < len(traj['values'])]
        aveFvec.append(np.mean(values_at_time))
        varFvec.append(np.var(values_at_time))
    
    return np.array(list(all_trajectory_data[next(iter(all_trajectory_data))]['times'])), np.array(aveFvec), np.array(varFvec)

# Example usage:


lp = 0

times1, aveFvec1, varFvec1 = FmeanVar(1.58, lp)
times2, aveFvec2, varFvec2 = FmeanVar(4.75, lp)
times1=times1*1000
times2=times2*1000

def func1(t, a):
    return a*t**0
def func2(t, a, b):
    return a*t+b

indices0=np.where(times1 > 100000)[0]
indices1 = np.where(times1 > 1000000)[0]
indices2 = np.where(times2 > 10000000)[0]
popt1, pcov1 = curve_fit(func1,times1[indices1],aveFvec1[indices1])
popt2, pcov2 = curve_fit(func2,np.log10(times2[indices2]),np.log10(aveFvec2[indices2]))

############# Compute the confidence interval for the fits ##############
# Function to calculate standard errors and confidence intervals
def calculate_confidence_intervals(X, Y, m, b):
    n = len(X)
    predicted = m * X + b
    residuals = Y - predicted
    SER = np.sqrt(np.sum(residuals**2) / (n - 2))
    SE_b = SER / np.sqrt(np.sum((X - np.mean(X))**2))
    SE_a = SER * np.sqrt(1/n + np.mean(X)**2 / np.sum((X - np.mean(X))**2))
    alpha = 0.05
    t_crit = stats.t.ppf(1 - alpha/2, n - 2)
    CI_m = m - t_crit * SE_b, m + t_crit * SE_b
    CI_b = b - t_crit * SE_a, b + t_crit * SE_a
    return CI_m, CI_b

#m1=popt1[1]
#b1=popt1[0]
m2=popt2[0]
b2=popt2[1]
# Calculate confidence intervals for (m1, b1)
#X1 = np.log10(realtwvec[:])
#CI_m1, CI_b1 = calculate_confidence_intervals(X1, np.log10(tauc1vec[:]), m1, b1)

# Calculate confidence intervals for (m2, b2)
X2 = np.log10(times2[indices2])
CI_m2, CI_b2 = calculate_confidence_intervals(X2, np.log10(aveFvec2[indices2]), m2, b2)

#deltaCI_m1=abs(np.round((CI_m1[0]-CI_m1[1])/2,2))
deltaCI_m2=abs(np.round((CI_m2[0]-CI_m2[1])/2,5))
#print(f'error for m1={deltaCI_m1}')
print(f'error for m2={deltaCI_m2}')

plt.figure(figsize=(8, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
plt.plot(times1[indices0],aveFvec1[indices0],'-o',color='blue')
plt.plot(times1[indices0],np.ones(len(times1[indices0]))*popt1[0],'--',c='black')
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$P_{u}(t)$')
plt.xlabel(r'$t$')
#plt.title('H=1.58')
mticker.LogFormatterSciNotation(base=10.0, labelOnlyBase=False, minor_thresholds=None, linthresh=None)
plt.show()

significand, exponent = f"{abs(np.round(popt2[0],5)):e}".split('e')
exponent = int(exponent)  # Convert exponent to integer to remove leading zeros
# Construct the label string in LaTeX format
label_string = r"${} \times 10^{{{}}}$".format(significand[:4], exponent)
plt.figure(figsize=(8, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
plt.plot(times2[indices0],aveFvec2[indices0],'-s',color='red')
formatted_label = f"{abs(np.round(popt2[0],3)):10}"
plt.plot(times2[indices0],10**popt2[1]*times2[indices0]**popt2[0],'--',c='black',label=label_string)
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$P_{u}(t)$')
plt.xlabel(r'$t$')
#plt.title('H=4.75')
plt.legend(fontsize=20)
plt.show()

plt.figure(figsize=(8, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
plt.plot(times1[indices0],aveFvec1[indices0],'-o',color='blue')
plt.plot(times1[indices0],np.ones(len(times1[indices0]))*popt1[0],'--',c='tab:blue',label=r'$\alpha=0$')
plt.plot(times2[indices0],aveFvec2[indices0],'-s',color='red')
plt.plot(times2[indices0],10**popt2[1]*times2[indices0]**popt2[0],'--',c='tab:red',label=r'%s' %np.round(popt2[0],3))
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$P_{u}(t)$')
plt.xlabel(r'$t$')
#plt.title('H=4.75')
plt.legend()
plt.show()