#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:38:04 2023

@author: ryota
"""

import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import os
import sys
from scipy.optimize import curve_fit 
import scipy.stats as stats
# Updating the file path with the provided file
file_path1 = '/Volumes/My Passport/Data/LLPS_Glass_Simulation/H=1.58/PotEnergy2_seq1_lp=0_eps=2.dat'
file_path2 = '/Volumes/My Passport/Data/LLPS_Glass_Simulation/H=4.75/PotEnergy2_seq1_lp=0_eps=2.dat'


def loaddata(file_path):
    # Load the file and create arrays for each column
    # Skipping the first line which contains headers
    data = np.loadtxt(file_path, skiprows=1)

    # Extracting each column into separate arrays
    step = data[:, 0]
    time = data[:, 1]
    potential_energy = data[:, 2]
    total_energy = data[:, 3]
    
    return(step,abs(potential_energy))


time1,potential_energy_1=loaddata(file_path1)
time2,potential_energy_2=loaddata(file_path2)

def smooth_with_box_kernel(data, window_size):
    # Create a box kernel (all ones) of the specified window size
    kernel = np.ones(window_size) / window_size

    # Apply convolution to smooth the data
    smoothed_data = np.convolve(data, kernel, mode='same')

    return smoothed_data

# Choose a window size, for example, 5
window_size = 10

# Smooth the Potential Energy array
smoothed_potential_energy_1 = smooth_with_box_kernel(potential_energy_1, window_size)
smoothed_potential_energy_2 = smooth_with_box_kernel(potential_energy_2, window_size)





def func1(t, a):
    return a*t**0
def func2(t, a, b):
    return a*t+b

indices0 = np.where(time1 > 100000)[0][:-10]
indices1 = np.where(time1 > 10000000)[0][:-10]
indices2 = np.where(time2 > 10000000)[0][:-10]
popt1, pcov1 = curve_fit(func1,time1[indices1],smoothed_potential_energy_1[indices1])
popt2, pcov2 = curve_fit(func2,np.log10(time2[indices2]),np.log10(smoothed_potential_energy_2[indices2]))
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
X2 = np.log10(time2[indices2])
CI_m2, CI_b2 = calculate_confidence_intervals(X2, np.log10(smoothed_potential_energy_2[indices2]), m2, b2)

#deltaCI_m1=abs(np.round((CI_m1[0]-CI_m1[1])/2,2))
deltaCI_m2=abs(np.round((CI_m2[0]-CI_m2[1])/2,5))
#print(f'error for m1={deltaCI_m1}')
print(f'error for m2={deltaCI_m2}')


plt.figure(figsize=(8, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
plt.plot(time1[indices0],smoothed_potential_energy_1[indices0],'-',color='blue')
plt.plot(time1[indices0],popt1*np.ones(len(time1[indices0])),'--',c='black')
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$-U(k_BT)$')
plt.xlabel(r'$t$')
#plt.plot('legend')
plt.show()


significand, exponent = f"{abs(np.round(popt2[0],5)):e}".split('e')
exponent = int(exponent)  # Convert exponent to integer to remove leading zeros

# Construct the label string in LaTeX format
label_string = r"${} \times 10^{{{}}}$".format(significand[:4], exponent)
plt.figure(figsize=(8, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
plt.plot(time2[indices0],smoothed_potential_energy_2[indices0],'-',color='red')
#formatted_label = f"{abs(np.round(popt2[0],4)):10}"
plt.plot(time2[indices0],10**popt2[1]*time1[indices0]**popt2[0],'--',c='black',label=label_string)
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$-U(k_BT)$')
plt.xlabel(r'$t$')
plt.legend(fontsize=20)
plt.show()