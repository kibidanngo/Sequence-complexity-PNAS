#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 21:19:15 2023

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

# Memo: For H=1.58 the trajectory is divided to 100 starting t=100, For H=4.75 the trajectory is divided to 3 starting t=100


L=84
N=64

#pairs = [[i, j] for i in range(L * N) for j in range(i + 2, L * N)]
npairs=14442625

def FmeanVar(H, lp):
    ensvec=[1,2]
    # Initialize an empty dictionary to hold all trajectory data
    all_trajectory_data = {}

    for ens in ensvec:
        #filepath =  f'/Volumes/My Passport/Data/LLPS_Glass_Simulation/H={H}/OverlapF_trajectroy_block_threshold_H={H}_lp={lp}_ens={ens}.dat'        
        filepath =  f'/Volumes/My Passport/Data/LLPS_Glass_Simulation/H={H}/OverlapF_trajectroy_block_H={H}_lp={lp}_ens={ens}.dat'
        #filepath =  f'/Volumes/My Passport/Data/LLPS_Glass_Simulation/H={H}/ISF/OverlapF_trajectroy_block_H={H}_lp={lp}_ens={ens}.dat'
        #filepath =  f'/Volumes/My Passport/Data/LLPS_Glass_Simulation/H={H}/ISF_H={H}_lp={lp}_ens={ens}.dat'
        
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


lp = 4

times1, aveFvec1, varFvec1 = FmeanVar(1.58, lp)
times2, aveFvec2, varFvec2 = FmeanVar(4.75, lp)
times1=times1*1000
times2=times2*1000

def func1(t, a):
    return a*t**0
def func2(t, a, b):
    return a*t+b

indices1 = np.where(times1 > 1000)[0]
indices2 = np.where(times2 > 100)[0]
popt1, pcov1 = curve_fit(func1,times1[indices1],aveFvec1[indices1])
popt2, pcov2 = curve_fit(func2,np.log10(times2[indices2]),np.log10(aveFvec2[indices2]))


plt.figure(figsize=(8, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
plt.plot(times1[1:],npairs*varFvec1[1:],'-o',color='blue',markersize=10)

#plt.yscale('log')
plt.ylabel(r'$\chi_4(\tau)$')
plt.xlabel(r'$\tau$')
#plt.title('H=1.58')
plt.ticklabel_format(style='sci', axis='y', scilimits=(3,3))
plt.gca().ticklabel_format(useMathText=True)
plt.xscale('log')
plt.show()

plt.figure(figsize=(8, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
plt.plot(times1[1:],aveFvec1[1:],'-o',color='blue',markersize=10)
plt.xscale('log')
#plt.yscale('log')
plt.ylabel(r'$F_{t_w}(\tau)$')
plt.xlabel(r'$\tau$')
#plt.title('H=1.58')
plt.show()



plt.figure(figsize=(8, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
plt.plot(times2[1:],npairs*varFvec2[1:],'-s',color='red',markersize=10)


#plt.yscale('log')
plt.ylabel(r'$\chi_4(\tau)$')
plt.xlabel(r'$\tau$')
#plt.title('H=4.75')
#plt.legend()
plt.ticklabel_format(style='sci', axis='y', scilimits=(3,3))
plt.gca().ticklabel_format(useMathText=True)
plt.xscale('log')
plt.show()


plt.figure(figsize=(8, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
plt.plot(times2[1:],aveFvec2[1:],'-s',color='red',markersize=10)

plt.xscale('log')
#plt.yscale('log')
plt.ylabel(r'$F_{t_w}(\tau)$')
plt.xlabel(r'$\tau$')
#plt.title('H=4.75')
#plt.legend()
plt.show()

plt.figure(figsize=(8, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
plt.plot(times1[1:],aveFvec1[1:],'-o',color='blue',label='H=1.58')
plt.plot(times2[1:],aveFvec2[1:],'-o',color='red',label='H=4.75')
plt.xscale('log')
#plt.yscale('log')
plt.ylabel(r'$F_{t_w}(\tau)$')
plt.xlabel(r'$\tau$')
#plt.title('H=4.75')
plt.legend()
plt.show()

