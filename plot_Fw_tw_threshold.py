#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 20:42:10 2023

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
import matplotlib.ticker as ticker

L=84
N=64

#pairs = [[i, j] for i in range(L * N) for j in range(i + 2, L * N)]
npairs=14442625

def FmeanVar(H, lp,tw):
    ensvec=[1,2]
    # Initialize an empty dictionary to hold all trajectory data
    all_trajectory_data = {}

    for ens in ensvec:
        filepath =  f'/Volumes/My Passport/Data/LLPS_Glass_Simulation/H={H}/OverlapF_Nothreshold_tw={tw}_H={H}_lp={lp}_ens={ens}.dat'
        
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


freq=10 # frequency for reading dcd file 
steps=1000
    
#twvec=[100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000]
twvec=[300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000]
lp = 4

realtwvec=steps*freq*np.array(twvec) # 1 step in the simulation is 10**3


### Stretched exp Fit ###
def func(t, a, b):
    return -(np.array(t) /a)**b

### Stretched exp Fit ###
def func_fixed_afa(t, a):
    #return -(np.array(t) /a)**0.2 # epb=0 H=1.58
    #return -(np.array(t) /a)**0.15 # epb=2 
    return -(np.array(t) /a)**0.12 # epb=4 
    
    

# def func1(t, a):
#     return a*t**0
# def func2(t, a, b):
#     return a*t+b



plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
fig,ax=plt.subplots(1, 1, figsize=(10, 6),dpi=400)
#plt.style.use('plot.mplstyle') 
c=np.round(np.array(realtwvec),1)
cmap = plt.get_cmap("jet", len(c))
#norm = matplotlib.colors.BoundaryNorm(c,len(c))
norm = LogNorm(vmin=min(realtwvec), vmax=max(realtwvec))
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
i=0
popt1vec=[]
std1vec=[]
for tw in twvec:
    times1, aveFvec1, varFvec1 = FmeanVar(1.58, lp,tw)
    times1=steps*freq*times1
    
    #popt1, pcov1 = curve_fit(func, times1, np.log(aveFvec1),maxfev=10000)
    popt1, pcov1 = curve_fit(func_fixed_afa, times1, np.log(aveFvec1))
    std1vec=std1vec+[np.sqrt(pcov1[0])[0]]
    
    popt1vec=popt1vec+[popt1]#popt2vec=popt2vec +[popt2]
    #print(popt1vec)
    plt.plot(times1,aveFvec1,"-o",label='H=%s'%(1.58),c=cmap(i),markersize=13)
    #plt.plot(times1,np.exp(func_fixed_afa(np.array(times1), *popt1)),'-.', color="black")
    i+=1
plt.xlabel(r"$\tau$")
plt.ylabel(r"$F_{t_w}(\tau)$")
#plt.xscale('log')
#plt.title(r'$\epsilon_b=$'+str(lp)+r'$(k_BT)$')
#plt.title(r'$H=1.58$')
#plt.title('Low')
cbar=fig.colorbar(sm,ax=ax)
cbar.set_label(r'$t_w$',rotation=0,labelpad=20)
plt.grid(False)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.gca().ticklabel_format(useMathText=True)
plt.show()

plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
fig,ax=plt.subplots(1, 1, figsize=(10, 6),dpi=400)
#plt.style.use('plot.mplstyle') 
c=np.round(np.array(realtwvec),1)
cmap = plt.get_cmap("jet", len(c))
#norm = matplotlib.colors.BoundaryNorm(c,len(c))
norm = LogNorm(vmin=min(realtwvec), vmax=max(realtwvec))
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
i=0
popt2vec=[]
std2vec=[]
for tw in twvec:
    times2, aveFvec2, varFvec2 = FmeanVar(4.75, lp,tw)
    times2=steps*freq*times2
    #popt2, pcov2 = curve_fit(func, times2, np.log(aveFvec2),maxfev=10000)
    popt2, pcov2 = curve_fit(func_fixed_afa, times2, np.log(aveFvec2))
    popt2vec=popt2vec +[popt2]
    std2vec=std2vec+[np.sqrt(pcov2[0])[0]]
    #print(popt1vec)
    plt.plot(times2,aveFvec2,'-s',label='H=%s'%(4.75),c=cmap(i),markersize=13)
    #plt.plot(times2,np.exp(func_fixed_afa(np.array(times2), *popt2)),'-.', color="black")
    i+=1
plt.xlabel(r"$\tau$")
plt.ylabel(r"$F_{t_w}(\tau)$")
#plt.xscale('log')
#plt.title(r'$\epsilon_b=$'+str(lp)+r'$(k_BT)$')
#plt.title(r'$H=1.58$')
#plt.title('Low')
cbar=fig.colorbar(sm,ax=ax)
cbar.set_label(r'$t_w$',rotation=0,labelpad=20)
plt.grid(False)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.gca().ticklabel_format(useMathText=True)
plt.show()

#################### Fitting ################################################

tauc1vec=np.array([popt1vec[i][0] for i in range(len(realtwvec))])
tauc2vec=np.array([popt2vec[i][0] for i in range(len(realtwvec))])
#afa1vec=np.array([popt1vec[i][1] for i in range(len(realtwvec))])
#afa2vec=np.array([popt2vec[i][1] for i in range(len(realtwvec))])





m1, b1 = np.polyfit(np.log10(realtwvec[:]), np.log10(tauc1vec[:]), 1)
print(m1,b1)
m2, b2 = np.polyfit(np.log10(realtwvec[:]), np.log10(tauc2vec[:]), 1)
print(m2,b2)

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


# Calculate confidence intervals for (m1, b1)
X1 = np.log10(realtwvec[:])
CI_m1, CI_b1 = calculate_confidence_intervals(X1, np.log10(tauc1vec[:]), m1, b1)

# Calculate confidence intervals for (m2, b2)
X2 = np.log10(realtwvec[:])
CI_m2, CI_b2 = calculate_confidence_intervals(X2, np.log10(tauc2vec[:]), m2, b2)

deltaCI_m1=abs(np.round((CI_m1[0]-CI_m1[1])/2,2))
deltaCI_m2=abs(np.round((CI_m2[0]-CI_m2[1])/2,2))
print(f'error for m1={deltaCI_m1}')
print(f'error for m2={deltaCI_m2}')
############################ relaxation time ################################
def int_decimal_formatter(x, pos):
    #return f"{x:.1f}"  # Adjust '.2f' to change the number of decimal places
    return f"{int(x)}"
multiple=0
#multiple=4 #lp=0

plt.figure(figsize=(10, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
#line=np.linspace(twcom[0],twcom[-1],100)
#plt.plot(twcom,taucvec0,'o',c='black')
#plt.plot(realtwvec,tauc1vec,'o',c='blue',markersize=13)
plt.errorbar(realtwvec, 10**(-multiple)*tauc1vec, 10**(-multiple)*np.array(std1vec)/2,marker='o',color='blue',markersize=13,ls='none')
#plt.plot(realtwvec[:],10**(-multiple)*np.mean(tauc1vec)*np.ones(len(realtwvec)),'--',c='black')
plt.plot(realtwvec[:],10**b1*realtwvec[:]**m1,'--',c='black',label=r"$\mu=$ %s"  %(round(m1,2)))
#plt.xscale('log')
plt.xlabel(r"$t_w$")
plt.ylabel(r"$\tau_c$")
#plt.ylim([10**(-multiple)*3*10**4,10**(-multiple)*5.5*10**4]) # epb=0
plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize=30)

#plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(int_decimal_formatter))
#plt.gca().yaxis.set_minor_formatter(ticker.FuncFormatter(int_decimal_formatter))
#ax.xaxis.set_minor_formatter(ticker.NullFormatter())  # Hide labels for minor ticks
# Annotate the multiplier
#ax = plt.gca()
#ax.annotate('x$10^{%s}$'%multiple, xy=(0.01, 1.02), xycoords='axes fraction', fontsize=30, verticalalignment='bottom')

plt.show()

multiple=0
#multiple=8 #lp=0
plt.figure(figsize=(10, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
#line=np.linspace(twcom[0],twcom[-1],100)
#plt.plot(twcom,taucvec0,'o',c='black')
#plt.plot(realtwvec,tauc2vec,'s',c='red',markersize=13)
plt.errorbar(realtwvec, 10**(-multiple)*tauc2vec,  10**(-multiple)*np.array(std2vec)/2,marker='s',color='red',markersize=13,ls='none')
plt.plot(realtwvec[:], 10**(-multiple)*10**b2*realtwvec[:]**m2,'--',c='black',label=r"$\mu=$ %s"  %(round(m2,2)))
#plt.xscale('log')
plt.xlabel(r"$t_w$")
plt.ylabel(r"$\tau_c$")
plt.yscale('log')
plt.xscale('log')
#plt.ylim([10**4,10**5])
#plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(int_decimal_formatter))
#plt.gca().yaxis.set_minor_formatter(ticker.FuncFormatter(int_decimal_formatter))
#ax.xaxis.set_minor_formatter(ticker.NullFormatter())  # Hide labels for minor ticks
# Annotate the multiplier
#ax = plt.gca()
#ax.annotate('x$10^{%s}$'%multiple, xy=(0.01, 1.02), xycoords='axes fraction', fontsize=30, verticalalignment='bottom')

plt.legend(fontsize=30)
plt.show()

######################## Collpase of courves ###########################################
c=np.round(np.array(realtwvec[:]),1)
cmap = plt.get_cmap("jet", len(c))
#norm = matplotlib.colors.BoundaryNorm(c,len(c))
norm = LogNorm(vmin=min(realtwvec[:]), vmax=max(realtwvec[:]))
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
fig,ax=plt.subplots(figsize=(10, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')

for i,tw in enumerate(twvec):
    times1, aveFvec1, varFvec1 = FmeanVar(1.58, lp,tw)
    plt.plot(np.array(times1)/tauc1vec[i],aveFvec1,label='H=%s'%(1.58),marker="o",markersize=13,c=cmap(i))
    #plt.plot(space*np.array(deltimevec1)/taucvec0[i],aveFvec1/func(deltimevec1, taucvec0[i], betavec0[i]),label='H=%s'%(1.58),marker="s",markersize=10,c=color)
plt.xlabel(r"$\tau/\tau_c$")
plt.ylabel(r"$F_{t_w}(\tau/\tau_c)$")
#plt.title(r'$\epsilon_b=$'+str(lp)+r'$(k_BT)$')
#plt.title(r'$H=1.58$')
#plt.title('Low')
cbar=plt.colorbar(sm,ax=ax)
cbar.set_label(r'$t_w$',rotation=0,labelpad=20)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.gca().ticklabel_format(useMathText=True)
plt.show()


c=np.round(np.array(realtwvec[:]),1)
cmap = plt.get_cmap("jet", len(c))
#norm = matplotlib.colors.BoundaryNorm(c,len(c))
norm = LogNorm(vmin=min(realtwvec[:]), vmax=max(realtwvec[:]))
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
fig,ax=plt.subplots(figsize=(10, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')

for i,tw in enumerate(twvec):
    times2, aveFvec2, varFvec2 = FmeanVar(4.75, lp,tw)
    plt.plot(np.array(times2)/tauc2vec[i],aveFvec2,label='H=%s'%(4.75),marker="s",markersize=13,c=cmap(i))
    #plt.plot(space*np.array(deltimevec1)/taucvec0[i],aveFvec1/func(deltimevec1, taucvec0[i], betavec0[i]),label='H=%s'%(1.58),marker="s",markersize=10,c=color)
plt.xlabel(r"$\tau/\tau_c$")
plt.ylabel(r"$F_{t_w}(\tau/\tau_c)$")
#plt.title(r'$\epsilon_b=$'+str(lp)+r'$(k_BT)$')
#plt.title(r'$H=1.58$')
#plt.title('Low')
cbar=plt.colorbar(sm,ax=ax)
cbar.set_label(r'$t_w$',rotation=0,labelpad=20)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.gca().ticklabel_format(useMathText=True)
plt.show()


