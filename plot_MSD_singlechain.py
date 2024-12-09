#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 14:41:58 2024

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

L=84
N=1

#pairs = [[i, j] for i in range(L * N) for j in range(i + 2, L * N)]
npairs=3403

def FmeanVar(H, lp,tw):
    ensvec=[1,2]
    # Initialize an empty dictionary to hold all trajectory data
    all_trajectory_data = {}

    for ens in ensvec:
        filepath =  f'/Volumes/My Passport/Data/LLPS_Glass_Simulation/H={H}/MSD_pairs_SingleChain_tw={tw}_H={H}_lp={lp}_ens={ens}.dat'
        
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
freq=1
steps=1000
    
#twvec=[100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000]
twvec=[300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000]
lp = 4

realtwvec=steps*freq*np.array(twvec)


### Stretched exp Fit ###
def func(t, a, b):
    return -(np.array(t) /a)**b

### Stretched exp Fit ###
def func_fixed_afa(t, a):
    #return -(np.array(t) /a)**0.2 # epb=0 H=1.58
    #return -(np.array(t) /a)**0.15 # epb=2 
    return -(np.array(t) /a)**0.12 # epb=4 
    
def func_linear(t, a, b):
    return a*t+b
    
def func_linear_fixed_1(t,b):
    
    #return 0.43*t+b # H=1.58 lp=0
    #return 0.22*t+b # H=1.58 lp=2
    return 0.19*t+b # H=1.58 lp=4



def func_linear_fixed_2(t,b):
    
    #return 0.21*t+b # H=4.75 lp=0
    #return 0.20*t+b # H=4.75 lp=2
    return 0.16*t+b # H=4.75 lp=4




plt.style.use("/Users/ryota/Documents/Glass/Analysis/plot.mplstyle")
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
    times1=times1*steps*freq
    #popt1, pcov1 = curve_fit(func_linear, np.log10(times1[1:]), np.log10(aveFvec1[1:]))
    popt1, pcov1 = curve_fit(func_linear_fixed_1, np.log10(times1[1:]), np.log10(aveFvec1[1:]))
    #print(popt1,pcov1)
    std1vec=std1vec+[np.sqrt(pcov1[0])[0]]
    popt1vec=popt1vec+[popt1]
    #print(popt1vec)
    plt.plot(times1,aveFvec1,"-o",label='H=%s'%(1.58),c=cmap(i),markersize=13)
    #plt.plot(times1,10**(popt1[0])*times1**0.43,'-.', color="black")
    i+=1
plt.xlabel(r"$\tau$")
plt.ylabel(r"$\langle \Delta R^2_{t_w}(\tau)\rangle$")
plt.xscale('log')
plt.yscale('log')
#plt.title(r'$\epsilon_b=$'+str(lp)+r'$(k_BT)$')
#plt.title(r'$H=1.58$')
#plt.title('Low')
cbar=fig.colorbar(sm,ax=ax)
cbar.set_label(r'$t_w$',rotation=0,labelpad=20)
plt.grid(False)
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
    times2=times2*steps*freq
    #popt2, pcov2 = curve_fit(func_linear, np.log10(times2[1:]), np.log10(aveFvec2[1:]))
    popt2, pcov2 = curve_fit(func_linear_fixed_2, np.log10(times2[1:]), np.log10(aveFvec2[1:]))
    popt2vec=popt2vec +[popt2]
    std2vec=std2vec+[np.sqrt(pcov2[0])[0]]
    #print(popt1vec)
    plt.plot(times2,aveFvec2,'-s',label='H=%s'%(4.75),c=cmap(i),markersize=13)
    #plt.plot(times2,10**(popt2[0])*times2**0.21,'-.', color="black")
    i+=1
plt.xlabel(r"$\tau$")
plt.ylabel(r"$\langle \Delta R^2_{t_w}(\tau)\rangle$")
plt.xscale('log')
plt.yscale('log')
#plt.title(r'$\epsilon_b=$'+str(lp)+r'$(k_BT)$')
#plt.title(r'$H=1.58$')
#plt.title('Low')
cbar=fig.colorbar(sm,ax=ax)
cbar.set_label(r'$t_w$',rotation=0,labelpad=20)
plt.grid(False)
plt.show()

#################### Fitting ################################################

#slope1vec=np.array([popt1vec[i][0] for i in range(len(realtwvec))])
#slope2vec=np.array([popt2vec[i][0] for i in range(len(realtwvec))])
IC1vec=np.array([popt1vec[i][0] for i in range(len(realtwvec))])
IC2vec=np.array([popt2vec[i][0] for i in range(len(realtwvec))])

#print('mean1=',np.mean(slope1vec))
#print('mean2=',np.mean(slope2vec))




m1, b1  = np.polyfit(np.log10(realtwvec[:]), IC1vec[:], 1) #IC is in log scale alreaady 
print(m1,b1)
m2, b2 = np.polyfit(np.log10(realtwvec[:]), IC2vec[:], 1)
print(m2,b2)

############################ relaxation time ################################

plt.figure(figsize=(10, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
#plt.plot(realtwvec[:],np.mean(tauc1vec)*np.ones(len(realtwvec)),'--',c='black',label=r"$\alpha=$ %s"  %(round(m1,2)))
#plt.plot(realtwvec[:],np.ones(len(realtwvec))*np.mean(10**IC1vec),'--',c='black',label=r"$\alpha=$ %s"  %(round(m1,2)))
plt.plot(realtwvec[:],realtwvec[:]**m1*10**b1,'--',c='black',label=r"$\alpha=$ %s"  %(abs(round(m1,2))))
plt.errorbar(realtwvec[:], 10**IC1vec,np.array(std1vec)/2,marker='o',color='b',markersize=13,ls='none')
#plt.xscale('log')
plt.xlabel(r"$t_w$")
plt.ylabel(r"$D({t_w})$")
#plt.ylim([6*10**(-2),7.5*10**(-2)])
plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize=20)

plt.show()


plt.figure(figsize=(10, 6),dpi=400)
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
#line=np.linspace(twcom[0],twcom[-1],100)
#plt.plot(twcom,taucvec0,'o',c='black')
#plt.plot(realtwvec,10**IC2vec,'s',c='red',markersize=13)
plt.plot(realtwvec[:],realtwvec[:]**m2*10**b2,'--',c='black',label=r"$\alpha=$ %s"  %(abs(round(m2,2))))

plt.errorbar(realtwvec[:], 10**IC2vec,np.array(std2vec)/2,marker='s',color='r',markersize=13,ls='none')
#plt.xscale('log')
plt.xlabel(r"$t_w$")
plt.ylabel(r"$D({t_w})$")
plt.yscale('log')
plt.xscale('log')
#plt.ylim([10**4,10**5])
plt.legend(fontsize=20)
plt.show()

######################## Collpase of courves ###########################################
# c=np.round(np.array(realtwvec[:]),1)
# cmap = plt.get_cmap("jet", len(c))
# #norm = matplotlib.colors.BoundaryNorm(c,len(c))
# norm = LogNorm(vmin=min(realtwvec[:]), vmax=max(realtwvec[:]))
# sm = cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# plt.figure(figsize=(10, 6),dpi=400)
# plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')

# for i,tw in enumerate(twvec):
#     times1, aveFvec1, varFvec1 = FmeanVar(1.58, lp,tw)
#     plt.plot(np.array(times1)/tauc1vec[i],aveFvec1,label='H=%s'%(1.58),marker="o",markersize=13,c=cmap(i))
#     #plt.plot(space*np.array(deltimevec1)/taucvec0[i],aveFvec1/func(deltimevec1, taucvec0[i], betavec0[i]),label='H=%s'%(1.58),marker="s",markersize=10,c=color)
# plt.xlabel(r"$\tau/\tau_c$")
# plt.ylabel(r"$F_{t_w}(\tau/\tau_c)$")
# #plt.title(r'$\epsilon_b=$'+str(lp)+r'$(k_BT)$')
# #plt.title(r'$H=1.58$')
# #plt.title('Low')
# cbar=plt.colorbar(sm)
# cbar.set_label(r'$t_w$',rotation=0,labelpad=20)
# plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
# plt.show()


# c=np.round(np.array(realtwvec[:]),1)
# cmap = plt.get_cmap("jet", len(c))
# #norm = matplotlib.colors.BoundaryNorm(c,len(c))
# norm = LogNorm(vmin=min(realtwvec[:]), vmax=max(realtwvec[:]))
# sm = cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# plt.figure(figsize=(10, 6),dpi=400)
# plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')

# for i,tw in enumerate(twvec):
#     times2, aveFvec2, varFvec2 = FmeanVar(4.75, lp,tw)
#     plt.plot(np.array(times2)/tauc2vec[i],aveFvec2,label='H=%s'%(4.75),marker="s",markersize=13,c=cmap(i))
#     #plt.plot(space*np.array(deltimevec1)/taucvec0[i],aveFvec1/func(deltimevec1, taucvec0[i], betavec0[i]),label='H=%s'%(1.58),marker="s",markersize=10,c=color)
# plt.xlabel(r"$\tau/\tau_c$")
# plt.ylabel(r"$F_{t_w}(\tau/\tau_c)$")
# #plt.title(r'$\epsilon_b=$'+str(lp)+r'$(k_BT)$')
# #plt.title(r'$H=1.58$')
# #plt.title('Low')
# cbar=plt.colorbar(sm)
# cbar.set_label(r'$t_w$',rotation=0,labelpad=20)
# plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
# plt.show()


