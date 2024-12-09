#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  8 14:40:37 2021

@author: ryota
"""


import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import random
import math
L=84
N=64
box=100

def trajectory(H,ens,seq,lp):
    
    if H==1.58:
        seq=1
    #file='/home/rt24494/Research/IDP/Simulation/Glass/H=%s/output%s_seq%s_lp=%s_eps=2.dcd' %(H,ens,seq,lp)
    #topo=md.load('/home/rt24494/Research/IDP/Simulation/Glass/H=%s/topology.pdb' %H).topology 
    file='/Volumes/My Passport/Data/LLPS_Glass_Simulation/H=%s/output%s_seq%s_lp=%s_eps=2.dcd' %(H,ens,seq,lp)
    topo=md.load('/Volumes/My Passport/Data/LLPS_Glass_Simulation/H=1.58/initial_lp=0.pdb').topology 
    
    if H==1.58 or H==4.75:
        st=1000
    else:
        st=10
    traj=md.load(file,stride=st,top=topo)
    #traj=traj[200:]
    traj=traj[-1]
    time = traj.time
    
    return traj,time

def Rg(traj):   
    Rg=md.compute_rg(traj)
    return Rg


def gytensor(traj):
    gy=md.compute_gyration_tensor(traj)
    return gy




def nematic(traj):          
    chain_indices = [[n+x for x in range(L)] for n in range(0, int(L*N), L)]
    #rigid_ini=40
    #rigid_end=50
    #chain_indices = [[n+x for x in range(rigid_ini,rigid_end+1)] for n in range(0, int(L*N), L)]
    S2 = md.compute_nematic_order(traj, indices=chain_indices)
    return S2

def nemave(H,ens,seq,lp):
    #i=0
    vec=[]
    #for ens in ensvec:
    #    for seq in seqvec:    
    (traj,time)= trajectory(H,ens,seq,lp)
    nemavec=nematic(traj)            
    avenem=np.mean(nemavec)            
    vec=vec+[avenem]
    #i=i+1
            
    nemave=np.mean(vec)
    nemstd=np.std(vec)
    
    return (nemave,nemstd)





#Hvec=[1.58,2.0,2.5,3.0,3.5,4.0,4.5,4.75]
Hvec=[1.58,4.75,2.0,3.0]
lpvec=[0,1,2,3,4,5]
seqvec=[1]
ensvec=[1]
#Tvec=[2.0,1.8,0.2]
#sigstvec=[1.0,8.0]
#ensvec=[1,2]


#Rg_ensmean=np.zeros(len(ensvec)*len(seqvec))
nema_ensmean=np.zeros(len(ensvec)*len(seqvec))
#Rg_ensstd=np.zeros(len(ensvec)*len(seqvec))
nema_ensstd=np.zeros(len(ensvec)*len(seqvec))
#RgMat_mean=np.zeros((len(lpvec),len(Hvec)))
nemaMat_mean=np.zeros((len(lpvec),len(Hvec)))
#RgMat_ste=np.zeros((len(lpvec),len(Hvec)))
nemaMat_ste=np.zeros((len(lpvec),len(Hvec)))


j=0
for lp in lpvec:
    k=0
    for H in Hvec:
        i=0
        for ens in ensvec:
            for seq in seqvec:
                
                ave,std=nemave(H,ens,seq,lp)                              
                nema_ensmean[i]=ave
                #Rg_ensstd[i]=stdRg
                #k2_ensstd[i]=stdk2
            
                i=i+1
        
        #RgMat_mean[j,k]=np.mean(Rg_ensmean)
        nemaMat_mean[j,k]=np.mean(nema_ensmean)
        #RgMat_ste[j,k]=np.std(Rg_ensmean)/np.sqrt(len(ensvec)*len(seqvec))
        nemaMat_ste[j,k]=np.std(nema_ensmean)/np.sqrt(len(ensvec)*len(seqvec))
        k=k+1
    j=j+1
    

#psivec= [i/l for i in lpvec]

# for ii in range(len(lpvec)):
#     plt.errorbar(Tvec,RgMat_mean[ii,:],yerr=RgMat_ste[ii,:],marker="o",label=r"$\psi_{st}=$"+str(psivec[ii]))
#     plt.legend()
# plt.xlabel(r"$T^*$")
# plt.ylabel(r"$R_g$")
# plt.show()


# ### Rg normalized by sigst######
# for ii in range(len(lpvec)):
#     plt.errorbar(Tvec,RgMat_mean[ii,:]/lpvec[ii],yerr=RgMat_ste[ii,:],marker="o",label=r"$\psi_{st}=$"+str(psivec[ii]))
#     plt.legend()
# plt.xlabel(r"$T^*$")
# plt.ylabel(r"$R_g$")
# plt.show()
# #################################

# for ii in range(len(Tvec)):
#     plt.errorbar(psivec,RgMat_mean[:,ii],yerr=RgMat_ste[:,ii],marker="o",label=r"$T^*=$"+str(Tvec[ii]))
#     plt.legend()
# plt.xlabel(r"$\psi_{st}$")
# plt.ylabel(r"$R_g$")
# plt.show()

# for ii in range(len(lpvec)):
#     plt.errorbar(Tvec,k2Mat_mean[ii,:],yerr=k2Mat_ste[ii,:],marker="o",label=r"$\psi_{st}=$"+str(psivec[ii]))
#     plt.legend()
# plt.xlabel(r"$T^*$")
# plt.ylabel(r"$\kappa^2$")
# plt.show()

# for ii in range(len(Tvec)):
#     plt.errorbar(psivec,k2Mat_mean[:,ii],yerr=k2Mat_ste[:,ii],marker="o",label=r"$T^*=$"+str(Tvec[ii]))
#     plt.legend()
# plt.xlabel(r"$\psi_{st}$")
# plt.ylabel(r"$\kappa^2$")
# plt.show()

#plt.plot(time,k2vec)
#plt.show()

import matplotlib.colors
plt.style.use('plot.mplstyle')
# c=lpvec
# cmap = plt.get_cmap("jet", len(c))
# norm = matplotlib.colors.BoundaryNorm(c,len(c))
# sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
# sm.set_array([])  # this line may be ommitted for matplotlib >= 3.1
# fig, ax = plt.subplots(dpi=400)
# for ii in range(len(lpvec)):
#     ax.errorbar(Hvec,nemaMat_mean[ii,:],yerr=nemaMat_ste[ii,:],c=cmap(ii),marker="o",label=r"$\epsilon_b=$"+str(lpvec[ii])+r"$(k_BT)$")
# #cbar=fig.colorbar(sm, ticks=c)
# #cbar.set_label(r"$\epsilon_b(k_BT)$",rotation=0,labelpad=15)
# plt.xlabel(r"$H$(bit)")
# plt.ylabel(r"$S$")
# plt.legend(fontsize=20,markerscale=2)
# plt.savefig('SvsH_maxmin.png', bbox_inches='tight',dpi = 450)
# plt.show()

#c=Hvec
#cmap = plt.get_cmap("jet", len(c))
#norm = matplotlib.colors.BoundaryNorm(c,len(c))
#sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
#sm.set_array([])  # this line may be ommitted for matplotlib >= 3.1
fig, ax = plt.subplots(dpi=400)
for ii in range(len(Hvec)):
    if Hvec[ii]==1.58:
        mk="o"
    elif Hvec[ii]==4.75:
    	mk="s"
    elif Hvec[ii]==2.0:
        mk="^"
    elif Hvec[ii]==3.0:
        mk="v"
    ax.errorbar(lpvec,nemaMat_mean[:,ii],yerr=nemaMat_ste[:,ii],label=r'$H_3=$%s'%Hvec[ii],marker=mk)
#cbar=fig.colorbar(sm, ticks=c)
#cbar.set_label(r'$H$(bit)',rotation=0,labelpad=15)
plt.xlabel(r"$\epsilon_b(k_BT)_revision$")
plt.ylabel(r"$S$")
plt.legend(fontsize=20,markerscale=2)
#plt.savefig('Svslp_maxmin_revise.png', bbox_inches='tight',dpi = 450)
plt.show()



# ##########Density Plot ############
# import seaborn as sns


# fig, ax = plt.subplots(dpi=400)
# im = ax.imshow(nemaMat_mean,cmap='jet')
# cbar = fig.colorbar(im)
# cbar.set_label(r'$S$',rotation=0,labelpad=15)

# # We want to show all ticks...
# ax.set_xticks(np.arange(len(Hvec)))
# ax.set_yticks(np.arange(len(lpvec)))
# # ... and label them with the respective list entries

# ax.xaxis.set_ticks_position('bottom')
# ax.xaxis.set_label_position('bottom')

# ax.set_xticklabels(Hvec)
# ax.set_yticklabels(lpvec)

# # Rotate the tick labels and set their alignment.
# plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
#          rotation_mode="anchor")

# plt.xlabel(r"$H$(bit)")
# plt.ylabel(r"$\epsilon_b(k_BT)$")


# fig.tight_layout()
# plt.savefig('S_H_lp_maxmin.png', bbox_inches='tight',dpi = 450)
# plt.show()





