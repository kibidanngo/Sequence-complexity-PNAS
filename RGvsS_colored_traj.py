

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 08:49:37 2020

@author: ryota
"""
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import csv

L=84
N=64
box=100


def trajectory(H,ens,seq,lp):
    
    if H==1.58:
        seq=1
    file='/home/rt24494/Research/IDP/Simulation/Glass/H=%s/output%s_seq%s_lp=%s_eps=2.dcd' %(H,ens,seq,lp)
    topo=md.load('/home/rt24494/Research/IDP/Simulation/Glass/H=%s/initial_lp=0.pdb' %H).topology 
    
    if H==1.58 or H==4.75:
        st=100
    else:
        st=1
    traj=md.load(file,stride=st,top=topo)
    traj=traj[-10000:]
    time = traj.time
    
    return traj,time

def nematic(traj):
      
    #dist=md.compute_rg(traj) ## Rg
    
    
    chain_indices = [[n+x for x in range(L)] for n in range(0, int(L*N), L)]
    #rigid_ini=40
    #rigid_end=50
    #chain_indices = [[n+x for x in range(rigid_ini,rigid_end+1)] for n in range(0, int(L*N), L)]
    S2 = md.compute_nematic_order(traj, indices=chain_indices)
    return S2



def Rg(traj):   
    Rg=md.compute_rg(traj)
    return Rg


def Rgave(H,ens,seq,lp):
    #i=0
    vec=[]
    #for ens in ensvec:
    #    for seq in seqvec:    
    (traj,time)= trajectory(H,ens,seq,lp)
    Rgavevec=Rg(traj)            
    aveRg=np.mean(Rgavevec)            
    vec=vec+[aveRg]
    #i=i+1
            
    Rgave=np.mean(vec)
    Rgstd=np.std(vec)
    
    return (Rgave,Rgstd)     

#nemavec=nematic(traj)
#gamavec=[gamma(traj,i) for i in range(len(time))]


def RgNemavec(H,ens,seq,lp):
    (traj,time)= trajectory(H,ens,seq,lp)
    nemavec=nematic(traj)
    Rgvec=Rg(traj)
    
    return (time,nemavec,Rgvec)
             


# for H in Hvec:
    
#     if H==1.58:
#         (traj1,time1)= trajectory(H,lp,ens,seq)
#         nemavec=nematic(traj1)
#         Rgvec=Rg(traj1)
#         Rgvec1=Rgvec
#         nemavec1=nemavec
        
#     if H==4.75:
#         (traj2,time2)= trajectory(H,lp,ens,seq)
#         nemavec=nematic(traj2)
#         Rgvec=Rgvec(traj2)
#         Rgvec2=Rgvec
#         nemavec2=nemavec
        

#H=1.58
#lpvec=[0]
lp=2
ens=1
seq=1
#colors=['Blues','Purples','Greens','Reds']
#ii=0
#for lp in lpvec:
#    (time,nemavec,Rgvec)=RgNemavec(H,ens,seq,lp)
#    #fig, ax = plt.subplots()
#    plt.scatter(Rgvec/box, nemavec, c=time,cmap=colors[ii],label=r'$\epsilon_b=$' + str(lp)+r'$(k_BT)$')
#    ii+=1

(time1,nemavec1,Rgvec1)=RgNemavec(1.58,ens,seq,lp)
(time2,nemavec2,Rgvec2)=RgNemavec(4.75,ens,seq,lp)


def write_data_to_file(filename, time, nemavec, Rgvec):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Time', 'Nemavec', 'Rgvec'])
        for t, nema, rg in zip(time, nemavec, Rgvec):
            writer.writerow([t, nema, rg])

# Assuming 'box' is defined earlier in your code
write_data_to_file(f'phase_H=1.58_lp={lp}_ens={ens}_seq={seq}.csv', time1, nemavec1, Rgvec1/box)
write_data_to_file(f'phase_H=4.75_lp={lp}_ens={ens}_seq={seq}.csv', time2, nemavec2, Rgvec2/box)

# fig, ax = plt.subplots()
# plt.figure()
# plt.style.use('plot.mplstyle')
# plt.scatter(Rgvec2/box, nemavec2, c=time2,cmap = 'Reds',label=r'$H_3=$' + str(4.75))
# plt.scatter(Rgvec1/box, nemavec1, c=time1,cmap = 'Blues',label=r'$H_3=$' + str(1.58))
# # #plt.xlim([8,11])
# plt.xlabel(r"$R_G/L$")
# plt.ylabel("S")
# plt.xlim([0.09,0.115])
# #ax.legend(loc="best")
# #cbar = plt.colorbar()
# #plt.legend()
# #plt.show()

# plt.savefig('SvsRG_lp=%s_seq=%s_ens=%s.png'%(lp,seq,ens), bbox_inches='tight',dpi = 450)

