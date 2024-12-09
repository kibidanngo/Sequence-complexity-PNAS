#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 13:26:01 2023

@author: ryota
"""





import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import os
import sys 

L=84
N=64

twval=int(sys.argv[1])
ens=int(sys.argv[2])
H=float(sys.argv[3])
lp=int(sys.argv[4])
#twval=100
#ens=1
#H=1.58
#lp=0
seq=1


def trajectory(H,lp,ens,seq):
    if H==1.58:
        seq=1
    
   
    file='/home/rt24494/Research/IDP/Simulation/Glass/H=%s/output%s_seq%s_lp=%s_eps=2.dcd' %(H,ens,seq,lp)  
    #file='/Volumes/My Passport/Data/LLPS_Glass_Simulation/H=1.58/output1_seq1_lp=0_eps=2.dcd'  
    #topo=md.load('/Volumes/My Passport/Data/LLPS_Glass_Simulation/H=1.58/topology.pdb' )
    topo=md.load('/home/rt24494/Research/IDP/Simulation/Glass/H=%s/topology.pdb' % H )
    
    
    traj=md.load(file,stride=10,top=topo)
    #traj=traj[100:]
    time = traj.time
    
    return traj,time



#cutoff=1.0
deltaT=10




def file(tw,H,lp,ens,seq):
    #file= 'OverlapFfullshort_H=%s_lp=%s.dat'%(H,lp)
    file= 'MSD_pairs_tw=%s_H=%s_lp=%s_ens=%s.dat'%(tw,H,lp,ens)
    return file 



def compute_MSD(trajectory,pairslist):
    
    num_frames = trajectory.n_frames 

    msd_vec=np.zeros(num_frames)

    for t in range(num_frames):
        #current_positions = trajectory.xyz[t]
        dist0tvec= (pairslist[t] - pairslist[0])**2  
        msd=np.mean(dist0tvec)
        msd_vec[t]=msd
        

   
        

    return msd_vec  # Return the real part of F_s(q, t)




def calculateF(deltaT,ini,las,H,lp,ens,seq):
    traj,time=trajectory(H,lp,ens,seq)
    pairs=[[i,j] for i in range(L*N) for j in range (i+2,L*N)]
    for i in range(N):
        pairs = pairs + [[((L-1)+L*i)%(N*L),((L-1)+L*i+1)%(N*L)]]
            
        
    pairs=np.array(pairs)
    
    
    filename = file(ini,H,lp,ens,seq)
    if os.path.exists(filename):
            os.remove(filename)
    ntraj=1
    for ii in range(ini,las,deltaT):
    
        traji=traj[ii:ii+deltaT]
        timei=time[ii:ii+deltaT]
        #dist0tvec=np.zeros((len(timei),npairs))
        #dist0tvec=np.zeros(len(timei),dtype=object)
        #pairdistless=np.zeros(len(timei))
        #pairslist=md.squeucompute_distances(traji, pairs)
        #all_pairs_distances = md.compute_distances(traji, pairs)
        #threshold = 5.0 # for near particles 
        #initial_pairs_within_threshold = all_pairs_distances[0] < threshold
        
        pairslist=md.compute_distances(traji, pairs)
        
        #Fw_data = np.zeros(len(traji))
        
        
        
        #for ii in range(len(timei)):            
        #    pairslist.append(np.where(initial_pairs_within_threshold, all_pairs_distances[ii], np.nan))
        pairslist=np.array(pairslist) 
        
        MSD_vec=compute_MSD(traji,pairslist)
        
        #for i in range(len(timei)):            
            #dist0tvec =pairslist[i] - pairslist[0]
            #print(np.isnan(dist0tvec))
           
            #Fw_data[i]=np.mean(np.abs(dist0tvec) < cutoff)
            #print(Fw_data[i])
        

        with open(filename, 'a') as f:
            #counttime = 1
            for idx,item in enumerate(MSD_vec):
                f.write(f"{ntraj}-th {timei[idx]} {item}\n")
        ntraj+=1



tau=100

calculateF(deltaT,twval,twval+tau,H,lp,ens,seq) # calculate F and wirte data in text file 






