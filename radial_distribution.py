#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 13:26:23 2020

@author: ryota
"""


import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import loadtxt

L=84
N=64
box=100



def trajectory(H,ens,seq,lp):
    if H==1.58:
        seq=1
    
    #file='/Users/ryota/mnt/lakers/Research/IDP/Simulation/L=84_N=64_GPU/H=%s/output%s_seq%s_lp=%s_eps=2_.dcd' %(H,ens,seq,lp)    
    #topo=md.load('/Users/ryota/mnt/lakers/Research/IDP/Simulation/L=84_N=64_GPU/H=%s/initial.pdb' % H )
    
    #file='/home/rt24494/Research/IDP/Simulation/Glass/H=%s/output%s_seq%s_lp=%s_eps=2.dcd' %(H,ens,seq,lp)
    #topo=md.load('/home/rt24494/Research/IDP/Simulation/Glass/H=%s/topology.pdb'%H).topology 
    
    
      
    file='/Volumes/My Passport/Data/LLPS_Glass_Simulation/H=%s/output%s_seq%s_lp=%s_eps=2.dcd' %(H,ens,seq,lp)
    topo=md.load('/Volumes/My Passport/Data/LLPS_Glass_Simulation/H=%s/initial_lp=0.pdb'%H).topology 
    
    if H==1.58 or H==4.75:
        st=1000
    else:
        st=10
    traj=md.load(file,stride=st,top=topo)
    traj=traj[200:]
    time = traj.time
    
    return traj,time

# def trajectory(H,ens,seq,lp):
#     if H==1.58:
#         seq=1
    
#     #file='/Users/ryota/mnt/lakers/Research/IDP/Simulation/L=84_N=64_GPU/H=%s/output%s_seq%s_lp=%s_eps=2_.dcd' %(H,ens,seq,lp)    
#     #topo=md.load('/Users/ryota/mnt/lakers/Research/IDP/Simulation/L=84_N=64_GPU/H=%s/initial.pdb' % H )
    
#     file='/Users/ryota/mnt/lakers/Research/IDP/Simulation/Glass/H=%s/output%s_seq%s_lp=%s_eps=2.dcd' %(H,ens,seq,lp)
#     topo=md.load('/Users/ryota/mnt/lakers/Research/IDP/Simulation//Glass/H=%s/topology.pdb'%H).topology 
    
#     if H==1.58 or H==4.75:
#         st=100000
#     else:
#         st=1000
#     traj=md.load(file,stride=st,top=topo)
#     #traj=traj[20:]
#     time = traj.time
    
#     return traj,time

def radial(H,lp,ens,seq,pairs,traj):
    #pairs=np.array([[i,j] for j in range(L*N) for i in range(j+1,L*N)])
    #traj,time=trajectory(H,lp,ens,seq)
    rini=0.05
    rlas=5.05
    rbin=0.05
    rvec=np.linspace(rini,rlas,int((rlas-rini)/rbin))
    dist=md.compute_distances(traj,pairs)
    sumvec=np.zeros(len(rvec))
    for i in range(len(dist)):
        freqvec=np.histogram(dist[i],range=(rini,rlas),bins=int((rlas-rini)/rbin))[0]
        sumvec=sumvec+freqvec
    avesumvec=sumvec/len(dist)
    nparticle=np.sqrt(2*np.sum(avesumvec))
    ndist=avesumvec/(np.sum(avesumvec))*nparticle
    #ndensity=nparticle*avesumvec
    Vr=4*np.pi*(rlas-rini)**3/3
    rho=nparticle/Vr
    #Nr=L*N
    rdist=ndist/(4*np.pi*rvec**2*rbin*rho)
    #normadist=dist/dist[-1]
    return rvec[1:],rdist[1:]
        


def seqtype(H,seq):
    file='/Volumes/My Passport/Data/LLPS_Glass_Simulation/H=%s/seq%s.out' %(H,seq) 
    #file='/Users/ryota/mnt/lakers/Research/IDP/Simulation/Glass/H=%s/seq%s.out' %(H,seq) 
    
    fparticles = open(file, "r")
    particlesDat = fparticles.read()
    particlesDat = particlesDat.split('\n')
    particleListPerChain = [[particlesDat[i].split()[0], particlesDat[i].split()[1]] for i in range(len(particlesDat) - 1)]
    nchains=N
    particleList = particleListPerChain*nchains
    fparticles.close()
    return particleList

def pairs(H,seq,stype1,stype2):
    seqtypestr1=str(stype1)
    seqtypestr2=str(stype2)
    seqtypevec=seqtype(H,seq)
    typevec1=[]
    typevec2=[]
    
    if seqtypestr1 == seqtypestr2:
        for i in range(len(seqtypevec)):
            if seqtypevec[i][1]==seqtypestr1:
                typevec1=typevec1+[i]
        typepairs=np.array([[typevec1[i],typevec1[j]] for i in range(len(typevec1)) for j in range (i+1,len(typevec1))])
    
    else:
    
        for i in range(len(seqtypevec)):
            if seqtypevec[i][1]==seqtypestr1:
                typevec1=typevec1+[i]
            elif seqtypevec[i][1]==seqtypestr2:
                typevec2=typevec2+[i]
                
        typepairs=np.array([[typevec1[i],typevec2[j]] for i in range(len(typevec1)) for j in range (len(typevec2))])       
    
                
    return np.array(typepairs)
    

###### different H at same lp #########
lp=4
ens=1
seq=1

(traj1,time1)= trajectory(1.58,ens,seq,lp)
(traj2,time2)= trajectory(4.75,ens,seq,lp)

pairvec1=[[i,j] for i in range(L*N) for j in range (i+2,L*N)]
for i in range(N):
    pairvec1 = pairvec1 +  [[((L-1)+L*i)%(N*L),((L-1)+L*i+1)%(N*L)]]
pairvec1=np.array(pairvec1)

#pairvec1AA=np.array([[i,j] for i in range(L*N) for j in range (i,L*N)])
#rdf1=md.compute_rdf(traj1[1:],pairvec1,periodic=True,r_range=[0,5],bin_width=0.05)
rvec1,rdf1=radial(1.58,lp,ens,seq,pairvec1,traj1)

pairvec2=[[i,j] for i in range(L*N) for j in range (i+2,L*N)]
for i in range(N):
    pairvec2 = pairvec2 + [[((L-1)+L*i)%(N*L),((L-1)+L*i+1)%(N*L)]]
pairvec2=np.array(pairvec2)
#rdf2=md.compute_rdf(traj2[1:],pairvec2,periodic=True,r_range=[0,5],bin_width=0.05)
rvec2,rdf2=radial(4.75,lp,ens,seq,pairvec2,traj2)

#plt.plot(rdf1[0][1:],rdf1[1][1:]/rdf1[1][-1],label='H=%s'%(1.58),marker="o",markersize=5)
#plt.plot(rdf2[0][1:],rdf2[1][1:]/rdf2[1][-1],label='H=%s'%(4.75),marker="s",markersize=5)
plt.figure()
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
plt.plot(rvec2,rdf2,label=r'$H_3=$%s'%(4.75),marker="s",markersize=5,color="tab:orange")
plt.plot(rvec1,rdf1,label=r'$H_3=$%s'%(1.58),marker="o",markersize=5,color="tab:blue")

#plt.title(r'$l_p=$'+str(lp))
plt.xlabel(r"$r(\sigma)$")
plt.ylabel(r"$g(r)$")
plt.legend(fontsize=20,markerscale=2)
#plt.xlim([0,15])

plt.savefig('Pr_lp='+str(lp)+'.png', bbox_inches='tight',dpi = 450)
plt.show()
###### different H at same lp for specific pair#########
# #lp=0
# #ens=1
# #seq=1
# #typeseq1=0
# #typeseq2=0

# for i in ([0,0],[0,1],[0,2]):
#     typeseq1=i[0]
#     typeseq2=i[1]

#     if typeseq1==0:
#         typstr1='A'
#     elif typeseq1==1:
#         typstr1='B'
#     elif typeseq1==2:
#         typstr1='C'
        
#     if typeseq2==0:
#         typstr2='A'
#     elif typeseq2==1:
#         typstr2='B'
#     elif typeseq2==2:
#         typstr2='C'
    
#     #(traj1,time1)= trajectory(1.58,ens,seq,lp)
#     #(traj2,time2)= trajectory(4.75,ens,seq,lp)
    
#     pairvec1=pairs(1.58,seq,typeseq1,typeseq2)
#     #pairvec1AA=np.array([[i,j] for i in range(L*N) for j in range (i,L*N)])
#     #rdf1=md.compute_rdf(traj1[1:],pairvec1,periodic=True,r_range=[0,5],bin_width=0.05)
#     rvec1,rdf1=radial(1.58,lp,ens,seq,pairvec1,traj1)
    
#     pairvec2=pairs(4.75,seq,typeseq1,typeseq2)
#     #rdf2=md.compute_rdf(traj2[1:],pairvec2,periodic=True,r_range=[0,5],bin_width=0.05)
#     rvec2,rdf2=radial(4.75,lp,ens,seq,pairvec2,traj2)
    
#     #plt.plot(rdf1[0][1:],rdf1[1][1:]/rdf1[1][-1],label='H=%s'%(1.58),marker="o",markersize=5)
#     #plt.plot(rdf2[0][1:],rdf2[1][1:]/rdf2[1][-1],label='H=%s'%(4.75),marker="s",markersize=5)
#     plt.figure()
#     plt.style.use('plot.mplstyle')
#     plt.plot(rvec2,rdf2,label=r'$H_3=$%s'%(4.75),marker="s",markersize=5,color="tab:orange")
#     plt.plot(rvec1,rdf1,label=r'$H_3=$%s'%(1.58),marker="o",markersize=5,color="tab:blue")
    
#     #plt.title(r'$l_p=$'+str(lp))
#     plt.xlabel(r"$r($\sigma)$")
#     plt.ylabel(r"$g(r)$")
#     plt.title(typstr1+'-'+typstr2)
#     plt.legend()
#     #plt.xlim([0,15])
    
#     plt.savefig('Pr'+typstr1 + typstr2+'_lp='+str(lp)+'.png', bbox_inches='tight',dpi = 450)
#     plt.show()


####### different lp at same H #########
# ens=1
# seq=1
# H=4.75
# lpvec=[0,2,4]
# rdfmat=np.zeros((len(lpvec),99))
# ii=0
# for lp in lpvec:
#     (traj,time)= trajectory(H,lp,ens,seq)
    
#     pairvec=[[i,j] for i in range(L*N) for j in range (i,L*N)]
#     #pair=np.array([0,50]).reshape((1,2))
#     rdf=md.compute_rdf(traj[1:],pairvec,periodic=True,r_range=[0,5],bin_width=0.05)
    
#     maxi=max(rdf[1][1:])
   
#     rdfmat[ii]=rdf[1][1:]/maxi
#     ii=ii+1


# marker=["o","^","s"]
# for i in range(len(lpvec)):
#     plt.plot(rdf[0][1:],rdfmat[i],label=r'$l_p=$%s'%(lpvec[i]),marker=marker[i],markersize=5)


# #plt.title(r'$l_p=$'+str(lp))
# plt.xlabel(r"$r$(nm)")
# plt.ylabel(r"$\hat{g}$(r)")
# plt.legend()
# #plt.show()
# plt.savefig('Pr_H='+str(H)+'.png', bbox_inches='tight',dpi = 450)
