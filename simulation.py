#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from simtk import openmm
from simtk import unit
from simtk.openmm import app
from simtk.openmm import vec3
import numpy as np
from numpy import diag
# import matplotlib.pyplot as plt
# from mdtraj.reporters import DCDReporter
# from mdtraj.reporters import NetCDFReporter
from simtk.openmm.app.dcdreporter import DCDReporter
import mdtraj as md
import os
import time
import sys 


        
start = time.time()
np.random.seed()


# Input parameters in openmm  units#

kBoltzmann = 1.380648 * 10 ** (-23) * 6.023 * 10 ** (23) / 1000 * unit.kilojoule_per_mole / unit.kelvin
viscosity = 8.90 * 10 ** -4 * 1000 * 6.023 * 10 ** 23 / (10 ** 9 * 10 ** 12) * \
            unit.amu / unit.nanometer / unit.picoseconds

kBT = 1 * unit.kilojoule_per_mole
input_file="initial_lp=0.pdb"
input_coords = md.load(input_file)
nparticles = input_coords.n_atoms
# Basic units
sigmaBase = 1 * unit.nanometer  # 1 Angstrom
temperature = kBT/unit.MOLAR_GAS_CONSTANT_R

# Set box size
edge = 1 * sigmaBase 
#Creat topology
topology=input_coords.topology.to_openmm()

#topology=input_coords.xyz.input_coords.xyz
xperiod=100 * sigmaBase
yperiod=100 * sigmaBase
zperiod=100 * sigmaBase
topology.setPeriodicBoxVectors([[xperiod.value_in_unit(unit.nanometers),0,0], [0,yperiod.value_in_unit(unit.nanometers),0], [0,0,zperiod.value_in_unit(unit.nanometers)]])


# No need to write angstrom as mdtraj import already converted the units
positions = unit.Quantity(np.array(input_coords.xyz.tolist()[0][0:nparticles]), unit.nanometer)


lp = int(sys.argv[1])
seqnum = int(sys.argv[2])
ensnum = int(sys.argv[3])


print("Loaded everything")

unitLen = sigmaBase
unitEn = 1 * kBT
unitTimeHighFric = 1 * unit.picosecond
diffCoeff =  1 * unitLen**2/unitTimeHighFric
unitTimeLowFric = 1 * unit.picosecond

tau_L = unitTimeLowFric.value_in_unit(unit.picosecond)
low_fric = 0.01 / tau_L
timeStepLD = 0.01 * tau_L



# Create an empty system object
system = openmm.System()
# set periodic Box in system
system.setDefaultPeriodicBoxVectors([xperiod, 0, 0], [0, yperiod, 0], [0, 0, zperiod])



filename="seq"+str(seqnum)+".out"


nchains=64
fparticles = open(filename, "r")
particlesDat = fparticles.read()
particlesDat = particlesDat.split('\n')
particleListPerChain = [[particlesDat[i].split()[0], particlesDat[i].split()[1]] for i in range(len(particlesDat) - 1)]
particleList = particleListPerChain*nchains
fparticles.close()

   
## =============================================================================
#k_H = 10 * unit.kilojoules_per_mole/unit.nanometer**2
#HarmonicBonds=openmm.HarmonicBondForce()
#
#
#lseq=len(particleListPerChain)
#for n in range(nchains):
#    for i in range(lseq-1):
#        HarmonicBonds.addBond(n*lseq+i, n*lseq+(i+1) ,0,k_H)
#        #print(n*lseq+i, n*lseq+(i+1))
#        
#
#system.addForce(HarmonicBonds)
## =============================================================================




# =============================================================================
k_bend = lp * kBT/sigmaBase**2 # 150 is the literature value
#bendAngle = openmm.CustomCompoundBondForce(3, "kbend*(1+(cos(angle(p1,p2,p3))))")
bendAngle = openmm.CustomAngleForce("kbend*(1+(cos(theta)))")
bendAngle.addGlobalParameter('kbend',k_bend)
bendAngle.setUsesPeriodicBoundaryConditions(True)
print(bendAngle.usesPeriodicBoundaryConditions())

lseq=len(particleListPerChain)
#lrig = 40 # start of rigid part in sequence
#nrig = 8 # number of rigid part in sequence
lrig=0
nrig=lseq
for n in range(nchains):
    for i in range(lseq-2):
    #for i in range(int(lseq/2)-2):
        if lrig <= i and i <= lrig+nrig:
            bendAngle.addAngle(n*lseq+i,n*lseq+(i+1),n*lseq+(i+2))
        #print(n*lseq+i, n*lseq+(i+1),n*lseq+(i+2))

    
system.addForce(bendAngle)    
# =============================================================================



# =============================================================================
#delta = [
#     #  A     B     C     (internal unit - nm)
#        0,    0,    0,   #  A
#        0,    0,   0,   #  B
#        0,    0,   0,   #  C
#                        ]  



delta = [
     #  A     B     C     (internal unit - nm)
        1,    0,    0,   #  A
        0,    1,   0,   #  B
        0,    0,   1,   #  C                        ]  
			]

eps=2
epsilon = eps * unit.kilojoules_per_mole
sigma= 1* sigmaBase
cutoff = 4 * sigma


LJNonBound = openmm.CustomNonbondedForce("4*ep*(sig/r)^12-del*4*ep*(sig/r)^6 ; del=delta(type1,type2)")



LJNonBound.addGlobalParameter('ep', epsilon)
LJNonBound.addGlobalParameter('sig', sigma)
LJNonBound.addTabulatedFunction('delta', openmm.Discrete2DFunction(3, 3, delta))
LJNonBound.addPerParticleParameter('type')
LJNonBound.setCutoffDistance(cutoff)
LJNonBound.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
print(LJNonBound.usesPeriodicBoundaryConditions())

for i in range(nparticles):
         if particleList[i][1] == '0':
             LJNonBound.addParticle([0])
         elif particleList[i][1] == '1':
             LJNonBound.addParticle([1])
         elif particleList[i][1] == '2':
             LJNonBound.addParticle([2])
         
for n in range(nchains):
    for i in range(lseq-1):
        LJNonBound.addExclusion(n*lseq+i,n*lseq+(i+1))
                


system.addForce(LJNonBound)    
# =============================================================================

# =============================================================================
R0_fene = 1.5 * sigma
#R0_fene_other = 3 * sigmaBase
k_fene = 30 * epsilon/sigma**2
diameter = 1 * sigmaBase

feneBonds = openmm.CustomBondForce("-k/2.0*R0^2*log(1-((r-r0)/R0)^2)")
feneBonds.addGlobalParameter('k', k_fene)
feneBonds.addGlobalParameter('R0', R0_fene)
feneBonds.addGlobalParameter('r0',diameter)
feneBonds.setUsesPeriodicBoundaryConditions(True)
print(feneBonds.usesPeriodicBoundaryConditions())

lseq=len(particleListPerChain)
for n in range(nchains):
    for i in range(lseq-1):
        feneBonds.addBond(n*lseq+i, n*lseq+(i+1))
        
 
system.addForce(feneBonds)    
# =============================================================================
    
cm = openmm.CMMotionRemover()
#system.addForce(cm)
print("Setup forces")

Marray = [1] * nparticles

for index in range(nparticles):
    #nonBondedTotalForce.addParticle([index])
    #print (index)
    system.addParticle(Marray[index])

   
# =============================================================================
from simtk.openmm import LangevinIntegrator, VerletIntegrator, BrownianIntegrator
from simtk.openmm.app import Simulation, PDBReporter, StateDataReporter, PDBFile
from sys import stdout
from openmmtools.integrators import BAOABIntegrator
 
platform = openmm.Platform.getPlatformByName('CUDA')
 
integratorBAOAB = BAOABIntegrator(temperature, low_fric, timeStepLD) 
  
simulation = Simulation(topology, system, integratorBAOAB, platform)



simulation.context.setPositions(positions)

boxvector = [[xperiod/unit.nanometer,0,0], [0,yperiod/unit.nanometer,0], [0,0,zperiod/unit.nanometers]]*unit.nanometer


simulation.context.setPeriodicBoxVectors(*boxvector)



simulation.minimizeEnergy(1*unit.kilocalorie_per_mole, 10000)

print (system.usesPeriodicBoundaryConditions())

nsteps = 1000000000
freq = 100000

simulation.reporters.append(DCDReporter('output'+str(ensnum)+'_'+filename[0:-4]+'_lp='+str(lp)+'_'+'eps='+str(eps)+'.dcd', freq,append=False, enforcePeriodicBox=True))

simulation.reporters.append(StateDataReporter('PotEnergy'+str(ensnum)+'_'+filename[0:-4]+'_lp='+str(lp)+'_'+'eps='+str(eps)+'.dat',freq,potentialEnergy=True,\
                                              totalEnergy=True,step=True,time=True,separator='   '))
with open('topology.pdb', 'w') as f:
    PDBFile.writeFile(topology, positions, f)
 

simulation.step(nsteps)
     




 
end = time.time()
print("Time Elapsed : ", end - start)
# =============================================================================
