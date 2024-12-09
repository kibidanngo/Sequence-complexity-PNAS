#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 21:15:29 2024

@author: ryota
"""

import pandas as pd
import matplotlib.pyplot as plt

lp=0

def read_data_from_file(filename):
    df = pd.read_csv(filename)
    return df['Time'].tolist(), df['Nemavec'].tolist(), df['Rgvec'].tolist()

path='/Volumes/My Passport/Data/LLPS_Glass_Simulation/'
time1, nemavec1, Rgvec1 = read_data_from_file(path+f'H=1.58/phase_H=1.58_lp={lp}_ens=1_seq=1.csv')
time2, nemavec2, Rgvec2 = read_data_from_file(path+f'H=4.75/phase_H=4.75_lp={lp}_ens=1_seq=1.csv')

plt.figure(dpi=400,figsize=(6,5))
plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
plt.scatter(Rgvec2, nemavec2, c=time2, cmap='Reds', label=r'$H_3=$' + str(4.75),s=40)
plt.scatter(Rgvec1, nemavec1, c=time1, cmap='Blues', label=r'$H_3=$' + str(1.58),s=40)
plt.xlabel(r"$R_G/L$")
plt.ylabel("S")
plt.xlim([0.09,0.115]) # lp=0
#plt.xlim([0.1,0.113]) # lp=2
#plt.xlim([0.22,0.5])
#plt.colorbar()
#plt.title(f'$\epsilon_b=${lp}')
#plt.legend()
plt.show()