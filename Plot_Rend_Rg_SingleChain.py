#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 11:43:32 2023

@author: ryota
"""
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit




def t(L, lp):
    return L / ((2 * lp) / 3)

def afa(L, lp):
    return 3 * t(L, lp) / 4

def Nc(L, lp):
    afa_val = afa(L, lp)
    return (4 * (afa_val)**(3/2) * np.exp(afa_val)) / (
        np.pi**(3/2) * (4 + 12 * (afa_val)**-1 + 15 * (afa_val)**-2))

def PR(R, L, lp):
    return Nc(L, lp) * (4 * np.pi * (R / L)**2) / (1 - (R / L)**2)**(9/2) * \
           np.exp((-3 * t(L, lp)) / (4 * (1 - (R / L)**2))) / L



def compute_distributions(seq, lp, paths, ens_values):
    # Set the aesthetic style of the plots
    #sns.set_style("whitegrid")
    plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
    plt.figure(figsize=(20, 5),dpi=400)

    # Define colors for each path
    colors = {
        'H=1.58': 'blue',
        'H=4.75': 'red'
    }

    for path in paths:
        end_to_end_distances_all = []
        radii_of_gyration_all = []

        for ens in ens_values:
            # Load the trajectory
            topology_path = f"{path}/topology_singlechain.pdb"
            trajectory_path = f"{path}/outputSingleChain{ens}_seq{seq}_lp={lp}_eps=2.dcd"
            trajectory = md.load(trajectory_path, top=topology_path)

            # Compute end-to-end distance for each frame
            end_to_end_distances = np.linalg.norm(trajectory.xyz[:, 0, :] - trajectory.xyz[:, -1, :], axis=1)
            end_to_end_distances_all.extend(end_to_end_distances)

            # Compute radius of gyration for each frame
            radii_of_gyration = md.compute_rg(trajectory)
            radii_of_gyration_all.extend(radii_of_gyration)

        # Probability distribution of end-to-end distance
        
        plt.subplot(1, 2, 1)
        plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
        sns.histplot(end_to_end_distances_all, stat="density",bins=50, kde=False, color=colors[path[-6:]], label=f'$H_3=${path[-4:]}')
        plt.title(r'$\epsilon_b$=%s'%lp)
        plt.xlabel(r'$R_{ee}\ (\sigma)$ ')
        plt.ylabel(r'$P(R_{ee})$')
        plt.legend(fontsize=25)
        plt.grid(False)
        
        
            
        # Fitting PR to the density distribution
        L = 84
        hist, bin_edges = np.histogram(end_to_end_distances_all, bins=50, density=True)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        popt, _ = curve_fit(lambda R, lp: PR(R, L, lp), bin_centers, hist, p0=lp+1)
        lp_fit = popt[0]
    
        # Plot fitted PR function
        R_values = np.linspace(min(end_to_end_distances_all), max(end_to_end_distances_all), 100)
        PR_values = PR(R_values, L, lp_fit)
        plt.plot(R_values, PR_values, color=colors[path[-6:]], label=r'$l_p=$%s'%np.round(lp_fit,2))
    
        plt.legend(fontsize=25)
        #plt.show()
    
            

        # Probability distribution of radius of gyration
        
        plt.subplot(1, 2, 2)
        sns.histplot(radii_of_gyration_all, stat="density",bins=50, kde=True, color=colors[path[-6:]], label=f'$H_3=${path[-4:]}')
        plt.title(r'$\epsilon_b$=%s'%lp)
        plt.xlabel(r'$R_g\ (\sigma)$')
        plt.ylabel(r'$P(R_g)$')
        plt.legend(fontsize=25)
        plt.grid(False)

    plt.subplot(1, 2, 1)
    #lt.legend()
    plt.subplot(1, 2, 2)
    #plt.legend(fontsize=25)
    plt.tight_layout()
    plt.show()


def compute_means(seq, lp, paths, ens_values):
    end_to_end_means = []
    rg_means = []
    
    for path in paths:
        end_to_end_distances_all = []
        radii_of_gyration_all = []

        for ens in ens_values:
            # Load the trajectory
            topology_path = f"{path}/topology_singlechain.pdb"
            trajectory_path = f"{path}/outputSingleChain{ens}_seq{seq}_lp={lp}_eps=2.dcd"
            trajectory = md.load(trajectory_path, top=topology_path)

            # Compute end-to-end distance for each frame
            end_to_end_distances = np.linalg.norm(trajectory.xyz[:, 0, :] - trajectory.xyz[:, -1, :], axis=1)
            end_to_end_distances_all.extend(end_to_end_distances)

            # Compute radius of gyration for each frame
            radii_of_gyration = md.compute_rg(trajectory)
            radii_of_gyration_all.extend(radii_of_gyration)

        end_to_end_means.append(np.mean(end_to_end_distances_all))
        rg_means.append(np.mean(radii_of_gyration_all))
    
    return end_to_end_means, rg_means


def plot_mean_values(seq_values, lp_values, paths, ens_values):
    # Set the plot style
    plt.style.use('/Users/ryota/Documents/Glass/Analysis/plot.mplstyle')
    
    # Define colors for each path
    colors = {
        'H=1.58': 'blue',
        'H=4.75': 'red'
    }
    
    markers={'H=1.58': 'o',
        'H=4.75': 's'}

    for seq in seq_values:
        end_to_end_means_all = []
        rg_means_all = []
        for lp in lp_values:
            end_to_end_means, rg_means = compute_means(seq, lp, paths, ens_values)
            end_to_end_means_all.append(end_to_end_means)
            rg_means_all.append(rg_means)
        
        end_to_end_means_all = np.array(end_to_end_means_all).T
        rg_means_all = np.array(rg_means_all).T

        # Plotting the mean values
        plt.figure(figsize=(15, 5), dpi=400)
        
        # End-to-End Distance
        plt.subplot(1, 2, 1)
        for i, path in enumerate(paths):
            plt.plot(lp_values, end_to_end_means_all[i], marker=markers[path[-6:]], label=f'$H_3=${path[-4:]}', color=colors[path[-6:]])
        #plt.title(f'Mean End-to-End Distance\nseq={seq}')
        plt.xlabel(r'$\epsilon_b\ (k_BT)$')
        plt.ylabel(r'$\langle R_{ee}\rangle \ (\sigma)$')
        plt.legend(fontsize=25)
        plt.grid(False)
        
        # Radius of Gyration
        plt.subplot(1, 2, 2)
        for i, path in enumerate(paths):
            plt.plot(lp_values, rg_means_all[i], marker=markers[path[-6:]], label=f'$H_3=${path[-4:]}', color=colors[path[-6:]])
        #plt.title(f'Mean Radius of Gyration\nseq={seq}')
        plt.xlabel(r'$\epsilon_b\ (k_BT)$')
        plt.ylabel(r'$\langle R_{g} \rangle \ (\sigma)$')
        plt.legend(fontsize=25)
        plt.grid(False)

        plt.tight_layout()
        plt.show()







# Base paths to the data files
paths = ['/Users/ryota/Documents/Glass/H=1.58', '/Users/ryota/Documents/Glass/H=4.75']

# Define the parameter ranges
ens_values = [1, 2]
seq_values = [1]
lp_values = [0, 1,2,3, 4,5,]

# Loop over all combinations of parameters (except ens) and compute/plot the distributions
for seq in seq_values:
    for lp in lp_values:
        compute_distributions(seq, lp, paths, ens_values)

# Plot mean values as functions of lp
plot_mean_values(seq_values, lp_values, paths, ens_values)