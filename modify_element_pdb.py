#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 16:04:56 2023

@author: ryota
"""

H=1.58
lp=4
N=64

# Define a mapping from sequence numbers to atom types
seq_to_atom_type = {
    0: 'C',
    1: 'H',
    2: 'S'
}

# Parse sequence information from seq1.out into a Python list
seq_info = []
path='/Volumes/My Passport/Data/LLPS_Glass_Simulation/PNAS_revise/H=%s' %H
seq1_out_path = path+'/seq1.out'  # Replace with your file path
with open(seq1_out_path, 'r') as f:
    for line in f.readlines():
        _, bead_type = line.strip().split()
        seq_info.append(int(bead_type))

# Read the initial_lp=0.pdb file and modify both the 3rd column and the atom types based on the sequence information
#pdb_path = path+f'/initial_N={N}.pdb' 
pdb_path = path+f'/very_initial_N={N}.pdb'  
fully_corrected_pdb_lines = []
with open(pdb_path, 'r') as f:
    for i, line in enumerate(f.readlines()):
        if line.startswith('ATOM') or line.startswith('HETATM'):  # Check for both 'ATOM' and 'HETATM'
            # Determine which bead in the polymer this atom corresponds to
            bead_index = i % 84  # 84 beads in a polymer
            # Determine the new atom type based on the sequence information
            new_atom_type = seq_to_atom_type[seq_info[bead_index]]
            # Modify the line to update the 3rd column (columns 13-16) and the atom type (columns 77-78)
            fully_corrected_line = line[:12] + new_atom_type.ljust(4) + line[16:76] + new_atom_type.ljust(2) + line[78:]
            fully_corrected_pdb_lines.append(fully_corrected_line)
        else:
            fully_corrected_pdb_lines.append(line)

# Write the fully corrected PDB lines back to a new file
#fully_corrected_pdb_path = path+'SingleChain_H=%s_lp=%s.pdb' %(H,lp)  # Replace with your desired output path
fully_corrected_pdb_path =path+f'/modified_very_initial_N={N}.pdb'  # Replace with your desired output path
with open(fully_corrected_pdb_path, 'w') as f:
    f.writelines(fully_corrected_pdb_lines)
