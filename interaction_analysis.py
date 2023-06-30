#!/usr/bin/env python
"""
analysis interaction between ligand and receptor

usage:
    interaction_analysis.py  -il input_ligand -ir input_receptor [-o output_dir]
output:
    output_dir/
        ligand_receptor_residues_dist_arr_residue.txt (an array of distance between ligand and receptor residues)
        ligand_receptor_sidechain_dist_arr_residue.txt (an array of distance between ligand and receptor sidechain)
load data:
    data = np.loadtxt('output/a.csv', delimiter=',')
"""

import argparse
import matplotlib.pyplot as plt
from math import sqrt, ceil
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import os
import json

# 0. varibales and utilities
# set argparser
parser = argparse.ArgumentParser()
parser.add_argument('-il', '--input_ligand', help='input ligand (single) pdbqt file ')
parser.add_argument('-ir', '--input_receptor', help='input receptor pdbfile')
parser.add_argument('-o', '--output_dir', default='output/dist_output', help='output directory')

# arg passed from command line
args = parser.parse_args()
ligand_file = args.input_ligand
receptor_file = args.input_receptor
output_dir = args.output_dir

ligand_name = os.path.basename(ligand_file).split('.')[0]
receptor_name = os.path.basename(receptor_file).split('.')[0]

# varibales

# ligand_file = "/Users/muwang/PycharmProjects/PFleet/data/docking_result/AF-Q8NH18-F1-model_v4/2-BUTANONE_regen_ligand_1.pdbqt"
# receptor_file = "data/receptor/AF-Q8NH18-F1-model_v4.pdb"
# output_dir = './dist_output'


# check if output_dir exsits and mkdir
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
assert os.path.exists(ligand_file), f'{ligand_file} not found'
assert os.path.exists(receptor_file), f'{receptor_file} not found'
assert ligand_file.endswith('pdbqt'), f'{ligand_file} is not pdbqt file'
assert receptor_file.endswith('pdb') or receptor_file.endswith('cif'), f'{receptor_file} is not pdb&cif file'


# utilities
def report_receptor(u: mda.Universe):
    protein = u.select_atoms('protein')
    report = f"""
    receptor file: {u.filename}
    number of chain: {len(protein.segments)}
    number of residue: {len(protein.residues)}
    number of atoms: {len(protein.atoms)}
    """
    print(report)
    return len(protein.segments), len(protein.residues), len(protein.atoms)


def report_ligand(u: mda.Universe):
    report = f"""
    receptor file: {u.filename}
    number of atoms: {len(u.atoms)}
    """
    print(report)
    return len(u.atoms)


def get_interact_residue_idx(data, threshold):
    """
    This function takes a numpy array called 'data' with a shape of (A, B) and a list called 'threshold' containing N float elements. It iterates over each element 'i' in the 'threshold' list, and for each row in the 'data' array, it checks if there is an element 'j' that is smaller than 'i'. If such an element exists, it records the index of that row.
    :param data: numpy array
    :param threshold: list
    :return: row_indices: dict[threshold[i]: residue idx of mda university, ...]
    """
    print('residue idx start from 1!')
    assert len(threshold) >= 1, "threshold list must contain at least one element"
    print(f'threshold is: {threshold}')

    row_in_threshold = {}
    for i, element in enumerate(threshold):
        row_in_threshold[i] = [row_num for row_num, row in enumerate(data) if
                               any(element <= threshold[i] for j, element in enumerate(row))]
    return row_in_threshold


def write_residue_idx(data_dict: dict, output_path):
    """
    This function takes a dict called 'row_in_threshold' and write it to a txt file
    :param row_in_threshold: dict
    :param output_dir: str
    :return: None
    """

    # 写入Json文件
    assert output_path.endswith('.json'), 'output path must be csv file'
    with open(output_path, 'w') as f:
        json.dump(data_dict, f)


# 1.get ligand atom
# ligand pdbqt contains one ligand
# use vina_split to split output into single ligand (9 ligands in total)
pdbqt_file = ligand_file
u = mda.Universe(pdbqt_file)
ligand_atoms_nr = report_ligand(u)
ligand_atoms = u.atoms
ligand_atoms_coords = u.atoms.positions

# 2.get protein residues & atoms
# receptor use pdb file, because pdbqt has duplicated info and cause bug when parsing by MDAnalysis
pdb_file = receptor_file
u = mda.Universe(pdb_file)
receptor_chain_nr, receptor_resi_nr, receptor_atom_nr = report_receptor(u)
select_cmd = "protein"
protein = u.select_atoms(select_cmd)
receptor_chains = protein.segments  # should be only one Chain in this Universe
receptor_residues = protein.residues
receptor_atoms = protein.atoms
receptor_atoms_coords = protein.atoms.positions

# 3.calcualte distance between ligand atoms vs receptor center_of_mass
# get positions of center of mass of each residue
residue_center_of_mass = []
for residue in receptor_residues:
    residue_center_of_mass.append(residue.atoms.center_of_mass())
residue_center_of_mass = np.array(residue_center_of_mass)

# get postions of center of mass of the side chain atom of each residue

sidechain_center_of_mass = []
for res in receptor_residues:
    sidechain = res.atoms.select_atoms('not name BB and not name H*')
    if len(sidechain) > 0:
        com = sidechain.center_of_mass()
        sidechain_center_of_mass.append(com)
    else:
        sidechain_center_of_mass.append(np.nan)
sidechain_center_of_mass = np.array(sidechain_center_of_mass)

# distance anaylse
dist_arr_residue = distances.distance_array(ligand_atoms.positions,  # reference
                                            # receptor_atoms.positions, # configuration
                                            residue_center_of_mass,
                                            box=u.dimensions)
print(f'dist_arr_residue.shape.shape(ligand atoms vs resi center of mass):\n {dist_arr_residue.shape}')
dist_arr_sidechain = distances.distance_array(ligand_atoms.positions,  # reference
                                              # receptor_atoms.positions, # configuration
                                              sidechain_center_of_mass,
                                              box=u.dimensions)
print(f'dist_arr_sidechain.shape(ligand atoms vs resi center of mass):\n {dist_arr_sidechain.shape}')

# extract residue name with distance under threshold
resids_residue = get_interact_residue_idx(dist_arr_residue, threshold=[7.0, 6.5, 5.5, 4.5, 4.0])
resids_sidechain = get_interact_residue_idx(dist_arr_sidechain, threshold=[7.0, 6.5, 5.5, 4.5, 4.0])

## prepare for plot
# ligand_atoms_nr = dist_arr_residue.shape[0]
# receptor_atoms_nr = dist_arr_residue.shape[1]
#
## Define the number of subplots you want to create
# num_subplots = int(sqrt(dist_arr_residue.shape[0] * dist_arr_residue.shape[1]) // ligand_atoms_nr) + 1
# num_col_subplot = ceil(sqrt(dist_arr_residue.shape[0] * dist_arr_residue.shape[1]))
# subplots_per_row = 1
#
# target_shape = [ligand_atoms_nr, int(num_col_subplot * num_subplots)]
# assert target_shape[1] * num_subplots >= receptor_atoms_nr
# new_array = np.ones(target_shape, dtype=dist_arr_residue.dtype) * 9999
# new_array[:, :dist_arr_residue.shape[1]] = dist_arr_residue
# dist_arr_residue = new_array
# dist_arr_residue = dist_arr_residue.reshape([num_subplots, target_shape[0], num_col_subplot])
# print(f'new dist_arr.shape: {dist_arr_residue.shape}')

## Create the subplots
# fig, axes = plt.subplots(num_subplots, subplots_per_row, figsize=(6, 10))
#
## Flatten the axes array if necessary
# if num_subplots == 1:
#    axes = np.array([axes])
#
#
## Iterate over each subplot and plot the corresponding data
# def visualization():
#    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=2.9, hspace=0.9)
#    for i, ax in enumerate(axes.flat):
#        if i < num_subplots:
#            # add residue ID labels to axes
#            tick_interval_ligand = 1
#            tick_interval_receptor = 5
#
#            block_nr, ligand_atoms_nr, receptor_atoms_nr = dist_arr_residue.shape
#
#            ax.set_yticks(np.arange(ligand_atoms_nr)[::tick_interval_ligand])
#            xtickvalue = np.arange(i * receptor_atoms_nr, (i + 1) * receptor_atoms_nr)[::tick_interval_receptor]
#            ax.set_xticks(xtickvalue)
#            ax.set_yticklabels(ligand_atoms.ids[::tick_interval_ligand])
#            print(np.arange(1, receptor_atoms_nr + 1)[::tick_interval_receptor])
#            print(ligand_atoms.ids[::tick_interval_ligand])
#            # print(receptor_residues.residues.resnums[::tick_interval_receptor])
#            # print(receptor_residues.residues.resnums)
#            xtickvalue_new = np.zeros(block_nr * ligand_atoms_nr * receptor_atoms_nr, dtype=int)
#            xtickvalue_new[:len(receptor_residues.residues.resnums)] += receptor_residues.residues.resnums
#            xticklabels = xtickvalue_new[i * receptor_atoms_nr:(i + 1) * receptor_atoms_nr:tick_interval_receptor]
#            assert len(xtickvalue) == len(xticklabels)
#            ax.set_xticklabels(xticklabels)
#            print(">>>")
#            print(f'xtickvalue:{xtickvalue}')
#            print(f"xticklabels:{xticklabels}")
#
#            ax.set_xlabel('Receptor Residues')
#            ax.set_ylabel('Ligand')
#
#            im = ax.imshow(dist_arr_residue[i], label=f'Ligand Atom {i + 1}')
#            fig.colorbar(im, ax=ax)
#
#    plt.show()
#
#
## visualization() # visualization is not good enough and meaningful enough

# save distance array to csv file
# a1, a2, a3 = dist_arr_residue.shape
# np.savetxt('output/a.csv', dist_arr_residue.reshape((a2, a1 * a3)).T, delimiter=',')

residue_out = f'{output_dir}/{ligand_name}_{receptor_name}_dist_arr_residue.csv'
sidechain_out = f'{output_dir}/{ligand_name}_{receptor_name}_dist_arr_sidechain.csv'
np.savetxt(residue_out, dist_arr_residue.T, delimiter=',')
print(f'save residue distance array to {residue_out}')
np.savetxt(sidechain_out, dist_arr_sidechain.T, delimiter=',')
print(f'save sidechain distance array to {sidechain_out}')

# write residue idx to json file
residue_out = f'{output_dir}/{ligand_name}_{receptor_name}_resids_residue.json'
sidechain_out = f'{output_dir}/{ligand_name}_{receptor_name}_resids_sidechain.json'
write_residue_idx(resids_residue, output_path=residue_out)
write_residue_idx(resids_sidechain, output_path=sidechain_out)
print(f'save residue idx to {residue_out}')
print(f'save sidechain idx to {sidechain_out}')

# read distance
# data = np.loadtxt('output/a.csv', delimiter=',')
