"""
A script for batch vina_split ligands pdbqt files from docking results.

usage:
    python batch_vina_split.py  -i output/
input:
    /path/to/docking_results (path)
output:
    batch_vina_split_results.csv

file structure:
    default output directory of autodock vina should look like this:
    project/
        receptor/
        ligand/
        output/
            receptor_1/
                ligand_1.pdbqt
            receptor_2/
                ligand_1.pdbqt
"""

import os
from glob import glob
import multiprocessing
from multiprocessing import Pool
import argparse

# set argparser
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_dir', help='path to docking results')

# arg passed from command line
args = parser.parse_args()
input_dir = args.input_dir

# get path of ligand files
ligand_path = glob(input_dir + '/*/*.pdbqt')
print()
print(f'> {len(ligand_path)} ligand files found')


def vina_split(in_file: str):
    cmd = f'vina_split --input {in_file}'
    print(cmd)
    os.system(cmd)


# write ligand results to csv
def write_split_ligands(out_path: str = 'output/batch_vina_split_results.csv'):
    with open(out_path, 'w') as f:
        # /Users/muwang/PycharmProjects/PFleet/data/docking_result/AF-Q8NH18-F1-model_v4/2-BUTANONE_regen_ligand_5.pdbqt

        f.write('receptor_path,ligand_name,ligand_path\n')
        split_ligand_path = glob(input_dir + '/*/*_ligand_?.pdbqt')

        def get_names(path):
            receptor_name = 'receptor/' + path.split('/')[-2]+'.pdb'
            ligand_name = 'ligand/' + os.path.basename(path).split('_regen_')[0]+'_regen.sdf'
            return receptor_name, ligand_name, path

        for ligand in split_ligand_path:
            receptor, ligand, result_lig = get_names(ligand)
            f.write(f'{receptor},{ligand},{result_lig}\n')
    print(f'> save results to {out_path}')


# use multiprocessing to speed up conversion
def mulit_process():
    cpu_nr = multiprocessing.cpu_count()
    if cpu_nr >= 32:
        cpu_nr = 32
    print(f'cpu number: {cpu_nr}')

    with Pool(cpu_nr) as p:
        p.map(vina_split, ligand_path)


def main():
    mulit_process()
    write_split_ligands()


if __name__ == '__main__':
    main()
