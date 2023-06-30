"""
This script is used to run the interaction analysis in batch mode.

usage:
    python batch_interaction_analysis.py [-ic /path/to/dock_results]  [-o /path/to/output]

"""
import os
from glob import glob
import multiprocessing
from multiprocessing import Pool
import argparse
import pandas as pd

# set argparser
parser = argparse.ArgumentParser()
parser.add_argument('-ic', '--vina_split_csv', default='output/batch_vina_split_results.csv',
                    help='path to docking results')
# parser.add_argument('-ir', '--receptor_dir', help='path to receptor files')
parser.add_argument('-o', '--output_dir', default='output/dist_output', help='path to output directory')

# arg passed from command line
args = parser.parse_args()
vina_split_csv = args.vina_split_csv
# receptor_dir = args.receptor_dir
output_dir = args.output_dir

# get data from vina_split_csv by pandas
df = pd.read_csv(vina_split_csv)
# get receptor_path, ligand_name, ligand_path from df
receptor_path = df['receptor_path'].tolist()
ligand_path = df['ligand_path'].tolist()
data_path = zip(receptor_path, ligand_path)


# use multiprocessing to run the interaction analysis
def interaction_analysis(data_temp):
    receptor_temp = data_temp[0]
    ligand_temp = data_temp[1]
    cmd = f'python interaction_analysis.py -il {ligand_temp} -ir {receptor_temp} -o {output_dir}'
    print(cmd)
    os.system(cmd)


# use multiprocessing to run the interaction analysis
def batch_interaction_analysis(data_path_temp):
    cpu_nr = multiprocessing.cpu_count()
    if cpu_nr >= 32:
        cpu_nr = 32
    with Pool(multiprocessing.cpu_count()) as p:
        p.map(interaction_analysis, data_path_temp)


# run the batch_interaction_analysis
if __name__ == '__main__':
    batch_interaction_analysis(data_path)
