#!/usr/bin/env python
"""
Convert pdbqt files (of ligands or receptors) to pdb files by openbabel.
This script will serach pdbqt files under /path/to/pdbqt, convert and save to the same directory.

requirement:
    openbabel 2.4.90 -> obabel

input:
    /path/to/pdbqt (path)

usage:
    pdbqt2pdb.py /path/to/pdbqt

"""

from glob import glob
import os
import argparse

# set argparser
parser = argparse.ArgumentParser()
parser.add_argument('pdbqt_path', help='path to pdbqt files')

# arg passed from command line
args = parser.parse_args()
pdbqt_path = args.pdbqt_path


def convert(pdbqt_file):
    pdb_file = pdbqt_file.replace('.pdbqt', '.pdb')
    cmd = f'obabel -ipdbqt {pdbqt_file} -opdb -O {pdb_file}'
    print(cmd)
    os.system(cmd)


def main():
    if not os.path.exists(pdbqt_path):
        raise FileNotFoundError(f'{pdbqt_path} not found')

    # use multiprocessing to speed up conversion
    import multiprocessing
    from multiprocessing import Pool
    # check number of cpu
    cpu_nr = multiprocessing.cpu_count()
    print(f'cpu number: {cpu_nr}')
    if cpu_nr >= 8:
        cpu_nr = 8

    with Pool(cpu_nr) as p:
        p.map(convert, glob(f'{pdbqt_path}/*.pdbqt'))


if __name__ == '__main__':
    main()
