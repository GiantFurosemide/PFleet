#!/usr/bin/env python
"""
Convert pdbqt files (of ligands or receptors) to pdb files by openbabel.
This script will serach pdbqt files under /path/to/pdbqt, convert and save to the same directory.

requirement:
    openbabel 2.4.90 -> obabel

input:
    /path/to/pdbqt (path)

usage:
    batch_format_convert.py /path/to/pdbqt [-i input_format] [-o output_format]

"""

from glob import glob
import os
import argparse
from functools import partial

# set argparser
parser = argparse.ArgumentParser()
parser.add_argument('pdbqt_path', help='path to pdbqt files')
parser.add_argument('-i', '--input_format', default='pdbqt', help='input format')
parser.add_argument('-o', '--output_format', default='pdb', help='output format')

# arg passed from command line
args = parser.parse_args()
pdbqt_path = args.pdbqt_path
in_format = args.input_format
out_format = args.output_format


def convert(in_file: str, in_format: str = 'pdbqt', out_format: str = 'pdb'):
    in_format = in_format.lower()
    out_format = out_format.lower()
    out_file = in_file.replace(f'.{in_format}', f'.{out_format}')

    cmd = f'obabel -i{in_format} {in_file} -o{out_format} -O {out_file}'
    print(cmd)
    os.system(cmd)


# make partial function
convert = partial(convert, in_format=in_format, out_format=out_format)


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
        p.map(convert, glob(f'{pdbqt_path}/*.{in_format}'))


if __name__ == '__main__':
    main()
