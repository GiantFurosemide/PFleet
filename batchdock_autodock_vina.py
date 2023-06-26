"""
This script is used to batch dock ligands to a receptor using autodock vina.

requirements:
    openbabel -> obabel
    autodock vina -> vina

input:
    receptor file (pdb) (already aligned and keep only the chain of interest, no ligands)
    ligand files (sdf)
    config file (conf.txt) (set grid box size and center)

output:
    docked ligand files (pdbqt)
    log files (txt)
    out_pdbqt_path.txt (a list of path of out pdbqt files)
    out_pdbqt_log.txt (a list of path of out log files)

file structure:
    batchdock_autodock_vina.py
    receptor/
        receptor_1.pdb
        receptor_2.pdb
    ligand/
        ligand_1.sdf
        ligand_2.sdf
    output/
        receptor_1/
            ligand_1.pdbqt
            ligand_1_log.txt
            ligand_2.pdbqt
            ligand_2_log.txt
        receptor_2/
            ligand_1.pdbqt
            ligand_1_log.txt
            ligand_2.pdbqt
            ligand_2_log.txt

"""

import os
from glob import glob
# 0
# set output directory
if not os.path.exists('output'):
    os.mkdir('output')
# 1
# get receptor and ligand file paths
receptor_path = glob('receptor/*.pdb')
ligand_path = glob('ligand/*.sdf')

# extract basename
receptor_name = [os.path.basename(i).split('.')[0] for i in receptor_path]
ligand_name = [os.path.basename(i).split('.')[0] for i in ligand_path]

# set output directory for each receptor
for i in receptor_name:
    if not os.path.exists('output/' + i):
        os.mkdir('output/' + i)
        print('output/' + i + ' created')

# 2
# convert receptors and ligands to pdbqt
# convert receptors

receptor_pdbqt_path = []
for i in receptor_path:
    out_pdbqt_temp = 'output/' + os.path.basename(i).split('.')[0] + '/' + os.path.basename(i).split('.')[0] + '.pdbqt'
    cmd = 'obabel -ipdb ' + i + ' -opdbqt -O ' + out_pdbqt_temp + ' -p 7.5 -r -xr'
    print(cmd)
    os.system(cmd)

    # save converted file path to a list
    receptor_pdbqt_path.append(out_pdbqt_temp)

    print(i + ' converted to '+ out_pdbqt_temp)
    del out_pdbqt_temp

# convert ligands
ligand_pdbqt_path = []
for ligand in ligand_path:
    out_smi_temp = 'ligand/' + os.path.basename(ligand).split('.')[0].strip() + '.smi'
    cmd = 'obabel -isdf ' + '\'' + ligand + '\'' + ' -osmi -O ' + out_smi_temp + ' -p 7.5'
    print(cmd)
    os.system(cmd)

    out_sdf_temp = 'ligand/' + os.path.basename(ligand).split('.')[0].strip() + '_regen' + '.sdf'
    cmd = 'obabel -ismi ' + out_smi_temp + ' -osdf -O ' + out_sdf_temp + ' --gen3d --conformer --nconf 50 --systematic'
    print(cmd)
    os.system(cmd)

    out_pdbqt_temp = 'ligand/' + os.path.basename(ligand).split('.')[0].strip() + '_regen' + '.pdbqt'
    cmd = 'obabel -isdf ' + out_sdf_temp + ' -opdbqt -O ' + out_pdbqt_temp + ' -p 7.5 '
    print(cmd)
    os.system(cmd)

    # save converted file path to a list
    ligand_pdbqt_path.append(out_pdbqt_temp)

    print(ligand + ' converted to ' + out_pdbqt_temp)
    del out_pdbqt_temp
    del out_sdf_temp
    del out_smi_temp

# 3
# run autodock vina
print(f" please check the config file (conf.txt) before running the script")
# create a list to store path of out pdbqt files
out_pdbqt_path = []
out_log_path = []
for receptor in receptor_pdbqt_path:
    for ligand in ligand_pdbqt_path:
        out_pdbqt_temp = 'output/' + os.path.basename(receptor).split('.')[0] + '/' + os.path.basename(ligand).split('.')[0] + '.pdbqt'
        out_log_temp = 'output/' + os.path.basename(receptor).split('.')[0] + '/' + os.path.basename(ligand).split('.')[0] + '_log.txt'
        cmd = 'vina --config conf.txt --receptor ' + receptor + ' --ligand ' + ligand + ' --out ' + out_pdbqt_temp + ' --log ' + out_log_temp
        print(cmd)
        os.system(cmd)
        # save converted file path to a list
        out_pdbqt_path.append(out_pdbqt_temp)
        out_log_path.append(out_log_temp)
        print(ligand + ' docked to ' + receptor)
        del out_pdbqt_temp
        del out_log_temp
    print(receptor + ' finished')

# 4
# save output path to a out_pdbqt_path.txt file
with open('out_pdbqt_path.txt', 'w') as f:
    for i in out_pdbqt_path:
        f.write(i + '\n')

# save output path to a out_log_path.txt file
with open('out_log_path.txt', 'w') as f:
    for i in out_log_path:
        f.write(i + '\n')

print('You can search out pdbqt path by "out_pdbqt_path.txt"')
