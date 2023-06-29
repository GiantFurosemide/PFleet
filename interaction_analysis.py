import matplotlib.pyplot as plt
from math import sqrt,ceil
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import distances

# get ligand atom
pdbqt_file = "/Users/muwang/PycharmProjects/PFleet/data/docking_result/AF-Q8NH18-F1-model_v4/2-BUTANONE_regen_ligand_1.pdbqt"
u = mda.Universe(pdbqt_file)
ligand_atoms = u.atoms
ligand_atoms_coords = u.atoms.positions

# get protein residues & atoms
pdbqt_file = "data/receptor/AF-Q8NH18-F1-model_v4.pdb"
u = mda.Universe(pdbqt_file)
receptor_chains = u.segments # should be only one Chain in this Universe
receptor_residues = u.residues
receptor_atoms = u.atoms
receptor_atoms_coords = u.atoms.positions

# get positions of center of mass of each residue
com_positions = []
for residue in receptor_residues:
    com_positions.append(residue.atoms.center_of_mass())
com_positions = np.array(com_positions)

# distance anaylse
dist_arr = distances.distance_array(ligand_atoms.positions, # reference
                                    #receptor_atoms.positions, # configuration
                                    com_positions,
                                    box= u.dimensions)
print(f'dist_arr.shape: {dist_arr.shape}')

# Assuming you have the array of distances named 'distances' with shape (5, 618)


ligand_atoms_nr = dist_arr.shape[0]
receptor_atoms_nr = dist_arr.shape[1]

# Define the number of subplots you want to create
num_subplots = int(sqrt(dist_arr.shape[0]*dist_arr.shape[1]) // ligand_atoms_nr) + 1
num_col_subplot = ceil(sqrt(dist_arr.shape[0]*dist_arr.shape[1]))
subplots_per_row = 1

target_shape = [ligand_atoms_nr,int(num_col_subplot*num_subplots)]
assert target_shape[1] * num_subplots >= receptor_atoms_nr
new_array = np.ones(target_shape,dtype = dist_arr.dtype)*9999
new_array[:,:dist_arr.shape[1]] = dist_arr
dist_arr = new_array
dist_arr = dist_arr.reshape([num_subplots,target_shape[0],num_col_subplot])
print(f'new dist_arr.shape: {dist_arr.shape}')

## Calculate the number of rows and columns for the subplots
#num_rows = (num_subplots + subplots_per_row - 1) // subplots_per_row
#num_cols = min(num_subplots, subplots_per_row)

## Create the subplots
fig, axes = plt.subplots(num_subplots,subplots_per_row,figsize=(6,10))

# Flatten the axes array if necessary
if num_subplots == 1:
    axes = np.array([axes])

#plt.imshow(dist_arr[0])


# Iterate over each subplot and plot the corresponding data
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=2.9, hspace=0.9)
for i, ax in enumerate(axes.flat):
    if i < num_subplots:
        ligand_index = i % ligand_atoms_nr
        receptor_index = i % receptor_atoms_nr

        # add residue ID labels to axes
        tick_interval_ligand = 1
        tick_interval_receptor = 5

        block_nr,ligand_atoms_nr, receptor_atoms_nr = dist_arr.shape

        ax.set_yticks(np.arange(ligand_atoms_nr)[::tick_interval_ligand])
        xtickvalue = np.arange(i*receptor_atoms_nr,(i+1)*receptor_atoms_nr)[::tick_interval_receptor]
        ax.set_xticks(xtickvalue)
        ax.set_yticklabels(ligand_atoms.ids[::tick_interval_ligand])
        print(np.arange(1,receptor_atoms_nr+1)[::tick_interval_receptor])
        print(ligand_atoms.ids[::tick_interval_ligand])
        #print(receptor_residues.residues.resnums[::tick_interval_receptor])
        #print(receptor_residues.residues.resnums)
        xtickvalue_new = np.zeros(block_nr*ligand_atoms_nr*receptor_atoms_nr,dtype=int)
        xtickvalue_new[:len(receptor_residues.residues.resnums)] += receptor_residues.residues.resnums
        xticklabels = xtickvalue_new[i * receptor_atoms_nr:(i + 1) * receptor_atoms_nr:tick_interval_receptor]
        assert len(xtickvalue) == len(xticklabels)
        ax.set_xticklabels(xticklabels)
        print(">>>")
        print(f'xtickvalue:{xtickvalue}')
        print(f"xticklabels:{xticklabels}")


        ax.set_xlabel('Receptor Residues')
        ax.set_ylabel('Ligand')

        im = ax.imshow(dist_arr[ligand_index], label=f'Ligand Atom {ligand_index + 1}')
        fig.colorbar(im, ax=ax)

        #ax.set_title(f'Distance between ligand-(receptor:center of mass of residue){i}')
        #ax.legend()

## Remove any unused subplots
#if num_subplots < num_rows * num_cols:
#    for i in range(num_subplots, num_rows * num_cols):
#        fig.delaxes(axes.flatten()[i])

# Adjust the layout and spacing
#fig.tight_layout()

# Show the plot
#cbar = fig.colorbar(fig,ax)
#cbar.ax.set_ylabel('Distance (Angstrom)')
plt.show()

# save dist_arr
a1,a2,a3 = dist_arr.shape
np.savetxt('output/a.csv',dist_arr.reshape((a2,a1*a3)).T,delimiter=',')
#data = np.loadtxt('output/a.csv', delimiter=',')
