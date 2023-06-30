# PFleet
A python package to annotate, manipulate(select/extract/combine), visualize multiple groups of related structures and ligands.

# protocol
* batch autodock: 
  ```bast
  python batchdock_autodock_vina.py
  ```
* analysis :
  ```bash 
  python interaction_analysis.py -il input_ligand -ir input_receptor [-o output_dir]
  ```