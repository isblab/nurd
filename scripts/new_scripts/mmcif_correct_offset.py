import os,sys
from Bio.PDB import *

in_file = sys.argv[1]
out_file = 'offset_corrected_'+in_file.split('/')[-1]

chain_in = input("Write a list of chain names for offset correction (Please give a comma separated list):\n")
chains = chain_in.split(',')

offset_in = input("Write a list of offsets for each chain provided (Please give a comma separated list):\n")                                        # Positive offsets increment the residue numbers
offsets = offset_in.split(',')

chain_off = {}

if len(offsets) != len(chains):
    print('Number of entries in chains and offsets should be equal')
    exit(1)
else:
    for chain in chains:
        chain_off[chain] = int(offsets[chains.index(chain)])


mmcif_dict = MMCIF2Dict.MMCIF2Dict(in_file)

chain_labels = mmcif_dict['_atom_site.pdbx_auth_asym_id']

for chain in chains:
    residue_num_change_indices = []
    for i, label in enumerate(chain_labels):
        if label==chain:
            residue_num_change_indices.append(i)

    for atom_id in residue_num_change_indices:
        mmcif_dict['_atom_site.pdbx_auth_seq_id'][atom_id] = str(int(mmcif_dict['_atom_site.pdbx_auth_seq_id'][atom_id]) + chain_off[chain])
        mmcif_dict['_atom_site.auth_seq_id'][atom_id] = str(int(mmcif_dict['_atom_site.auth_seq_id'][atom_id]) + chain_off[chain])
        # mmcif_dict['_atom_site.label_seq_id'][atom_id] = str(int(mmcif_dict['_atom_site.label_seq_id'][atom_id]) + chain_off[chain])

io = MMCIFIO()
io.set_dict(mmcif_dict)
io.save(out_file)
