from Bio import PDB
import glob


########################################################################################
######################################## Inputs ########################################
########################################################################################

pdb_files = ['pdbs/HDAC1-MTA1_ext.pdb',\
            'pdbs/MTA1.0_1-164_BAH.pdb',\
            'pdbs/MTA1.0_334-353_ns.pdb',\
            'pdbs/MTA1.0_389-431_ZF.pdb',\
            'pdbs/MTA1.1_1-164_BAH.pdb',\
            'pdbs/MTA1.1_334-353_ns.pdb',\
            'pdbs/MTA1.1_389-431_ZF.pdb',\
            'pdbs/RBBP4.0-MTA1_468-546_5FXY.pdb',\
            'pdbs/RBBP4.1-MTA1_468-546_5FXY.pdb',\
            'pdbs/RBBP4.2-MTA1_670-691_4PBZ.pdb',\
            'pdbs/RBBP4.3-MTA1_670-691_4PBZ.pdb',\
            'pdbs/MBD3.0_1-71.pdb',\
            'pdbs/GATAD2B-Anumbers-MBD3-cc.pdb']

offsets = [{'A':0,'B':-1000,'C':-2000,'D':-3000},\
            {'J':0},\
            {'A':0},\
            {'A':0},\
            {'J':0},\
            {'A':0},\
            {'A':0},\
            {'A':0,'B':0},\
            {'A':0,'B':0},\
            {'A':0,'B':0},\
            {'A':0,'B':0},\
            {'A':0},\
            {'A':0,'B':5}]


########################################################################################
###################################### Actual work #####################################
########################################################################################

for file_index in range(len(pdb_files)):
    pdbfile = pdb_files[file_index]
    offset = offsets[file_index]
    print(file_index,pdbfile,offsets[file_index])
    pdb_io = PDB.PDBIO()
    pdb_parser = PDB.PDBParser()
    structure = pdb_parser.get_structure(" ", pdbfile)

    for model in structure:
        for chain in model:
            print(chain.id)
            for res in chain.get_residues():
                res_id_list = list(res.id)
                res_id_list[1] = res_id_list[1] + offset[chain.id]
                res.id = tuple(res_id_list)

    pdb_io.set_structure(structure)
    pdb_io.save('pdbs/offset_corrected/offset_corrected_'+pdbfile.split('/')[-1])
