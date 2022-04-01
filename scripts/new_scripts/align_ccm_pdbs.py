import os, sys
import IMP
import RMF
import IMP.core
import IMP.rmf
import IMP.pmi.analysis
import glob
from Bio.PDB import *


###################################################################################################
############################################# Inputs ##############################################
###################################################################################################

input_file = '../cluster_center_model.rmf3'
# pdb_dir = 'pdb_files/pdbs/offset_corrected/offset_corrected_*'

pdb_files = ['offset_corrected/offset_corrected_HDAC1-MTA1_ext.pdb',\
            'offset_corrected/offset_corrected_MTA1.0_1-164_BAH.pdb',\
            'offset_corrected/offset_corrected_MTA1.0_334-353_ns.pdb',\
            'offset_corrected/offset_corrected_MTA1.0_389-431_ZF.pdb',\
            'offset_corrected/offset_corrected_MTA1.1_1-164_BAH.pdb',\
            'offset_corrected/offset_corrected_MTA1.1_334-353_ns.pdb',\
            'offset_corrected/offset_corrected_MTA1.1_389-431_ZF.pdb',\
            # 'offset_corrected/offset_corrected_RBBP4.0-MTA1_468-546_5FXY.pdb',\
            # 'offset_corrected/offset_corrected_RBBP4.1-MTA1_468-546_5FXY.pdb',\
            # 'offset_corrected/offset_corrected_RBBP4.2-MTA1_670-691_4PBZ.pdb',\
            # 'offset_corrected/offset_corrected_RBBP4.3-MTA1_670-691_4PBZ.pdb',\
            'offset_corrected/offset_corrected_MBD3.0_1-71.pdb',\
            'offset_corrected/offset_corrected_GATAD2B.0-Anumbers-MBD3-cc.pdb',\
            'offset_corrected/offset_corrected_MBD3.1_1-71.pdb',\
            'offset_corrected/offset_corrected_GATAD2B.1-Anumbers-MBD3-cc.pdb']


# The all proteins list has the following architecture:
# [{protein:{chain_id,residue range}}, {protein:{chain_id,residue range}]
# The order of entries in the offset list must be the same as that in the pdb_files list

all_proteins = [{'MTA1':{'A':[0,range(165,334)]},'MTA1':{'C':[1,range(165,334)]}, 'HDAC1':{'B':[0,range(8,377)]},'HDAC1':{'D':[1,range(8,377)]}},\
            {'MTA1':{'J':[0,range(1,165)]}},\
            {'MTA1':{'A':[0,range(334,354)]}},\
            {'MTA1':{'A':[0,range(389,432)]}},\
            {'MTA1':{'J':[1,range(1,165)]}},\
            {'MTA1':{'A':[1,range(334,354)]}},\
            {'MTA1':{'A':[1,range(389,432)]}},\
            # {'RBBP4':{'A':[0,range(10,412)]},'MTA1':{'B':[0,range(468,547)]}},\
            # {'RBBP4':{'A':[1,range(10,412)]},'MTA1':{'B':[1,range(468,547)]}},\
            # {'RBBP4':{'A':[2,range(2,412)]},'MTA1':{'B':[0,range(670,692)]}},\
            # {'RBBP4':{'A':[3,range(2,412)]},'MTA1':{'B':[1,range(670,692)]}},\
            {'MBD3':{'A':[0,range(1,72)]}},\
            {'MBD3':{'B':[0,range(221,250)]},'P66A':{'A':[0,range(137,179)]}},\
            {'MBD3':{'A':[1,range(1,72)]}},\
            {'MBD3':{'B':[1,range(221,250)]},'P66A':{'A':[1,range(137,179)]}}]


# chains =
# residue_range = range(1,165)


###################################################################################################
##################################### Get transformations #########################################
###################################################################################################

for file_index in range(len(pdb_files)):
    pdb_file = pdb_files[file_index]
    proteins = all_proteins[file_index]
    # print(proteins)
    ccm_mdl = IMP.Model()
    ccm = RMF.open_rmf_file_read_only(input_file)
    hier = IMP.rmf.create_hierarchies(ccm, ccm_mdl)[0]
    IMP.rmf.load_frame(ccm, 0)
    ccm_mdl.update()
    pdb_ca_mdl = IMP.Model()
    pdb_ca = IMP.atom.read_pdb(pdb_file,pdb_ca_mdl,IMP.atom.CAlphaPDBSelector())
    pdb_ca_mdl.update()

    new_mdl = IMP.Model()
    reload = IMP.atom.read_pdb(pdb_file, new_mdl)

    coords_pdb_ca = {}
    coords_ccm = {}

    for prot in proteins.keys():
        for chain_id in proteins[prot]:
            protein_name = prot

            sel_ca_pdb = IMP.atom.Selection(pdb_ca,resolution=1,chain_id=chain_id,residue_indexes=[i for i in proteins[prot][chain_id][1]]).get_selected_particles()
            sel_ccm = IMP.atom.Selection(hier,resolution=1,molecule=protein_name,copy_index=proteins[prot][chain_id][0],residue_indexes=[i for i in proteins[prot][chain_id][1]]).get_selected_particles()

            # print(len(sel_ccm),len(sel_ca_pdb))
            # Remove coarse grained beads
            new_ccm_sel = []
            for selection in sel_ccm:
                if not IMP.atom.Fragment.get_is_setup(selection):
                    new_ccm_sel.append(selection)
                    # print(selection)

            # print(len(sel_ca_pdb),'\n\n', len(new_ccm_sel))
            # for i in range(len(new_ccm_sel)):
                # print(new_ccm_sel[i],sel_ca_pdb[i])


            coords_pdb_ca[protein_name] = [IMP.core.XYZ(i).get_coordinates() for i in sel_ca_pdb]
            coords_ccm[protein_name] = [IMP.core.XYZ(i).get_coordinates() for i in new_ccm_sel]
            # print(len(coords_pdb_ca[protein_name]),len(coords_ccm[protein_name]))
    _, transformation = IMP.pmi.analysis.Alignment(query=coords_pdb_ca, template=coords_ccm).align()
    print(transformation)


    ###################################################################################################
    #################################### Transform and write PDB ######################################
    ###################################################################################################

    IMP.atom.transform(reload, transformation)
    IMP.atom.write_pdb(reload, f"./aligned_{pdb_file.split('/')[-1]}.pdb")
