import os,sys
from Bio.PDB import *


########################################################################################
######################################## Inputs ########################################
########################################################################################

cif_files = ['mmcif/aligned_offset_corrected_HDAC1-MTA1_ext.pdb.pdb.cif',\
            'mmcif/aligned_offset_corrected_MTA1.0_1-164_BAH.pdb.pdb.cif',\
            'mmcif/aligned_offset_corrected_MTA1.0_334-353_ns.pdb.pdb.cif',\
            'mmcif/aligned_offset_corrected_MTA1.0_389-431_ZF.pdb.pdb.cif',\
            'mmcif/aligned_offset_corrected_MTA1.1_1-164_BAH.pdb.pdb.cif',\
            'mmcif/aligned_offset_corrected_MTA1.1_334-353_ns.pdb.pdb.cif',\
            'mmcif/aligned_offset_corrected_MTA1.1_389-431_ZF.pdb.pdb.cif',\
            # 'mmcif/aligned_offset_corrected_RBBP4.0-MTA1_468-546_5FXY.pdb.pdb.cif',\
            # 'mmcif/aligned_offset_corrected_RBBP4.1-MTA1_468-546_5FXY.pdb.pdb.cif',\
            # 'mmcif/aligned_offset_corrected_RBBP4.2-MTA1_670-691_4PBZ.pdb.pdb.cif',\
            # 'mmcif/aligned_offset_corrected_RBBP4.3-MTA1_670-691_4PBZ.pdb.pdb.cif',\
            'mmcif/aligned_offset_corrected_GATAD2B.0-Anumbers-MBD3-cc.pdb.pdb.cif',\
            'mmcif/aligned_offset_corrected_MBD3.0_1-71.pdb.pdb.cif',\
            'mmcif/aligned_offset_corrected_GATAD2B.1-Anumbers-MBD3-cc.pdb.pdb.cif',\
            'mmcif/aligned_offset_corrected_MBD3.1_1-71.pdb.pdb.cif']

chain_assgns = [{'A':'elmsant1','B':'hdac1','C':'elmsant2','D':'hdac2'},\
            {'J':'bah1'},\
            {'A':'h1'},\
            {'A':'zf1'},\
            {'J':'bah2'},\
            {'A':'h2'},\
            {'A':'zf2'},\
            # {'A':'rbbp0','B':'mta1r1'},\
            # {'A':'rbbp1','B':'mta2r1'},\
            # {'A':'rbbp2','B':'mta1r2'},\
            # {'A':'rbbp3','B':'mta2r2'},\
            {'A':'gatacc1','B':'mbdcc1'},\
            {'A':'mbdmbd1'},\
            {'A':'gatacc2','B':'mbdcc2'},\
            {'A':'mbdmbd2'}]


########################################################################################
###################################### Actual work #####################################
########################################################################################

for file_index in range(len(cif_files)):
    cif_file = cif_files[file_index]
    new_name_assn = chain_assgns[file_index]

    mmcif_dict = MMCIF2Dict.MMCIF2Dict(cif_file)

    old_chain_lables = mmcif_dict['_atom_site.auth_asym_id']
    new_chain_lables = []

    for label in old_chain_lables:
        if label in new_name_assn.keys():
            new_chain_lables.append(new_name_assn[label])
        else:
            new_chain_lables.append(label)

    mmcif_dict['_atom_site.pdbx_auth_asym_id'] = new_chain_lables
    mmcif_dict['_atom_site.auth_asym_id'] = mmcif_dict['_atom_site.pdbx_auth_asym_id']
    mmcif_dict['_atom_site.label_asym_id'] = mmcif_dict['_atom_site.pdbx_auth_asym_id']

    out_file = str('mmcif/chainname_changed_' + cif_file.split('/')[-1])

    mmcif_io2 = MMCIFIO()
    mmcif_io2.set_dict(mmcif_dict)
    mmcif_io2.save(out_file)
