import glob
import chimera
from chimera import openModels
from chimera import runCommand

mmcifs = ['mmcifs/chainname_changed_aligned_offset_corrected_HDAC1-MTA1_ext.pdb.pdb.cif', 'mmcifs/chainname_changed_aligned_offset_corrected_MTA1.0_1-164_BAH.pdb.pdb.cif', 'mmcifs/chainname_changed_aligned_offset_corrected_MTA1.0_334-353_ns.pdb.pdb.cif', 'mmcifs/chainname_changed_aligned_offset_corrected_MTA1.0_389-431_ZF.pdb.pdb.cif', 'mmcifs/chainname_changed_aligned_offset_corrected_MTA1.1_1-164_BAH.pdb.pdb.cif', 'mmcifs/chainname_changed_aligned_offset_corrected_MTA1.1_334-353_ns.pdb.pdb.cif', 'mmcifs/chainname_changed_aligned_offset_corrected_MTA1.1_389-431_ZF.pdb.pdb.cif', 'mmcifs/chainname_changed_aligned_offset_corrected_RBBP4.0-MTA1_468-546_5FXY.pdb.pdb.cif', 'mmcifs/chainname_changed_aligned_offset_corrected_RBBP4.1-MTA1_468-546_5FXY.pdb.pdb.cif', 'mmcifs/chainname_changed_aligned_offset_corrected_RBBP4.2-MTA1_670-691_4PBZ.pdb.pdb.cif', 'mmcifs/chainname_changed_aligned_offset_corrected_RBBP4.3-MTA1_670-691_4PBZ.pdb.pdb.cif']
runCommand('open '+'multiscaling_striped_cluster_center_model.rmf3')
for p in mmcifs:
	runCommand('open '+p)