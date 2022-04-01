from Bio.PDB import *
import glob
import sys

pdb_io = PDBIO()
mmcif_io = MMCIFIO()

pdb_dir = sys.argv[1] +'*'          # The path to pdb files with the first few characters for the glob
pdbfiles = glob.glob(pdb_dir)

for pdbfile in pdbfiles:
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(" ", pdbfile)

    mmcif_io.set_structure(structure)
    mmcif_io.save('mmcifs/'+pdbfile.split('/')[-1][0:-4]+'.cif')