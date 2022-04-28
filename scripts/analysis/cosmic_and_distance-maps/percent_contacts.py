import numpy as np
import os
import sys
import glob

dist_threshold = sys.argv[1]
dist_threshold = float(dist_threshold)
i = 0
interface_residues = []
distance_matrices = glob.glob(os.getcwd()+'/*Distance-matrix*')

for dm in distance_matrices:
    dm1 = dm.split('/')[-1]
    prot1 = dm1.rstrip("_Distance-matrix").split("-")[0]
    prot2 = dm1.rstrip("_Distance-matrix").split("-")[1].rstrip('_Distance') # "prot1.copy1-prot2.copy2_distance_matrix"
    print(f'{i}\t{prot1}\t{prot2}')
    i += 1
    with open(dm,'r') as inf:
        distances = np.loadtxt(inf, delimiter=',')
        for y in range(1,len(distances)):
            '''
            (y,x) will be be index of any given distance
            where y will be residue number of protein1 and x will be residue number of protein2
            ax is the individual distance
            '''
            ay = distances[y]
            for x in range(1,len(ay)):
                ax = ay[x]
                if ax <= dist_threshold:
                    res1 = f'{prot1} : {y}'
                    res2 = f'{prot2} : {x}'
                    if res1 not in interface_residues:
                        interface_residues.append(res1)
                    if res2 not in interface_residues:
                        interface_residues.append(res2)

total_residues = (2*715) + (2*482) + (4*425) + 295 + 43         # 2MTA + 2HDAC + 4RBBP + MBD + GATA
interface_residues = list(set(interface_residues))
res_at_interface = len(interface_residues)

percent = (res_at_interface/total_residues)*100
print(f'Residues at interface: {res_at_interface} \t Total residues: {total_residues} \t Percentage: {percent}')
