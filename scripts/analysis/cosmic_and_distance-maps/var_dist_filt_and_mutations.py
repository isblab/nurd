import numpy as np
import os
import sys
import glob

def parse_mutations_file(mutations_csv):
    # keys_lst = ['MTA1','MBD3','HDAC1','RBBP4']
    mutations = {"MBD3":[],"MTA1":[],"HDAC1":[],"RBBP4":[],'P66A':[]}
    with open(mutations_csv,'r') as mutf:
        for ln in mutf.readlines():
            if ln.startswith('>'):
                continue
            else:
                ln1 = ln.strip().split(',')
                k = ln1[0]
                for i in range(1,len(ln1)):
                    mutations[k].append(int(ln1[i]))
    return mutations


###############################################################################


mutations_csv = sys.argv[1]
dist_threshold = sys.argv[2]
dist_threshold = float(dist_threshold)

# STEP 1 Parse mutations
mutations = parse_mutations_file(mutations_csv) # {"MBD3":[],"MTA1":[],}

# STEP 2 Create interface mutated residues list
interface_mutated = {'MBD3':[],'MTA1':[],'HDAC1':[],'RBBP4':[],'P66A':[]}
all_pairs = []

# STEP 2 Go through each distance matrix and identify interface mutations for each protein
distance_matrices = glob.glob(os.getcwd()+'/*Distance-matrix*')

for dm in distance_matrices:
    dm1 = dm.split('/')[-1]
    prot1 = dm1.rstrip("_Distance-matrix").split("-")[0].split('.')[0]
    prot2 = dm1.rstrip("_Distance-matrix").split("-")[1].split('.')[0] # "prot1.copy1-prot2.copy2_distance_matrix"

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

                if float(ax) <= dist_threshold:
                    pair = f'{prot1} : {y} - {prot2} : {x}'
                    all_pairs.append(pair)
                    if y in mutations[prot1]:
                        interface_mutated[prot1].append(y)

                    if x in mutations[prot2] :
                        interface_mutated[prot2].append(x)


# TODO make list unique
unique_interface_mutated = {}
for p in interface_mutated.keys():
    unique_interface_mutated[p] = list(set(interface_mutated[p]))
    print(p,unique_interface_mutated[p])

with open('mapped_mutations.txt','w') as outf:
    for protein in unique_interface_mutated.keys():
        out_str1 = '>' + protein + '\n'
        out_str2 = ', '.join(str(element) for element in sorted(unique_interface_mutated[protein])) + '\n'
        outf.write(out_str1)
        outf.write(out_str2)


with open(f'pairs_{dist_threshold}a.txt','w') as pf:
    for i in all_pairs:
        i_to_write = i + '\n'
        pf.write(i_to_write)
