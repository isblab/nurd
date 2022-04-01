import argparse
import IMP
import RMF
import IMP.rmf
import os, sys
import random

__doc__ = "Given a mutations file, generate random mutations. \
This will generate 'n' number of random mutations for all the proteins in the mutations file (where 'n = number of true mutaions for that protein)"

################################################################################
#################################### Inputs ####################################
################################################################################

n_residues = {'MTA1':(1,716), 'HDAC1':(1,483), 'RBBP4':(1,426), 'MBD3':(1,296), 'P66A':(1,531)}


################################################################################
################################## Functions ###################################
################################################################################

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mutation_file','-mf',dest="mutation_file", required=True, type=str, help='location of mutations csv file')
    parser.add_argument('--output', '-o', dest="output", help='mutation-colored model in RMF format.', required=True)
    return parser.parse_args()



def get_mutations(input_file):
    mutations = {}
    with open(input_file,'r') as mf:
        for ln in mf.readlines():
            if ln.strip():
                ln = ln.strip()
                if (not ln.startswith('>')) and not ln.startswith('\n'):
                    protein = ln.split(',')[0]
                    if protein not in mutations.keys():
                        values1 = ln.split(',')[1:]
                        values2 = []
                        for val in values1:
                            values2.append(int(val))
                        mutations[protein] = values2
                    else:
                        for mut in ln.split(',')[1:]:
                            mutations[protein].append(int(mut))
    return mutations



def generate_random_mutations(mutations_dict,residues):
    random_mutations1 = {}
    random_mutations2 = {}

    for prot in mutations_dict.keys():
        n_true_mutations = len(mutations_dict[prot])
        if not prot in random_mutations1.keys():
            random_mutations1[prot] = []
        
        while len(random_mutations1[prot]) < n_true_mutations:
            new_mut = random.randint(residues[prot][0],residues[prot][1])
            if new_mut not in random_mutations1[prot]:
                random_mutations1[prot].append(new_mut)

    for protein in random_mutations1.keys():
        random_mutations2[protein] = sorted(random_mutations1[protein])

    return random_mutations2


################################################################################
##################################### Main #####################################
################################################################################

args = parse_args()

all_mutations = get_mutations(args.mutation_file)

random_mutations = generate_random_mutations(all_mutations,n_residues)
# print(random_mutations)

with open(args.output,'w') as outf:
    out_str = ''
    for protein in random_mutations.keys():
        out_str = f"{out_str}{protein},"
        for mutation in random_mutations[protein]:
            out_str = f"{out_str} {mutation},"
        out_str = f"{out_str}\n\n"
    # print(out_str)
    outf.write(out_str)


'''
Steps:
1. Get all the mutations
2. Get the number of mutations and the residue range for each protein
3. Generate n random numbers between 1 and the length of the protein. n should be equal to the number of true mutations on that protein
4. Export a file that has these n random mutations in the format similar to the input mutations file

'''