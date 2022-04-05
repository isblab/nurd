import argparse
import IMP
import RMF
import IMP.rmf
import os, sys, glob
from color_mutation import get_selected_particles, get_bead_name

__doc__ = "Generate input for the color_mutations script"


def parse_args():
    parser = argparse.ArgumentParser(description='Generate input for the color_mutations script')
    parser.add_argument('-m', '--mutations_file_prefix', dest='input', help='Input mutations csv file')
    parser.add_argument('-c', '--cl_model', dest='cl_model', help= 'Cluster center model')
    parser.add_argument('-r', '--raw_mut', dest='raw_inf', default=None, help='Mutations.txt file')
    parser.add_argument('-res', '--resolution', type=int, dest='res', default=1, help='Resolution')
    return parser.parse_args()


def get_all_beads(model, cl_model,res):
    sel0 = get_selected_particles(model,cl_model,res)
    sel_particles = sel0.get_selected_particles()                                   # Selects all the particles from the IMP.Selection class

    beads = []
    for i, leaf in enumerate(sel_particles):
        bead_name = get_bead_name(leaf)
        beads.append(bead_name)
    return beads


def write_bead_name_file(beads_list,penalty_list,res):
    print('Writting to all_bead.txt')
    with open(f"all_beads_penalty_{res}.txt",'w') as outf:
        for bi, pi in zip(beads_list,penalty_list):
            out_str = str(bi) + ',' + str(pi) + '\n'
            outf.write(out_str)



def get_mol_cp_dict(beads_list):
    mols = []
    cp_num = {}
    for bead in beads:                                                              # mol_name:cp_num:st_res:end_res
        mol_name = bead.split(':')[0]
        if mol_name not in mols:
            mols.append(mol_name)
    for mol in mols:
        mol_name = 'string'
        for bead in beads:
            if bead.split(':')[0] == mol:
                cp_num[mol] = int(bead.split(':')[1])
    return cp_num


def get_all_residues(beads_list):
    all_res = []
    for bead in beads:
        if bead.split(':')[2] == bead.split(':')[3]:
            ln = bead.split(':')[0] +':'+ bead.split(':')[1] +':'+ bead.split(':')[2]
            all_res.append(ln)
        else:
            for i in range(int(bead.split(':')[2]),int(bead.split(':')[2])+1):
                ln = bead.split(':')[0] +':'+ bead.split(':')[1] +':'+ str(i)
                all_res.append(ln)
    return all_res


def write_all_res_file(all_residues_list):
    print('Writting to all_residues.txt')
    with open('all_residues.txt','w') as outf:
        for residue in all_residues_list:
            outline = residue + '\n'
            outf.write(outline)




def parse_mutation_files(fname,d):
    with open(fname,'r') as inf:
        for ln in inf.readlines():
            if not len(ln)==0:
                if ln.startswith('>1to2'):
                    penalty = 1
                elif ln.startswith('>3to4'):
                    penalty = 2
                elif ln.startswith('>5plus'):
                    penalty = 3
                else:
                    prot = ln.strip().split(',')[0]
                    res_lst = ln.strip().split(',')[1:]
                    new_res_lst = []
                    for elem in res_lst:
                        if '-' in elem:
                            for i in range(int(elem.split('-')[0]),int(elem.split('-')[1])+1):
                                new_res_lst.append((i,penalty))
                        else:
                            new_res_lst.append((int(elem),penalty))
                    # print(new_res_lst)
                    d[prot].extend(new_res_lst)
    # print(d)
    return d



################################################################################
################################## Main ########################################
################################################################################
if __name__ == '__main__':
    args = parse_args()

    mdl = IMP.Model()
    beads = get_all_beads(mdl, args.cl_model,args.res)
    # write_bead_name_file(beads)

mutations_files = glob.glob(args.input+'*')
print(mutations_files)

cosmic_m = {"MTA1":[],"HDAC1":[],"RBBP4":[],"MBD3":[],"P66A":[]}
for file in mutations_files:
    cosmic_m = parse_mutation_files(file,cosmic_m)

# print(cosmic_m)
penalties = []
for bead in beads:
    prot = bead.split(':')[0]
    res_range = [r for r in range(int(bead.split(':')[2]),int(bead.split(':')[3])+1)]  
    # print(bead,res_range)
    penalty = 0

    for mut in cosmic_m[prot]:
        if mut[0] in res_range:
            if mut[1]>penalty:
                penalty = mut[1]
    penalties.append(penalty)
write_bead_name_file(beads,penalties,args.res)