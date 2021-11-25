import argparse
import numpy as np
from matplotlib import pyplot as plt
import IMP
import RMF
import IMP.rmf
import os, sys

__doc__ = "Get contact maps for a given protein pair"


def parse_args():
    parser = argparse.ArgumentParser(description="Get contact maps for a given protein pair")
    parser.add_argument('--inputA', '-ia', dest="ia", help='cluster list of sample A RMFs. (eg. cluster.0.sample_A.txt)', required=True)
    parser.add_argument('--inputB', '-ib', dest="ib", help='cluster list of sample B RMFs. (eg. cluster.0.sample_B.txt)', required=True)
    parser.add_argument('--rmfA', '-ra', dest="ra", help='rmfA file. (eg. A_gsm_clust0.rmf3)', required=True)
    parser.add_argument('--rmfB', '-rb', dest="rb", help='rmfB file. (eg. B_gsm_clust0.rmf3)', required=True)
    parser.add_argument('--clmdl', '-c', dest="cmdl", help='Cluster center model rmf3 file. (eg. cluster_center_model.rmf3)', required=True)
    parser.add_argument('--textA', '-ta', dest="ta", help='Text file associated with rmfA. (eg. A_gsm_clust0.txt)', required=True)

    return parser.parse_args()


def get_bead_name(particle):

    '''
    Input: particle
    Output: bead name in the format molecule_name:copy_number:start_residue:end_residue
    '''

    mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(particle))
    copy_number=IMP.atom.get_copy_index(IMP.atom.Hierarchy(particle))

    if IMP.atom.Fragment.get_is_setup(particle): # CG bead                      ###################### Did not understand get_is_setup
        residues_in_bead = IMP.atom.Fragment(particle).get_residue_indexes()
        bead_name = mol_name+":"+str(copy_number)+":"+str(min(residues_in_bead))+":"+str(max(residues_in_bead))
    else:
        residue_in_bead = str(IMP.atom.Residue(particle).get_index())
        bead_name = mol_name+":"+str(copy_number)+":"+residue_in_bead+":"+residue_in_bead

    return bead_name


def get_nmodels_in_A(ta_file):
    with open(ta_file,'r') as taf:
        ln_count = 0
        for ln in taf.readlines():
            ln_count += 1
    return ln_count

################################################################################
############################# Main #############################################
################################################################################
if __name__ == '__main__':
    args = parse_args()
    nA = get_nmodels_in_A(args.ta)
    RESOLUTION = 10

    # Step 1. Go through each RMF
    # Step 2. In each RMF, get the two selections corresponding to the 2 proteins.
    #          Create a matrix of appropriate size, get the res ranges per bead in a list, if reading first RMF.
    # Step 3. For each bead pair, get the distance between the bead centers. Add it to the distances in the appropriate matrix entry.
    # Step 4. After 1-3, get the average for all matrix entries by dividing by num of models

    # Step 5. Fix the matrix for plotting.

    sizes_dict = {"MTA1":716, "HDAC1":483, "MBD3":292, "RBBP4":426}
    proteins = ["MTA1.0", "MTA1.1", "HDAC1.0", "HDAC1.1", "MBD3.0", "RBBP4.0", "RBBP4.1", "RBBP4.2", "RBBP4.3"]

    # Create list of model indices for sampleA, sample_B
    sample_A_models = []
    sample_B_models = []

    with open(args.ia,'r') as iaf:
         for ln in iaf.readlines():
             sample_A_models.append(int(ln.strip()))
    with open(args.ib,'r') as ibf:
         for ln in ibf.readlines():
             sample_B_models.append(int(ln.strip()))

    sample_A_models.sort()
    sample_B_models.sort()


    nModels = len(sample_A_models)+len(sample_B_models)
    print(f'Total number of models:\t {nModels}')


    # Step 1. Go through each RMF
    sample_A_rmf = args.ra
    sample_B_rmf = args.rb


    mdl_a = IMP.Model()
    rmf_fh_a = RMF.open_rmf_file_read_only(sample_A_rmf)
    h_a = IMP.rmf.create_hierarchies(rmf_fh_a, mdl_a)[0]

    for index1 in range(len(proteins)):
        for index2 in range(index1+1,len(proteins)):
            prot1 = proteins[index1]
            prot2 = proteins[index2]

            p1name = prot1.split('.')[0]
            p2name = prot2.split('.')[0]

            copy1 = prot1.split('.')[1]
            copy2 = prot2.split('.')[1]

            size1 = sizes_dict[prot1.split('.')[0]]
            size2 = sizes_dict[prot2.split('.')[0]]
            # print(prot1, copy1, prot2, copy2)

            print(prot1,'\t',prot2)
            print('Reading rmfA file')
            i = 0

            distances = {}


            for sa in sample_A_models:
                IMP.rmf.load_frame(rmf_fh_a, int(sa))
                mdl_a.update()
                sel1 = IMP.atom.Selection(h_a, resolution=RESOLUTION, molecule=p1name, copy_index=int(copy1))
                sel2 = IMP.atom.Selection(h_a, resolution=RESOLUTION, molecule=p2name, copy_index=int(copy2))


                if i%10 == 0:
                    print(f'Reading frame {i} out of {len(sample_A_models)}')
                i+=1
                for bead1 in sel1.get_selected_particles():
                    for bead2 in sel2.get_selected_particles():
                        # print(bead1,bead2,'##############')
                        dist = IMP.core.get_distance(IMP.core.XYZR(bead1),IMP.core.XYZR(bead2))
                        if dist<0:
                            dist=0
                        key = get_bead_name(bead1) + "--" + get_bead_name(bead2)
                        if key in distances.keys():
                            distances[key] += dist
                        else:
                            distances[key] = dist


            mdl_b = IMP.Model()
            rmf_fh_b = RMF.open_rmf_file_read_only(sample_B_rmf)
            h_b = IMP.rmf.create_hierarchies(rmf_fh_b, mdl_b)[0]

            print('Reading rmfB file')
            i = 0
            for sb in sample_B_models:
                fname = sb - nA
                IMP.rmf.load_frame(rmf_fh_b, fname)
                mdl_b.update()
                sel1 = IMP.atom.Selection(h_b, resolution=RESOLUTION, molecule=p1name, copy_index=int(copy1))
                sel2 = IMP.atom.Selection(h_b, resolution=RESOLUTION, molecule=p2name, copy_index=int(copy2))

                if i%10 == 0:
                    print(f'Reading frame {i} out of {len(sample_B_models)}')
                i+=1
                for bead1 in sel1.get_selected_particles():
                    for bead2 in sel2.get_selected_particles():
                        dist = IMP.core.get_distance(IMP.core.XYZR(bead1),IMP.core.XYZR(bead2))
                        if dist<0:
                            dist=0
                        key = get_bead_name(bead1) + "--" + get_bead_name(bead2)
                        # print(key)
                        if key in distances.keys():
                            distances[key] += dist
                        else:
                            distances[key] = dist

            for k in distances.keys():
                distances[k] = distances[k]/nModels
            nBeadsA = len(sel1.get_selected_particles())+1
            nBeadsB = len(sel2.get_selected_particles())+1

            print('Distance Matrix Size: ',size1,size2,'\n')
            mat = np.zeros((size1,size2))

            for k in distances.keys():
                # print(k)
                bead1 = k.split('--')[0]
                bead2 = k.split('--')[1]
                # print(bead2)
                res_start_1 = int(bead1.split(':')[2])
                res_end_1 =   int(bead1.split(':')[3])
                res_start_2 = int(bead2.split(':')[2])
                res_end_2 = int(bead2.split(':')[3])
                for i in range(res_start_1,res_end_1+1):
                    for j in range(res_start_2,res_end_2+1):
                        mat[i,j] = distances[k]


            np.savetxt(f'{p1name}.{copy1}-{p2name}.{copy2}_Distance-matrix.csv',mat,delimiter=',')


            plt.figure(dpi=400)
            plt.imshow(mat, cmap='hot')
            plt.xlabel(prot2)
            plt.ylabel(prot1)
            plt.colorbar()
            plt.savefig(f'{p1name}.{copy1}-{p2name}.{copy2}_contact-map.png',dpi=600)

            # plt.show()
