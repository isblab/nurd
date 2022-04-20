import argparse
import IMP
import RMF
import IMP.rmf
import os, sys

__doc__ = "Get an RMF where beads of the cluster representative model are colored \
based on the mutation data provided"

def parse_args():
    parser = argparse.ArgumentParser(description="Get an RMF where beads of the cluster representative model are colored based on the mutation data provided")
    parser.add_argument('--input', '-i', dest="input", help='representative model in RMF format', default="cluster.0/cluster_center_model.rmf3")
    parser.add_argument('--resolution', '-r', dest="resolution", type=int, help='bead size (residues per bead) for annotating mutations. Must be same as the one used in input generation script', default=1)
    parser.add_argument('--mutation_file','-mf',dest="mutation_file", required=True, type=str, help='location of mutations file e.g. all_beads_penalty.txt')
    parser.add_argument('--output', '-o', dest="output", help='mutation-colored model in RMF format.', default="mutation_colored_cluster_center_model.rmf3")

    return parser.parse_args()


def get_selected_particles(m,input_file, resolution):
    sel0 = None
    inf = RMF.open_rmf_file_read_only(input_file)
    hier = IMP.rmf.create_hierarchies(inf, m)[0]

    IMP.rmf.load_frame(inf, 0)
    m.update()
    sel0 = IMP.atom.Selection(hier, resolution=resolution)

    del inf
    return sel0


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


################################################################################
############################# Main #############################################
################################################################################
if __name__ == '__main__':
    args = parse_args()

    penalties_lst = []
    with open(args.mutation_file,'r') as mf:
        for ln in mf.readlines():
            penalties_lst.append(int(ln.strip().split(',')[1]))

    white_color = (1.0,1.0,1.0)
    grey_color = (0.6,0.6,0.6)                        
    very_lightred_color = (0.99, 0.6, 0.6)
#     lightred_color = (1, 0.3, 0.3)
    red_color = (1,0,0)


    ### Open cluster center model
    if not os.path.exists(args.input):                                              # If cluster center model is not found, exit
        print(f"Cluster center file not found at {args.input}")
        exit(1)

    mdl = IMP.Model()                                                               # Lets begin by creating a model
    sel0 = get_selected_particles(mdl, args.input, args.resolution)                 # Returns an object of IMP.atom.Selection class, that can be used to select particles
    sel_particles = sel0.get_selected_particles()                                   # Selects all the particles from the IMP.Selection class


    ### Now, create a new model in which the beads will be colored according to the mutation data
    mdl_new = IMP.Model()
    p_root = mdl_new.add_particle('System')
    h_root = IMP.atom.Hierarchy.setup_particle(mdl_new,p_root)
    prev_prot = "Dummy.0"

    for i, leaf in enumerate(sel_particles):                                        # Leaf is every individual bead?
        bead_name = get_bead_name(leaf)
        prot_name = bead_name.split(':')[0]
        copy_num = bead_name.split(':')[1]
        start_res = bead_name.split(':')[2]
        end_res = bead_name.split(':')[3]
        curr_prot = prot_name+'.'+copy_num

        if start_res==end_res:                                                      # i.e one residue bead
            res_range = start_res
        else:
            res_range = start_res+'-'+end_res


        if curr_prot != prev_prot:                                                  # If it is a new protein, create a new protein in the new model and add it to the hierarchy of the new model
            p_curr_prot = mdl_new.add_particle(curr_prot)
            h_curr_prot = IMP.atom.Hierarchy.setup_particle(mdl_new,p_curr_prot)
            h_root.add_child(h_curr_prot)                                           # Dont forget to make it a child of root hierarchy in the new model
            prev_prot = curr_prot

        p_new = mdl_new.add_particle(res_range)                                     # Make a new particle in the new model for every entry in the sel_particles
        xyzr_new = IMP.core.XYZR.setup_particle(
                    mdl_new,p_new,IMP.core.XYZR(leaf).get_sphere())
        mass_new = IMP.atom.Mass.setup_particle(mdl_new, p_new, 1.0)

        if penalties_lst[i] == 0:
            color = white_color
        elif penalties_lst[i] == 1:
            color = grey_color
        elif penalties_lst[i] == 2:
            color = very_lightred_color
        elif penalties_lst[i] == 3:
            color = red_color

        c_new  = IMP.display.Colored.setup_particle(mdl_new,p_new, IMP.display.Color(color[0],color[1],color[2]))
        h_new = IMP.atom.Hierarchy.setup_particle(mdl_new,p_new)
        h_curr_prot.add_child(h_new)
    rmf_new = RMF.create_rmf_file(args.output)

    IMP.rmf.add_hierarchy(rmf_new,h_root)
    IMP.rmf.save_frame(rmf_new)

    del rmf_new
