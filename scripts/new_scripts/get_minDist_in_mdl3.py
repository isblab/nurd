import os, sys,glob
import pandas as pd
import argparse
import subprocess
import IMP
import IMP.core
import IMP.pmi.tools
import RMF
import IMP.rmf


'''
Steps:
1. Parse all the crosslinks
2. Read the stat file and fetch the lines that correspond to the relevant crosslinks
3.

'''
structured_regions = {'MTA1':[('bah1','bah2',range(1,165)),('elmsant1','elmsant2',range(165,334)),('h1','h2',range(334,354)),('zf1','zf2',range(389,432)),\
                            ('mta1r1','mta2r1',range(468,547)),('mta1r2','mta2r2',range(670,692))],\
                    'HDAC1':[('hdac1','hdac2',range(8,377))],\
                    'RBBP4':[('rbbp0','rbbp1','rbbp2','rbbp3',range(1,412))],\
                    'MBD3':[('mbdmbd1','mbdmbd2',range(1,72)),('mbdcc1','mbdcc2',range(221,250))],\
                    'P66A':[('gatacc1','gatacc2',range(137,179))]}


model_ids = {'elmsant1':1,'elmsant2':1,'hdac1':1,'hdac2':1,\
             'bah1':2,\
             'bah2':5,\
             'h1':3,\
             'h2':6,\
             'zf1':4,\
             'zf2':7,\
             # 'mta1r1':8,'rbbp0':8,\
             # 'mta2r1':9,'rbbp1':9,\
             # 'mta1r2':10,'rbbp2':10,\
             # 'mta2r2':11,'rbbp3':11,\
             'mbdcc1':8,'gatacc1':8,\
             'mbdcc2':9,'gatacc2':9,\
             'mbdmbd1':10,\
             'mbdmbd2':11}                                         # Chain name: model_id

__doc__ = "Get the minimum distance pair for each crosslink"

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', dest="input", help='Representative model in RMF format', default="cluster.0/cluster_center_model.rmf3")
    parser.add_argument('-ca', dest="csv_file_A", help='csvA file')
    parser.add_argument('-cb', dest="csv_file_B", help='csvB file')
    parser.add_argument('-m', dest='mdls', help='Path to models (Something like: "/sampcon/cluster.0.all.txt")')
    parser.add_argument('-r', dest='run_dirs', type=str, help='Path to the run directories')
    parser.add_argument('--type', dest='xltype', type=str, help='XL type. Eg. adh, dmtmm, etc. Enter it exactly as in the stat files')
    parser.add_argument('-x', dest='xlviol', type=str, help='Path to the file containing violated crosslinks')
    parser.add_argument('-t', dest='threshold', type=float, help='Distance threshold for the crosslinks')
    parser.add_argument('--mmcif', dest='cif_dir', type=str, help='Path to the structure .cif files')

    return parser.parse_args()


def strip_multiscaling(cluster_center_model_file):
    mdl = IMP.Model()
    ccm = RMF.open_rmf_file_read_only(cluster_center_model_file)
    hier = IMP.rmf.create_hierarchies(ccm, mdl)[0]
    IMP.rmf.load_frame(ccm, 0)
    mdl.update()

    sel_particles = IMP.atom.Selection(hierarchy=hier,resolution=1)
    particles = sel_particles.get_selected_particles()

    particle_details = {}
    cg_beads = {}
    for i, leaf in enumerate(particles):
        name = get_bead_name(leaf)
        min_res = int(name.split(':')[1].split('-')[0])
        max_res = int(name.split(':')[1].split('-')[1])
        if min_res == max_res:
            new_name = name.split('-')[0]
            particle_details[(new_name)] = i
        else:
            while min_res <= max_res:
                new_name = name.split(':')[0]+':'+str(min_res)
                particle_details[(new_name)] = i
                cg_beads[(new_name)] = i
                min_res+=1

    m_new = IMP.Model()
    p_root = m_new.add_particle("System")
    h_root = IMP.atom.Hierarchy.setup_particle(m_new,p_root)

    prev_prot ="DUMMY.0"

    for i,leaf in enumerate(particles):
        bead_name = get_bead_name(leaf)
        curr_prot = bead_name.split(':')[0]
        prot = curr_prot.split('.')[0]
        cp = curr_prot.split('.')[1]
        start_res = bead_name.split(':')[1].split('-')[0]
        end_res = bead_name.split(':')[1].split('-')[1]
        if start_res ==end_res:
            res_range =start_res
        else:
            res_range = start_res+"-"+end_res


        if curr_prot != prev_prot:
            # Create a new hierarchy particle for the protein
            p_curr_prot = m_new.add_particle(prot)
            # Add to the new model's hierarchy
            h_curr_prot = IMP.atom.Hierarchy.setup_particle(m_new,p_curr_prot)
            h_root.add_child(h_curr_prot)
            prev_prot = curr_prot

        p_new = m_new.add_particle(res_range)
        mass_new = IMP.atom.Mass.setup_particle(m_new,p_new,1.0)
        xyzr_new = IMP.core.XYZR.setup_particle(m_new,p_new,IMP.core.XYZR(leaf).get_sphere())
        r,g,b = str(IMP.display.Colored(leaf).get_color()).split(' ')
        c_new  = IMP.display.Colored.setup_particle(m_new,p_new,IMP.display.Color(float(r),float(g),float(b)))
        h_new = IMP.atom.Hierarchy.setup_particle(m_new,p_new)
        h_curr_prot.add_child(h_new)

    rmf_new = RMF.create_rmf_file('multiscaling_striped_cluster_center_model.rmf3')
    IMP.rmf.add_hierarchy(rmf_new,h_root)
    IMP.rmf.save_frame(rmf_new)
    return particle_details,cg_beads


def get_models_from_file(fil):
    lst = []
    fl = open(fil,'r')
    for ln in fl.readlines():
        lst.append(int(ln.strip()))
    fl.close()
    return lst



def get_model_identity_from_files(csv_file_A,csv_file_B):
    model_ids = {}      # key by model index. Value = (runid,replicaid,frameid) tuple for each model
    m_i = 0             # model index, first all A models followed by all B models

    for fl in [csv_file_A,csv_file_B]:
        df = pd.read_csv(fl)
        unique_replicas = df['rmf3_file'].unique() # the models are stored in the order of RMF3 files

        for rep in unique_replicas:
            rep_df = df[df['rmf3_file']==rep] # get all the rows where the RMF3 file is the current file

            for i,row in rep_df.iterrows():
                runid = (row['rmf3_file'].split('/')[0]).split('_')[1]
                replicaid = (row['rmf3_file'].split('//')[-1]).split('.')[0]
                frameid = row['rmf_frame_index']
                # print(runid,replicaid,frameid)
                model_ids[m_i] = (int(runid),int(replicaid),int(frameid))
                m_i +=1
    return model_ids



def get_violated_xls(xl_min_distances_file):
    violated_xls =[]
    with open(xl_min_distances_file,'r') as xlvf:
        for ln in xlvf.readlines():
            min_dist = ln.strip().split(',')[-1]
            if float(min_dist) > args.threshold:
                violated_xls.append(','.join(ln.split(',')[0:-1]))
    return violated_xls



def get_min_dist(stat_lines,xl_type):
    xl_lines = []
    min_dist_xl ={}

    for ln in stat_lines:
        if ln.startswith(f'CrossLinkingMassSpectrometryRestraint_Distance_|{xl_type}'):
            ln = ','.join(ln.split('|')[3:7])+':'+ln.split(' ')[-1]
            key = ln.split(',')[0].split('.')[0] + ',' +ln.split(',')[1] + ',' + ln.split(',')[2].split('.')[0] + ',' + ln.split(',')[3].split(':')[0]

            if key not in min_dist_xl.keys():
                min_dist_xl[key] = ln
            else:
                if float(ln.split(':')[-1]) < float(min_dist_xl[key].split(':')[-1]):
                    min_dist_xl[key] = ln
    return list(min_dist_xl.values())



def get_bead_name(particle):
    '''
    Input: particle
    Output: bead name in the format molecule_name:copy_number:start_residue:end_residue
    '''
    mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(particle))
    copy_number=IMP.atom.get_copy_index(IMP.atom.Hierarchy(particle))

    if IMP.atom.Fragment.get_is_setup(particle): # CG bead
        residues_in_bead = IMP.atom.Fragment(particle).get_residue_indexes()
        bead_name = mol_name+"."+str(copy_number)+":"+str(min(residues_in_bead))+"-"+str(max(residues_in_bead))
    else:
        residue_in_bead = str(IMP.atom.Residue(particle).get_index())
        bead_name = mol_name+"."+str(copy_number)+":"+residue_in_bead+"-"+residue_in_bead

    return bead_name



def get_chainid_modelid(protein,residue_id,structured_regions,model_ids,link,p1=False):
    domains = structured_regions[protein]
    for domain in domains:
        if residue_id in domain[-1]:
            if p1==True:
                if len(link.split(',')[0].split('.')) <= 1:
                    chainid = domain[0]
                else:
                    chainid = domain[int(link.split(',')[0].split('.')[1])]
            else:
                if len(link.split(',')[2].split('.')) <= 1:
                    chainid =domain[0]
                else:
                    chainid = domain[int(link.split(',')[2].split('.')[1])]

            modelid = model_ids[chainid]
    return modelid,chainid


############################################################################################################################
########################################################### Main ###########################################################
############################################################################################################################

args = parse_args()
violated_xls = get_violated_xls(args.xlviol)
# a = glob.glob(args.cif_dir+'/chainname_changed_aligned_offset_corrected_*.cif')
# for a1 in a:
#     print(a1)

# ------------------------------------------------------------------------------------
### 1. Get minimum distance crosslinks:
# ------------------------------------------------------------------------------------
ccm_model_index = get_models_from_file(args.mdls)[0]
print(f'Model index of cluster center model: {ccm_model_index}')

ccm_model_runid, ccm_model_replicaid, ccm_model_frameid = get_model_identity_from_files(args.csv_file_A, args.csv_file_B).get(ccm_model_index)
print(f'For the cluster center model: \nrunID = {ccm_model_runid} \treplicaID = {ccm_model_replicaid} \tframeID = {ccm_model_frameid}\n')

print('Reading the stat file')
stat_file_line_command = subprocess.Popen(["/home/shreyasarvindekar/imp-clean/build/setup_environment.sh","python3","/home/shreyasarvindekar/imp-clean/imp/modules/pmi/pyext/src/process_output.py",\
                                            "-f",args.run_dirs+'run_'+str(ccm_model_runid)+"/stat."+str(ccm_model_replicaid)+".out","-n",str(ccm_model_frameid)],\
                                            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
# For the stat file line command, remember to remove the output dir when running for subcomplexes other than MHR.
frame_out,frame_err = stat_file_line_command.communicate()
if len(str(frame_err))>3:
    print('Error in reading the stat file')
    exit(1)
stat_lines = str(frame_out)[2:-1].split('\\n')

min_dist_xl = get_min_dist(stat_lines,args.xltype)
if len(min_dist_xl)==0:
    print("Warning: No crosslinks of the specified type were found")

# all_particles = IMP.atom.Selection(hierarchy=hier).get_selected_particles()
particle_details,cg_beads = strip_multiscaling(args.input)
# for k in particle_details:
#     print(k,particle_details[k])

mdl = IMP.Model()
ccm = RMF.open_rmf_file_read_only('multiscaling_striped_cluster_center_model.rmf3')
hier = IMP.rmf.create_hierarchies(ccm, mdl)[0]
IMP.rmf.load_frame(ccm, 0)
mdl.update()


# all_particles_indices = IMP.atom.Selection(hierarchy=hier).get_selected_particle_indexes()



# ------------------------------------------------------------------------------------
### 2. Generate pseudobonds:
# ------------------------------------------------------------------------------------

pseudobonds = []
for xl in min_dist_xl:
    if xl.split(',')[0][-2] != '.':
        particle1 = f"{xl.split(',')[0]}.0:{xl.split(',')[1]}"
    else:
        particle1 = ':'.join(xl.split(',')[0:2])
    if xl.split(',')[2][-2] != '.':
        particle2 = f"{xl.split(',')[2]}.0:{xl.split(':')[0].split(',')[3]}"
    else:
        particle2 = ':'.join(xl.split(':')[0].split(',')[2:4])

    index1 = int(particle_details[particle1])
    index2 = int(particle_details[particle2])

    # Decide the color of the link
    crosslink = f"{xl.split(',')[0].split('.')[0]},{xl.split(',')[1]},{xl.split(',')[2].split('.')[0]},{xl.split(',')[3].split(':')[0]}"
    if crosslink in violated_xls:
        color = '#ff0000'
    else:
        color = color = '#adadad'

    # No loop bonds
    if index1 != index2:
        if particle1 in cg_beads.keys() and particle2 in cg_beads.keys():
            psb = f'#0:{index1}@b #0:{index2}@b {color}\n'
            pseudobonds.append(psb)

        elif particle1 in cg_beads.keys() and particle2 not in cg_beads.keys():
            prot = xl.split(',')[2].split('.')[0]
            res = int(xl.split(',')[3].split(':')[0])

            model, chain = get_chainid_modelid(prot,res,structured_regions,model_ids,xl)
            psb = f'#0:{index1}@b #{model}:{res}.{chain}@ca {color}\n'
            pseudobonds.append(psb)

        elif particle1 not in cg_beads and particle2 in cg_beads:
            prot = xl.split(',')[0].split('.')[0]
            res = int(xl.split(',')[1])

            model, chain = get_chainid_modelid(prot,res,structured_regions,model_ids,xl,p1=True)
            psb = f'#{model}:{res}.{chain}@ca #0:{index2}@b {color}\n'
            pseudobonds.append(psb)

        else:
            prot1 = xl.split(',')[0].split('.')[0]
            res1 = int(xl.split(',')[1])
            prot2 = xl.split(',')[2].split('.')[0]
            res2 = int(xl.split(',')[3].split(':')[0])
            model1, chain1 = get_chainid_modelid(prot1,res1,structured_regions,model_ids,xl,p1=True)
            model2, chain2 = get_chainid_modelid(prot2,res2,structured_regions,model_ids,xl)

            psb = f'#{model1}:{res1}.{chain1}@ca #{model2}:{res2}.{chain2}@ca {color}\n'
            pseudobonds.append(psb)

# for ln in pseudobonds:
#     print(ln.strip())


# ------------------------------------------------------------------------------------
### 3. Write the output files:
# ------------------------------------------------------------------------------------

pseudobond_fname = args.xltype + '_run_pseudobond'
with open(pseudobond_fname,'w') as pbf:
    for bond in pseudobonds:
        pbf.write(bond)

chimera_script_fname = 'pseudobond_chimera_script.py'
with open(chimera_script_fname,'w') as csf:
    csf.write(f"import glob\nimport chimera\nfrom chimera import openModels\nfrom chimera import runCommand\n\nmmcifs = {glob.glob(args.cif_dir+'/chainname_changed_aligned_offset_corrected_*.cif')}\n\
runCommand('open '+'multiscaling_striped_cluster_center_model.rmf3')\nfor p in mmcifs:\n\trunCommand('open '+p)")

print('-------------------------------------------------------------------------------------------------------------------------\
    \nPlease note that the total number of crosslinks in the pseudobond file will be less than the input crosslinks.\
    \nThis is because the crosslinks within each coarse-grained bead are removed from the pseudobond file.\
    \n-------------------------------------------------------------------------------------------------------------------------')
print(f"Output files are: \t{pseudobond_fname}\t{chimera_script_fname}")
