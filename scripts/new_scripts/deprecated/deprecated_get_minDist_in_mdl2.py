'''
1. Get min XLs
2. Select particle. Check if structured . imp.fragment.is setup
'''
import os, sys,glob
import pandas as pd
import argparse
import subprocess
import IMP
import IMP.core
import IMP.pmi.tools
import RMF
import IMP.rmf

__doc__ = "Get the minimum distance pair for each crosslink"


structured_regions = {'MTA1':[('bah1','bah2',range(1,165)),('elmsant1','elmsant2',range(165,334)),('h1','h2',range(334,354)),('zf1','zf2',range(389,432)),\
                            ('mta1r1','mta2r1',range(468,547)),('mta1r2','mta2r2',range(670,692))],\
                    'HDAC1':[('hdac1','hdac2',range(8,377))],\
                    'RBBP4':[('rbbp1','rbbp2','rbbp3','rbbp4',range(1,412))]}

model_ids = {'elmsant1':1,'elmsant2':1,'hdac1':1,'hdac2':1,\
             'bah1':2,\
             'bah2':5,\
             'h1':3,\
             'h2':6,\
             'zf1':4,\
             'zf2':7,\
             'mta1r1':8,'r1_rbbp0':8,\
             'mta2r1':9,'r1_rbbp1':9,\
             'mta1r2':10,'r2_rbbp2':10,\
             'mta2r2':11,'r2_rbbp3':11}                                         # Chain name: model_id



def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', dest="input", help='Representative model in RMF format', default="cluster.0/cluster_center_model.rmf3")
    parser.add_argument('-ca', dest="csv_file_A", help='csvA file')
    parser.add_argument('-cb', dest="csv_file_B", help='csvB file')
    parser.add_argument('-m', dest='mdls', help='Path to models (Something like: "/sampcon/cluster.0.all.txt")')
    parser.add_argument('-r', dest='run_dirs', type=str, help='Path to the run directories')
    parser.add_argument('--type', dest='xltype', type=str, help='XL type. Eg. adh, dmtmm, etc. Enter it exactly as in the stat files')
    parser.add_argument('-v', dest='xlviol', type=str, help='Path to the file containing violated crosslinks')
    parser.add_argument('-t', dest='threshold', type=float, help='Distance threshold for the crosslinks')
    parser.add_argument('--mmcif', dest='cif_dir', type=str, help='Path to the structure .cif files')
    # parser.add_argument('-x', dest="xl_file", help='Path to the XLs file')

    return parser.parse_args()


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


def get_bead_name(particle):
    '''
    Input: particle
    Output: bead name in the format molecule_name:copy_number:start_residue:end_residue
    '''
    mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(particle))
    copy_number=IMP.atom.get_copy_index(IMP.atom.Hierarchy(particle))

    if IMP.atom.Fragment.get_is_setup(particle): # CG bead                      ###################### Did not understand get_is_setup
        residues_in_bead = IMP.atom.Fragment(particle).get_residue_indexes()
        bead_name = mol_name+"."+str(copy_number)+":"+str(min(residues_in_bead))+"-"+str(max(residues_in_bead))
    else:
        residue_in_bead = str(IMP.atom.Residue(particle).get_index())
        bead_name = mol_name+"."+str(copy_number)+":"+residue_in_bead+"-"+residue_in_bead

    return bead_name


def get_min_dist_pair_for_xls(stat_lines):
# Takes lines from stat file as input and returns minimum pair xl dictionary
    xl_min_dist = {}                             # Here, name of the XL is the key and minimum distance entry "prot1.copy1,res1,prot2.copy2,res2 : distance" is the value
    for ln in stat_lines:
        if ln.startswith('CrossLinkingMassSpectrometryRestraint_Distance_|' + args.xltype):
            line = ','.join(ln.split('|')[3:7]) +' : '+ str(ln.split(' ')[-1])
            
            p1 = line.split(',')[0].split('.')[0]
            p2 = line.split(',')[2].split('.')[0]
            r1 = line.split(',')[1]
            r2 = line.split(',')[3].split(' : ')[0]
            dist = float(line.split(' ')[-1])

            k = p1+','+r1+','+p2+','+r2
            if k not in xl_min_dist.keys():
                # xl_min_dist[k] = line         # with distance
                xl_min_dist[k] = line.split(' : ')[0] + ':' + str(dist)
                
            else:
                if dist < float(xl_min_dist[k].split(':')[-1]):
                    # xl_min_dist[k] = line     # with distance
                    xl_min_dist[k] = line.split(' : ')[0] + ':' + str(dist)                    
    print(xl_min_dist)            
    return xl_min_dist


############################################################################################################################
########################################################### Main ########################################################### 
############################################################################################################################

args = parse_args()

violated_xls =[]
with open(args.xlviol,'r') as xlvf:
    for ln in xlvf.readlines():
        min_dist = ln.strip().split(',')[-1]
        if float(min_dist) > args.threshold:
            violated_xls.append(','.join(ln.split(',')[0:-1]))


# ------------------------------------------------------------------------------------
### 1. Get minimum distance crosslinks:
# ------------------------------------------------------------------------------------
ccm_model_index = get_models_from_file(args.mdls)[0]
print(f'Model index of cluster center model: {ccm_model_index}')

ccm_model_runid, ccm_model_replicaid, ccm_model_frameid = get_model_identity_from_files(args.csv_file_A, args.csv_file_B).get(ccm_model_index)
print(f'For the cluster center model: \nrunID = {ccm_model_runid} \treplicaID = {ccm_model_replicaid} \tframeID = {ccm_model_frameid}\n')

print('Reading the stat file')
stat_file_line_command = subprocess.Popen(["/home/shreyasarvindekar/imp-clean/build/setup_environment.sh","python3","/home/shreyasarvindekar/imp-clean/imp/modules/pmi/pyext/src/process_output.py",\
                                            "-f",args.run_dirs+'run_'+str(ccm_model_runid)+"/output/stat."+str(ccm_model_replicaid)+".out","-n",str(ccm_model_frameid)],\
                                            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
# For the stat file line command, remember to remove the output dir when running for subcomplexes other than MHR. 
frame_out,frame_err = stat_file_line_command.communicate()

stat_lines = str(frame_out)[2:-1].split('\\n')
xl_min_dist = get_min_dist_pair_for_xls(stat_lines)             # key = prot1,res1,prot2,res2    value = prot1.cp1,res1,prot2.cp2,res2
# for i in xl_min_dist:
#     print(i)


# ------------------------------------------------------------------------------------
### 2. Check if the particles are structured or unstructured
# ------------------------------------------------------------------------------------
mdl = IMP.Model()
ccm = RMF.open_rmf_file_read_only(args.input)
hier = IMP.rmf.create_hierarchies(ccm, mdl)[0]
IMP.rmf.load_frame(ccm, 0)
mdl.update()

struc_xl = []
semistuc_xl = []
unstruc_xl = []

for xl in xl_min_dist.values():
    p1 = xl.split(',')[0].split('.')[0]
    cp1 = xl.split(',')[0].split('.')[1]
    r1 = xl.split(',')[1]
    p2 = xl.split(',')[2].split('.')[0]
    cp2 = xl.split(',')[2].split('.')[1]
    r2 = xl.split(',')[3]

    particle1 = str(IMP.atom.Selection(hierarchy=hier, molecule=p1, copy_index=int(cp1), residue_index=int(r1)).get_selected_particles()[0])[1:-1]
    particle2 = str(IMP.atom.Selection(hierarchy=hier, molecule=p2, copy_index=int(cp2), residue_index=int(r2)).get_selected_particles()[0])[1:-1]
    
    if particle1.endswith('bead') and particle2.endswith('bead'):
        unstruc_xl.append(xl)
    elif particle1.endswith('bead') and not particle2.endswith('bead'):
        semistuc_xl.append(xl)
    elif not particle1.endswith('bead') and particle2.endswith('bead'):
        semistuc_xl.append(xl)
    elif not particle1.endswith('bead') and not particle2.endswith('bead'):
        struc_xl.append(xl)
    else:
        print('Something went wrong while splitting crosslinks into structured, semi-structured and unstructured')
        exit(1)


# ------------------------------------------------------------------------------------
### 3. Make pseudobonds
# ------------------------------------------------------------------------------------
# all_particles = IMP.atom.Selection(hierarchy=hier).get_selected_particles()
# all_particles_indices = IMP.atom.Selection(hierarchy=hier).get_selected_particle_indexes()
# print(all_particles)
all_particles_indices = []

particles = IMP.atom.Selection(hier, resolution=1).get_selected_particles()

particle_names = {}                                 # key = Bead name (prot.copy:residue) and value = particle index

for i, leaf in enumerate(particles):
    name = get_bead_name(leaf)
    min_res = int(name.split(':')[1].split('-')[0])
    max_res = int(name.split(':')[1].split('-')[1])
    if min_res == max_res:
        new_name = name.split('-')[0]
        particle_names[(new_name)] = i
        all_particles_indices.append(i)
        
    else:
        while min_res <= max_res:
            new_name = name.split(':')[0]+':'+str(min_res)
            particle_names[(new_name)] = i
            min_res+=1
            all_particles_indices.append(i)
# print(particle_names)



# for i in range(len(all_particles)):
#     name = get_bead_name(all_particles[i])
#     min_res = int(name.split(':')[1].split('-')[0])
#     max_res = int(name.split(':')[1].split('-')[1])
#     if min_res == max_res:
#         new_name = name.split('-')[0]
#         particle_names[(new_name)] = all_particles_indices[i]
#     else:
#         while min_res <= max_res:
#             new_name = name.split(':')[0]+':'+str(min_res)
#             particle_names[(new_name)] = all_particles_indices[i]
#             min_res+=1
            
# for name in particle_names.keys():
#     print(name, particle_names[name])

pseudobonds_list = []

for link in unstruc_xl:
    p1 = link.split(',')[0].split('.')[0]
    cp1 = link.split(',')[0].split('.')[1]
    p2 = link.split(',')[2].split('.')[0]
    cp2 = link.split(',')[2].split('.')[1]
    r1 = link.split(',')[1]
    r2 = link.split(',')[3]
    lnk = f'{p1},{r1},{p2},{r2}'
    if lnk in violated_xls:
        color = '#ff0000'
    else:
        color = '#87ceeb'

    particle1_name = f"{p1}.{cp1}:{r1}"
    particle2_name = f"{p2}.{cp2}:{r2}"
    particle1_index = particle_names[particle1_name]
    particle2_index = particle_names[particle2_name]
    if particle1_index!=particle2_index:
        pseudobond_ln = f'#0.2:{particle1_index} #0.2:{particle2_index} {color}\n'
        pseudobonds_list.append(pseudobond_ln)
        # print(f'Uns - {pseudobond_ln.strip()}')

for xlt in struc_xl, semistuc_xl:
    for link in xlt:
        p1 = link.split(',')[0].split('.')[0]
        cp1 = link.split(',')[0].split('.')[1]
        p2 = link.split(',')[2].split('.')[0]
        cp2 = link.split(',')[2].split('.')[1]
        r1 = link.split(',')[1]
        r2 = link.split(',')[3]
        lnk = f'{p1},{r1},{p2},{r2}'
        if lnk in violated_xls:
            color = '#ff0000'
        else:
            color = '#87ceeb'

        particle1_name = f"{p1}.{cp1}:{r1}"
        particle2_name = f"{p2}.{cp2}:{r2}"
    
        if IMP.atom.Fragment.get_is_setup(particles[particle_names[particle1_name]]) and not IMP.atom.Fragment.get_is_setup(particles[particle_names[particle2_name]]):            # If first particle is from unstructured regions
            particle1_index = particle_names[particle1_name]

            doms = structured_regions[p2]
            for dom in doms:
                if int(r2) in dom[-1]:
                    chain2id = dom[int(cp2)]
                    model2id = model_ids[chain2id]
                    break
            pseudobond_ln = f'#0.2:{particle1_index} #{model2id}:{r2}.{chain2id}@ca {color}\n'
            pseudobonds_list.append(pseudobond_ln)
            # print(f'Semi -  {pseudobond_ln.strip()}')

        elif not IMP.atom.Fragment.get_is_setup(particles[particle_names[particle1_name]]) and IMP.atom.Fragment.get_is_setup(particles[particle_names[particle2_name]]):          # If second particle is from unstructured regions
            particle2_index = particle_names[particle2_name]

            doms = structured_regions[p1]
            for dom in doms:
                if int(r1) in dom[-1]:
                    chain1id = dom[int(cp1)]
                    model1id = model_ids[chain1id]
                    break

            pseudobond_ln = f'#{model1id}:{r1}.{chain1id}@ca #0.2:{particle2_index} {color}\n'
            pseudobonds_list.append(pseudobond_ln)
            # print(f'Semi -  {pseudobond_ln.strip()}')

        elif not IMP.atom.Fragment.get_is_setup(particles[particle_names[particle1_name]]) and not IMP.atom.Fragment.get_is_setup(particles[particle_names[particle2_name]]):      # If neither of the particles from unstructured regions 
            # print(particle1_name,'\t', particle2_name)
            doms1 = structured_regions[p1]
            for dom in doms1:
                if int(r1) in dom[-1]:
                    chain1id = dom[int(cp1)]
                    model1id = model_ids[chain1id]
                    break
            doms2 = structured_regions[p2]
            # print(doms2)
            flag =0
            for dom in doms2:
                if int(r2) in dom[-1]:
                    chain2id = dom[int(cp2)]
                    # print(chain1id, '\t', chain2id)
                    model2id = model_ids[chain2id]
                    flag=1
                    break
            if flag==0:
                print('bug found')  
  
                # else:
                    
            pseudobond_ln = f'#{model1id}:{r1}.{chain1id}@ca #{model2id}:{r2}.{chain2id}@ca {color}\n'
            pseudobonds_list.append(pseudobond_ln)
            # print(f'Struc - {pseudobond_ln.strip()}\n--------')



pseudobond_fname = args.xltype + '_run_pseudobond'
with open(pseudobond_fname,'w') as pbf:
    for bond in pseudobonds_list:
        pbf.write(bond)

chimera_script_fname = '_'.join(pseudobond_fname.split('_')[0:2])+'_chimera_script.py'
with open(chimera_script_fname,'w') as csf:
    csf.write(f"import glob\nimport chimera\nfrom chimera import openModels\nfrom chimera import runCommand\n\nmmcifs = {glob.glob(args.cif_dir+'/chainname_changed_aligned_offset_corrected_*.cif')}\nrunCommand('open '+'{args.input}')\nfor p in mmcifs:\n\trunCommand('open '+p)")

print('-------------------------------------------------------------------------------------------------------------------------\
    \nPlease note that the total number of crosslinks in the pseudobond file will be less than the input crosslinks.\
    \nThis is because the crosslinks within each coarse-grained bead are removed from the pseudobond file.\
    \n-------------------------------------------------------------------------------------------------------------------------')
print(f"Output files are: \t{pseudobond_fname}\t{chimera_script_fname}")