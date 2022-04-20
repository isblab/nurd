# Make sure to pass offset corrected structure files to this script 

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


structured_regions = {'MTA1':[('bah1','bah2',range(1,165)),('elmsant1','elmsant2',range(165,333)),('h1','h2',range(334,353)),('zf1','zf2',range(389,431)),\
                            ('mta1r1','mta2r1',range(468,546)),('mta1r2','mta2r2',range(670,691))],\
                    'HDAC1':[('hdac1','hdac2',range(8,376))],\
                    'RBBP4':[('rbbp1','rbbp2','rbbp3','rbbp4',range(1,411))]}

model_ids = {'elmsant1':1,'elmsant2':1,'hdac1':1,'hdac2':1,\
             'bah1':2,\
             'bah2':3,\
             'h1':4,\
             'h2':5,\
             'zf1':6,\
             'zf2':7,\
             'mta1r1':8,'rbbp1':8,\
             'mta2r1':9,'rbbp2':9,\
             'mta1r2':10,'rbbp3':10,\
             'mta2r2':11,'rbbp4':11}                                         # Chain name: model_id



def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', dest="input", help='Representative model in RMF format', default="cluster.0/cluster_center_model.rmf3")
    parser.add_argument('-ca', dest="csv_file_A", help='csvA file')
    parser.add_argument('-cb', dest="csv_file_B", help='csvB file')
    parser.add_argument('-m', dest='mdls', help='Path to models (Something like: "/sampcon/cluster.0.all.txt")')
    parser.add_argument('-r', dest='run_dirs', type=str, help='Path to the run directories')
    parser.add_argument('-v', dest='xlviol', type=str, help='Path to the file containing violated crosslinks')
    parser.add_argument('-t', dest='threshold', type=float, help='Distance threshold for the crosslinks')
    parser.add_argument('--type', dest='xltype', type=str, help='XL type. Eg. adh, dmtmm, etc. Enter it exactly as in the stat files')
    parser.add_argument('--mmcif', dest='cif_dir', type=str, help='Path to the structure .cif files')
    
    
    # parser.add_argument('-x', dest="xl_file", help='Path to the XLs file')

    return parser.parse_args()



# def get_selected_particles(m,cluster_center_model):
#     selection = None
#     inf = RMF.open_rmf_file_read_only(cluster_center_model)
#     hier = IMP.rmf.create_hierarchies(inf, m)[0]

#     IMP.rmf.load_frame(inf, 0)
#     m.update()
#     selection = IMP.atom.Selection(hier)											   # I believe, if I dont specify the resolution, it will take the default resolution of the bead.
#     # , resolution=resolution)
#     del inf
#     return selection



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
            dist = ln.split(' ')[-1]

            k = p1+','+r1+','+p2+','+r2
            if k not in xl_min_dist.keys():
                xl_min_dist[k] = line
            else:
                if dist < xl_min_dist[k]:
                    xl_min_dist[k] = line
                else:
                    continue
    return xl_min_dist



def get_xls_split(xl_min_dist):
    unstruc_xls = []
    struc_xls = []
    semi_struc_xls = []

    for key in xl_min_dist.keys():
        prot1,res1,prot2,res2 = key.split(',')[0], key.split(',')[1], key.split(',')[2], key.split(',')[3]

        
        xl = xl_min_dist[key].split(' : ')[0]
        bead1 = ':'.join(xl.split(',')[0:2])
        bead2 = ':'.join(xl.split(',')[2:4])
        if bead1 in struc_beads and bead2 in struc_beads:
            struc_xls.append(xl_min_dist[key].split(' : ')[0])
            
        elif bead1 in struc_beads and not bead2 in struc_beads:
            semi_struc_xls.append(xl_min_dist[key].split(' : ')[0])

        elif not bead1 in struc_beads and bead2 in struc_beads:
            semi_struc_xls.append(xl_min_dist[key].split(' : ')[0])
        
        elif not bead1 in struc_beads and not bead2 in struc_beads:
            unstruc_xls.append(xl_min_dist[key].split(' : ')[0])


    return unstruc_xls, struc_xls, semi_struc_xls



def get_chainid(struc_xl_list,semi_struc_xl_list):
    out_xl = {}
    # out_ ={}
    
    for xllist in (struc_xl_list,semi_struc_xl_list):
        for xl in xllist:
            prot1 = xl.split(',')[0].split('.')[0]
            cp1 = int(xl.split(',')[0].split('.')[1])
            r1 = int(xl.split(',')[1])
            
            prot2 = xl.split(',')[2].split('.')[0]
            cp2 = int(xl.split(',')[2].split('.')[1])
            r2 = int(xl.split(',')[3])
            
            chain1, chain2 = 'Uns', 'Uns'

            for prot in structured_regions.keys():
                if prot1==prot:
                    domains = structured_regions[prot1]
                    for dom in domains:
                        if r1 in dom[-1]:
                            chain1 = dom[cp1]

            for prot in structured_regions.keys():
                if prot2==prot:
                    domains = structured_regions[prot2]
                    for dom in domains:
                        if r2 in dom[-1]:
                            chain2 = dom[cp2]
            
            out_xl[xl] = f"{chain1},{chain2}"

    return out_xl
################################################################################
############################# Main #############################################
################################################################################

args = parse_args()

violated_xls = []
# with open(args.xlviol,'r') as xlvf:
#     for ln in xlvf.readlines():
#         min_dist = ln.strip().split(',')[-1]
#         if float(min_dist) > args.threshold:
#             violated_xls.append(ln)


mdl = IMP.Model()
ccm = RMF.open_rmf_file_read_only(args.input)
hier = IMP.rmf.create_hierarchies(ccm, mdl)[0]
IMP.rmf.load_frame(ccm, 0)
mdl.update()
particles = IMP.atom.Selection(hier, resolution=1).get_selected_particles()
# print(particles)
particle_details = {}
struc_beads = []
unstruc_beads = []
for i, leaf in enumerate(particles):
    name = get_bead_name(leaf)
    min_res = int(name.split(':')[1].split('-')[0])
    max_res = int(name.split(':')[1].split('-')[1])
    if min_res == max_res:
        new_name = name.split('-')[0]
        particle_details[new_name] = i
        struc_beads.append(new_name)
    else:
        while min_res <= max_res:
            new_name = name.split(':')[0]+':'+str(min_res)
            particle_details[new_name] = i
            min_res+=1
            unstruc_beads.append(new_name)


for i in particle_details.keys():
    print(particle_details[i],i)


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
xl_min_dist = get_min_dist_pair_for_xls(stat_lines)
# print(xl_min_dist)
unstruc_xls, struc_xls, semi_struc_xls = get_xls_split(xl_min_dist)
not_uns_xl_mapped = get_chainid(struc_xls,semi_struc_xls)


pseudobonds_list = []
for link in unstruc_xls:
    color = '#87ceeb'
    particle1 = link.split(',')[0]+':'+link.split(',')[1]
    particle2 = link.split(',')[2]+':'+link.split(',')[3]

    p1 = particle1.split(':')[0].split('.')[0]
    p2 = particle2.split(':')[0].split('.')[0]
    r1 = particle1.split(':')[1]
    r2 = particle1.split(':')[1]
    lnk = f'{p1},{r1},{p2},{r2}'
    if lnk in violated_xls:
        color = '#ff0000'
    else:
        color = '#87ceeb'

    pseudobond_ln = f'#0.2:{str(particle_details[particle1])} #0.2:{str(particle_details[particle2])} {color}\n'
    pseudobonds_list.append(pseudobond_ln)

for link in not_uns_xl_mapped.keys():
    p1_id = link.split(',')[0].split('.')[0]
    cp1_id = int(link.split(',')[0].split('.')[1])
    res1id = link.split(',')[1]
    p2_id = link.split(',')[2].split('.')[0]
    cp2_id = int(link.split(',')[2].split('.')[1])
    res2id = link.split(',')[3]

    chain1id = not_uns_xl_mapped[link].split(',')[0]
    chain2id = not_uns_xl_mapped[link].split(',')[1]
    
    lnk = f'{p1_id},{res1id},{p2_id},{res2id}'
    if lnk in violated_xls:
        color = '#ff0000'
    else:
        color = '#87ceeb'
    
    pseudobond_ln = f'#{model_ids[chain1id]}:{res1id}.{chain1id}@ca #{model_ids[chain2id]}:{res2id}.{chain2id}@ca {color}\n'
    pseudobonds_list.append(pseudobond_ln)


pseudobond_fname = args.xltype + '_run_pseudobond'
with open(pseudobond_fname,'w') as pbf:
    for bond in pseudobonds_list:
        pbf.write(bond)

chimera_script_fname = '_'.join(pseudobond_fname.split('_')[0:2])+'_chimera_script.py'
with open(chimera_script_fname,'w') as csf:
    csf.write(f"import glob\nimport chimera\nfrom chimera import openModels\nfrom chimera import runCommand\n\nmmcifs = {glob.glob(args.cif_dir+'/chainname-changed_*.cif')}\nrunCommand('open '+'{args.input}')\nfor p in mmcifs:\n\trunCommand('open '+p)")



# pseudobond = #modelid:residue.chainid@ca #modelid:residue.chainid@ca


# structured_regions = {'MTA1':{'BAH':(1,164),'ELM2SANT':(165,333),'H':(334,353),'ZF':(389,431),'R1':(468,546),'R2':(670,691)},\
#                             'HDAC1':{'HDAC1':(8,376)},\
#                             'RBBP4':{'RBBP4':(1,411)}}

# chain_names = {'BAH':('bah1','bah2'),\
#             'ELM2SANT':('elmsant1','elmsant2'),\
#             'H':('h1','h2'),\
#             'ZF':('zf1','zf2'),\
#             'R1':('mta1r1','mta2r1'),\
#             'R2':('mta1r2','mta2r2'),\
#             'HDAC1':('hdac1','hdac2'),\
#             'RBBP4':('rbbp1','rbbp2','rbbp3','rbbp4')}                       # Domain name: (chain name 1, chain name2,...)

















