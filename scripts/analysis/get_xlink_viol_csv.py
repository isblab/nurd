#!/usr/bin/env python

import os,sys
import glob
import argparse
import subprocess
import string
import numpy
import pandas as pd
from matplotlib import pyplot as plt



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


class Violations(object):

    def __init__(self, threshold):

        self.violation_threshold  = threshold
        self.violation_counts = {}  # dictionary with a key per restraint storing number of models it is violated in
        self.global_minimum_xlink_distances = {} # minimun distance across all models

    def get_number_violated_restraints(self, frame_out, rst_keyword):
        num_violated = 0
        stat_lines = frame_out.strip().split("\n")
        minimum_xlink_distances={}
        for ln in stat_lines:
            if not ln.startswith(rst_keyword):
                continue

            [rst,value] = ln.strip().split()

            # need to handle ambiguous restraints
            items = rst.split('|')
            (prot1,pos1,prot2,pos2) = items[3:7]
            if '.' in prot1:
                protonly1 = prot1.split('.')[0]  # separate the protein name from copy number
            else:
                protonly1 = prot1
            if '.' in prot2:
                protonly2 = prot2.split('.')[0]   # separate the protein name from copy number
            else:
                protonly2=prot2
            if  (protonly1,pos1,protonly2,pos2) in minimum_xlink_distances: # update with current copy's distance
                minimum_xlink_distances[(protonly1,pos1,protonly2,pos2)]=min(minimum_xlink_distances[(protonly1,pos1,protonly2,pos2)],float(value))
            else: # add a new entry
                minimum_xlink_distances[(protonly1,pos1,protonly2,pos2)]=float(value)
            # print(protonly1,pos1,protonly2,pos2,minimum_xlink_distances[(protonly1,pos1,protonly2,pos2)])
        # print(len(minimum_xlink_distances))

        # finally update violation counts and the min_distance from the model
        for xl in minimum_xlink_distances:
            if minimum_xlink_distances[xl] > self.violation_threshold:
                num_violated += 1
                if not xl in self.violation_counts:
                    self.violation_counts[xl] = 1
                else:
                    self.violation_counts[xl] += 1

            if not xl in self.global_minimum_xlink_distances:
                self.global_minimum_xlink_distances[xl] = minimum_xlink_distances[xl]
            else:
                self.global_minimum_xlink_distances[xl] = min(minimum_xlink_distances[xl], self.global_minimum_xlink_distances[xl])
            # print(self.global_minimum_xlink_distances[xl])
        # print(self.violation_counts)
        return num_violated, len(minimum_xlink_distances)




parser = argparse.ArgumentParser(description='Get the number of crosslinks violated across all the models. This script will also return the global minimum distances for all crosslinks.')
parser.add_argument('-c', type=int, help='Cluster number (from sampcon)')
parser.add_argument('-m', type=str, help='Path to the model_analysis direcory (include "/model_analysis/")')
parser.add_argument('-k', type=str, help='Keyword from stat file. Eg. adh, bs3dss, intra_adh, etc.')
parser.add_argument('-t', type=float, help='Distance cutoff for the crosslink')
parser.add_argument('-r', type=str, help='Path to the run directories')
args = parser.parse_args()

cluster_models_file = 'cluster.' + str(args.c) + '.all.txt'  #e.g. cluster.0.all.txt
out_fname = 'cl' + str(args.c) + '_' + args.k + '.txt'
viol_out_fname = 'cl' + str(args.c) + '_xl_viol' + args.k + '.txt'
csv_file_path = str(args.m + 'good_scoring_models_*')
run_dirs = str(args.r) + 'run_'

csv_files = glob.glob(csv_file_path)
print('CSV files path:', csv_file_path)
print(csv_files)
print('Run directories path:', run_dirs)
# print(csv_files)

if 'A' in csv_files[0]:
    csv_file_A = csv_files[0] # e.g. good_scoring_models_A_cluster2_detailed.csv
    csv_file_B = csv_files[1]  # e.g. good_scoring_models_B_cluster2_detailed.csv
else:
    csv_file_A = csv_files[1] # e.g. good_scoring_models_A_cluster2_detailed.csv
    csv_file_B = csv_files[0]  # e.g. good_scoring_models_B_cluster2_detailed.csv

# print(csv_file_A,'\n',csv_file_B)
# print(cluster_models_file)


model_indices = get_models_from_file(cluster_models_file)
model_ids = get_model_identity_from_files(csv_file_A,csv_file_B)
# print(len(model_ids),len(model_indices))


keyword = 'CrossLinkingMassSpectrometryRestraint_Distance_|' + str(args.k)
Analysis = Violations(args.t)

# if xltype == "adh":
#     keyword = "CrossLinkingMassSpectrometryRestraint_Distance_|adh"
#     Analysis = Violations(35.0)
# elif xltype == "bs3dss":
#     keyword = "CrossLinkingMassSpectrometryRestraint_Distance_|bs3dss"
#     Analysis = Violations(35.0)
# elif xltype == "dmtmm":
#     keyword = "CrossLinkingMassSpectrometryRestraint_Distance_|dmtmm"
#     Analysis = Violations(25.0)

mdl_num = 0
# percent_violated_in_each_model = []
number_of_xls_violated_in_the_model = []

for mdl in model_indices:
    mdl_left = len(model_indices) - mdl_num
    mdl_num += 1
    # print(f'{mdl_left} models left')
    (run,rep,frame) = model_ids[mdl]

    stat_file_line_command = subprocess.Popen(["/home/shreyasarvindekar/imp-clean/build/setup_environment.sh","python","/home/shreyasarvindekar/imp-clean/imp/modules/pmi/pyext/src/process_output.py","-f",run_dirs+str(run)+"/stat."+str(rep)+".out","-n",str(frame+1)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    frame_out,frame_err = stat_file_line_command.communicate()

    # print(frame_out)
    # print(str(frame_err, 'utf-8'))
    num_violated_in_model, total_xls = Analysis.get_number_violated_restraints(str(frame_out, 'utf-8'),keyword)
    # percent = num_violated_in_model*100/total_xls
    # percent_violated_in_each_model.append(percent)
    number_of_xls_violated_in_the_model.append(num_violated_in_model)
    print(f'{mdl_left} models left \t\t Number of crosslinks violated in this model: {num_violated_in_model}\t\t Total Xls this model: {total_xls}')
    # print(Analysis.global_minimum_xlink_distances)


# print(Analysis.violation_counts)
num_violated_in_all_models = 0
print('Number of models is: ',len(model_indices))

for rst in Analysis.violation_counts:
    # print(f'Violated in {Analysis.violation_counts[rst]} models')
    if Analysis.violation_counts[rst] == len(model_indices): # violated in all models
        num_violated_in_all_models+=1
        # print ("Violated ",rst)
    else:
        # print('Not Violated')
        continue

num_violated_in_all_models = str(num_violated_in_all_models)
print("Number of crosslinks violated in all models",num_violated_in_all_models)
with open(viol_out_fname,'w') as xl_viol_file:
    xl_viol_file.write(f"For xltype = {args.k} \nNumber of crosslinks violated in all models: {num_violated_in_all_models} \n")
#
#
# with open(out_fname,'w') as xl_out:
#     for rst in Analysis.global_minimum_xlink_distances:
        # xl_out.write(str(rst[0])+','+str(rst[1])+','+str(rst[2])+','+str(rst[3])+','+str(Analysis.global_minimum_xlink_distances[rst])+'\n')
        # print(str(rst[0])+','+str(rst[1])+','+str(rst[2])+','+str(rst[3])+','+str(Analysis.global_minimum_xlink_distances[rst]))



with open(f'xlviol_percent_lst_{args.k}','w') as tf:
    for entry in number_of_xls_violated_in_the_model:
        tf.write(str(entry)+'\n')
    tf.write('>'+str(total_xls))

# plt.figure()
# plt.hist(percent_violated_in_each_model, bins=100, range=[0,50], color='#0095FF',histtype='step')
# plt.title(f'Percent crosslinks violations for {args.k}')
# plt.xlabel('Percentage')
# plt.ylabel('Number of models')
# plt.savefig(f'percent_viol_{args.k}.png')
