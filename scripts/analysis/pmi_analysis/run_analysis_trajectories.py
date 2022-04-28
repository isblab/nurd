import numpy as np
import pandas as pd
import math
import glob
import sys
import os

# Change this to the location of your PMI_analysis folder
sys.path.append('/home/shreyas/imp-clean/pmi_analysis/pyext/src')
from analysis_trajectories import *

#################################
########### MAIN ################
#################################

nproc = 120
top_dir =  sys.argv[1] 
analys_dir = os.getcwd()+'/model_analysis/'

# Check if analysis dir exists
if not os.path.isdir(analys_dir):
    os.makedirs(analys_dir)

# How are the trajectories dir names
dir_head = sys.argv[2]
out_dirs = glob.glob(top_dir+'/'+dir_head+'*/')
print(out_dirs)
################################
# Get and organize fields for
# analysis
################################
# Read the total score, plot
# and check for score convengence
XLs_cutoffs = {'adh':35.0, 'bs3dss':35.0, 'dmtmm':25.0}

# Load module
AT = AnalysisTrajectories(out_dirs,
                          dir_name=dir_head,
                          analysis_dir = analys_dir,
                          nproc=nproc)

# Define restraints to analyze
AT.set_analyze_XLs_restraint(XLs_cutoffs = XLs_cutoffs,
				ambiguous_XLs_restraint = True)
AT.set_analyze_Connectivity_restraint()
AT.set_analyze_Excluded_volume_restraint()
AT.set_analyze_EM_restraint()
AT.restraint_names['GaussianEMRestraint']='GaussianEMRestraint'

# Read stat files
AT.read_stat_files()
AT.write_models_info()
AT.get_psi_stats()

# What scores do we cluster on?
AT.hdbscan_clustering(['EV_sum', 'XLs_sum', 'GaussianEMRestraint_None'])
AT.summarize_XLs_info(ambiguous_XLs_restraint = True)
exit()


