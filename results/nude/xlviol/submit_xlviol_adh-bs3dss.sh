#! /bin/bash

for xlfile in ../../../input/xlms/new_xls/filtered_adh.dat ../../../input/xlms/new_xls/filtered_bs3dss.dat ; do $imp python /home/shreyas/Projects/from_github/scripts/new_scripts/get_xl_viol_validation_set.py -ia ../cluster.0.sample_A.txt -ib ../cluster.0.sample_B.txt -ra ../../model_analysis/A_gsm_clust3.rmf3 -rb ../../model_analysis/B_gsm_clust3.rmf3 -ta ../../model_analysis/A_gsm_clust3.txt -c ../cluster.0/cluster_center_model.rmf3 -x $xlfile -t 35.0 & done
