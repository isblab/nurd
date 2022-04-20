import os,sys

proteins = ["MTA1.0", "MTA1.1", "HDAC1.0", "HDAC1.1", "MBD3.0", "RBBP4.0", "RBBP4.1", "RBBP4.2", "RBBP4.3"]
prots = []

for index1 in range(len(proteins)):
	for index2 in range(index1+1, len(proteins)):
		prot1 = proteins[index1]
		prot2 = proteins[index2]
		prots.append(f"{prot1},{prot2}")

os.system(f"for protein in {prots}; do ~/imp-clean/build/setup_environment.sh python contact_maps_all_pairs_surface.py \
	-ia ../cluster.0.sample_A.txt -ib ../cluster.0.sample_B.txt -ra ../../model_analysis/A_gsm_clust1.rmf3 -rb ../../model_analysis/B_gsm_clust1.rmf3 \
	-ta ../../model_analysis/A_gsm_clust1.txt -c ../cluster.0/cluster_center_model.rmf3 -p $protein & done")



