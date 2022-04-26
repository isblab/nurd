import os,sys

proteins = ["MTA1.0", "MTA1.1", "HDAC1.0", "HDAC1.1", "MBD3.0", "RBBP4.0", "RBBP4.1", "RBBP4.2", "RBBP4.3"]
prots = []

for index1 in range(len(proteins)):
	if index1 != len(proteins)-1:
		for index2 in range(index1+1, len(proteins)):
			prot1 = proteins[index1]
			prot2 = proteins[index2]
			prots.append(f"{prot1},{prot2}")
print(prots)

# os.system(f"for protein in {prots}; do echo $protein ; done")

i=0
while i < len(prots):
	pair = prots[i]
	os.system(f"~/imp-clean/build/setup_environment.sh python contact_maps_all_pairs_surface.py \
		-ia ../cluster.0.sample_A.txt -ib ../cluster.0.sample_B.txt -ra ../../model_analysis/A_gsm_clust3.rmf3 -rb ../../model_analysis/B_gsm_clust3.rmf3 \
		-ta ../../model_analysis/A_gsm_clust3.txt -c ../cluster.0/cluster_center_model.rmf3 -p {pair}  &")
	i+=1
