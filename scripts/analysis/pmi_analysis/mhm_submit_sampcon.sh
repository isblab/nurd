#! /bin/bash

$imp  python  ~/imp-clean/imp/modules/sampcon/pyext/src/exhaust.py -n mhm -a -m cpu_omp -c 3 -d ../density.txt -gp -g 2.0  -sa ../model_analysis/A_gsm_clust1.txt -sb ../model_analysis/B_gsm_clust1.txt  -ra ../model_analysis/A_gsm_clust1.rmf3 -rb ../model_analysis/B_gsm_clust1.rmf3 

echo "Done!"
