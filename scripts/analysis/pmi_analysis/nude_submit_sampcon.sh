#! /bin/bash

~/imp-clean/build/setup_environment.sh  python  ~/imp-clean/imp/modules/sampcon/pyext/src/exhaust.py -n nude -a -m cpu_omp -c 7 -d ../density.txt -gp -g 2.0  -sa ../model_analysis/A_gsm_clust3.txt -sb ../model_analysis/B_gsm_clust3.txt  -ra ../model_analysis/A_gsm_clust3.rmf3 -rb ../model_analysis/B_gsm_clust3.rmf3

echo "Done!"
