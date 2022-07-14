[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6674232.svg)](https://doi.org/10.5281/zenodo.6674232)

# Integrative model of the NuRD subcomplexes

This repository is of the integrative model of the NuRD subcomplexes based on data from negative stain EM, chemical crosslinking, X-ray crystallography, DIA-MS, SEC-MALLS and [COSMIC](https://cancer.sanger.ac.uk/cosmic) (*Cancer Mutations Database*). It contains input data, scripts for modeling and results including bead models and localization probability density maps. The modeling was performed using [IMP](https://integrativemodeling.org) (*Integrative Modeling Platform*).

![four_stage_figure](https://user-images.githubusercontent.com/8314735/165240737-c960c153-6af9-4014-9d76-32bc8c01cd03.png)


## Directory structure
1. [inputs](inputs/) : contains the subdirectories for the input data used for modeling all the subcomplexes.
2. [scripts](scripts/) : contains all the scripts used for modeling and analysis of the models.
3. [results](results/) : contains the models and the localization probability densities of the top cluster of the subcomplexes .
4. [test](test/) : scripts for testing the sampling


### Simulations
These are the independent simulations:
1. Modeling of MHR subcomplex : `mhr`
2. Modeling of MHM subcomplex : `mhm`
3. Modeling of NuDe subcomplex : `nude`
4. Modeling of MHR without using the EM data : `mhr_xl_ctrl`

## Protocol
### Sampling
To run the sampling, run modeling scripts like this \
```
for runid in `seq 1 NRUNS` ; do mpirun -np NCORES $IMP python scripts/sample/SUBCOMPLEXNAME_modeling.py prod $runid ; done
```

where, \
`$IMP` is the setup script corresponding to the IMP installation directory (omit for binary installation), \
`SUBCOMPLEXNAME` is `mhr`, `mhm`, `nude`, `mhr_ctrl`   \
`NRUNS` is the number of runs, \
and `NCORES` is the number of cores on which replica exchange is to be carried out.

1. For MHR: ```SUBCOMPLEXNAME = mhr```, ```NCORES = 8``` and ```NRUNS = 50 ```
2. For MHM: ```SUBCOMPLEXNAME = mhm```, ```NCORES = 8``` and ```NRUNS = 30 ```
3. For NuDe: ```SUBCOMPLEXNAME = nude```, ```NCORES = 8``` and ```NRUNS = 50 ```
4. For MHR without using the EM data: ```SUBCOMPLEXNAME = mhr_xl_ctrl```, ```NCORES = 8``` and ```NRUNS = 30 ```


### Analysis
#### 1. Getting the good scoring models
  Good scoring models were selected using `pmi_analysis` (Please refer to [pmi_analysis tutorial](https://github.com/salilab/PMI_analysis) for more detailed explaination) along with our `variable_filter_v1.py` script. These scripts are run as described below:
  1. First, run `run_analysis_trajectories.py` as follows:\
      `$IMP run_analysis_trajectories.py modeling run_ `\
      where, `$IMP` is the setup script corresponding to the IMP installation directory (omit for binary installation), \
      `modeling` is the directory containing all the runs and \
      `run_` is the prefix for the names of individual run directories.\
      Alternatively, one can also run the `submit_run_analysis_trajectories.sh` script from the `scripts/analysis/pmi_analysis` directory

  2. Then run `variable_filter_v1.py` on the major cluster obtained as follows: \
      `$IMP variable_filter_v1.py -c N -g MODEL_ANALYSIS_DIR`
      where, `$IMP` is the setup script corresponding to the IMP installation directory (omit for binary installation), \
      `N` is the cluster number of the major cluster, \
      `MODEL_ANALYSIS_DIR` is the location of the directory containing the selected_models*.csv. \
      This can also be run using the `submit_variable_filter_v1.sh` script from the `scripts/analysis/pmi_analysis` directory.\
  _Please also refer to the comments in the `variable_filter_v1.py` for more details._

  3. The selected good scoring models were then extracted using `run_extract_good_scoring_models.py` as follows: \
      `$IMP python run_extract_good_scoring_models.py modeling run_ CLUSTER_NUMBER` \
      where, `$IMP` is the setup script corresponding to the IMP installation directory (omit for binary installation), \
      `modeling` is the path to the directory containing all the individual runs and \
      `CLUSTER_NUMBER` is the number of the major cluster to be extracted.\
      This can also be run using the script `submit_run_extract_models.sh` from the `scripts/analysis/pmi_analysis` directory.

#### 2. Running the sampling exhaustiveness tests (Sampcon)
A separate directory named `sampcon` was created and a `density.txt` file was added to it. This file contains the details of the domains to be split for plotting the localisation probability densities. Finally, sampling exhaustiveness tests were performed using `imp-sampcon` as shown in `scripts/analysis/pmi_analysis/*_sampcon.sh`. \
where, `*` is the name of the complex.

#### 3. Analysing the major cluster
  1. Crosslink violations were analyzed as follows: \
      `for xltype in adh bs3dss; do python get_xlink_viol_csv.py -c CLUSTER_NUMBER -m MODELANALYSIS_DIR -r modeling -k $xltype -t 35.0 & done` \
      and \
      `python get_xlink_viol_csv.py -c CLUSTER_NUMBER -m MODELANALYSIS_DIR -r modeling -k dmtmm -t 25.0` \
      One acn also use the `get_xl_viol_validation_set.py` script from the `scripts/analysis/xlviol` directory after changing the inputs section in the script as follows: \
      `$imp python get_xl_viol_validation_set.py -ia ../cluster.0.sample_A.txt -ib ../cluster.0.sample_B.txt -ra ../../MODEL_ANALYSIS_DIR/A_gsm_clust0.rmf3 -rb ../../MODEL_ANALYSIS_DIR/B_gsm_clust0.rmf3 -ta -ra ../../MODEL_ANALYSIS_DIR/A_gsm_clust0.txt -c ../cluster.0/cluster_center_model.rmf3 -x XL_FILE -t THRESHOLD` \
      where, `XL_FILE` is a file containing the crosslinks to be analysed.

  2. The above scripts generate files mentioning the minimum distance for each crosslink. These files were then passed to `xl_distance_hist_plotter.py` as follows: \
     `python xl_distance_hist_plotter.py FILE_NAME XL_NAME THRESHOLD` \
     where, `FILE_NAME` is the name of the file,\
     `XL_NAME` is the name of the linker used, and \
     `THRESHOLD` is the distance threshold for that linker. \
     This script will generate a histogram of the minimum distances spanned by the crosslinks.

  3. Then, the files obtained from scripts in point 1 were passed to `binner_cx-circos.py`as follows: \
     `python binner_cx-circos.py FILE_NAME` \
     where, `FILE_NAME` is the name of the file. \
     This script generates a binned version of the input file which can then be used to make the crosslink plots using [CIRCOS](http://cx-circos.net/).

  4. Contact maps were plotted for the NuDe models as follows:
      `scripts/analysis/cosmic_and_distance-maps/submit_contact_maps_all_pairs_surface.py` \
      This script calls the `scripts/analysis/cosmic_and_distance-maps/contact_maps_all_pairs_surface.py` script.
      _Please use `--help` for `contact_maps_all_pairs_surface.py` script for more details._

  3. Finally, COSMIC cancer mutations were annotated on the models as follows: \
      `python color_mutations/color_mutation.py -i cluster.0/cluster_center_model.rmf3 -r 10 -mf mutations.txt`


### Results

For each of the simulations, the following files are in the [results](results/) directory
* `cluster_center_model.rmf3` : representative bead model of the major cluster
* `chimera_densities.py` : to view the localization densities (.mrc files)
* `xl_violations.txt` : list of violated crosslinks

For the NuDe models, `mutation_colored_model.rmf` and `Distance_Maps` are also added.

### Information
**Author(s):** Shreyas Arvindekar, Shruthi Viswanath\
**Date**: May 19th, 2022\
**License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
International License.\
**Last known good IMP version:** [![build info](https://integrativemodeling.org/systems/35/badge.svg?branch=main)](https://integrativemodeling.org/systems/) \
**Testable:** Yes\
**Parallelizeable:** Yes\
**Publications:** Submitted for publication ([bioRxiv](https://doi.org/10.1101/2021.11.25.469965))
