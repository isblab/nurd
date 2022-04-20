
##########################################################################################
############### IMP Modeling Script for MHR-4RBBP (using ab-initio EM map) ###############
##########################################################################################


# Imports
from __future__ import print_function
import IMP
import RMF
import IMP.rmf
import IMP.pmi
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.topology
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
#import IMP.pmi.restraints.saxs
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.em
import IMP.pmi.dof
import IMP.atom
#import IMP.saxs

import sys

runType = sys.argv[1] # Specify test or prod
runID = sys.argv[2]   # Specify the number of runs
run_output_dir = 'run_' + str(runID)

if runType == "test":
    num_frames = 7500
elif runType == "prod":
    num_frames = 10000

rex_max_temp = 4.0

# Identify data files
adh_xl_data = "../../input/mhr/xlms/mhr_adh_xl.dat"
bs3dss_xl_data = "../../input/mhr/xlms/mhr_bs3dss_xl.dat"
dmtmm_xl_data = "../../input/mhr/xlms/mhr_dmtmm_xl.dat"
gmm_data = "../../input/mhr/gmm/Abinitio_largedataset_homogenousrefinement_MHR.gmm.50.txt"

# Restraint weights
xl_weight = 10.0
em_weight = 100.0
#saxs_weight = 0.01

# Topology File
topology_file = "../../input/mhr/topology.txt"

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is where the work begins
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# All IMP systems start out with a Model
mdl = IMP.Model()

# Read the topology file for a given state
t = IMP.pmi.topology.TopologyReader(topology_file)

# Create a BuildSystem macro to and add a state from a topology file
bs = IMP.pmi.macros.BuildSystem(mdl)
bs.add_state(t)

# executing the macro will return the root hierarchy and degrees of freedom (dof) objects
root_hier, dof = bs.execute_macro(max_rb_trans= 1.3,
                                  max_rb_rot= 0.1,
                                  max_bead_trans= 4.0,
                                  max_srb_trans= 0.3,
                                  max_srb_rot=0.1)


# Fix all rigid bodies but not Rpb4 and Rpb7 (the stalk)
# First select and gather all particles to fix.
fixed_particles=[]
for prot in ["MTA1"]:
    for cp in [0,1]:
        fixed_particles+=IMP.atom.Selection(root_hier,molecule=prot,copy_index=cp,residue_indexes=range(165,334)).get_selected_particles()
for prot in ["HDAC1"]:
    for cp in [0,1]:
        fixed_particles+=IMP.atom.Selection(root_hier,molecule=prot,copy_index=cp,residue_indexes=range(8,377)).get_selected_particles()

#print(fixed_particles)
# Fix the Corresponding Rigid movers and Super Rigid Body movers using dof
# The flexible beads will still be flexible (fixed_beads is an empty list)!
fixed_beads,fixed_rbs=dof.disable_movers(fixed_particles,
                                         [IMP.core.RigidBodyMover,
                                          IMP.pmi.TransformMover])


# It's useful to have a list of the molecules.
molecules = t.get_components()


##### Uncomment the following lines to get test.rmf file to visualise the system representation
#
# # Uncomment this line for verbose output of the representation
# IMP.atom.show_with_representations(root_hier)
# # output to RMF
# fname = 'test.rmf'
# rh = RMF.create_rmf_file(fname)
# IMP.rmf.add_hierarchy(rh, root_hier)
# IMP.rmf.save_frame(rh)




#####################################################
##################### RESTRAINTS ####################
#####################################################

# Restraints define functions that score the model based on
# input information.
#
# Restraint objects are first created in the definition.
# To be evaluated, the restraint object must be add_to_model().
#
# In some cases, sampled parameters for restraints must be added to the DOF
# object

# The output_objects list is used to collect all restraints
# where we want to log the output in the STAT file.
# Each restraint should be appended to this list.
output_objects = []

# -----------------------------
# %%%%% CONNECTIVITY RESTRAINT
#
# Restrains residues/particles that are connected in sequence
# This should be used for any system without an atomic force field (e.g. CHARMM)
# We apply the restraint to each molecule

for m in root_hier.get_children()[0].get_children():
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
    cr.add_to_model()
    output_objects.append(cr)

print("Connectivity restraint applied")


# -----------------------------
# %%%%% EXCLUDED VOLUME RESTRAINT
#
# Keeps particles from occupying the same area in space.
# Here, we pass a list of all molecule chains to included_objects to apply this to every residue.
# We could also have passed root_hier to obtain the same behavior.
#
# resolution=1000 applies this expensive restraint to the lowest resolution for each particle.
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                            included_objects=[root_hier],
                                            resolution=1000)
output_objects.append(evr)

print("Excluded volume restraint applied")



# # -----------------------------
# # %%%%% MINIMUM PAIR RESTRAINT
# mpr1 = IMP.pmi.restraints.basic.MinimumPairRestraint(tuple_selection1=(468,546,"MTA1",0),tuple_selection2=(1,425,"RBBP4",2),distmax=10.0,root_hier=root_hier,label="MTA1-RBBP4.2_mpr1")
# mpr2 = IMP.pmi.restraints.basic.MinimumPairRestraint(tuple_selection1=(468,546,"MTA1",1),tuple_selection2=(1,425,"RBBP4",3),distmax=10.0,root_hier=root_hier,label="MTA1.1-RBBP4.3_mpr2")
#
# output_objects.append(mpr1)
# output_objects.append(mpr2)


# -------------------------
# %%%%% CROSSLINKING RESTRAINT
#
# Restrains two particles via a distance restraint based on
# an observed crosslink.
#
# First, create the crosslinking database from the input file
# The "standard keys" correspond to a crosslink csv file of the form:
#
# Protein1,Residue1,Protein2,Residue2
# A,18,G,24
# A,18,G,146
# A,50,G,146
# A,50,G,171
# A,50,G,189
#
# This restraint allows for ambiguity in the crosslinked residues,
# a confidence metric for each crosslink and multiple states.
# See the PMI documentation or the MMB book chapter for a
# full discussion of implementing crosslinking restraints.

# This first step is used to translate the crosslinking data file.
# The KeywordsConverter maps a column label from the xl data file
# to the value that PMI understands.
# Here, we just use the standard keys.
# One can define custom keywords using the syntax below.
# For example if the Protein1 column header is "prot_1"
# xldbkc["Protein1"]="prot_1"

# The CrossLinkDataBase translates and stores the crosslink information
# from the file "xl_data" using the KeywordsConverter.

xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc.set_standard_keys()


xldb_adh = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_adh.create_set_from_file(file_name=adh_xl_data,
                              converter=xldbkc)
xlr_adh = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                CrossLinkDataBase=xldb_adh, # The crosslink database.
                length=25,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.0001,          # This adds a linear term to the scoring function
                label="adh",                        #   to bias crosslinks towards each other
                weight=xl_weight)       # Scaling factor for the restraint score.

xldb_bs3dss = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_bs3dss.create_set_from_file(file_name=bs3dss_xl_data,
                                 converter=xldbkc)
xlr_bs3dss = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                CrossLinkDataBase=xldb_bs3dss, # The crosslink database.
                length=25,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.0001,           # This adds a linear term to the scoring function
                label="bs3dss",                        #   to bias crosslinks towards each other
                weight=xl_weight)       # Scaling factor for the restraint score.


xldb_dmtmm = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_dmtmm.create_set_from_file(file_name=dmtmm_xl_data,
                                converter=xldbkc)
xlr_dmtmm = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                CrossLinkDataBase=xldb_dmtmm, # The crosslink database.
                length=16,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.0001,           # This adds a linear term to the scoring function
                label="dmtmm",                        #   to bias crosslinks towards each other
                weight=xl_weight)       # Scaling factor for the restraint score.



output_objects.append(xlr_adh)
output_objects.append(xlr_bs3dss)
output_objects.append(xlr_dmtmm)

print("Cross-linking restraint applied")

# # -------------------------
 # %%%%% EM RESTRAINT
 #
 # Scores a model based on its cross-correlation to an EM density.
 # Since cross-sorrelation is very expensive, we approximate both
 # the EM map and model as a set of 3D Gaussians (done in Representation).
 #
 # First, collect all density particles from the model.
densities = IMP.atom.Selection(root_hier,representation_type=IMP.atom.DENSITIES).get_selected_particles()

emr = IMP.pmi.restraints.em.GaussianEMRestraint(
             densities,                 # Evaluate the restraint using these model densities
             target_fn=gmm_data,        # The EM map, approximated as a gaussian mixture model (GMM)
             slope=0.00000001,          # a small linear restraint to pull objects towards the EM map center
             scale_target_to_mass=True, # Normalizes the total density of the model wrs: EM map. Only set to true
                                        #   if the EM map and "densities" contain the same objects.
             weight=em_weight)          # the scaling factor for the EM score

output_objects.append(emr)

print("EM Restraint Applied")


#####################################################
###################### SAMPLING #####################
#####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.
print("The type of run is: " + str(runType))
print("Number of sampling frames: " + str(num_frames))
# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    max_translation=50,
                                    excluded_rigid_bodies=fixed_rbs,
                                    hierarchies_included_in_collision=fixed_particles)

# Shuffling randomizes the bead positions. It's good to
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
dof.optimize_flexible_beads(500)

# Now, add all of the other restraints to the scoring function to start sampling
evr.add_to_model()
emr.add_to_model()
# mpr1.add_to_model()
# mpr2.add_to_model()
xlr_adh.add_to_model()
xlr_bs3dss.add_to_model()
xlr_dmtmm.add_to_model()

print("Replica Exchange Maximum Temperature : " + str(rex_max_temp))

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
        root_hier=root_hier,                    # pass the root hierarchy
        crosslink_restraints=[xlr_adh,xlr_bs3dss,xlr_dmtmm],         # This allows viewing the crosslinks in Chimera
        monte_carlo_temperature = 1.0,
        replica_exchange_minimum_temperature = 1.0,
        replica_exchange_maximum_temperature = rex_max_temp,
	monte_carlo_sample_objects=dof.get_movers(),  # pass all objects to be moved ( almost always dof.get_movers() )
        global_output_directory=run_output_dir,       # The output directory for this sampling run.
        output_objects=output_objects,          # Items in output_objects write information to the stat file.
        monte_carlo_steps=10,                   # Number of MC steps between writing frames
        number_of_best_scoring_models=0,        # set >0 to store best PDB files (but this is slow)
        number_of_frames=num_frames)            # Total number of frames to run / write to the RMF file.
        #test_mode=test_mode)                    # (Ignore this) Run in test mode (don't write anything)

# Ok, now we finally do the sampling!
rex.execute_macro()

#print(sr.evaluate())
#print(emr.evaluate())
# print(xlr_adh.evaluate())
# print(xlr_bs3dss.evaluate())
# print(xlr_dmtmm.evaluate())
# print(mpr1.evaluate())
# print(mpr2.evaluate())

# Outputs are then analyzed in a separate analysis script.
