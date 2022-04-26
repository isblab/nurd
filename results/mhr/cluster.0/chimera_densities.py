import chimera
from chimera import openModels
from chimera import runCommand

em_map_name = 'Abinitio_largedataset_homogenousrefinement_MHR.mrc'
# Names of proteins/domains for which we have created densities
prots = ['HDAC1.0_C-term_377-482',
'HDAC1.1_C-term_377-482',\
'MTA1-HDAC1',\
'MTA1.0_BAH_1-164',\
'MTA1.1_BAH_1-164',\
'MTA1.0_mid_334-467',\
'MTA1.1_mid_334-467',\
'MTA1.0_R1_RBBP4.0',\
'MTA1.1_R1_RBBP4.1',\
'MTA1.0_between_R1_and_R2',\
'MTA1.1_between_R1_and_R2',\
'MTA1.0_R2_RBBP4.2',\
'MTA1.1_R2_RBBP4.3',\
'MTA1.0_C-term_692-715',\
'MTA1.1_C-term_692-715']

# Set visualization thresholds
threshold = {'HDAC1.0_C-term_377-482':0.0207,
'HDAC1.1_C-term_377-482':0.019,\
'MTA1-HDAC1':0.0736,\
'MTA1.0_BAH_1-164':0.0323,\
'MTA1.1_BAH_1-164':0.0261,\
'MTA1.0_mid_334-467':0.0428,\
'MTA1.1_mid_334-467':0.0368,\
'MTA1.0_R1_RBBP4.0':0.0412,\
'MTA1.1_R1_RBBP4.1':0.0549,\
'MTA1.0_between_R1_and_R2':0.0547,\
'MTA1.1_between_R1_and_R2':0.037,\
'MTA1.0_R2_RBBP4.2':0.0533,\
'MTA1.1_R2_RBBP4.3':0.0622,\
'MTA1.0_C-term_692-715':0.0129,\
'MTA1.1_C-term_692-715':0.00548}

# transp = {'HDAC1.0_C-term_377-482':,
# 'HDAC1.1_C-term_377-482':,\
# 'MTA1-HDAC1':,\
# 'MTA1.0_BAH_1-164':,\
# 'MTA1.1_BAH_1-164':,\
# 'MTA1.0_mid_334-467':,\
# 'MTA1.0_mid_334-467':,\
# 'MTA1.0_R1_RBBP4.0':,\
# 'MTA1.1_R1_RBBP4.1':,\
# 'MTA1.0_between_R1_and_R2':,\
# 'MTA1.1_between_R1_and_R2':,\
# 'MTA1.0_R2_RBBP4.2':,\
# 'MTA1.1_R2_RBBP4.3':,\
# 'MTA1.0_C-term_692-715':,\
# 'MTA1.1_C-term_692-715':}

# Color of each protein/domain
col = {'HDAC1.0_C-term_377-482':'yellow',
'HDAC1.1_C-term_377-482':'yellow',\
'MTA1-HDAC1':'orange',\
'MTA1.0_BAH_1-164':'sandy brown',\
'MTA1.1_BAH_1-164':'sandy brown',\
'MTA1.0_mid_334-467':'orange red',\
'MTA1.1_mid_334-467':'orange red',\
'MTA1.0_R1_RBBP4.0':'blue',\
'MTA1.1_R1_RBBP4.1':'blue',\
'MTA1.0_between_R1_and_R2':'cyan',\
'MTA1.1_between_R1_and_R2':'cyan',\
'MTA1.0_R2_RBBP4.2':'cornflower blue',\
'MTA1.1_R2_RBBP4.3':'cornflower blue',\
'MTA1.0_C-term_692-715':'red',\
'MTA1.1_C-term_692-715':'red'}


runCommand('set bgcolor white')
i=0

#Read localization density by component, both samples together
for p in prots:
    runCommand('open LPD_'+p+'.mrc')
    runCommand('volume #'+str(i)+' step 1 ')
    runCommand('volume #'+str(i)+' level '+str(threshold[p]))
    # runCommand('volume #'+str(i)+' transparency '+str(transp[p]))
    runCommand('color '+col[p]+' #'+str(i))
    runCommand('2dlabels create ' + str(p) + '_lab text ' + str(p) + ' color '+col[p]+' size 30 xpos .1 ypos ' + str( 0.85 - i / 20.0))
    i += 1

runCommand('open '+str(em_map_name))
runCommand('volume #'+str(i)+' step 1 level 0.211 style mesh color "dim grey"')
