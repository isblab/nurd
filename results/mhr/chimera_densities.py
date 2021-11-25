import chimera
from chimera import openModels
from chimera import runCommand

em_map_name = 'Abinitio_largedataset_homogenousrefinement_MHR.mrc'
# Names of proteins/domains for which we have created densities
prots = ['HDAC1.0_C-term_377-482',\
'HDAC1.1_C-term_377-482',\
'MTA1.0_BAH_1-164',\
'MTA1.1_BAH_1-164',\
'MTA1.0_ELM2-SANT_165-333',\
'MTA1.1_ELM2-SANT_165-333',\
'HDAC1.0_1-376',\
'HDAC1.1_1-376',\
'MTA1.0_mid_334-467',\
'MTA1.1_mid_334-467',\
'MTA1.0_R1_468-546',\
'MTA1.1_R1_468-546',\
'RBBP4.0_1-425',\
'RBBP4.1_1-425',\
'MTA1.0_between_R1_and_R2',\
'MTA1.1_between_R1_and_R2',\
'MTA1.0_R2_670-691',\
'MTA1.1_R2_670-691',\
'RBBP4.2_1-425',\
'RBBP4.3_1-425',\
'MTA1.0_C-term_692-715',\
'MTA1.1_C-term_692-715']

# Set visualization thresholds
threshold = {'HDAC1.0_C-term_377-482':0.0207,           # 10%
'HDAC1.1_C-term_377-482':0.019,                         # 10%
'MTA1.0_BAH_1-164':0.0323,                              # 10%
'MTA1.1_BAH_1-164':0.0261,                              # 10%
'MTA1.0_ELM2-SANT_165-333':0.0542,                      # 10%
'MTA1.1_ELM2-SANT_165-333':0.0451,                      # 10%
'HDAC1.0_1-376':0.0724,                                 # 10%
'HDAC1.1_1-376':0.07,                                   # 10%
'MTA1.0_mid_334-467':0.0428,                            # 10%
'MTA1.1_mid_334-467':0.0368,                            # 10%
'MTA1.0_R1_468-546':0.0217,                             # 10%
'MTA1.1_R1_468-546':0.0271,                             # 10%
'RBBP4.0_1-425':0.0345,                                 # 10%
'RBBP4.1_1-425':0.0476,                                 # 10%
'MTA1.0_between_R1_and_R2':0.0465,                      # 10%
'MTA1.1_between_R1_and_R2':0.037,                       # 10%
'MTA1.0_R2_670-691':0.0108,                             # 10%
'MTA1.1_R2_670-691':0.0133,                             # 10%
'RBBP4.2_1-425':0.0499,                                 # 10%
'RBBP4.3_1-425':0.0602,                                 # 10%
'MTA1.0_C-term_692-715':0.00756,                        # 10%
'MTA1.1_C-term_692-715':0.00548}                        # 10%


# Color of each protein/domain
col = {'HDAC1.0_C-term_377-482':'goldenrod',
'HDAC1.1_C-term_377-482':'goldenrod',\
'MTA1.0_ELM2-SANT_165-333':'orange',\
'MTA1.1_ELM2-SANT_165-333':'orange',\
'HDAC1.0_1-376':'yellow',\
'HDAC1.1_1-376':'yellow',\
'MTA1.0_BAH_1-164':'sandy brown',\
'MTA1.1_BAH_1-164':'sandy brown',\
'MTA1.0_mid_334-467':'orange red',\
'MTA1.1_mid_334-467':'orange red',\
'MTA1.0_R1_468-546':'red',\
'MTA1.1_R1_468-546':'red',\
'RBBP4.0_1-425':'blue',\
'RBBP4.1_1-425':'blue',\
'MTA1.0_between_R1_and_R2':'brown',\
'MTA1.1_between_R1_and_R2':'brown',\
'MTA1.0_R2_670-691':'red',\
'MTA1.1_R2_670-691':'red',\
'RBBP4.2_1-425':'cornflower blue',\
'RBBP4.3_1-425':'cornflower blue',\
'MTA1.0_C-term_692-715':'dark red',\
'MTA1.1_C-term_692-715':'dark red'}


runCommand('set bgcolor white')
i=0

#Read localization density by component, both samples together
for p in prots:
    runCommand('open LPD_'+p+'.mrc')
    runCommand('volume #'+str(i)+' step 1 ')
    runCommand('volume #'+str(i)+' level '+str(threshold[p]))
    # runCommand('volume #'+str(i)+' transparency '+str(transp[p]))
    runCommand('color '+col[p]+' #'+str(i))
    # runCommand('2dlabels create ' + str(p) + '_lab text ' + str(p) + ' color '+col[p]+' size 30 xpos .1 ypos ' + str( 0.85 - i / 20.0))
    i += 1

runCommand('open '+str(em_map_name))
runCommand('volume #'+str(i)+' step 1 level 0.211 style mesh color "slate grey"')
