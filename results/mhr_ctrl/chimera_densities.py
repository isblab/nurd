import chimera
from chimera import openModels
from chimera import runCommand

# Names of proteins/domains for which we have created densities
prots = ['HDAC1.0_C-term_377-482',
'HDAC1.1_C-term_377-482',\
'MTA1.0_BAH_1-164',\
'MTA1.1_BAH_1-164',\
'MTA1-HDAC1',\
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

# Color of each protein/domain
col = {'HDAC1.0_C-term_377-482':'yellow',
'HDAC1.1_C-term_377-482':'yellow',\
'MTA1.0_BAH_1-164':'sandy brown',\
'MTA1.1_BAH_1-164':'sandy brown',\
'MTA1-HDAC1':'orange',\
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

# Set thresholds for each protein
thresholds = {'HDAC1.0_C-term_377-482':0.0105,
'HDAC1.1_C-term_377-482':0.02,\
'MTA1.0_BAH_1-164':0.00744,\
'MTA1.1_BAH_1-164':0.00471,\
'MTA1-HDAC1':0.129,\
'MTA1.0_mid_334-467':0.005,\
'MTA1.1_mid_334-467':0.015,\
'MTA1.0_R1_RBBP4.0':0.08,\
'MTA1.1_R1_RBBP4.1':0.09,\
'MTA1.0_between_R1_and_R2':0.028,\
'MTA1.1_between_R1_and_R2':0.0115,\
'MTA1.0_R2_RBBP4.2':0.04,\
'MTA1.1_R2_RBBP4.3':0.059,\
'MTA1.0_C-term_692-715':0.00405,\
'MTA1.1_C-term_692-715':0.0025}

# transp = {'HDAC1.0_C-term_377-482':0.3,
# 'HDAC1.1_C-term_377-482':0.3,\
# 'MTA1-HDAC1':0.3,\
# 'MTA1.0_BAH_1-164':0.3,\
# 'MTA1.1_BAH_1-164':0.3,\
# 'MTA1.0_mid_334-467':0.5,\
# 'MTA1.1_mid_334-467':0.5,\
# 'MTA1.0_R1_RBBP4.0':0.2,\
# 'MTA1.1_R1_RBBP4.1':0.2,\
# 'MTA1.0_between_R1_and_R2':0.2,\
# 'MTA1.1_between_R1_and_R2':0.2,\
# 'MTA1.0_R2_RBBP4.2':0.2,\
# 'MTA1.1_R2_RBBP4.3':0.2,\
# 'MTA1.0_C-term_692-715':0.2,\
# 'MTA1.1_C-term_692-715':0.2}


runCommand('set bgcolor white')
i=0

#Read localization density by component, both samples together
for p in prots:
    runCommand('open LPD_'+p+'.mrc')
    runCommand('volume #'+str(i)+' step 1 ')
    #runCommand('volume #'+str(i)+' transparency '+str(transp[p]))
    runCommand('color '+col[p]+' #'+str(i))
    runCommand('volume #'+str(i)+' level '+str(thresholds[p]))
    #runCommand('2dlabels create ' + str(p) + '_lab text ' + str(p) + ' color '+col[p]+' size 30 xpos .1 ypos ' + str( 0.85 - i / 20.0))
    i += 1
