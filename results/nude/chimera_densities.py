import chimera
from chimera import openModels
from chimera import runCommand

# Names of proteins/domains for which we have created densities
prots = ['HDAC1.0_C-term_377-482',
        'HDAC1.1_C-term_377-482',
        'MTA1.0_BAH_1-164',
        'MTA1.1_BAH_1-164',
        'HDAC1.0_1-376',
        'HDAC1.1_1-376',
        'MTA1.0_ELM2-SANT_165-333',
        'MTA1.1_ELM2-SANT_165-333',
        'MTA1.0_mid_334-467',
        'MTA1.1_mid_334-467',
        'MBD3_C-term',
        'MBD3_mid_uns',
        'MBD3_mid_struc',
        'MBD3_N-term',
        'P66A',
        'MTA1.0_R1_468-546',
        'MTA1.1_R1_468-546',
        'MTA1.0_between_R1_and_R2',
        'MTA1.1_between_R1_and_R2',
        'MTA1.0_R2_670-691',
        'MTA1.0_R2_670-691',
        'MTA1.0_C-term_692-715',
        'MTA1.1_C-term_692-715',
        'RBBP4.0_1-425',
        'RBBP4.1_1-425',
        'RBBP4.2_1-425',
        'RBBP4.3_1-425']

# Set visualization thresholds
threshold = {'HDAC1.0_C-term_377-482':0.0285,           # 7%
        'HDAC1.1_C-term_377-482':0.0394,                # 7%
        'MTA1.0_BAH_1-164':0.0458,                      # ~10%
        'MTA1.1_BAH_1-164':0.0307,                      # ~10%
        'HDAC1.0_1-376':0.0684,                         # ~10%
        'HDAC1.1_1-376':0.0735,                         # ~10%
        'MTA1.0_ELM2-SANT_165-333':0.0567,              # ~10%
        'MTA1.1_ELM2-SANT_165-333':0.0594,              # ~10%
        'MTA1.0_mid_334-467':0.0234,                    # ~10%
        'MTA1.1_mid_334-467':0.0325,                    # ~10%
        'MBD3_C-term':0.0259,                           # ~10%
        'MBD3_mid_uns':0.0442,                          # ~10%
        'MBD3_mid_struc':0.0108,                        # ~10%
        'MBD3_N-term':0.0439,                           # ~10%
        'P66A':0.0182,                                  # ~10%
        'MTA1.0_R1_468-546':0.0178,                     # ~10%
        'MTA1.1_R1_468-546':0.0208,                     # ~10%
        'MTA1.0_between_R1_and_R2':0.049,               # ~10%
        'MTA1.1_between_R1_and_R2':0.051,               # ~10%
        'MTA1.0_R2_670-691':0.011,                      # ~10%
        'MTA1.0_R2_670-691':0.011,                      # ~10%
        'MTA1.0_C-term_692-715':0.011,                  # ~10%
        'MTA1.1_C-term_692-715':0.011,                  # ~10%
        'RBBP4.0_1-425':0.037,                          # ~10%
        'RBBP4.1_1-425':0.050,                          # ~10%
        'RBBP4.2_1-425':0.040,                          # ~10%
        'RBBP4.3_1-425':0.049}                          # ~10%



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
    'MTA1.1_C-term_692-715':'dark red',\
    'MBD3_C-term':'dark green',
    'MBD3_mid_uns':'chartreuse',
    'MBD3_mid_struc':'lime green',
    'MBD3_N-term':'spring green',
    'P66A':'dim gray'}


runCommand('set bgcolor white')
i=0

#Read localization density by component, both samples together
for p in prots:
    runCommand('open LPD_'+p+'.mrc')
    runCommand('volume #'+str(i)+' step 1 ')
    runCommand('volume #'+str(i)+' level '+str(threshold[p]))
    # runCommand('volume #'+str(i)+' transparency '+str(transp[p]))
    runCommand('color '+col[p]+' #'+str(i))
    # runCommand('2dlabels create ' + str(p) + '_lab text ' + str(p) + ' color '+col[p]+' size 30 xpos .1 ypos ' + str(0.9 - i / 25.0))
    i += 1

runCommand('open emd_22904.mrc')
runCommand('volume #'+str(i)+' step 2 level 0.294 style mesh color "dark slate gray"')
