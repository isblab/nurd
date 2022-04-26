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
threshold = {'HDAC1.0_C-term_377-482':0.0424,           # ~10%
        'HDAC1.1_C-term_377-482':0.0577,                # ~10%
        'MTA1.0_BAH_1-164':0.0601,                      # ~10%
        'MTA1.1_BAH_1-164':0.0349,                      # ~10%
        'HDAC1.0_1-376':0.0755,                         # ~10%
        'HDAC1.1_1-376':0.0774,                         # ~10%
        'MTA1.0_ELM2-SANT_165-333':0.0569,              # ~10%
        'MTA1.1_ELM2-SANT_165-333':0.06,                # ~10%
        'MTA1.0_mid_334-467':0.0259,                    # ~10%
        'MTA1.1_mid_334-467':0.033,                     # ~10%
        'MBD3_C-term':0.0241,                           # ~10%
        'MBD3_mid_uns':0.0376,                          # ~10%
        'MBD3_mid_struc':0.0125,                        # ~10%
        'MBD3_N-term':0.0321,                           # ~10%
        'P66A':0.0264,                                  # ~20%
        'MTA1.0_R1_468-546':0.0114,                     # ~10%
        'MTA1.1_R1_468-546':0.0139,                     # ~10%
        'MTA1.0_between_R1_and_R2':0.0515,              # ~10%
        'MTA1.1_between_R1_and_R2':0.0348,              # ~10%
        'MTA1.0_R2_670-691':0.0115,                     # ~10%
        'MTA1.0_R2_670-691':0.0115,                     # ~10%
        'MTA1.0_C-term_692-715':0.0144,                 # ~10%
        'MTA1.1_C-term_692-715':0.00998,                # ~10%
        'RBBP4.0_1-425':0.0418,                         # ~10%
        'RBBP4.1_1-425':0.0594,                         # ~10%
        'RBBP4.2_1-425':0.0403,                         # ~10%
        'RBBP4.3_1-425':0.0642}                         # ~10%



# Color of each protein/domain
col = {'HDAC1.0_1-376':'yellow',\
    'HDAC1.1_1-376':'yellow',\
    'HDAC1.0_C-term_377-482':'#c2a300',\
    'HDAC1.1_C-term_377-482':'#c2a300',\
    'MTA1.0_BAH_1-164':'#FFA500',\
    'MTA1.1_BAH_1-164':'#FFA500',\
    'MTA1.0_ELM2-SANT_165-333':'#FF8C00',\
    'MTA1.1_ELM2-SANT_165-333':'#FF8C00',\
    'MTA1.0_mid_334-467':'#f27130',\
    'MTA1.1_mid_334-467':'#f27130',\
    'MTA1.0_R1_468-546':'#ff5f24',\
    'MTA1.1_R1_468-546':'#ff5f24',\
    'MTA1.0_between_R1_and_R2':'#FF6347',\
    'MTA1.1_between_R1_and_R2':'#FF6347',\
    'MTA1.0_R2_670-691':'#ff4f38',\
    'MTA1.1_R2_670-691':'#ff4f38',\
    'MTA1.0_C-term_692-715':'#ff4040',\
    'MTA1.1_C-term_692-715':'#ff4040',\
    'RBBP4.0_1-425':'#0984ad',\
    'RBBP4.1_1-425':'#0984ad',\
    'RBBP4.2_1-425':'#88d7f7',\
    'RBBP4.3_1-425':'#88d7f7',\
    'MBD3_C-term':'dark green',
    'MBD3_mid_uns':'#84ff6b',
    'MBD3_mid_struc':'#4ca84c',
    'MBD3_N-term':'#4affa4',
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

# runCommand('open emd_22904.mrc')
# runCommand('volume #'+str(i)+' step 1 level 0.294 style surface color "#aeaeae" transparency 0.8')
