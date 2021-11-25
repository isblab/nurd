import chimera
from chimera import openModels
from chimera import runCommand

# Names of proteins/domains for which we have created densities
prots = ['HDAC1.0_C-term_377-482',
        'HDAC1.1_C-term_377-482',
        'HDAC1.0_1-376',
        'HDAC1.1_1-376',
        'MTA1.0_BAH_1-164',
        'MTA1.1_BAH_1-164',
        'MTA1.0_ELM2_SANT_165-333',
        'MTA1.1_ELM2_SANT_165-333',
        'MTA1.0_mid_334-431',
        'MTA1.1_mid_334-431',
        'MBD3.0_C-term',
        'MBD3.0_mid_struc',
        'MBD3.0_mid_uns',
        'MBD3.0_N-term',
        'P66A.0']

# Set visualization thresholds
# threshold = {'HDAC1.0_C-term_377-482':0.242,       # 30%
#             'HDAC1.1_C-term_377-482':0.205,        # 30%
#             'HDAC1.0_1-376':0.235,                 # 30%
#             'HDAC1.1_1-376':0.229,                 # 30%
#             'MTA1.0_BAH_1-164':0.221,              # 30%
#             'MTA1.1_BAH_1-164':0.202,              # 30%
#             'MTA1.0_ELM2_SANT_165-333':0.131,      # 20%
#             'MTA1.1_ELM2_SANT_165-333':0.121,      # 20%
#             'MTA1.0_mid_334-431':0.054,            # ~20%
#             'MTA1.1_mid_334-431':0.0788,           # ~20%
#             'MBD3.0_C-term':0.0937,                # ~20%
#             'MBD3.0_mid_struc':0.0388,             # ~20%
#             'MBD3.0_mid_uns':0.125,                # ~20%
#             'MBD3.0_N-term':0.0734,                # ~20%
#             'P66A.0':0.0374}                       # ~20%
            
mostly_old_threshold = {'HDAC1.0_C-term_377-482':0.163,           # 5%
                    'HDAC1.1_C-term_377-482':0.134,           # 5%
                    'HDAC1.0_1-376':0.235,                 # 30%
                    'HDAC1.1_1-376':0.229,                 # 30%
                    'MTA1.0_BAH_1-164':0.1472,                  # 20%
                    'MTA1.1_BAH_1-164':0.135,                   # 20%
                    'MTA1.0_ELM2_SANT_165-333':0.131,      # 20%
                    'MTA1.1_ELM2_SANT_165-333':0.121,      # 20%
                    'MTA1.0_mid_334-431':0.0605,                # 10%
                    'MTA1.1_mid_334-431':0.0863,                # 10%
                    'MBD3.0_C-term':0.02755,                      # 5%
                    'MBD3.0_mid_struc':0.0111,                   # 5%
                    'MBD3.0_mid_uns':0.0347,                     # 5%
                    'MBD3.0_N-term':0.0186,                      # 5%
                    'P66A.0':0.0154}                              # 5%

# Color of each protein/domain
col = {'HDAC1.0_C-term_377-482':'goldenrod',
    'HDAC1.1_C-term_377-482':'goldenrod',
    'HDAC1.0_1-376':'yellow',
    'HDAC1.1_1-376':'yellow',
    'MTA1.0_BAH_1-164':'sandy brown',
    'MTA1.1_BAH_1-164':'sandy brown',
    'MTA1.0_ELM2_SANT_165-333':'orange',
    'MTA1.1_ELM2_SANT_165-333':'orange',
    'MTA1.0_mid_334-431':'orange red',
    'MTA1.1_mid_334-431':'orange red',
    'MBD3.0_C-term':'lime green',
    'MBD3.0_mid_struc':'lime green',
    'MBD3.0_mid_uns':'lime green',
    'MBD3.0_N-term':'lime green',
    'P66A.0':'dim gray'}


runCommand('set bgcolor white')
i=0

#Read localization density by component, both samples together
for p in prots:
    runCommand('open LPD_'+p+'.mrc')
    runCommand('volume #'+str(i)+' step 1 ')
    runCommand('volume #'+str(i)+' level '+str(mostly_old_threshold[p]))
    # runCommand('volume #'+str(i)+' transparency '+str(transp[p]))
    runCommand('color '+col[p]+' #'+str(i))
    # runCommand('2dlabels create ' + str(p) + '_lab text ' + str(p) + ' color '+col[p]+' size 30 xpos .1 ypos ' + str( 0.85 - i / 20.0))
    i += 1

runCommand('open emd_21382.mrc')
runCommand('volume #'+str(i)+' step 1 level 0.12 style mesh color "dark slate gray"')
