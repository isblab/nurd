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
        'P66A.0',
        'MBD3.1_C-term',
        'MBD3.1_mid_struc',
        'MBD3.1_mid_uns',
        'MBD3.1_N-term',
        'P66A.1']

# Set visualization thresholds
# threshold = {'HDAC1.0_C-term_377-482':0.133,       # 20%
#             'HDAC1.1_C-term_377-482':0.111,        # 20%
#             'HDAC1.0_1-376':0.153,                 # 20%
#             'HDAC1.1_1-376':0.153,                 # 20%
#             'MTA1.0_BAH_1-164':0.126,              # 20%
#             'MTA1.1_BAH_1-164':0.125,              # 20%
#             'MTA1.0_ELM2_SANT_165-333':0.121,      # 20%
#             'MTA1.1_ELM2_SANT_165-333':0.127,      # 20%
#             'MTA1.0_mid_334-431':0.0876,           # 20%
#             'MTA1.1_mid_334-431':0.0658,           # 20%
#             'MBD3.0_C-term':0.0475,                # 10%
#             'MBD3.0_mid_struc':0.0148,             # 10%
#             'MBD3.0_mid_uns':0.0763,               # 10%
#             'MBD3.0_N-term':0.0372,                # 10%
#             'P66A.0':0.0191,                       # 10%
#             'MBD3.1_C-term':0.0554,                # 10%
#             'MBD3.1_mid_struc':0.0202,             # 10%
#             'MBD3.1_mid_uns':0.0581,               # 10%
#             'MBD3.1_N-term':0.0462,                # 10%
#             'P66A.1':0.02}                         # 10%

mostly_old_threshold = {'HDAC1.0_C-term_377-482':0.133,       # 20%
            'HDAC1.1_C-term_377-482':0.111,        # 20%
            'HDAC1.0_1-376':0.153,                 # 20%
            'HDAC1.1_1-376':0.153,                 # 20%
            'MTA1.0_BAH_1-164':0.126,              # 20%
            'MTA1.1_BAH_1-164':0.125,              # 20%
            'MTA1.0_ELM2_SANT_165-333':0.121,      # 20%
            'MTA1.1_ELM2_SANT_165-333':0.127,      # 20%
            'MTA1.0_mid_334-431':0.0876,           # 20%
            'MTA1.1_mid_334-431':0.0658,           # 20%
            'MBD3.0_C-term':0.0784,                # 20%
            'MBD3.0_mid_struc':0.0296,             # 20%
            'MBD3.0_mid_uns':0.0763,               # 10%
            'MBD3.0_N-term':0.0744,                # 20%
            'P66A.0':0.0517,                       # 27%
            'MBD3.1_C-term':0.111,                 # 20%
            'MBD3.1_mid_struc':0.0404,             # 20%
            'MBD3.1_mid_uns':0.0581,               # 10%
            'MBD3.1_N-term':0.0924,                # 20%
            'P66A.1':0.0482}                       # 26.69%



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
    'MBD3.0_C-term':'dark green',
    'MBD3.0_mid_struc':'lime green',
    'MBD3.0_mid_uns':'chartreuse',
    'MBD3.0_N-term':'spring green',
    'P66A.0':'dim gray',
    'MBD3.1_C-term':'#f62681',
    'MBD3.1_mid_struc':'#f85781',
    'MBD3.1_mid_uns':'#ff9eb6',
    'MBD3.1_N-term':'#ffd1e3',
    'P66A.1':'slate gray'}


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
