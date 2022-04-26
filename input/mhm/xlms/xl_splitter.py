import sys
import os

xl_master_path = os.getcwd() +'/' + str(sys.argv[1])

adh = []
bs3dss = []
dmtmm = []

with open(xl_master_path,'r') as xlf:
    for ln in xlf.readlines():
        print(ln)
        ln1 = ln.strip()
        if ln1.split(',')[0].startswith('ADH'):
            adh.append(ln1)
        elif ln1.split(',')[0].startswith('BS3_DSS'):
            bs3dss.append(ln1)
        elif ln1.split(',')[0].startswith('DMTMM'):
            dmtmm.append(ln1)

nADH = len(adh)
nBS3DSS = len(bs3dss)
nDMTMM = len(dmtmm)

print(nADH,nBS3DSS,nDMTMM)
print("Total: ", nADH + nBS3DSS + nDMTMM)

with open('adh.dat', 'w') as ad:
    for ia in adh:
        ad.write(ia+'\n')
with open('bs3dss.dat', 'w') as bs:
    for ib in bs3dss:
        bs.write(ib+'\n')
with open('dmtmm.dat', 'w') as dm:
    for ic in dmtmm:
        dm.write(ic+'\n')
