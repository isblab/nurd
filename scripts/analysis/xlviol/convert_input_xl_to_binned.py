import sys,os
import glob

in_dir = '/home/shreyasarvindekar/Projects/NuRD/revision/input/nude/xlms/*_master.dat'
files = glob.glob(in_dir)

outxls = []
for file in files:
    with open(file,'r') as inf:
        for ln in inf.readlines():
            f = file.split('/')[-1]
            if (not ln.startswith('Linker')) and (not ln.startswith('Protein')):
                if f.startswith('adh'):
                    bin_num = 1
                elif f.startswith('bs3dss'):
                    bin_num = 2
                elif f.startswith('dmtmm'):
                    bin_num = 3
                else:
                    print('I dont know this XL type. I am quitting')
                    exit(1)

                if len(ln.split(',')) > 4:
                    xl = f"{','.join(ln.strip().split(',')[1:])},{bin_num}\n"
                    outxls.append(xl)

with open('circos-ready_allxls.dat','w') as outf:
    for xl in outxls:
        outf.write(xl)
