import sys
import os

in_file = sys.argv[1]
intra_out_file = 'intra_' + in_file
inter_out_file = 'inter_' + in_file

with open(in_file,'r') as inf:
    with open(intra_out_file,'w') as intra_outf:
        with open(inter_out_file,'w') as inter_outf:
            for ln in inf.readlines():
                if ln.startswith('Protein1'):
                    intra_outf.write(ln)
                    inter_outf.write(ln)
                else:
                    ln = ln.strip().split(',')
                    if ln[0]==ln[2]:
                        ln = ','.join(ln) + '\n'
                        intra_outf.write(ln)
                    else:
                        ln = ','.join(ln) + '\n'
                        inter_outf.write(ln)

print('Done!')
