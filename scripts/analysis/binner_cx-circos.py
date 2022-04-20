import sys
import os

infile = sys.argv[1]
outfile = 'binned_'+infile
threshold1 = 10.0
threshold2 = 20.0
threshold3 = 30.0
threshold4 = 40.0
all_lengths = []
all_lines = []
with open(infile,'r') as f:
    with open(outfile,'w') as w:
        for ln in f.readlines():
            all_lines.append(ln)
            line = ln.strip().split(',')
            length = float(line[4])

            if length <= threshold1:
                len_bin = 1
            elif length <= threshold2:
                len_bin = 2
            elif length <= threshold3:
                len_bin = 3
            elif length <= threshold4:
                len_bin = 4
            else:
                len_bin = 5
            line[4] = str(len_bin)
            new_line = ','.join(line)
            print(new_line)
            w.write(new_line+'\n')
print('Done')
