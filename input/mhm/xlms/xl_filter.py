import os,sys

in_file = sys.argv[1]
first_line_in_file = sys.argv[2]
out_file = 'filtered_' + in_file
in_file_lines= []
out_file_lines = []
with open(in_file,'r') as inf:
    for ln in inf.readlines():
        if first_line_in_file == 'Linker':
            ln = ln.strip().split(',')
            ln = ','.join(ln[1:])
        else:
            ln = ln.strip()

        in_file_lines.append(ln)
print(f'Number of XLs in input file: {len(in_file_lines)}')
print('################################')
for ln in in_file_lines:
    ln = ln.split(',')
    if ln[0]=='MTA1' and int(ln[1])<=431:
        if (ln[2]=='MTA1' and int(ln[3])<=431) or (ln[2]=='HDAC1') or ln[2]=='MBD3' or (ln[2]=='P66A' and int(ln[3])>=136 and int(ln[3])<=178):
            ln = ','.join(ln)
            out_file_lines.append(ln)

    if ln[0]=='HDAC1':
        if (ln[2]=='MTA1' and int(ln[3])<=431) or ln[2]=='HDAC1' or ln[2]=='MBD3' or (ln[2]=='P66A' and int(ln[3])>=136 and int(ln[3])<=178):
            ln = ','.join(ln)
            out_file_lines.append(ln)

    if ln[0]=='MBD3':
        if (ln[2]=='MTA1' and int(ln[3])<=431) or ln[2]=='HDAC1' or ln[2]=='MBD3' or (ln[2]=='P66A' and int(ln[3])>=136 and int(ln[3])<=178):
            ln = ','.join(ln)
            out_file_lines.append(ln)

    if ln[0]=='P66A' and int(ln[1])>=136 and int(ln[1])<=178:
        if (ln[2]=='MTA1' and int(ln[3])<=431) or ln[2]=='HDAC1' or ln[2]=='MBD3' or (ln[2]=='P66A' and int(ln[3])>=136 and int(ln[3])<=178):
            ln = ','.join(ln)
            out_file_lines.append(ln)

print(f'Number of XLs in output file: {len(out_file_lines)}')

with open(out_file,'w') as outf:
    for i in range(len(out_file_lines)):
        outf.write(out_file_lines[i]+'\n')
