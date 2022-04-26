import os,sys

xl_file = sys.argv[1]
split_at = sys.argv[2]
outfile = 'binned_' + xl_file

threshold1 = 10.0
threshold2 = 20.0
threshold3 = 30.0
threshold4 = 40.0

with open(xl_file,'r') as inf:
	with open(outfile,'w') as outf:
		for ln in inf.readlines():
			key, distance = ln.strip().split(split_at)[0], float(ln.strip().split(split_at)[1]) 
			flag = None
			if distance <= threshold1:
				flag = 1
			elif distance <= threshold2:
				flag = 2
			elif distance <= threshold3:
				flag = 3
			elif distance <= threshold4:
				flag = 4
			else:
				flag = 5

			outline = f"{key},{str(flag)}\n"
			outf.write(outline)