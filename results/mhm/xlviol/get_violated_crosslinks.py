import os,sys

__doc__ = 'Given all the crosslink distances, this script will return two files. One will have all the XLs while the other will contain only the violated crosslinks. \
			Both these files will be in p1,r1,p2,r2,dist format. '

xl_dist_file = sys.argv[1]
separator_for_distance = sys.argv[2]
threshold = sys.argv[3]
violated_xls = []
out_fname = xl_dist_file.split('.')[0]+"_format_corrected.txt"
viol_f = xl_dist_file.split('.')[0]+"_violated_xls.txt"

with open(out_fname,'w') as outf:
	with open(viol_f,'w') as viof:
		with open(xl_dist_file,'r') as inf:
			for ln in inf.readlines():
				link = ln.strip().split(separator_for_distance)[0]
				distance = ln.strip().split(separator_for_distance)[1]

				outline = f"{link},{distance}\n"
				outf.write(outline)

				if float(distance) > float(threshold):
					violated_xls.append(outline)
					viof.write(outline)
				

