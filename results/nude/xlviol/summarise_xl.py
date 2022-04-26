import os,sys

input_xls = (sys.argv[1]).split(',')
threshold = float(sys.argv[2])

with open('summary_xls.dat','a') as outf:
	for file in input_xls:
		with open(file,'r') as inf:
			for ln in inf.readlines():
				dist = float(ln.strip().split(',')[-1])
				if dist > threshold:
					category = 2
				else:
					category = 1

				outf.write(f"{','.join(ln.strip().split(',')[:-1])},{category}\n")

print("Done!")