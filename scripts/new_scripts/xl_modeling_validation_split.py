import os
import sys
import random

xl_file = sys.argv[1]         # Set path to the parent crosslink.dat or crosslink.csv files
percent = 80				  # Percentage of crosslinks that go to modeling subset
training_xl_file = xl_file.split('.')[0] + '_training.dat'
validation_xl_file = xl_file.split('.')[0] + '_validation.dat' 


input_xls = []
with open(xl_file,'r') as xl_inf:
	for ln in xl_inf.readlines():
		if not ln.startswith('Protein1'):
			input_xls.append(ln.strip())

# input_xls = list(set(input_xls))

num_xl_in_train = (percent/100) * len(input_xls)


print(f'The total number of input crosslinks: {len(input_xls)}\nTotal number of crosslinks in training set: {int(num_xl_in_train)}')


random.shuffle(input_xls)
training_xls = []
validation_xls = []

i = 1
while i<= num_xl_in_train:
	training_xls.append(input_xls[i])
	i += 1

for index in range(i-1,len(input_xls)):
	validation_xls.append(input_xls[index])

n=0
for s in training_xls:
	if s in validation_xls:
		print(s)
		n+=1
print(len(training_xls),len(validation_xls),n)

# for xlink in input_xls:
# 	if xlink not in training_xls:
# 		validation_xls.append(xlink)

print(f'Number of crosslinks in training set: {len(training_xls)}\nNumber of crosslinks in validation set: {len(validation_xls)}')


with open(training_xl_file,'w') as txlf:
	txlf.write('Protein1,Residue1,Protein2,Residue2\n')
	for j in training_xls:
		tstr = j + '\n'
		txlf.write(tstr)

with open(validation_xl_file,'w') as vxlf:
	vxlf.write('Protein1,Residue1,Protein2,Residue2\n')
	for k in validation_xls:
		vstr = k + '\n'
		vxlf.write(vstr)
