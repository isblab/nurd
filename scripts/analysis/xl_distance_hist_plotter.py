import sys
import os
from matplotlib import pyplot as plt

infile = sys.argv[1]
xltype = sys.argv[2]
threshold = sys.argv[3]
distances = []
with open(infile,'r') as inf:
    for ln in inf.readlines():
        ln = ln.strip().split(',')
        dist = int(round(float(ln[4]),0))
        distances.append(dist)

x = []
y = []
for i in range(len(distances)+1):
    x.append(int(threshold))
    y.append(i)

plt.figure()
plt.plot(x,y, color='red')
plt.hist(distances, bins=10, range=[0,50], color='#0095FF')
plt.title(f'XL Distances for {xltype.upper()}')
plt.xlabel('Distance')
plt.ylabel('Number of XLs')
plt.savefig(f'{xltype}_hist.png')
plt.show()
