# Preprocessor script to export all maximum non-zero differences between sets
# of TEME vectors to CSV files

import sys

with open(sys.argv[1]) as f:
  data1 = f.readlines()
with open(sys.argv[2]) as f:
  data2 = f.readlines()

dataset1 = sys.argv[1].split('.')[0]
dataset2 = sys.argv[2].split('.')[0]
  
dotsp = {}
dotsv = {}

for i, line in enumerate(data1):
  if (line[:3] != '---'):
    asfc_vals = list(map(float, line.split()))
    impr_vals = list(map(float, data2[i].split()))
    if asfc_vals[0] not in dotsp.keys():
      dotsp[asfc_vals[0]] = []
    if asfc_vals[0] not in dotsv.keys():
      dotsv[asfc_vals[0]] = []
    diff = max(abs(asfc_vals[1] - impr_vals[1]), abs(asfc_vals[2] - impr_vals[2]), abs(asfc_vals[3] - impr_vals[3]))
    if (diff != 0):
      dotsp[asfc_vals[0]].append(diff)
    #diff = max(abs(asfc_vals[4] - impr_vals[4]), abs(asfc_vals[5] - impr_vals[5]), abs(asfc_vals[6] - impr_vals[6]))
    #if (diff != 0):
    #  dotsv[asfc_vals[0]].append(diff)
    
with open('{}-{}pos.csv'.format(dataset1, dataset2), 'w') as f:
  for timestamp in dotsp.keys():
    for diff in dotsp[timestamp]:
      f.write('{},{}\n'.format(timestamp, diff * 1000))

#with open('{}-{}vel.csv'.format(dataset1, dataset2), 'w') as f:
#  for timestamp in dotsp.keys():
#    for diff in dotsv[timestamp]:
#      f.write('{},{}\n'.format(timestamp, diff * 1000))