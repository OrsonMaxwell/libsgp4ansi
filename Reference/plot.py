#!/bin/env python
# Plot differences between sets of TEME vectors

import matplotlib.pyplot as plt
import re
import os, sys

with open(sys.argv[1] + '/a721.out') as f:
  ref = f.readlines()
with open(sys.argv[1] + '/i72.out') as f:
  diff1 = f.readlines()

current_sat = ''
orbital_period = 0.0
diff = 0.0
dotsp = {}

for i, line in enumerate(ref):
  m = re.fullmatch('([0-9]{1,5})\ \(([0-9\.]+)\)', line.strip())
  if m != None:
    if len(dotsp.keys()) > 0:
      dotsp[current_sat][0].append(orbital_period)
      dotsp[current_sat][1].append(1000 * diff) # get back to meters
      print('[{0}] Sat {1}, period: {2:8.3f}, '.format(i, current_sat, orbital_period), end='')
      print('max deviation: {0:10.9f}'.format(1000 * diff))
    current_sat = m.group(1)
    orbital_period = float(m.group(2))
    dotsp[current_sat] = ([], [])
    diff = 0.0
  else:
    ref_vals = list(map(float, line.split()))
    diff1_vals = list(map(float, diff1[i].split()))
    cur_line_max_diff = max(abs(ref_vals[1] - diff1_vals[1]), abs(ref_vals[2] - diff1_vals[2]), abs(ref_vals[3] - diff1_vals[3]))
    if diff < cur_line_max_diff:
      diff = cur_line_max_diff

dotsp[current_sat][0].append(orbital_period)
dotsp[current_sat][1].append(1000 * diff) # get back to meters
print('[   EOF] Sat {0}, period: {1:8.3f}, '.format(current_sat, orbital_period), end='')
print('max deviation: {0:10.9f}'.format(1000 * diff))
      
print('Loaded {} data points'.format(len(dotsp.keys())))

xmin = 0
xmax = 2000
xticksmin = range(xmin, xmax + 1, 100)
xticksmaj = range(xmin, xmax + 1, 500)
    
plt.subplot(111)
for st in dotsp.keys():  
  plt.plot(dotsp[st][0], dotsp[st][1], c='b', marker='D', ls="", mew=0.0, markersize=3, alpha=0.55)
plt.subplot(111).set_xlim(xmin, xmax)
plt.subplot(111).set_xticks(xticksmin, True);
plt.subplot(111).set_xticks(xticksmaj, False);
plt.subplot(111).set_yscale("log")
plt.subplot(111).autoscale(True, 'y', None)
plt.subplot(111).grid(which='major', axis='both', ls='dashed', alpha=0.7)
plt.subplot(111).grid(which='minor', axis='x', ls='dotted', alpha=0.3)
plt.title('Max difference in ephemerides, m: STR#3 vs AIAA (blue)')
plt.xlabel('Orbital period, min')
plt.ylabel('Difference, m')
plt.tight_layout(pad=0.0, w_pad=0.1, h_pad=0.1)

plt.show()