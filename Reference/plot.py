# Plot differences between sets of TEME vectors

import matplotlib.pyplot as plt

with open('STR3.txt') as f:
  ref = f.readlines()
with open('AIAA.txt') as f:
  diff1 = f.readlines()

sats = []
dotsp = {}

for i, line in enumerate(ref):
  if (line[:3] == '---'):
    current_sat = line[4:].strip()
    dotsp[current_sat] = ([], [])
  else:
    ref_vals = list(map(float, line.split()))
    diff1_vals = list(map(float, diff1[i].split()))
    dotsp[current_sat][0].append(ref_vals[0])
    dotsp[current_sat][1].append(1000 * max(abs(ref_vals[1] - diff1_vals[1]), abs(ref_vals[2] - diff1_vals[2]), abs(ref_vals[3] - diff1_vals[3])))

xmin = -1440
xmax = 1440
xticks = range(xmin, xmax + 1, 60)
xticksmaj = range(xmin, xmax + 1, 240)
    
plt.subplot(111)
for st in dotsp.keys():  
  plt.plot(dotsp[st][0], dotsp[st][1], c='b', marker='D', ls="", mew=0.0, markersize=3, alpha=0.55)
plt.subplot(111).set_xlim(xmin, xmax)
plt.subplot(111).set_xticks(xticks, True);
plt.subplot(111).set_xticks(xticksmaj, False);
plt.subplot(111).set_yscale("log")
plt.subplot(111).autoscale(True, 'y', None)
plt.subplot(111).grid(which='major', axis='both', ls='dashed', alpha=0.7)
plt.subplot(111).grid(which='minor', axis='x', ls='dotted', alpha=0.3)
plt.title('Position vector: STR#3 vs AIAA (blue)')

plt.show()