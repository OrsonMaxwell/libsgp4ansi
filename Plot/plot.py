#!/bin/env python
# Plot differences between sets of TEME vectors

import matplotlib.pyplot as plt
import re
import os, sys
import subprocess
import datetime

# Run reference executable
print('--- Running Sattrack Report #3 verification mode...')
start = datetime.datetime.now()
print(subprocess.Popen('sgp4.exe a v SGP4-VER.TLE', stdout=subprocess.PIPE, creationflags=0x08000000).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_av = stop - start
print('--- Done in {}s'.format(tdiff_av.total_seconds()))

print('--- Running Sattrack Report #3 full catalogue mode...')
start = datetime.datetime.now()
print(subprocess.Popen('sgp4.exe a c full.tle', stdout=subprocess.PIPE, creationflags=0x08000000).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_ac = stop - start
print('--- Done in {}s'.format(tdiff_ac.total_seconds()))

print('--- Running AIAA-2006-6753 verification mode...')
start = datetime.datetime.now()
print(subprocess.Popen('sgp4.exe i v SGP4-VER.TLE', stdout=subprocess.PIPE, creationflags=0x08000000).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_iv = stop - start
print('--- Done in {}s'.format(tdiff_iv.total_seconds()))

print('--- Running AIAA-2006-6753 full catalogue mode...')
start = datetime.datetime.now()
print(subprocess.Popen('sgp4.exe i c full.tle', stdout=subprocess.PIPE, creationflags=0x08000000).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_ic = stop - start
print('--- Done in {}s'.format(tdiff_ic.total_seconds()))

print('--- Running libsgp4ansi full catalogue mode...')
start = datetime.datetime.now()
#print(subprocess.Popen('ansi.exe v full.tle', stdout=subprocess.PIPE, creationflags=0x08000000).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_ic = stop - start
print('--- Done in {}s'.format(tdiff_ic.total_seconds()))

print('--- Loading full catalogue data...')

with open('a721.out') as f:
  ref = f.readlines()
with open('i72.out') as f:
  diff1 = f.readlines()
with open('ansi.out') as f:
  diff2 = f.readlines()

current_sat = ''
orbital_period = 0.0
delta1 = 0.0
delta2 = 0.0
dotsp1 = {}
dotsp2 = {}

for i, line in enumerate(ref):
  m = re.fullmatch('([0-9]{1,5})\ \(([0-9\.]+)\)', line.strip())
  if m != None:
    if len(dotsp1.keys()) > 0:
      dotsp1[current_sat][0].append(orbital_period)
      dotsp1[current_sat][1].append(1000 * delta1) # get back to meters
      dotsp2[current_sat][0].append(orbital_period)
      dotsp2[current_sat][1].append(1000 * delta2)
      print('[{0}] Sat {1}, per:{2:8.3f}, '.format(i, current_sat, orbital_period), end='')
      print('d1: {0:10.9f}, d2: {1:10.9f}\r'.format(1000 * delta1, 1000 * delta2), end='')
    current_sat = m.group(1)
    orbital_period = float(m.group(2))
    dotsp1[current_sat] = ([], [])
    dotsp2[current_sat] = ([], [])
    delta1 = 0.0
    delta2 = 0.0
  else:
    ref_vals = list(map(float, line.split()))
    diff1_vals = list(map(float, diff1[i].split()))
    diff2_vals = list(map(float, diff2[i].split()))
    tempdelta1 = max(abs(ref_vals[1] - diff1_vals[1]), abs(ref_vals[2] - diff1_vals[2]), abs(ref_vals[3] - diff1_vals[3]))
    if delta1 < tempdelta1:
      delta1 = tempdelta1
    tempdelta2 = max(abs(ref_vals[1] - diff2_vals[1]), abs(ref_vals[2] - diff2_vals[2]), abs(ref_vals[3] - diff2_vals[3]))
    if delta2 < tempdelta2:
      delta2 = tempdelta2

# Repeat for last sat
dotsp1[current_sat][0].append(orbital_period)
dotsp1[current_sat][1].append(1000 * delta1) # get back to meters
dotsp2[current_sat][0].append(orbital_period)
dotsp2[current_sat][1].append(1000 * delta2)
print('[   EOF] Sat {0}, per:{1:8.3f}, '.format(current_sat, orbital_period), end='')
print('d1: {0:10.9f}, d2: {1:10.9f}'.format(1000 * delta1, 1000 * delta2))
print('--- Loaded {} data points per set'.format(len(dotsp1.keys())))

print('--- Plotting data - please wait...')
# Plotting
xmin = 0
xmax = 2000
xticksmin = range(xmin, xmax + 1, 100)
xticksmaj = range(xmin, xmax + 1, 500)
    
plt.subplot(111)
for st in dotsp1.keys():  
  plt.plot(dotsp1[st][0], dotsp1[st][1], c='b', marker='D', ls='', mew=0.0, markersize=3, alpha=0.55)
for st in dotsp2.keys():  
  plt.plot(dotsp2[st][0], dotsp2[st][1], c='r', marker='D', ls='', mew=0.0, markersize=3, alpha=0.55)
plt.subplot(111).set_xlim(xmin, xmax)
plt.subplot(111).set_xticks(xticksmin, True);
plt.subplot(111).set_xticks(xticksmaj, False);
plt.subplot(111).set_yscale('log')
plt.subplot(111).autoscale(True, 'y', None)
plt.subplot(111).grid(which='major', axis='both', ls='dashed', alpha=0.7)
plt.subplot(111).grid(which='minor', axis='x', ls='dotted', alpha=0.3)
plt.title('Max difference in ephemerides: STR#3 vs AIAA (blue), STR#3 vs ANSI (red)')
plt.xlabel('Orbital period, min')
plt.ylabel('Difference, m')
plt.tight_layout(pad=0.0, w_pad=0.1, h_pad=0.1)

plt.show()