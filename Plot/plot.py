#!/bin/env python
# Plot differences between sets of TEME vectors

import matplotlib.pyplot as plt
import re
import os, sys
import subprocess
import datetime

sgp4bin = 'sgp4'

ver_filename = 'SGP4-VER.TLE'
fullcat_filename = 'full.tle'

if os.name == 'nt':
  print("Assuming Windows environment")
  sgp4bin += '.exe'
if os.name == 'posix':
  print("Assuming Linux (posix) environment")
  
# Run reference executable
print('Running Sattrack Report #3 verification mode...')
start = datetime.datetime.now()
print(subprocess.Popen([sgp4bin, 'a', 'v', ver_filename], stdout=subprocess.PIPE).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_av = stop - start

print('Running AIAA-2006-6753 verification mode...')
start = datetime.datetime.now()
print(subprocess.Popen([sgp4bin, 'i', 'v', ver_filename], stdout=subprocess.PIPE).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_iv = stop - start

print('Running libsgp4ansi verification mode...')
start = datetime.datetime.now()
print(subprocess.Popen(['ansi.exe', 'v', ver_filename], stdout=subprocess.PIPE).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_ansic = stop - start

print('Running Sattrack Report #3 full catalogue mode...')
start = datetime.datetime.now()
print(subprocess.Popen([sgp4bin, 'a', 'c', fullcat_filename], stdout=subprocess.PIPE).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_ac = stop - start

print('Running AIAA-2006-6753 full catalogue mode...')
start = datetime.datetime.now()
print(subprocess.Popen([sgp4bin, 'i', 'c', fullcat_filename], stdout=subprocess.PIPE).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_ic = stop - start

print('Running libsgp4ansi full catalogue mode...')
start = datetime.datetime.now()
print(subprocess.Popen(['ansi.exe', 'c', fullcat_filename], stdout=subprocess.PIPE).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_ansic = stop - start

print('Running Sattrack Report #3 timing mode...')
start = datetime.datetime.now()
print(subprocess.Popen([sgp4bin, 'a', 't'], stdout=subprocess.PIPE).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_at = stop - start
print('Done in {}s'.format(tdiff_at.total_seconds()))

print('Running AIAA-2006-6753 timing mode...')
start = datetime.datetime.now()
print(subprocess.Popen([sgp4bin, 'i', 't'], stdout=subprocess.PIPE).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_it = stop - start
print('Done in {}s'.format(tdiff_it.total_seconds()))

print('Running libsgp4ansi timing mode...')
start = datetime.datetime.now()
print(subprocess.Popen(['ansi.exe', 't'], stdout=subprocess.PIPE).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_ansit = stop - start
print('Done in {}s'.format(tdiff_ansit.total_seconds()))

print('Loading full catalogue data...')
start = datetime.datetime.now()
with open('a721.out') as f:
  ref = f.readlines()
with open('i72.out') as f:
  diff1 = f.readlines()
with open('ansi.out') as f:
  diff2 = f.readlines()

if ((len(ref) != len(diff1)) or (len(ref) != len(diff2))):
  print('[ERROR] Data files differ in length! - Please check:')
  print('        a721.out is {} lines'.format(len(ref)))
  print('        i72.out is  {} lines'.format(len(diff1)))
  print('        ansi.out is {} lines'.format(len(diff2)))
#  sys.exit(0)

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
      dotsp1[current_sat + str(orbital_period)][0].append(orbital_period)
      dotsp1[current_sat + str(orbital_period)][1].append(1000 * delta1) # get back to meters
      dotsp2[current_sat + str(orbital_period)][0].append(orbital_period)
      dotsp2[current_sat + str(orbital_period)][1].append(1000 * delta2)
      print('[{0}] Sat {1}, per:{2:8.3f}, '.format(i, current_sat, orbital_period), end='')
      print('d1: {0:10.9f}, d2: {1:10.9f}\r'.format(1000 * delta1, 1000 * delta2), end='')
    current_sat = m.group(1)
    orbital_period = float(m.group(2))
    dotsp1[current_sat + str(orbital_period)] = ([], [])
    dotsp2[current_sat + str(orbital_period)] = ([], [])
    delta1 = 0.0
    delta2 = 0.0
  else:
    try:
      ref_vals = list(map(float, line.split()))
      diff1_vals = list(map(float, diff1[i].split()))
      diff2_vals = list(map(float, diff2[i].split()))
    except ValueError as e:
      print('\n[ERROR]', i, current_sat, str(orbital_period) + '\n' + str(e))
      sys.exit(0)
    tempdelta1 = max(abs(ref_vals[1] - diff1_vals[1]), abs(ref_vals[2] - diff1_vals[2]), abs(ref_vals[3] - diff1_vals[3]))
    if delta1 < tempdelta1:
      delta1 = tempdelta1
    tempdelta2 = max(abs(ref_vals[1] - diff2_vals[1]), abs(ref_vals[2] - diff2_vals[2]), abs(ref_vals[3] - diff2_vals[3]))
    if delta2 < tempdelta2:
      delta2 = tempdelta2

# Repeat for last sat
dotsp1[current_sat + str(orbital_period)][0].append(orbital_period)
dotsp1[current_sat + str(orbital_period)][1].append(1000 * delta1) # get back to meters
dotsp2[current_sat + str(orbital_period)][0].append(orbital_period)
dotsp2[current_sat + str(orbital_period)][1].append(1000 * delta2)

stop = datetime.datetime.now()
tdiff_read = stop - start

print('[EOF] Sat {0}, per:{1:8.3f}, '.format(current_sat, orbital_period), end='')
print('d1: {0:10.9f}, d2: {1:10.9f}'.format(1000 * delta1, 1000 * delta2))
print('Loaded {} data points per set'.format(len(dotsp1.keys())))

print('Plotting data - please wait...')
start = datetime.datetime.now()
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

stop = datetime.datetime.now()
tdiff_plot = stop - start

plt.show()