#!/bin/env python
# Plot differences between sets of TEME vectors

import matplotlib.pyplot as plt
import re
import os, sys
import subprocess
import datetime

######################################################

str3_binary = 'str3'
aiaa_binary = 'sgp4'

if os.name == 'nt':
  print("Assuming Windows environment")
  aiaa_binary += '.exe'
if os.name == 'posix':
  print("Assuming Linux (posix) environment")

ver_filename = 'SGP4-VER.TLE'
fullcat_filename = 'full.tle'

sat_re = re.compile('([0-9]{1,5})\ \(([0-9. ]+)\)')

######################################################

# print('Running AIAA-2006-6753 timing mode...')
# start = datetime.datetime.now()
# print(subprocess.Popen([str3_binary, 't'], stdout=subprocess.PIPE).communicate()[0].decode())
# stop = datetime.datetime.now()
# tdiff_it = stop - start

print('Running AIAA-2006-6753 timing mode...')
start = datetime.datetime.now()
# print(subprocess.Popen([aiaa_binary, 'i', 't'], stdout=subprocess.PIPE).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_it = stop - start

print('Running libsgp4ansi timing mode...')
start = datetime.datetime.now()
# print(subprocess.Popen(['ansi.exe', 't'], stdout=subprocess.PIPE).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_ansit = stop - start

# print('Running SatTrack Report #3 verification mode...')
# print(subprocess.Popen([str3_binary, 'v', ver_filename], stdout=subprocess.PIPE).communicate()[0].decode())
print('Running AIAA-2006-6753 verification mode...')
print(subprocess.Popen([aiaa_binary, 'i', 'v', ver_filename], stdout=subprocess.PIPE).communicate()[0].decode())
print('Running libsgp4ansi verification mode...')
print(subprocess.Popen(['ansi.exe', 'v', ver_filename], stdout=subprocess.PIPE).communicate()[0].decode())

# print('Running SatTrack Report #3 full catalogue mode...')
# start = datetime.datetime.now()
# print(subprocess.Popen([aiaa_binary, 'a', 'c', fullcat_filename], stdout=subprocess.PIPE).communicate()[0].decode())
# stop = datetime.datetime.now()
# tdiff_ac = stop - start

print('Running AIAA-2006-6753 full catalogue mode...')
start = datetime.datetime.now()
print(subprocess.Popen([aiaa_binary, 'i', 'c', fullcat_filename], stdout=subprocess.PIPE).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_ic = stop - start

print('Running libsgp4ansi full catalogue mode...')
start = datetime.datetime.now()
print(subprocess.Popen(['ansi.exe', 'c', fullcat_filename], stdout=subprocess.PIPE).communicate()[0].decode())
stop = datetime.datetime.now()
tdiff_ansic = stop - start

######################################################

print('Loading full catalogue data...')

with open('i72.out') as f:
  aiaa_file = f.readlines()
with open('ansi.out') as f:
  ansi_file = f.readlines()

aiaa_sats = {}
ansi_sats = {}

for i, line in enumerate(aiaa_file):
  m = sat_re.findall(line.strip())
  if m != []:
    current_sat = m[0][0]
    aiaa_sats[current_sat] = [float(m[0][1]), 0.0, {}]
  else:
    try:
      vals = list(map(float, line.split()))
    except ValueError as e:
      print('\n[{0}] {1}'.format(i, str(e)))
    aiaa_sats[current_sat][2][vals[0]] = vals[1:]
  if ((i % 1000 == 0) | (i == len(aiaa_file) - 1)):
    progress = int(i / (len(aiaa_file) - 1) * 20)
    print('AIAA: [{0}] {1:7}/{2:7} lines, {3:5} satellites\r'.format('=' * progress + ' ' * (20 - progress), i, len(aiaa_file) - 1, len(aiaa_sats.keys())), end='') 

print('')
    
for i, line in enumerate(ansi_file):
  m = sat_re.findall(line.strip())
  if m != []:
    current_sat = m[0][0]
    ansi_sats[current_sat] = [float(m[0][1]), 0.0, {}]
  else:
    try:
      vals = list(map(float, line.split()))
    except ValueError as e:
      print('\n[{0}] {1}'.format(i, str(e)))
    ansi_sats[current_sat][2][vals[0]] = vals[1:]
  if ((i % 1000 == 0) | (i == len(ansi_file) - 1)):
    progress = int(i / (len(ansi_file) - 1) * 20)
    print('ANSI: [{0}] {1:7}/{2:7} lines, {3:5} satellites\r'.format('=' * progress + ' ' * (20 - progress), i, len(ansi_file) - 1, len(ansi_sats.keys())), end='') 

print('')

######################################################

if len(aiaa_sats.keys()) != len(ansi_sats.keys()):
  print('[WARNG] Satellite count discrepancy in full cat. data!')
  for sat in list(aiaa_sats.keys()):
    if sat not in list(ansi_sats.keys()):
      print('[   ->] ANSI list missing satellite {}'.format(sat))
  for sat in list(ansi_sats.keys()):
    if sat not in list(aiaa_sats.keys()):
      print('[   ->] AIAA list missing satellite {}'.format(sat))

if len(aiaa_file) != len(ansi_file):
  print('[WARNG] Data points discrepancy in full cat. data!')

fullcat_points = {}
  
for sat in list(aiaa_sats.keys()):
  maxdelta = 0
  for time in aiaa_sats[sat][2].keys():
    if time not in ansi_sats[sat][2].keys():
      print('[   ->] ANSI list missing data at t={0:7} for satellite {1}'.format(time, sat))
    else:
      delta = max(abs(aiaa_sats[sat][2][time][0] - ansi_sats[sat][2][time][0]), abs(aiaa_sats[sat][2][time][1] - ansi_sats[sat][2][time][1]), abs(aiaa_sats[sat][2][time][2] - ansi_sats[sat][2][time][2]))
      if delta > maxdelta:
        maxdelta = delta
  fullcat_points[sat] = maxdelta * 1000 # Convert to meters
for sat in list(ansi_sats.keys()):
  maxdelta = 0
  for time in ansi_sats[sat][2].keys():
    if time not in aiaa_sats[sat][2].keys():
      print('[   ->] AIAA list missing data at t={0:7} for satellite {1}'.format(time, sat))

######################################################

print('Loading verification data...')

with open('i72_ver.out') as f:
  aiaa_ver_file = f.readlines()
with open('ansi_ver.out') as f:
  ansi_ver_file = f.readlines()

aiaa_ver_sats = {}
ansi_ver_sats = {}

for i, line in enumerate(aiaa_ver_file):
  m = sat_re.findall(line.strip())
  if m != []:
    current_sat = m[0][0]
    aiaa_ver_sats[current_sat] = [float(m[0][1]), 0.0, {}]
  else:
    try:
      vals = list(map(float, line.split()))
    except ValueError as e:
      print('\n[{0}] {1}'.format(i, str(e)))
    aiaa_ver_sats[current_sat][2][vals[0]] = vals[1:]
  if ((i % 1000 == 0) | (i == len(aiaa_ver_file) - 1)):
    progress = int(i / (len(aiaa_ver_file) - 1) * 20)
    print('AIAA: [{0}] {1:7}/{2:7} lines, {3:5} satellites\r'.format('=' * progress + ' ' * (20 - progress), i, len(aiaa_ver_file) - 1, len(aiaa_ver_sats.keys())), end='') 

print('')
    
for i, line in enumerate(ansi_ver_file):
  m = sat_re.findall(line.strip())
  if m != []:
    current_sat = m[0][0]
    ansi_ver_sats[current_sat] = [float(m[0][1]), 0.0, {}]
  else:
    try:
      vals = list(map(float, line.split()))
    except ValueError as e:
      print('\n[{0}] {1}'.format(i, str(e)))
    ansi_ver_sats[current_sat][2][vals[0]] = vals[1:]
  if ((i % 1000 == 0) | (i == len(ansi_ver_file) - 1)):
    progress = int(i / (len(ansi_ver_file) - 1) * 20)
    print('ANSI: [{0}] {1:7}/{2:7} lines, {3:5} satellites\r'.format('=' * progress + ' ' * (20 - progress), i, len(ansi_ver_file) - 1, len(ansi_ver_sats.keys())), end='') 

print('')

######################################################

aiaa_pos23599_x1 = []
aiaa_pos23599_yx = []
aiaa_pos23599_yy = []
aiaa_pos23599_yz = []

ansi_pos23599_yx = []
ansi_pos23599_yy = []
ansi_pos23599_yz = []

for time in aiaa_ver_sats['23599'][2].keys():
  aiaa_pos23599_x1.append(time)
  aiaa_pos23599_yx.append(aiaa_ver_sats['23599'][2][time][0])
  aiaa_pos23599_yy.append(aiaa_ver_sats['23599'][2][time][1])
  aiaa_pos23599_yz.append(aiaa_ver_sats['23599'][2][time][2])
  ansi_pos23599_yx.append(ansi_ver_sats['23599'][2][time][0])
  ansi_pos23599_yy.append(ansi_ver_sats['23599'][2][time][1])
  ansi_pos23599_yz.append(ansi_ver_sats['23599'][2][time][2])

######################################################

# Plotting
print('Plotting data - please wait...')

# Max difference in ephemerides between maths
fig = plt.figure(1)
fig.canvas.set_window_title('Max. diff in ephems') 

xmin = 0
xmax = 2000
xticksmin = range(xmin, xmax + 1, 100)
xticksmaj = range(xmin, xmax + 1, 500)

plt.subplot(111)
for sat in list(aiaa_sats.keys()):
  plt.plot(aiaa_sats[sat][0], fullcat_points[sat], c='b', marker='D', ls='', mew=0.0, markersize=3, alpha=0.55)
  
plt.subplot(111).set_xlim(xmin, xmax)
plt.subplot(111).set_xticks(xticksmin, True);
plt.subplot(111).set_xticks(xticksmaj, False);
plt.subplot(111).set_yscale('log')
plt.subplot(111).autoscale(True, 'y', None)
plt.subplot(111).grid(which='major', axis='both', ls='dashed', alpha=0.7)
plt.subplot(111).grid(which='minor', axis='x', ls='dotted', alpha=0.3)
plt.title('Max difference in ephemerides: AIAA (zero) vs ANSI (blue)')
plt.xlabel('Orbital period, min')
plt.ylabel('Difference, m')
plt.tight_layout(pad=0.0, w_pad=0.1, h_pad=0.1)

# Runtime difference between static reference and library code
fig = plt.figure(2)
fig.canvas.set_window_title('Timing results') 

plt.subplot(111)
plt.bar([1,2], [tdiff_it.total_seconds(), tdiff_ansit.total_seconds()])
plt.subplot(111).autoscale(True, 'y', None)
plt.title('Timing')
plt.xlabel('')
plt.ylabel('Run time, s')
plt.tight_layout(pad=0.0, w_pad=0.1, h_pad=0.1)

# Lunar-solar modifications for Satellite 23599
fig = plt.figure(4)
fig.canvas.set_window_title('Lunar-solar modifications for Satellite 23599') 

xmin = int(min(aiaa_pos23599_x1))
xmax = int(max(aiaa_pos23599_x1))
xticksmin = range(xmin, xmax + 1, 20)
xticksmaj = range(xmin, xmax + 1, 60)

plt.plot(aiaa_pos23599_x1, aiaa_pos23599_yx, c='b', marker='', ls=':')
plt.plot(aiaa_pos23599_x1, aiaa_pos23599_yy, c='g', marker='', ls=':')
plt.plot(aiaa_pos23599_x1, aiaa_pos23599_yz, c='m', marker='', ls=':')

plt.plot(aiaa_pos23599_x1, ansi_pos23599_yx, c='b', marker='', ls='-', alpha=0.6)
plt.plot(aiaa_pos23599_x1, ansi_pos23599_yy, c='g', marker='', ls='-', alpha=0.6)
plt.plot(aiaa_pos23599_x1, ansi_pos23599_yz, c='m', marker='', ls='-', alpha=0.6)

plt.subplot(111).set_xlim(xmin, xmax)
plt.subplot(111).set_xticks(xticksmin, True);
plt.subplot(111).set_xticks(xticksmaj, False);
plt.subplot(111).autoscale(True, 'y', None)
plt.subplot(111).grid(which='major', axis='both', ls='dashed', alpha=0.7)
plt.subplot(111).grid(which='minor', axis='x', ls='dotted', alpha=0.3)
plt.title('Lunar-solar modifications for Satellite 23599')
plt.xlabel('Time from epoch, min')
plt.ylabel('Position components, km')
plt.tight_layout(pad=0.0, w_pad=0.1, h_pad=0.1)

plt.show()