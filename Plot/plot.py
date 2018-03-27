#!/bin/env python
# Plot differences between sets of TEME vectors

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import re
import os, sys
import subprocess
import datetime
import csv
from math import degrees, isnan

# Configuration ###############################################################
dry_run = False
elev_only = False
if (len(sys.argv) > 1):
  if (sys.argv[1] == 'd'):
    dry_run = True
  elif (sys.argv[1] == 'e'):
    elev_only = True

aiaa_binary = 'aiaa'
ansi_binary = 'ansi'

if os.name == 'nt':
  print("Assuming Windows environment")
  aiaa_binary += '.exe'
  ansi_binary += '.exe'
if os.name == 'posix':
  print("Assuming Linux (posix) environment")

ver_filename = 'SGP4-VER.TLE'
fullcat_filename = 'full.tle'

aiaa_out_filename = 'i72.out'
ansi_out_filename = 'ansi.out'

sat_re = re.compile(r'([0-9]{1,5})\ \(([0-9. e+]+)\)')

tdiff_it    = datetime.timedelta()
tdiff_ansit = datetime.timedelta()
tdiff_ic    = datetime.timedelta()
tdiff_ansic = datetime.timedelta()

# Test runs ###################################################################

def elev():
  print('Running libsgp4ansi in elevations generation mode...')
  print(subprocess.Popen([os.path.join(os.getcwd(), ansi_binary), 'p'],
        stdout=subprocess.PIPE).communicate()[0].decode())

  el_x = []
  el_y = []
  # el_a = []

  with open('elevations.out') as f:
    elevations_data = f.readlines()

  for line in elevations_data:
    el_x.append(int(line.split(',')[0]))
    el_y.append(float(line.split(',')[1]))
    # el_a.append(float(line.split(',')[2]))
    
  fig, ax = plt.subplots(1, 1)
  fig.canvas.set_window_title('Elevation') 

  xmin = 0
  xmax = el_x[len(el_x) - 1]
  xticksmin = range(xmin, xmax + 1, 1440 * 60)
  xticksmaj = range(xmin, xmax + 1, 360 * 60)

  yticksmin = [-75, -60, -45, -30, -15, 15, 30, 45, 60, 75]
  yticksmaj = [0]

  ax.plot(el_x, el_y, c='b', marker='', ls='-', lw=0.5)
  # ax.plot(el_x, el_a, c='r', marker='', ls='-', lw=0.5)
    
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(-90, 90)
  ax.set_xticks(xticksmin, True);
  ax.set_xticks(xticksmaj, False);
  ax.grid(which='major', axis='x', ls='dashed', alpha=0.7)
  ax.grid(which='minor', axis='x', ls='dotted', alpha=0.3)
  ax.set_yticks(yticksmin, True);
  ax.set_yticks(yticksmaj, False);
  ax.grid(which='major', axis='y', ls='-', alpha=1)
  ax.grid(which='minor', axis='y', ls='dotted', alpha=0.5)
  ax.set_title('Elevation of test satellite-observer pair')
  ax.set_xlabel('Time, min')
  ax.set_ylabel('Elevation, deg')
  fig.set_tight_layout({'pad':0.0, 'w_pad':0.1, 'h_pad':0.1})

  plt.show()
  
  exit(0)

def verif():
  print('Running AIAA-2006-6753 verification mode...')
  print(subprocess.Popen([os.path.join(os.getcwd(), aiaa_binary), 'i', 'v', ver_filename],
        stdout=subprocess.PIPE).communicate()[0].decode())
  print('Running libsgp4ansi verification mode...')
  print(subprocess.Popen([os.path.join(os.getcwd(), ansi_binary), 'v', ver_filename],
        stdout=subprocess.PIPE).communicate()[0].decode())
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
      print('AIAA: [{0}] {1:7}/{2:7} lines, {3:5} satellites\r'.format(
            '=' * progress + ' ' * (20 - progress), i, len(aiaa_ver_file) - 1, 
            len(aiaa_ver_sats.keys())), end='') 

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
      print('ANSI: [{0}] {1:7}/{2:7} lines, {3:5} satellites\r'.format(
            '=' * progress + ' ' * (20 - progress), i, len(ansi_ver_file) - 1,
            len(ansi_ver_sats.keys())), end='') 

  print('')

  verif_points = {}
  verif_annot = []

  for sat in list(aiaa_ver_sats.keys()):
    maxdelta = 0
    for time in aiaa_ver_sats[sat][2].keys():
      try:
        delta = max(abs(aiaa_ver_sats[sat][2][time][0] - ansi_ver_sats[sat][2][time][0]),
        abs(aiaa_ver_sats[sat][2][time][1] - ansi_ver_sats[sat][2][time][1]),
        abs(aiaa_ver_sats[sat][2][time][2] - ansi_ver_sats[sat][2][time][2]))
      except:
        delta = 0
      if delta > maxdelta:
        maxdelta = delta
    verif_points[sat] = maxdelta * 1000 # Convert to meters
    verif_annot.append(sat)

  # Max difference in verification ephemerides ##################################
  fig, axo = plt.subplots(1, 1)
  fig.canvas.set_window_title('Verif. ephemeris difference') 

  for i, sat in enumerate(list(aiaa_ver_sats.keys())):
    axo.plot(i, verif_points[sat], c='r', marker='o', ls='', mew=0.0, markersize=8, alpha=0.7)
    axo.annotate(sat, (i, verif_points[sat]))

  axo.set_xlim(-1, len(aiaa_ver_sats.keys()) + 5)
  axo.set_yscale('log')
  axo.autoscale(True, 'y', None)
  axo.grid(which='major', axis='both', ls='dashed', alpha=0.7)
  axo.grid(which='minor', axis='x', ls='dotted', alpha=0.3)
  axo.set_title('Max difference in ephemerides: AIAA (zero) vs ANSI (red)')
  axo.set_xlabel('Satellite, number')
  axo.set_ylabel('Difference, m')
  fig.set_tight_layout({'pad':0.0, 'w_pad':0.1, 'h_pad':0.1})

  # Solving Kepler's Equation for satellite 23333 ###############################
  aiaa_23333_t = []
  aiaa_23333_i = []
  ansi_23333_i = []
  aiaa_23333_m = []
  ansi_23333_m = []

  for time in aiaa_ver_sats['23333'][2].keys():
    aiaa_23333_t.append(time)
    aiaa_23333_i.append(aiaa_ver_sats['23333'][2][time][5])
    aiaa_23333_m.append(aiaa_ver_sats['23333'][2][time][9])
    ansi_23333_i.append(ansi_ver_sats['23333'][2][time][5])
    ansi_23333_m.append(ansi_ver_sats['23333'][2][time][9])

  fig, (axi, axm) = plt.subplots(2, 1, sharex=True)
  fig.canvas.set_window_title('Kepler\'s Equation for 23333') 

  xmin = 240
  xmax = 720
  xticksmin = range(xmin, xmax + 1, 60)
  xticksmaj = range(xmin, xmax + 1, 120)

  y2min = 0
  y2max = 360
  y2ticksmin = [30, 60, 210, 240, 300, 330]
  y2ticksmaj = [90, 180, 270, 360]

  axi.plot(aiaa_23333_t, aiaa_23333_i, c='k', marker='', ls=':')
  axi.plot(aiaa_23333_t, ansi_23333_i, c='b', marker='', ls='-', alpha=0.6)

  axm.plot(aiaa_23333_t, aiaa_23333_m, c='k', marker='', ls=':')
  axm.plot(aiaa_23333_t, ansi_23333_m, c='g', marker='', ls='-', alpha=0.6)

  axi.set_xlim(xmin, xmax)
  axi.set_xticks(xticksmin, True);
  axi.set_xticks(xticksmaj, False);
  axi.autoscale(True, 'y', None)
  axi.grid(which='major', axis='both', ls='dashed', alpha=0.7)
  axi.grid(which='minor', axis='x', ls='dotted', alpha=0.3)
  axi.set_title('Kepler\'s for 23333, AIAA (black dashed) vs ANSI (color)')
  axi.set_ylabel('Inclination, deg')

  axm.set_xticks(xticksmin, True);
  axm.set_xticks(xticksmaj, False);
  axm.set_ylim(y2min, y2max)
  axm.set_yticks(y2ticksmin, True);
  axm.set_yticks(y2ticksmaj, False);
  axm.grid(which='major', axis='both', ls='dashed', alpha=0.7)
  axm.grid(which='minor', axis='x', ls='dotted', alpha=0.3)
  axm.set_xlabel('Time from epoch, min')
  axm.set_ylabel('Mean anomaly, deg')
  fig.set_tight_layout({'pad':0.0, 'w_pad':0.1, 'h_pad':0.1})
    
  # Lunar-solar modifications for Satellite 23599 ###############################
  aiaa_argper = []
  ansi_argper = []

  aiaa_pos23599_x1 = []
  aiaa_pos23599_yx = []
  aiaa_pos23599_yy = []
  aiaa_pos23599_yz = []

  ansi_pos23599_yx = []
  ansi_pos23599_yy = []
  ansi_pos23599_yz = []

  for time in aiaa_ver_sats['23599'][2].keys():
    aiaa_pos23599_x1.append(time)
    aiaa_argper.append(aiaa_ver_sats['23599'][2][time][7])
    ansi_argper.append(ansi_ver_sats['23599'][2][time][7])
    aiaa_pos23599_yx.append(aiaa_ver_sats['23599'][2][time][0])
    aiaa_pos23599_yy.append(aiaa_ver_sats['23599'][2][time][1])
    aiaa_pos23599_yz.append(aiaa_ver_sats['23599'][2][time][2])
    ansi_pos23599_yx.append(ansi_ver_sats['23599'][2][time][0])
    ansi_pos23599_yy.append(ansi_ver_sats['23599'][2][time][1])
    ansi_pos23599_yz.append(ansi_ver_sats['23599'][2][time][2])

  fig, (axa, axp) = plt.subplots(2, 1, sharex=True)
  fig.canvas.set_window_title('L-S mods for 23599') 

  xmin = int(min(aiaa_pos23599_x1))
  xmax = int(max(aiaa_pos23599_x1))
  xticksmin = range(xmin, xmax + 1, 20)
  xticksmaj = range(xmin, xmax + 1, 60)

  axa.plot(aiaa_pos23599_x1, aiaa_argper, c='k', marker='', ls=':')
  axa.plot(aiaa_pos23599_x1, ansi_argper, c='r', marker='', ls='-', alpha=0.6)

  axa.set_xlim(xmin, xmax)
  axa.set_xticks(xticksmin, True);
  axa.set_xticks(xticksmaj, False);
  axa.autoscale(True, 'y', None)
  axa.grid(which='major', axis='both', ls='dashed', alpha=0.7)
  axa.grid(which='minor', axis='x', ls='dotted', alpha=0.3)
  axa.set_title('L-S mod. for 23599, AIAA (black dashed) vs ANSI (color)')
  axa.set_ylabel('Argument of perigee, deg')

  axp.plot(aiaa_pos23599_x1, aiaa_pos23599_yx, c='k', marker='', ls=':')
  axp.plot(aiaa_pos23599_x1, aiaa_pos23599_yy, c='k', marker='', ls=':')
  axp.plot(aiaa_pos23599_x1, aiaa_pos23599_yz, c='k', marker='', ls=':')

  axp.plot(aiaa_pos23599_x1, ansi_pos23599_yx, c='b', marker='', ls='-', alpha=0.6)
  axp.plot(aiaa_pos23599_x1, ansi_pos23599_yy, c='g', marker='', ls='-', alpha=0.6)
  axp.plot(aiaa_pos23599_x1, ansi_pos23599_yz, c='m', marker='', ls='-', alpha=0.6)

  axp.set_xticks(xticksmin, True);
  axp.set_xticks(xticksmaj, False);
  axp.autoscale(True, 'y', None)
  axp.grid(which='major', axis='both', ls='dashed', alpha=0.7)
  axp.grid(which='minor', axis='x', ls='dotted', alpha=0.3)
  axp.set_xlabel('Time from epoch, min')
  axp.set_ylabel('Position components, km')
  fig.set_tight_layout({'pad':0.0, 'w_pad':0.1, 'h_pad':0.1})

  # Lyddane Choice Modification for Satellite 14128 #############################

  incl_bar = []
  aiaa_14128_t = []
  aiaa_14128_i = []
  ansi_14128_i = []
  dpy_14128 = []
  dpz_14128 = []

  for time in aiaa_ver_sats['14128'][2].keys():
    incl_bar.append(11.4592)
    aiaa_14128_t.append(time)
    aiaa_14128_i.append(aiaa_ver_sats['14128'][2][time][5])
    ansi_14128_i.append(ansi_ver_sats['14128'][2][time][5])
    dpy_14128.append(aiaa_ver_sats['14128'][2][time][1] - ansi_ver_sats['14128'][2][time][1])
    dpz_14128.append(aiaa_ver_sats['14128'][2][time][2] - ansi_ver_sats['14128'][2][time][2])

  fig, (axi, axp) = plt.subplots(2, 1, sharex=True)
  fig.canvas.set_window_title('Lyddane Choice for 14128') 

  xmin = 0
  xmax = 2880
  xticksmin = range(xmin, xmax + 1, 60)
  xticksmaj = range(xmin, xmax + 1, 240)

  axi.plot(aiaa_14128_t, incl_bar, c='k', marker='', ls='-', lw=0.5)

  axi.plot(aiaa_14128_t, aiaa_14128_i, c='k', marker='', ls=':')
  axi.plot(aiaa_14128_t, ansi_14128_i, c='b', marker='', ls='-', alpha=0.6)

  axp.plot(aiaa_14128_t, dpy_14128, c='g', marker='', ls='-')
  axp.plot(aiaa_14128_t, dpz_14128, c='m', marker='', ls='-')

  axi.set_xlim(xmin, xmax)
  axi.set_xticks(xticksmin, True);
  axi.set_xticks(xticksmaj, False);
  axi.autoscale(True, 'y', None)
  axi.grid(which='major', axis='both', ls='dashed', alpha=0.7)
  axi.grid(which='minor', axis='x', ls='dotted', alpha=0.3)
  axi.set_title('Lyddane Choice for 14128, AIAA (black dashed) vs ANSI (color)')
  axi.set_ylabel('Inclination, deg')

  axp.set_xticks(xticksmin, True);
  axp.set_xticks(xticksmaj, False);
  axp.set_ylim(y2min, y2max)
  axp.set_yticks(y2ticksmin, True);
  axp.set_yticks(y2ticksmaj, False);
  axp.grid(which='major', axis='both', ls='dashed', alpha=0.7)
  axp.grid(which='minor', axis='x', ls='dotted', alpha=0.3)
  axp.set_xlabel('Time from epoch, min')
  axp.set_ylabel('Position difference, m')
  fig.set_tight_layout({'pad':0.0, 'w_pad':0.1, 'h_pad':0.1}) 

def fullcat():
  print('Running AIAA-2006-6753 full catalogue mode...')
  start = datetime.datetime.now()
  print(subprocess.Popen([os.path.join(os.getcwd(), aiaa_binary), 'i', 'c', fullcat_filename],
        stdout=subprocess.PIPE).communicate()[0].decode())
  stop = datetime.datetime.now()
  global tdiff_ic
  tdiff_ic = stop - start
  print('Running libsgp4ansi full catalogue mode...')
  start = datetime.datetime.now()
  print(subprocess.Popen([os.path.join(os.getcwd(), ansi_binary), 'c', fullcat_filename],
        stdout=subprocess.PIPE).communicate()[0].decode())
  stop = datetime.datetime.now()
  global tdiff_ansic
  tdiff_ansic = stop - start
  print('Loading full catalogue data...')
  with open(aiaa_out_filename) as f:
    aiaa_file = f.readlines()
  with open(ansi_out_filename) as f:
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
      print('AIAA: [{0}] {1:7}/{2:7} lines, {3:5} satellites\r'.format(
            '=' * progress + ' ' * (20 - progress), i, len(aiaa_file) - 1,
            len(aiaa_sats.keys())), end='') 

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
      print('ANSI: [{0}] {1:7}/{2:7} lines, {3:5} satellites\r'.format(
            '=' * progress + ' ' * (20 - progress), i, len(ansi_file) - 1,
            len(ansi_sats.keys())), end='') 

  print('')

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
        delta = max(abs(aiaa_sats[sat][2][time][0] - ansi_sats[sat][2][time][0]),
        abs(aiaa_sats[sat][2][time][1] - ansi_sats[sat][2][time][1]),
        abs(aiaa_sats[sat][2][time][2] - ansi_sats[sat][2][time][2]))
        if delta > maxdelta:
          maxdelta = delta
    fullcat_points[sat] = maxdelta * 1000 # Convert to meters
  for sat in list(ansi_sats.keys()):
    for time in ansi_sats[sat][2].keys():
      if time not in aiaa_sats[sat][2].keys():
        print('[   ->] AIAA list missing data at t={0:7} for satellite {1}'.format(time, sat))
        
  # Max difference in ephemerides between maths #################################
  fig, ax = plt.subplots(1, 1)
  fig.canvas.set_window_title('Max. diff in ephems') 

  xmin = 0
  xmax = 2000
  xticksmin = range(xmin, xmax + 1, 100)
  xticksmaj = range(xmin, xmax + 1, 500)

  for sat in list(aiaa_sats.keys()):
    ax.plot(aiaa_sats[sat][0], fullcat_points[sat], c='b', marker='D', ls='',
            mew=0.0, markersize=3, alpha=0.55)
    
  ax.set_xlim(xmin, xmax)
  ax.set_xticks(xticksmin, True);
  ax.set_xticks(xticksmaj, False);
  ax.set_yscale('log')
  ax.autoscale(True, 'y', None)
  ax.grid(which='major', axis='both', ls='dashed', alpha=0.7)
  ax.grid(which='minor', axis='x', ls='dotted', alpha=0.3)
  ax.set_title('Max difference in ephemerides: AIAA (zero) vs ANSI (blue)')
  ax.set_xlabel('Orbital period, min')
  ax.set_ylabel('Difference, m')
  fig.set_tight_layout({'pad':0.0, 'w_pad':0.1, 'h_pad':0.1})
  
def timing():
  print('Running AIAA-2006-6753 timing mode...')
  start = datetime.datetime.now()
  print(subprocess.Popen([os.path.join(os.getcwd(), aiaa_binary), 'i', 't'],
        stdout=subprocess.PIPE).communicate()[0].decode())
  stop = datetime.datetime.now()
  global tdiff_it
  tdiff_it = stop - start
  print('Running libsgp4ansi timing mode...')
  start = datetime.datetime.now()
  print(subprocess.Popen([os.path.join(os.getcwd(), ansi_binary), 't'],
        stdout=subprocess.PIPE).communicate()[0].decode())
  stop = datetime.datetime.now()
  global tdiff_ansit
  tdiff_ansit  = stop - start
  
  # Runtime difference between static reference and library code ################
  fig, (axtt, axtc) = plt.subplots(2, 1)
  fig.canvas.set_window_title('Timing results') 
  
  axtt.barh([1.9], [tdiff_it.total_seconds()], height=0.5, color='g')
  axtt.barh([1.1], [tdiff_ansit.total_seconds()], height=0.5, color='b')
  axtt.autoscale(True, 'x', None)
  axtt.set_yticks([], True)
  axtt.set_yticks([], False)
  axtt.set_title('Timing')
  axtt.set_ylabel('Synth. test')
  axtt.text(tdiff_it.total_seconds() / 2, 1.85, 'AIAA', ha='center', color='w', weight='bold')
  axtt.text(tdiff_ansit.total_seconds() / 2, 1.05, 'ANSI', ha='center', color='w', weight='bold')
  
  axtc.barh([1.9], [tdiff_ic.total_seconds()], height=0.5, color='g')
  axtc.barh([1.1], [tdiff_ansic.total_seconds()], height=0.5, color='b')
  axtc.autoscale(True, 'x', None)
  axtc.set_yticks([], True)
  axtc.set_yticks([], False)
  axtc.set_xlabel('Run time, s')
  axtc.set_ylabel('Real run')
  axtc.text(tdiff_ic.total_seconds() / 2, 1.85, 'AIAA', ha='center', color='w', weight='bold')
  axtc.text(tdiff_ansic.total_seconds() / 2, 1.05, 'ANSI', ha='center', color='w', weight='bold')
  
  fig.set_tight_layout({'pad':0.0, 'w_pad':0.1, 'h_pad':0.1})
  
def shadow():
  print('Running libsgp4ansi shadow plotting mode...')
  print(subprocess.Popen([os.path.join(os.getcwd(), ansi_binary), 's'],
        stdout=subprocess.PIPE).communicate()[0].decode())
  print('Loading shadow data...')

  sh_sun_az = []
  sh_sun_el = []
  sh_moon_az = []
  sh_moon_el = []
  sh_sat_az = []
  sh_sat_el = []
  sh_sun_lat = []
  sh_sun_lon = []
  sh_moon_lat = []
  sh_moon_lon = []
  sh_sat_lat = []
  sh_sat_lon = []
  sun_lat = []
  sun_lon = []
  moon_lat = []
  moon_lon = []
  with open('ansi_shadows.csv', newline='') as csvfile:
      csvreader = csv.reader(csvfile, lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC)
      for row in csvreader:
          sh_sun_az.append(row[1])
          sh_sun_el.append(90 - degrees(row[2]))
          sh_moon_az.append(row[3])
          sh_moon_el.append(90 - degrees(row[4]))
          sh_sat_az.append(row[5])
          sh_sat_el.append(90 - degrees(row[6]))
          sh_sun_lat.append(degrees(row[7]))
          sh_sun_lon.append(degrees(row[8]))
          sh_moon_lat.append(degrees(row[9] if not isnan(row[9]) else -100))
          sh_moon_lon.append(degrees(row[10] if not isnan(row[10]) else 0))
          sh_sat_lat.append(degrees(row[11]))
          sh_sat_lon.append(degrees(row[12]))
          sun_lat.append(degrees(row[13]))
          sun_lon.append(degrees(row[14]))
          moon_lat.append(degrees(row[15]))
          moon_lon.append(degrees(row[16]))
      print('Loaded {} shadow points'.format(len(sh_sat_az)))
  # Shadow map plot #############################################################
  fig, axs = plt.subplots(1, 1)
  fig.canvas.set_window_title('Solar and lunar shadows')

  img = plt.imread('map.png')
  axs.imshow(img, extent=[-180, 180, -90, 90])

  xmin = -180
  xmax = 180
  xticksmin = range(xmin, xmax + 1, 10)
  xticksmaj = range(xmin, xmax + 1, 30)
  ymin = -90
  ymax = 90
  yticksmin = range(ymin, ymax + 1, 10)
  yticksmaj = range(ymin, ymax + 1, 30)

  sat_track, = axs.plot(sh_sat_lon, sh_sat_lat, c='r', marker='', ls='-', lw=2, alpha=0.8)
  sh_moon_l, = axs.plot(sh_moon_lon, sh_moon_lat, c='k', marker='', ls='-', lw=4, alpha=0.3)
  sh_moon_m  = axs.scatter(sh_moon_lon, sh_moon_lat, edgecolors='b', s=16, facecolors='none', alpha=1)
  sh_sun_l,  = axs.plot(sh_sun_lon, sh_sun_lat, c='k', marker='', ls='-', lw=4, alpha=0.3)
  sh_sun_m   = axs.scatter(sh_sun_lon, sh_sun_lat, edgecolors='#FFFF00', s=16, facecolors='none', alpha=1)
  moon_m     = axs.scatter(moon_lon, moon_lat, edgecolors='none', facecolors='b', s=49, alpha=1)
  sun_m      = axs.scatter(sun_lon, sun_lat, edgecolors='none', facecolors='#FFFF00', s=49, alpha=1)

  axs.legend([sat_track, (sh_moon_l, sh_moon_m), (sh_sun_l, sh_sun_m), sun_m, moon_m],
             ['ISS ground track', 'Satellite lunar shadow', 'Satellite solar shadow', 'Sun', 'Moon'])

  # Banner, WY
  coords = (-106.865443, 44.601882)
  axs.plot(coords[0], coords[1], c='w', marker='P', ls='', mew=0.0, markersize=5, alpha=1)
  axs.annotate('Banner, WY', xy=coords, color='white', horizontalalignment='right', verticalalignment='bottom',)
    
  axs.set_xticks(xticksmin, True);
  axs.set_xticks(xticksmaj, False);
  axs.set_yticks(yticksmin, True);
  axs.set_yticks(yticksmaj, False);
  axs.set_xlim(xmin, xmax)
  axs.set_ylim(ymin, ymax)
  axs.grid(which='major', axis='both', ls='dashed', alpha=0.7)
  axs.grid(which='minor', axis='both', ls='dotted', alpha=0.3)
  axs.set_title('Solar eclipse 2017 scenario')
  axs.set_xlabel('Longitude')
  axs.set_ylabel('Latitude')
  fig.set_tight_layout({'pad':0.0, 'w_pad':0.1, 'h_pad':0.1})

  # Sky radar plot ##############################################################
  fig, axr = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})
  fig.canvas.set_window_title('Sky in Banner, WY')

  axr.set_theta_zero_location('N')
  axr.set_rlim(0, 90, 1)
  axr.set_yticklabels([])
  axr.set_facecolor('#ddf3ff')

  axr.xaxis.grid(True,color='k',linestyle=':')
  axr.yaxis.grid(True,color='k',linestyle=':')

  sat_track, = axr.plot(sh_sat_az, sh_sat_el, c='r', marker='', ls='-', lw=1, alpha=0.8)
  moon_m     = axr.scatter(sh_moon_az, sh_moon_el, edgecolors='b', s=24, facecolors='none', alpha=1)
  sun_m      = axr.scatter(sh_sun_az, sh_sun_el, edgecolors='#FFFF00', s=24, facecolors='none', alpha=1)

  axr.legend([sat_track, sun_m, moon_m], ['ISS', 'The Sun', 'The Moon'])
  axr.set_title('Sky radar as seen from Banner, WY')


if __name__ == '__main__':
  if elev_only:
    elev()
  if not dry_run:
    # verif()
    # fullcat()
    # timing()
    shadow()
    print('Plotting data - please wait...')
    plt.show()
