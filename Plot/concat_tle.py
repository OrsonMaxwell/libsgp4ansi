#!/bin/env python
from os import listdir
from re import compile, M
from sys import argv

main_re = compile(r'''^(.{1,24})$
^1\ ([0-9]{5})([A-Z\ ])\ ([0-9A-Z ]{8})\ ([0-9]{2})([0-9\.]{12})\ ([-|\ ]\.[0-9]{8})\ ([\ |-][0-9]{5}[+|-][0-9])\ ([\ |-][0-9]{5}[+|-][0-9])\ ([0-9])\ ([\ 0-9]{4})([0-9])$
^2\ ([0-9]{5})\ ([ 0-9.]{8})\ ([ 0-9.]{8})\ ([0-9]{7})\ ([ 0-9.]{8})\ ([ 0-9.]{8})\ ([ 0-9.]{11})([ 0-9]{5})([0-9])$''', M)

sats = []

#outfile = open(argv[1] + '/' + argv[1] + '.tle', 'w')

with open('full.tle', 'r') as f:
  for match in main_re.finditer(f.read()):
    if match.group(2) not in sats:
      sats.append(match.group(2))
      print(match.group(0))
     #outfile.write(match.group(0))

#outfile.close()
