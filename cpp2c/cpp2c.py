import re

regexp = re.compile("double\*\ ([\w]+),*")
id_re = "[_a-zA-Z][_a-zA-Z0-9]{0,30}"

with open("func.cpp", 'r') as infile:
  indata = infile.read()

vars = regexp.findall(indata)

blocks = indata.split('{', 1)

added_length = 0

for id in re.finditer(id_re, blocks[1]):
  if (id[0] in vars):
    blocks[1] = blocks[1][:id.start() + added_length] + '*' + id[0] + blocks[1][id.end() + added_length:]
    added_length += 1
    print("Replaced {} by *{}\tat\t{}, {}".format(id[0], id[0], id.start() + added_length, id.end() + added_length))

with open("func.c", 'w') as outfile:
  outfile.write(blocks[0])
  outfile.write('{')
  outfile.write(blocks[1])