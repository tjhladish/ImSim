#!/usr/bin/python
from sys import argv, exit

infile = argv[1]
outfile = argv[2]

data = dict()

for line in file(infile):
    parts = line.strip().split('\t')
    if len(parts) != 3:
        print "Something's wrong with input! Expected 3 fields per line, found", len(parts), ":", line
        exit()
    year, loc, incidence = parts
    if loc not in data:
        data[loc] = dict()
    data[loc][year] = str(float(incidence) / 100.0) #### Currently, incidence data is counts per 100

fo = open(outfile, 'w')

for loc in sorted(data.keys()):
    for year in sorted(data[loc].keys()):
        fo.write( ','.join([loc, year, data[loc][year]]) + '\n' )

fo.close()
