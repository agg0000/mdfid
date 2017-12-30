#!/usr/bin/env python
'''
as the main program to run the molecular dynamics
'''

import os
import re
import sys
import copy
import datetime
import numpy as np

import ic
import softapi
import writefile

start = datetime.datetime.now()
pyfile, infile = sys.argv

autoan = 0.529177211
mas = {'H' : 1823.0 , 'Li': 12761.0, 'Be': 16407.0, 'B' : 19688.4, 'I' :231521.0,
       'C' : 21876.0, 'N' : 25522.0, 'O' : 29168.0, 'F' : 34637.0, 'Ne': 36460.0,
       'Na': 41929.0, 'Mg': 44298.9, 'Al': 49221.0, 'Si': 51044.0, 'P' : 56513.0,
       'S' : 58336.0, 'Cl': 64716.5, 'Ca': 72920.0, 'K' : 71097.0, 'Br':145840.0,}

'''
the parameter is from the input file inp
the 1 row is time step
the 2 row is number of all atoms
the 3 row is the inital energy
the 4 row is the file name 
the 5 row is loop cycle, if not zero
the 6 row is convergence energy, if not zero
the 7 row is the quantum chemistry soft you needed
'''

inp = open(infile).read().splitlines()
fd={}
for keyword in inp:
	ka = keyword.split('=')
	fd[ka[0]] = ka[1]

dt = float(fd['dt'])
numele = int(fd['numele'])
numcycle = int(fd['numcycle'])
scfenry = float(fd['energy'])
filename = fd['filename']
softuse = fd['softuse']
keyname = fd['keyname']

totalmom = []
for word in fd:
	if word[0] == 'p':
		totalmom.append(fd[word])

quansoft = softapi.softname(softuse)
out0 = "%s.zzop" %filename

if keyname == 'i':
	totalv, pos0, symb, far0, mom0, intcycle = ic.initp(numele, totalmom, softuse, filename, out0, dt, start)
elif keyname == 'c':
	totalv, pos0, symb, far0, mom0, intcycle = ic.continuep(numele, out0, start)
else:
	exit('no keyname for i or c')

while True:
	pos0 = pos0 / autoan
	for i in range(numele):
		mass = float(mas[symb[i]])
		x1 = pos0[i] + mom0[i] * dt / mass + far0[i] * dt ** 2 / (2 * mass)
		pos0[i] = x1
	pos0 = pos0 * autoan

	quansoft.reinpotfile(filename, pos0, symb, intcycle, numele, filename)
	runfile = quansoft.runsoft(filename, intcycle)

	far1 = copy.deepcopy(far0)
	far0 = quansoft.getforce(filename, numele, out0, totalv, symb, start, intcycle)
	for i in range(numele):
		p1 = mom0[i] + dt * (far0[i] + far1[i]) / 2
		mom0[i] = p1

	outfile = open(out0, 'a+')
	outfile.write('\n' + '='*84 + '\n')
	outfile.write(' '*41 + 'Cycle %d' %intcycle + '\n')
	outfile.write('='*84 + '\n'*2)
	outfile.close()

	writefile.writecon(out0, symb, pos0, 'coordinate', numele)
	writefile.writecon(out0, symb, far0, 'Forces', numele)
	totalv = quansoft.getenergy(filename, out0, totalv, start, intcycle)
	writefile.writecon(out0, symb, mom0, 'momentum', numele)
	writefile.kine(mom0, symb, out0)

	if numcycle:
		if intcycle == numcycle:
			writefile.endtime(out0, 'cycle', intcycle, start)
			exit('cycle normal exit')

	if scfenry:
		if abs(totalv[-1] - totalv[-2]) < scfenry:
			writefile.endtime(out0, 'energy', intcycle, start)
			exit('energy normal exit')

	intcycle += 1
