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

import softapi
import writefile

start = datetime.datetime.now()
pyfile, inputfile = sys.argv

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

inp = open(inputfile).readlines()
dt = float(inp[0])
numele = int(inp[1])
totalmom = re.split(';|\n', inp[2])
filename = inp[3].strip()
numcycle = int(inp[4])
scfenry = float(inp[5])
softuse = inp[6].strip()

quansoft = softapi.softname(softuse)
inputfile = quansoft.getfile(filename)
out0 = "%s.out" %filename

totalv = []
totalt = []

while '' in totalmom:
	totalmom.remove('')

numbermom = {}
for mom in totalmom:
	everymom = mom.split(',')
	symbolm = int(everymom[0]) - 1
	numbermom[symbolm] = everymom[1:]

for i in range(numele):
	if i not in numbermom:
		numbermom[i] = [0, 0, 0]

originmom = []
for i in range(numele):
	originmom.append(numbermom[i])

outfile = open(out0, 'w')
outfile.write('*'*99 + '\n')
outfile.write('*'*36 + ' Welcome to zongzi program ' + '*'*36 + '\n')
outfile.write('*'*99 + '\n'*2)

outfile.write('input file name is <%s>' %inputfile + '\n')
outfile.write('step(time) is %d a.u.' %dt + '\n')
outfile.write('number of atoms need to calculate is %d' %numele + '\n'*2)

outfile.write('='*84 + '\n')
outfile.write(' '*41 + 'cycle 0' + '\n')
outfile.write('='*84 + '\n'*2)
outfile.close()

pos0, symb = quansoft.getpos(inputfile, numele)
writefile.writecon(out0, symb, pos0, 'coordinate', numele)

runfile = quansoft.runsoft(inputfile)
far0 = quansoft.getforce(runfile, numele, out0, totalv, symb, start)

mom0 = np.array(originmom, dtype = float)
writefile.writecon(out0, symb, mom0, 'momentum', numele)
writefile.kine(mom0, symb, totalt, out0)

intcycle = 1
while True:
	outfile = open(out0, 'a+')
	outfile.write('\n' + '='*84 + '\n')
	outfile.write(' '*41 + 'Cycle %d' %intcycle + '\n')
	outfile.write('='*84 + '\n'*2)
	outfile.close()

	for i in range(numele):
		mass = float(mas[symb[i]])
		x1 = pos0[i] + mom0[i] * dt / mass + far0[i] * dt ** 2 / (2 * mass)
		pos0[i] = x1

	writefile.writecon(out0, symb, pos0, 'coordinate', numele)

	inputfile = quansoft.reinpotfile(inputfile, pos0, intcycle, numele, filename)
	runfile = quansoft.runsoft(inputfile)

	far1 = copy.deepcopy(far0)
	far0 = quansoft.getforce(runfile, numele, out0, totalv, symb, start)

	for i in range(numele):
		p1 = mom0[i] + dt * (far0[i] + far1[i]) / 2
		mom0[i] = p1

	writefile.writecon(out0, symb, mom0, 'momentum', numele)
	writefile.kine(mom0, symb, totalt, out0)

	if numcycle:
		if intcycle == numcycle:
			writefile.endtime(out0, 'cycle', start)
			exit('cycle normal exit')

	if scfenry:
		if abs(totalv[-1] - totalv[-2]) < scfenry:
			writefile.endtime(out0, 'energy', start)
			exit('energy normal exit')

	intcycle += 1
