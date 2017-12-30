#!/usr/bin/env python
'''
for choose init or continue
'''

import numpy as np

import softapi
import writefile

def initp(numele, totalmom, softuse, filename, out0, dt, start):
	'''
	init function
	'''
	intcycle = 0
	totalv = []
	quansoft = softapi.softname(softuse)

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
	outfile.write(' Welcome to zongzi program '.center(99, '*') + '\n')
	outfile.write('*'*99 + '\n'*2)

	outfile.write('input file name is <%s> \n' %filename)
	outfile.write('step(time) is %d a.u. \n' %dt)
	outfile.write('number of atoms need to calculate is %d \n\n' %numele)

	outfile.write('='*84 + '\n')
	outfile.write('cycle 0'.center(84) + '\n')
	outfile.write('='*84 + '\n'*2)
	outfile.close()

	pos0, symb = quansoft.getpos(filename, numele)

	quansoft.runsoft(filename, intcycle)
	far0 = quansoft.getforce(filename, numele, out0, totalv, symb, start, intcycle)
	mom0 = np.array(originmom, dtype = float)

	writefile.writecon(out0, symb, pos0, 'coordinate', numele)
	writefile.writecon(out0, symb, far0, 'Forces', numele)
	totalv = quansoft.getenergy(filename, out0, totalv, start, intcycle)
	writefile.writecon(out0, symb, mom0, 'momentum', numele)
	writefile.kine(mom0, symb, out0)

	intcycle += 1

	return totalv, pos0, symb, far0, mom0, intcycle

#------------------------------------------------------------------------------------

def continuep(numele, out0, start):
	'''
	continue function
	'''
	outfile = open(out0)
	txt = outfile.readlines()
	outfile.close()

	cyclinearray = []
	for i, line in enumerate(txt):
		if 'Cycle' in line:
			cyclinearray.append(i)

	newtxt = open(out0, 'a+')
	newtxt.write('\n\n' + '#'*99 + '\n')
	newtxt.write(('restart at %s' %str(start)).center(99))
	newtxt.write('\n' + '#'*99 + '\n')
	newtxt.close()

	lastline = cyclinearray[-1]
	intcycle = int(txt[lastline].split()[1]) + 1

	for i, line in enumerate(txt[lastline:]):
		if 'coordinate' in line:
			coornum = i + lastline

		if 'momentum' in line:
			momenum = i + lastline

		if 'Forces' in line:
			forcnum = i + lastline

	pos0 = []
	symb = []
	for i in range(numele):
		fs = txt[coornum + 4 + i].split()
		xq = fs[1:]
		symb.append(fs[0])
		pos0.append(xq)
	
	pos0 = np.array(pos0, dtype = float)

	far0 = []
	for i in range(numele):
		f = txt[forcnum + 4 + i].split()[1:]
		far0.append(f)

	far0 = np.array(far0, dtype = float)

	ev = float(txt[forcnum + 6 + numele].split()[4])
	totalv = [ev]

	mom0 = []
	for i in range(numele):
		p = txt[momenum + 4 + i].split()[1:]
		mom0.append(p)

	mom0 = np.array(mom0, dtype = float)

	return totalv, pos0, symb, far0, mom0, intcycle
