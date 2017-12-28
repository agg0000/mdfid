#!/usr/bin/env python
'''
set up the API for use quantum chemistry soft
'''

import os
import re
import datetime
import numpy as np

import writefile

def softname(softuse):
	'''
	convince the quantum chemistry from the keyword
	'''
	softall = {"gaussian":gaussian(), "molpro":molpro()}

	if softuse not in softall:
		exit("please use useful key words")

	return softall.get(softuse)

#####################################################################################
#                            class gaussian                                         #
#####################################################################################

class gaussian():
	'''
	use the gaussian soft
	'''
	def __init__(self):
		self = self

#-------------------------------------------------------------------------------------

	def getfile(self, getname):
		return "%s.gjf" %getname

#-------------------------------------------------------------------------------------
	'''
	def getfreq(self, filename, starttime, outname):
	
		origincon = open('%s.gjf' %filename).readlines()
		lencon = len(origincon)

		setnumber = lencon + 1
		for i in range(setnumber):
			if '#' in origincon[i]:
				setnumber = i

		basisset = 'mark'
		if setnumber < lencon:
			setarray = origincon[setnumber].split(' ')

			for keyele in setarray:
				if '\\' in keyele:
					basisset = keyele
			else:
				exit('error input gjf')
				writefile.werror(outname, 'error input gjf', starttime)
				
		if basisset == 'mark':
			exit('error basis set')
			writefile.werror(outname, 'error basis set', starttime)

		origincon[setnumber] = '# %s opt freq\n' %basisset
		newgjf = filename + 'freq.gjf'

		newgjfcon = open(newgjf, 'a+')
		for lin in origincon:
			newgjfcon.write(lin)

		newgjfcon.close()

		newlog = self.runsoft(newgjf)
		logcon = open(newlog).readlines()

		linenum = []
		for i ,line in enumerate(logcon):
			if 'Frequencies' in line:
				linenum.append(i + 5)

		if not len(linenum):
			exit('error freq')
			writefile.werror(outname, 'error freq', starttime)

		freq = {}
		for l in linenum:
			freqele = logcon[l].split()

			for ii, ele in enumerate(freqele):
				if float(ele) < 10:
					continue
	'''
#-------------------------------------------------------------------------------------

	def reinpotfile(self, oldname, newpos, symb, intcycle, numele, getname):
		'''
		use the new position to create a new gjf file
		'''
		strcycle = str(intcycle) # mark the number of gjf file
		newfname = getname + strcycle
		newname = "%s.gjf" %newfname
		oldfile = open(oldname).readlines()
		newfile = open(newname, 'w')

		linenu = 4 # ensure the positine line number.
		for i, line in enumerate(oldfile):
			if line[0] in ['%', '#']:
				linenu += 1
		
			if getname in line:
				oldfile[i] = re.sub(r"%s\d*" %getname, "%s" %newfname, line)

		for i in range(linenu): 
			newfile.write(oldfile[i])

		for i in range(linenu, linenu + numele): # the new gjf file is the same as old gjf file except position
			newfile.write(' ')
			newfile.write(symb[i - linenu])
			newfile.write('                  ')

			for ii in range(3):
				newfile.write('    ')
				newfile.write(str(newpos[i -linenu][ii]))
		
			newfile.write('\n')

		for i in range(linenu + numele, len(oldfile)):
			newfile.write(oldfile[i])

		newfile.close()
		return newname

#-------------------------------------------------------------------------------------

	def getpos(self, oldname, numele):
		'''
		from the first gjf file to get the first position
		'''
		oldfile = open(oldname).readlines()
		filelen = len(oldfile)

		linenu = 4 # ensure the positine line number.
		for i in range(filelen):
			if oldfile[i][0] in ['%', '#']:
				linenu += 1

		numpos  = oldfile[linenu: linenu + numele]

		pos = []
		sym = []
		for i in range(numele):
			posarray = numpos[i].split()

			pos.append(posarray[1:])
			sym.append(posarray[0])

		pos = np.array(pos, dtype = float)

		return pos, sym
	
#-------------------------------------------------------------------------------------

	def runsoft(self, inputfile):
		'''
		for run the gaussian progam, and return the log file
		'''
		inputname = inputfile.split('.')[0]
		os.system('g09 %s' %inputfile)
	
		return '%s.log' %inputname

#------------------------------------------------------------------------------------

	def getforce(self, resultname, numele, outname, potential, symb, starttime):
		'''
		get the force and position from the soft output file
		only use to calculate the Hartree-Fork theory
		'''
		resultfile = open(resultname).readlines()

		signf = 0
		for i, line in enumerate(resultfile):
			if 'Forces ' in line:
				signf = i

		if 'Forces ' in resultfile[signf]:
			numf = resultfile[signf + 3: signf + 3 + numele]
		else:
			writefile.werror(outname, 'forces', starttime)
			exit('no force in outfile')

		sfar = []
		for i in range(numele):
			eachforce = numf[i].split()
			sfar.append(eachforce[2:])

		far = np.array(sfar, dtype = float)
		
		return far

#------------------------------------------------------------------------------------

	def getenergy(self, resultname, outname, potential, starttime):
		'''
		get the energy
		'''
		resultfile = open(resultname).readlines()

		signe = []
		for i, line in enumerate(resultfile):
			if 'SCF Done' in line:
				signe.append(i)

			if '\\MP2=' in line:
				signe.append(i)

			if '\\CISD=' in line:
				signe.append(i)

			if 'ITN=' in line:
				signe.append(i)

		ii = signe[-1]
		if 'SCF Done' in resultfile[ii]:
			energyrow = resultfile[ii].split()

			outenergy = float(energyrow[4])
			potential.append(outenergy)

			outfile = open(outname, 'a+')
			outfile.write('the total potential is %f a.u. \n' %outenergy)
			outfile.close()

		elif '\\MP2=' in resultfile[ii]:
			energyrow = resultfile[ii].split('\\')

			for linemp in energyrow:
				if 'MP2' in linemp:
					outenergy = linemp.split('=')
					outenergy = float(outenergy[1])

			potential.append(outenergy)

			outfile = open(outname, 'a+')
			outfile.write('the total potential is %f a.u. \n' %outenergy)
			outfile.close()

		elif '\\CISD=' in resultfile[ii]:
			energyrow = resultfile[ii].split('\\')

			for linemp in energyrow:
				if 'CISD' in linemp:
					outenergy = linemp.split('=')
					outenergy = float(outenergy[1])

			potential.append(outenergy)

			outfile = open(outname, 'a+')
			outfile.write('the total potential is %f a.u. \n' %outenergy)
			outfile.close()

		elif 'ITN' in resultfile[ii]:
			energyrow = resultfile[ii].split()

			outenergy = float(energyrow[5])
			potential.append(outenergy)

			outfile = open(outname, 'a+')
			outfile.write('the total potential is %f a.u. \n' %outenergy)
			outfile.close()

		else:
			writefile.werror(outname, 'energy', starttime)
			exit('no energy in outfile')

		return potential

#####################################################################################
#                             class molpro                                          #
#####################################################################################

class molpro():
	'''
	use the molpro soft
	'''
	def __init__(self):
		self = self

#-------------------------------------------------------------------------------------

	def getfile(self, getname):
		return "%s.in" %getname

#-------------------------------------------------------------------------------------

	def reinpotfile(self, oldname, newpos, symb, intcycle, numele, getname):
		'''
		use the new position to create a new in file
		'''
		strcycle = str(intcycle)
		newfname = getname + strcycle
		newname = "%s.in" %newfname
		oldfile = open(oldname).readlines()
		newfile = open(newname, 'w')

		geoline = len(oldfile) + 1
		for i, line in enumerate(oldfile):
			if 'geometry' in line:
				geoline = i + 3

		if 'geometry' not in oldfile[geoline - 3]:
			writefile.werror(out0, 'geometry', starttime)

		for i in range(geoline):
			newfile.write(oldfile[i])

		for i in range(geoline, geoline + numele):
			newfile.write(symb[i - geoline])	
			newfile.write('                  ')

			for ii in range(3):
				newfile.write('    ')
				newfile.write(str(newpos[i - geoline][ii]))
		
			newfile.write('\n')

		for i in range(geoline + numele, len(oldfile)):
			newfile.write(oldfile[i])

		newfile.close()
		return newname

#-------------------------------------------------------------------------------------

	def getpos(self, oldname, numele):
		'''
		from the first gjf file to get the first position
		'''
		oldfile = open(oldname).readlines()
		geoline = len(oldfile) + 1
		for i, line in enumerate(oldfile):
			if 'geometry' in line:
				geoline = i + 3

		numpos = oldfile[geoline: geoline + numele]

		pos = []
		sym = []
		for i in range(numele):
			posarray = numpos[i].split()

			pos.append(posarray[1:])
			sym.append(posarray[0])

		pos = np.array(pos, dtype = float)

		return pos, sym

#-------------------------------------------------------------------------------------

	def runsoft(self, inputfile):
		'''
		for run the gaussian progam, and return the log file
		'''
		inputname = inputfile.split('.')[0]
		os.system('molpro %s' %inputfile)
	
		return '%s.out' %inputname

#------------------------------------------------------------------------------------

	def getforce(self, resultname, numele, outname, potential, symb, starttime):
		'''
		get the force and position from the soft output file
		only use to calculate the Hartree-Fork theory
		'''
		resultfile = open(resultname).readlines()

		signf = 0
		for i, line in enumerate(resultfile):
			if 'GRADIENT FOR ' in line:
				signf = i

		if 'GRADIENT FOR ' in resultfile[signf]:
			numf = resultfile[signf + 4: signf + 4 + numele]
		else:
			writefile.werror(outname, 'forces', starttime)
			exit('no force in outfile')

		sfar = []
		for i in range(numele):
			eachforce = numf[i].split()
			sfar.append(eachforce[1:])

		far = np.array(sfar, dtype = float)
		far = far * -1.0
		
		return far

#------------------------------------------------------------------------------------

	def getenergy(self, resultname, outname, potential, starttime):
		'''
		get the energy
		'''
		resultfile = open(resultname).readlines()

		try:
			ev = resultfile[-3].split()[0]
			potential.append(ev)
			ev = float(ev)		
			
			outfile = open(outname, 'a+')
			outfile.write('SCF Done energy is %f a.u.' %ev)
			outfile.close()

		except:
			writefile.werror(outname, 'energy', starttime)
			exit('no energy in outfile')

		return
