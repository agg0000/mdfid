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
	softall = {"gaussian":gaussian()}

	if softuse not in softall:
		exit("please use 'gaussian'")

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

	def reinpotfile(self, oldname, newpos, intcycle, numele, getname):
		'''
		use the new position to create a new gjf file
		'''
		strcycle= str(intcycle) # mark the number of gjf file
		newname = "%s%s.gjf" %(getname, strcycle)
		oldfile = open(oldname).readlines()
		newfile = open(newname, 'w')
		filelen = len(oldfile)

		linenu = 4 # ensure the positine line number.
		chkline = filelen + 1 #ensure chkline exit 
		for i in range(filelen):
			if oldfile[i][0] in ['%', '#']:
				linenu += 1
		
			if "%chk=" in oldfile[i]:
				chkline = i

		elearray = [] # elearray means element array
		elementa = oldfile[linenu: numele + linenu] # elementa means all elements of rows
		for i in range(numele):# to get the symbols of element
			element = elementa[i].split(' ')

			while '' in element:
				element.remove('')

			elearray.append(element[0])
	
		if chkline < filelen:
			oldchkrow = re.split('\.|\d', oldfile[chkline]) # to rewrite the chk file name

			while '' in oldchkrow:
				oldchkrow.remove('')

			chkrow = oldchkrow[0] + strcycle + '.' + oldchkrow[1]

		for i in range(linenu): 
			if i == chkline:
				newfile.write(chkrow)
				continue

			newfile.write(oldfile[i])

		for i in range(linenu, linenu + numele): # the new gjf file is the same as old gjf file except position
			newfile.write(' ')
			newfile.write(elearray[i - linenu])
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

#	def getinitfile():

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
		signe = 0
		for i in range(len(resultfile)):
			if 'Forces 'in resultfile[i]:
				signf = i

			if 'SCF Done' in resultfile[i]:
				signe = i

		if 'Forces ' in resultfile[signf]:
			numf = resultfile[signf + 3: signf + 3 + numele]
		else:
			writefile.werror(outname, 'forces', starttime)
			exit('no force in outfile')

		sfar = []
		for i in range(numele):
			eachforce = numf[i].split(' ')

			while '' in eachforce:
				eachforce.remove('')

			sfar.append(eachforce[2:])

		far = np.array(sfar, dtype = float)
		writefile.writecon(outname, symb, far, 'Forces', numele)

		if 'SCF Done' in resultfile[signe]:
			energyrow = resultfile[signe].split()

			outenergy = float(energyrow[4])
			potential.append(outenergy)

			outfile = open(outname, 'a+')
			outfile.write(resultfile[signe])
			outfile.close()
		else:
			writefile.werror(outname, 'SCF Done', starttime)
			exit('no energy in outfile')

		return far

#####################################################################################
