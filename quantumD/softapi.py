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

	def getfreq(self, filename, starttime):
		'''
		get the opt and freq infomation for origin configuration
		the detail of wigner sampling is from sharc(Surface Hopping including ARbitrary Couplings)
		http://sharc-md.org
		'''
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
				writefile.werroe(outname, 'error input gjf', starttime)
				
		if basisset == 'mark':
			exit('error basis set')
			writefile.werroe(outname, 'error basis set', starttime)

		origincon[setnumber] = '# %s opt freq\n' %basisset
		newgjf = filename + 'freq.gjf'

		newgjfcon = open(newgjf, 'a+')
		for lin in origincon:
			newgjfcon.write(lin)

		newgjfcon.close()

		newlog = self.runsoft(newgjf)
		logcon = open(newlog).readlines()

#-------------------------------------------------------------------------------------

	def reinpotfile(self, oldname, newpos, intcycle, numele, getname):
		'''
		use the new position to create a new gjf file
		'''
		strcycle= str(intcycle) # mark the number of gjf file
		newname = "%s%s.gjf" %(getname, strcycle)
		oldfile = open(oldname).readlines()
		newfile = open(newname, 'a+')
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

	def runsoft(self, inputfile):
		'''
		for run the gaussian progam, and return the log file
		'''
		inputname = inputfile.split('.')[0]
		os.system('g09 %s' %inputfile)
	
		return '%s.log' %inputname

#------------------------------------------------------------------------------------

	def getinfo(self, resultname, numele, outname, potential, starttime):
		'''
		get the force and position from the soft output file
		only use to calculate the Hartree-Fork theory
		'''
		resultfile = open(resultname).readlines()

		signf = 0
		signe = 0
		signq = 0
		for i in range(len(resultname)):
			if 'Forces 'in resultfile[i]:
				signf = i

			if 'SCF Done' in resultfile[i]:
				signe = i

			if 'Input orientation' in resultfile[i]:
				signq = i

		if 'Force ' in resultfile[signf]:
			numf = resultfile[signf + 3: signf + 3 + numele]
		else:
			writefile.werroe(outname, 'forces', starttime)
			exit('no force in outfile')

		sfar = []
		for i in range(numele):
			eachforce = numf[i].split(' ')

			while '' in eachforce:
				eachforce.remove('')

			sfar.append(eachforce[2:])

		far = np.array(sfar, dtype = float)

		if 'SCF Done' in resultfile[signe]:
			energyrow = resultfile[signe].split(' ')

			while '' in energyrow:
				energyrow.remove('')

			outenergy = float(energyrow[5])
			potential.append(outenergy)
		else:
			writefile.werroe(outname, 'SCF Done', starttime)
			exit('no energy in outfile')

		if 'Input orientation' in resultfile[signq]:
			numq = resultfile[signq + 5: sign + 5 + numele]
		else:
			writefile.werroe(outname, 'position', starttime)
			exit('no position in outfile')

		symb = []
		spos = []
		for i in range(numele):
			eachpos = numq.split(' ')

			while '' in eachpos:
				eachpos.remove('')

			spos.append(eachpos[3:])
			symb.append(eachpos[1])

		pos = np.array(spos, dtype = float)

		return pos, far

#####################################################################################
