#!/usr/bin/env python
'''
set up the API for use quantum chemistry soft
'''

import os
import re
import datetime
import numpy as np

def softname(softuse):
	'''
	convince the quantum chemistry from the keyword
	'''
	softall = {"gaussian":gaussian()}
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

	def regjf(self, oldgjf, newpos, intcycle, numele, gjfname):
		'''
		use the new position to create a new gjf file
		'''
		strcycle  = str(intcycle) # mark the number of gjf file
		newgjf  = "%s%s.gjf" %(gjfname, strcycle)
		oldfile = open(oldgjf).readlines()
		newfile = open(newgjf, 'a+')
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
			element = re.split(' |\n', elementa[i])

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
		return newgjf

#------------------------------------------------------------------------------------

	def getpos(self, gjfname, numele, outname):
		'''
		from the first gjf file to get the first position
		'''
		gjffile = open(gjfname).readlines()
		filelen = len(gjffile)
	
		linenu = 4 # ensure the positine line number.
		for i in range(filelen):
			if gjffile[i][0] in ['%', '#']:
				linenu += 1

		numpos  = gjffile[linenu: linenu + numele]

		pos = []
		sym = []
		for i in range(numele):
			posarray = re.split(' |\n', numpos[i])

			while '' in posarray:
				posarray.remove('')

			if posarray[0] == 1:
				posarray[0] = ' ' + posarray[0]

			pos.append(posarray[1:])
			sym.append(posarray[0])

		pos = np.array(pos, dtype = float)

		return pos, sym

#-------------------------------------------------------------------------------------

	def rungau(newgjffile):
		'''
		for run the gaussian progam, and return the log file
		'''
		newfilename = newgjffile.split('.')[0]
		os.system('g09 %s' %newgjffile)
	
		return '%s.log' %newfilename

#------------------------------------------------------------------------------------

	def getforce(logname, sym, numele, outname, potential, start):
		'''
		from the log file getting the force and potential energy
		and write down to the out file
		'''
		logfile = open(logname).readlines()

		signforce = 0
		signdone  = 0
		itnenergy = []
		for i in range(len(logfile)):
			if 'Forces ' in logfile[i]:
				signforce = i

			if 'SCF Done' in logfile[i]:
				signdone = i
			elif 'ITN=' in logfile[i]:
				itnenergy.append(i)
	
		if 'Forces ' in logfile[signforce]:
			numforce = logfile[signforce + 3: signforce + 3 + numele]
		else:
			endtime = datetime.datetime.now()
			usetime = endtime - start

			outfile = open(outname, 'a+')
			outfile.write('*'*99)
			outfile.write('\n')
			outfile.write('ERROR FORCES')
			outfile.write('\n')
			outfile.write('Gd'*49)
			outfile.write('\n')
			outfile.write('*'*99)
			outfile.write('\n')
			outfile.write('no force in log')
			outfile.write('\n')
			outfile.write('Start program at  ')
			outfile.write(str(start))
			outfile.write('\n')
			outfile.write('Use time  ')
			outfile.write(str(usetime))
			outfile.close()
			exit('no force in the log')

		far = []
		for i in range(numele):
			eachforce = re.split(' |\n', numforce[i])

			while '' in eachforce:
				eachforce.remove('')

			far.append(eachforce[2:])

		far = np.array(far, dtype = float)

		writecon(outname, sym, far, 'Forces')
	
		'''
		casscf write down the energy using "ITN="
		'''
		outfile = open(outname, 'a+')
		if 'SCF Done' in logfile[signdone]:
			outfile.write(logfile[signdone])
			outfile.write('\n')
			outfile.close()

			energyrow = logfile[signdone].split(' ')

			while '' in energyrow:
				energyrow.remove('')

			outenergy = float(energyrow[4])
			potential.append(outenergy)
	
		elif 'ITN=' in logfile[itnenergy[-1]]:
			outfile.write(logfile[itnenergy[-1]])
			outfile.write('\n')
			outfile.close()

			energyrow = logfile[itnenergy[-1]].split(' ')

			while '' in energyrow:
				energyrow.remove('')

			outenergy = float(energyrow[5])
			potential.append(outenergy)

		else:
			endtime = datetime.datetime.now()
			usetime = endtime - start

			outfile.write('*'*99)
			outfile.write('\n')
			outfile.write('ERROR ENERGY')
			outfile.write('\n')
			outfile.write('Gd'*49)
			outfile.write('\n')
			outfile.write('*'*99)
			outfile.write('\n')
			outfile.write('no SCF Done in log')
			outfile.write('\n')
			outfile.write('Start program at  ')
			outfile.write(str(start))
			outfile.write('\n')
			outfile.write('Use time  ')
			outfile.write(str(usetime))
			outfile.close()
			exit('no SCF Done in the log')

		return far

#####################################################################################
