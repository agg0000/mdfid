#!/usr/bin/env python
'''
set up the API for use quantum chemistry soft
'''

import os
import re
import datetime
import numpy as np

import writefile

autoan = 0.529177211
def softname(softuse):
	'''
	convince the quantum chemistry from the keyword
	'''
	softall = {"gaussian": gaussian(),  "molpro": molpro(),
                   "turbomole":turbomole()}

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

	def reinpotfile(self, filename, newpos, symb, intcycle, numele, getname):
		'''
		use the new position to create a new gjf file
		'''
		oldcycle = intcycle - 1

		if oldcycle:
			oldname = filename + str(oldcycle)

		else:
			oldname = filename

		oldfile = open("%s.gjf" %oldname).readlines()
		newfile = open("%s%d.gjf" %(filename, intcycle), 'w')

		linenu = 4 # ensure the positine line number.
		for i, line in enumerate(oldfile):
			if line[0] in ['%', '#']:
				linenu += 1
		
			if filename in line:
				oldfile[i] = re.sub(r"%s\d*" %filename, "%s%d" %(filename, intcycle), line)

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
		return

#-------------------------------------------------------------------------------------

	def getpos(self, filename, numele):
		'''
		from the first gjf file to get the first position
		'''
		oldfile = open('%s.gjf' %filename).readlines()
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

	def runsoft(self, filename, intcycle):
		'''
		for run the gaussian progam, and return the log file
		'''
		if intcycle:
			runfile = filename + str(intcycle)

		else:
			runfile = filename

		os.system('g09 %s.gjf' %runfile)
		return

#------------------------------------------------------------------------------------

	def getforce(self, filename, numele, outname, potential, symb, starttime, intcycle):
		'''
		get the force and position from the soft output file
		only use to calculate the Hartree-Fork theory
		'''
		if intcycle:
			resultname = filename + str(intcycle)

		else:
			resultname = filename

		resultfile = open("%s.log" %resultname).readlines()

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

	def getenergy(self, filename, outname, potential, starttime, intcycle):
		'''
		get the energy
		'''
		if intcycle:
			resultname = filename + str(intcycle)

		else:
			resultname = filename

		resultfile = open("%s.log" %resultname).readlines()

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

	def reinpotfile(self, filename, newpos, symb, intcycle, numele, getname):
		'''
		use the new position to create a new in file
		'''
		oldcycle = intcycle - 1

		if oldcycle:
			oldname = filename + str(oldcycle)

		else:
			oldname = filename

		oldfile = open("%s.in" %oldname).readlines()
		newfile = open("%s%d.in" %(filename, intcycle), 'w')

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
		return

#-------------------------------------------------------------------------------------

	def getpos(self, filename, numele):
		'''
		from the first gjf file to get the first position
		'''
		oldfile = open('%s.in' %filename).readlines()
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

	def runsoft(self, filename, intcycle):
		'''
		for run the gaussian progam, and return the log file
		'''
		if intcycle:
			runfile = filename + str(intcycle)

		else:
			runfile = filename
	
		os.system('molpro %s.in' %runfile)
		return

#------------------------------------------------------------------------------------

	def getforce(self, filename, numele, outname, potential, symb, starttime, intcycle):
		'''
		get the force and position from the soft output file
		only use to calculate the Hartree-Fork theory
		'''
		if intcycle:
			resultname = filename + str(intcycle)

		else:
			resultname = filename

		resultfile = open("%s.out" %resultname).readlines()

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

	def getenergy(self, filename, outname, potential, starttime, intcycle):
		'''
		get the energy
		'''
		if intcycle:
			resultname = filename + str(intcycle)

		else:
			resultname = filename

		resultfile = open("%s.out" %resultname).readlines()

		try:
			ev = resultfile[-3].split()[0]
			outenergy = float(ev)		
			potential.append(outenergy)
			
			outfile = open(outname, 'a+')
			outfile.write('the total potential is %f a.u. \n' %outenergy)
			outfile.close()

		except:
			writefile.werror(outname, 'energy', starttime)
			exit('no energy in outfile')

		return potential

#####################################################################################
#                             class turbomole                                       #
#####################################################################################

class turbomole():
	'''
	use the molpro soft
	'''
	def __init__(self):
		self = self

#-------------------------------------------------------------------------------------

	def reinpotfile(self, filename, newpos, symb, intcycle, numele, getname):
		oldcycle = intcycle - 1

		if oldcycle:
			oldname = filename + str(oldcycle)

		else:
			oldname = filename

		os.system("mkdir -p %s%d" %(filename, intcycle))
		oldfile = open("%s//coord" %oldname).readlines()
		newfile = open("%s%d//coord" %(filename, intcycle), 'w')

		newpos = newpos / autoan
		newfile.write(oldfile[0])
		for i in range(1, numele + 1):
			for ii in range(3):
				newfile.write('    ')
				newfile.write(str(newpos[i - 1][ii]))

			newfile.write('  ')
			newfile.write(symb[i - 1])

			newfile.write('\n')

		for i in range(1 + numele, len(oldfile)):
			newfile.write(oldfile[i])

		newfile.close()
		return 

#-------------------------------------------------------------------------------------

	def getpos(self, filename, numele):
		oldfile = open("%s//coord" %filename).readlines()
		
		pos = []
		sym = []
		for i in range(1, numele + 1):
			posarray = oldfile[i].split()
			
			pos.append(posarray[:3])
			sym.append(posarray[-1].upper())

		pos = np.array(pos, dtype = float)
		pos = pos * autoan

		return pos, sym

#-------------------------------------------------------------------------------------

	def runsoft(self, filename, intcycle):
		oldcycle = intcycle - 1
		intname = filename + str(intcycle)
		oldname = filename + str(oldcycle)

		if not intcycle:
			os.chdir(filename)
			os.system("dscf")
			os.system("grad")
			os.chdir("..")

		elif not oldcycle:
			os.system("cp %s//control %s" %(filename, intname))
			os.system("cp %s//basis %s" %(filename, intname))
			os.system("cp %s//mos %s" %(filename, intname))
			os.chdir(intname)
			os.system("dscf")
			os.system("grad")
			os.chdir("..")

		else:
			os.system("cp %s//control %s" %(oldname, intname))
			os.system("cp %s//basis %s" %(oldname, intname))
			os.system("cp %s//mos %s" %(oldname, intname))
			os.chdir(intname)
			os.system("dscf")
			os.system("grad")
			os.chdir("..")
		
		return

#-------------------------------------------------------------------------------------

	def getforce(self, filename, numele, outname, potential, symb, starttime, intcycle):
		if intcycle:
			resultname = filename + str(intcycle)

		else:
			resultname = filename

		if not os.path.isfile("%s//gradient" %resultname):
			writefile.werror(outname, 'forces', starttime)
			exit('no force in outfile')
			
		resultfile = open("%s//gradient" %resultname).readlines()

		numf = resultfile[numele + 2: numele + 2 + numele]
		
		sfar = []
		for i in range(numele):
			nuf = numf[i].replace("D", "e")
			eachforce = nuf.split()
			sfar.append(eachforce)

		far = np.array(sfar, dtype = float)
		far = far * -1

		return far

#------------------------------------------------------------------------------------

	def getenergy(self, filename, outname, potential, starttime, intcycle):
		if intcycle:
			resultname = filename + str(intcycle)

		else:
			resultname = filename

		if not os.path.isfile("%s//energy" %resultname):
			writefile.werror(outname, 'energy', starttime)
			exit('no energy in outfile')

		resultfile = open("%s//energy" %resultname).readlines()

		ev = resultfile[1].split()[1]
		outenergy = float(ev)		
		potential.append(outenergy)

		outfile = open(outname, 'a+')
		outfile.write('the total potential is %f a.u. \n' %outenergy)
		outfile.close()

		return potential


