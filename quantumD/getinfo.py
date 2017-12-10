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
