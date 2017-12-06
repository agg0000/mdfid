#!/usr/bin/env python
'''
for create output file
'''

def writecon(outname, sym, con, group):
	'''
	writecon means write content.
	because the out file has to many similarity content about force, momentum and position
	so use this function to write down the similarity content in every cycle
	'''
	outfile = open(outname, 'a+')
	outfile.write('\n')
	outfile.write(' '*35)
	outfile.write('The %s'%group)
	outfile.write('\n')
	outfile.write(' Symbol')
	outfile.write(' '*24)
	outfile.write('x')
	outfile.write(' '*16)
	outfile.write('y')
	outfile.write(' '*16)
	outfile.write('z')
	outfile.write('\n')
	outfile.write('-'*77)
	outfile.write('\n')

	copysym = copy.deepcopy(sym)
	for i in range(numele):
		outfile.write('   ')

		if len(copysym[i]) == 1:
			copysym[i] = copysym[i] + ' '

		outfile.write(copysym[i])
		outfile.write(' '*15)

		'''
		this loop is for aline of the number
		'''
		for ii in range(3):
			floatcon = float(con[i][ii])
			if floatcon == 0:
				sign = str(abs(floatcon))
				size = 14 - len(sign)
				sign = sign + size * '0'
				len0 = 4

			elif floatcon > 0:
				sign = str(con[i][ii])
				size = 14 - len(sign)
				sign = sign + size * '0'
				len0 = 4

			elif float(con[i][ii]) < 0:
				sign = str(con[i][ii])
				size = 15 - len(sign)
				sign = sign + size * '0'
				len0 = 3

			outfile.write(' '*len0)
			outfile.write(sign)

		outfile.write('\n')

	outfile.write('-'*77)
	outfile.write('\n'*2)
	outfile.close()

def kine(mom, sym, kinetic, outname):
	'''
	get the kinetic energy from the momentum
	'''
	kine = 0

	for i in range(len(sym)):
		mass = float(mas[sym[i]])
		momu = mom[i] ** 2 / mass
		kine = kine + momu.sum()

	kinetic.append(kine)
	
	outfile = open(outname, 'a+')
	outfile.write('the total kinetic is %f a.u.' %kinetic[-1])
	outfile.write('\n'*2)
	outfile.close()


