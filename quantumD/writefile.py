#!/usr/bin/env python
'''
for create output file
'''

mas = {'H' : 1823.0 , 'Li': 12761.0, 'Be': 16407.0, 'B' : 19688.4, 'I' :231521.0,
       'C' : 21876.0, 'N' : 25522.0, 'O' : 29168.0, 'F' : 34637.0, 'Ne': 36460.0,
       'Na': 41929.0, 'Mg': 44298.9, 'Al': 49221.0, 'Si': 51044.0, 'P' : 56513.0,
       'S' : 58336.0, 'Cl': 64716.5, 'Ca': 72920.0, 'K' : 71097.0, 'Br':145840.0,}


import copy
import datetime

def writecon(outname, sym, con, group, numele):
	'''
	writecon means write content.
	because the out file has to many similarity content about force, momentum and position
	so use this function to write down the similarity content in every cycle
	'''
	outfile = open(outname, 'a+')
	outfile.write('\n' + ' '*35 + 'The %s'%group + '\n')
	outfile.write(' Symbol' + ' '*24 + 'x' + ' '*16 + 'y' + ' '*16 + 'z\n')
	outfile.write('-'*77 + '\n')

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

	outfile.write('-'*77 + '\n'*2)
	outfile.close()

	return

#------------------------------------------------------------------------------------

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

	return

#------------------------------------------------------------------------------------

def werror(outname, keyword, starttime):
	'''
	for write down the error infomation
	'''

	endtime = datetime.datetime.now()
	usetime = endtime - starttime

	outfile = open(outname, 'a+')
	outfile.write('*'*99 + '\n')
	outfile.write('ERROR %s' %keyword + '\n')
	outfile.write('Grad'*24 + '\n')
	outfile.write('*'*99 + '\n')
	outfile.write('no %s in outfile\n' %keyword)
	outfile.write('start program at ' + str(starttime) + '\n')
	outfile.write('use time ' + str(usetime))
	outfile.close()

	return

#------------------------------------------------------------------------------------

def endtime(outname, keyword, starttime):
	endtime = datetime.datetime.now()
	usetime = endtime - starttime

	outfile = open(outname, 'a+')
	outfile.write('\n' + '*'*99 + '\n')
	outfile.write('cycle %d normal exit' %intcycle + '\n')
	outfile.write('Start program at  ' + str(starttime) + '\n')
	outfile.write('Use time  ' + str(usetime) + '\n')
	outfile.wrrite('%s normal exit' %keyword)
	outfile.close()
	
	return
