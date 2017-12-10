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
filename = inp[3],strip()
numcycle = int(inp[4])
scfenry = float(inp[5])
softuse = inp[6].strip()

quansoft = softapi.softname(softuse)
inputfile = quansoft.getfile(filename)
outfile = "%s.out" %filename

totalv = []
totalt = []

