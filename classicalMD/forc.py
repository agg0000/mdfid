#!/usr/bin/env python

import os
import re
import copy
import math
import time
import random
import numpy as np
from sys import argv

ff, nm, inpt = argv
name = nm.split('.')[0]

mas = {'1' : 1823.0 , '3' : 12761.0, '4' : 16407.0, '5' : 19688.4, '53':231521.0,
       '6' : 21876.0, '7' : 25522.0, '8' : 29168.0, '9' : 34637.0, '10': 36460.0,
	   '11': 41929.0, '12': 44298.9, '13': 49221.0, '14': 51044.0, '15': 56513.0,
	   '16': 58336.0, '17': 64716.5, '20': 72920.0, '19': 71097.0, '35':145840.0,}

def refor(efg, x1, nnn):
	rfg  = "%s%s" %(name, str(nnn))
	sefg = efg.split('.')[0]
	fgjf = open("%s.gjf" %sefg).readlines()
	agjf = open("%s.gjf" %rfg, 'w')

	elea = []
	for i in range(nat):
		elem = fgjf[6: nat + 6]
		
		for i in range(nat):
			elet = re.split(' |\n', elem[i])

			while '' in elet:
				elet.remove('')

			elea.append(elet[0])

	for i in range(6):
		agjf.write(fgjf[i])

	for i in range(6, 6 + nat):
		agjf.write(' ')
		agjf.write(elea[i - 6])
		agjf.write('                  ')

		for ii in range(3):
			agjf.write('    ')
			agjf.write(str(x1[i -6][ii]))
		
		agjf.write('\n')

	for i in range(6 + nat, len(fgjf)):
		agjf.write(fgjf[i])

	agjf.close()
	return "%s.gjf" %rfg

def gauss(rfg):
	erfg = rfg.split('.')[0]
	os.system("subg09 %s.gjf" %erfg)

	while not os.path.exists("%s.log" %erfg):
		time.sleep(3)
	
	return "%s.log" %erfg

def deltx(x0, p0, f0, m):
	x1 = x0 + p0 * dt / m + f0 * dt ** 2 / (2 * m )
	return x1

def deltp(p0, f0, f1):
	p1 = p0 + dt * (f0 + f1) / 2
	return p1

def momen(en):
	phi = random.uniform(0, math.pi)
	tha = random.uniform(0, 2 * math.pi)
	
	pzz = en * math.sin(phi)
	pyy = en * math.cos(phi) * math.cos(tha)
	pxx = en * math.cos(phi) * math.sin(tha)
	
	mom = [pxx, pyy, pzz]
	return mom

def poforc(log, ele, far, pos, sig):
	flog = open(log).readlines()

	row = 0
	raw = 0
	for i in range(len(flog)):
		if 'Forces ' in flog[i]:
			row = i

		if 'Standard o' in flog[i]:
			raw = i

	if row:
		nlog = flog[row + 3: row + 3 + nat]
	else:
		print("no 'Forces ' in the list")
		exit()
	
	if raw:
		plog = flog[raw + 5: raw + 5 + nat]
	else:
		print("no positon")
		exit()

	deo = copy.deepcopy(sig)
	det = copy.deepcopy(ele)
	for i in range(nat):
		anlog = re.split(' |\n', nlog[i])
		aplog = re.split(' |\n', plog[i])
		
		while '' in anlog:
			anlog.remove('')
		
		while '' in aplog:
			aplog.remove('')

		far.append(anlog[2:])
		pos.append(aplog[3:])
		if not len(det):
			ele.append(anlog[0])
		
		if not len(deo):
			sig.append(anlog[1])

	pos = np.array(pos, dtype = float)
	far = np.array(far, dtype = float)

	if not len(det) or not len(deo):
		return ele, far, pos, sig
	else:
		return far, pos

'''
输入文件inp为所需参数
第一行为间隔时间
第二行为总原子个数
第三行为初始原子所需给定能量
若不指定则默认为0,
所给定原子个数以分号（;）隔开
给定原子和给定动量以（,）隔开
'''
inp = open(inpt).readlines()
dt  = float(inp[0])
nat = int(inp[1])
eum = re.split(';|\n', inp[2])

while '' in eum:
	eum.remove('')

num = len(eum)
nem = []
for i in range(num):
	neme = eum[i].split(',')
	neme[0] = str(int(neme[0]) - 1)
	nem.append(neme)

pem = []
for i in range(len(nem)):
	pem.append(nem[0])

pnm = len(pem)
for i in range(nat):
	a = str(i)
	if a not in pem:
		nem.append([i, 0])

dum = {}
for i in range(nat):
	dum[str(nem[i][0])] = int(nem[i][1])

gjf0 = "%s.gjf" %name
ele0 = []
far0 = []
pos0 = []
sig0 = []
mom1 = []
for i in range(nat):
	a = str(i)
	moma = momen(dum[a])
	mom1.append(moma)

ele1, far1, pos1, sig1 = poforc(nm, ele0, far0, pos0, sig0)

nnn = 1
while 1:
	for i in range(nat):
		x1 = deltx(mom1[i], pos1[i], far1[i], mas[sig1[i]])
		pos1[i] = x1

	gjf0 = refor(gjf0, pos1, nnn)
	newlog = gauss(gjf0)
	
	far2 = copy.deepcopy(far1)
	far1, pos1 = poforc(newlog, ele1, far1, pos1, sig1)

	for i in range(nat):
		p1 = deltp(mom1[i], far2[i], far1[i])
		mom1[i] = p1

	deltf = far1 - far2
	if abs(deltf.sum()) < 0.001:
		break
	
	nnn += 1
	print('cycle %d' %nnn)
