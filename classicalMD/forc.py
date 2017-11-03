#!/usr/bin/env python

import os
import re
import copy
import numpy as np
from sys import argv

ff, nm, inpt = argv
name = inpt.split('.')[0]

def deltx(x0, p0, f0, m):
	x1 = x0 + p0 * dt / m + f0 * dt ** 2 / (2 * m )
	return x1

def deltp(p0, f0, f1):
	p1 = p0 + dt * (f0 + f1) / 2
	return p1

def momen(en):
	phi = random.uniform(0, math.pi)
	tha = random.uniform(0, 2 * math.pi)

def poforc(log, ele, far, pos):
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
	
	pos = np.array(pos, dtype = float)
	far = np.array(far, dtype = float)
	if not len(det):
		return ele, far, pos
	else:
		return far, pos

def md():
	pass

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
	nem.append(eum[i].split(','))

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

ele0 = []
far0 = []
pos0 = []
ele1, far1, pos1 = poforc(nm, ele0, far0, pos0)
'''
nnn = 1
while 1:
'''	
