#!/usr/bin/env python

import os
import re
import numpy as np
from sys import argv

ff, nm, inpt = argv

def deltx(x0, p0, f0, m):
	x1 = x0 + p0 * dt / m + f0 * dt ** 2 / (2 * m )
	return x1

def deltp(p0, f0, f1):
	p1 = p0 + dt * (f0 + f1) / 2
	return p1

def aforc(log, ele, far):
	flog = open(log).readlines()

	row = 0
	for i in range(len(flog)):
		if 'Forces ' in flog[i]:
			row = i

	if row:
		nlog = flog[row + 3: row + 3 + natom]
	else:
		print "no 'Forces ' in the list"
		exit()

	det = ele
	for i in range(len(nlog)):
		anlog = re.split(' |\n', nlog[i])
		
		while '' in anlog:
			anlog.remove('')

		far.append(np.array(anlog[1:], dtype=float))
		if not len(det):
			ele.append(anlog[0])

	if not len(det):
		return ele, far

	return far

def md():
	pass

'''
输入文件inp为所需参数
第一行为间隔时间
第二行为总原子个数
第三行为初始原子所需给定动量
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
nem = np.array([])
for i in range(num):
	nem.append(num[i].split('.'))

pem = nem[0, 0:]
pnm = len(pem)
for i in range(nat):
	if i not in pem:
		eum.append([i, 0])

dum = {}
for i in range(nat):
	dum[eum[i, 0]] = eum[i, 1]

ele0 = []
far0 = []
ele1, far1 = aforc(nm, ele0, far0)

nnn = 1
while 1:
	far1 = aforc(nm, ele, far)
