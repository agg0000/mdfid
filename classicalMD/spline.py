#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 11:53:59 2017
2017.11.12
中午，完成下三角行列式
下午，完成通过下三角求解方程组
晚上，完成矩阵部分
剩下系数部分未解决
@author: zongzi
"""

import copy
import numpy as np
from sys import argv

'''
将行列式化为下三角行列式
若直接计算下三角行列式，
则令n=0
若计算AX=B
将B的转置放于A末尾
令n != 0
'''
def ut(det, n):
	if n:
		qet = det[-1]
		det = det[:-1]

	for i in range(1, len(det)):
		for ii in range(i, len(det)):
			if det[-ii-1][-i]:
				rat = det[-i][-i] / det[-ii-1][-i]
				det[-ii-1] = det[-ii-1] * rat - det[-i]
				
				if n:
					qet[-ii-1] = qet[-ii-1] * rat - qet[-i]

	if n:
		det = list(det)
		det.append(qet)
		det = np.array(det)

	return det

'''
若计算AX=B
调用sol函数
'''
def sol(x):
	qet = x[-1]
	det = x[:-1]

	soa = np.zeros(len(det))
	for i in range(len(soa)):
		rao = 0
		for ii in range(i):
			rao = rao + det[i][ii] * soa[ii] 

		soa[i] = (qet[i] - rao) / det[i][i]

	return soa

'''
计算三次样条中m值的行列式
若k = 'a'，则采用自由边界条件，
若k = 'b'，则采用固定边界条件，
若k = 'c'，则采用非节点边界
'''
def spl(x, y, k):
	x1 = x[1:]
	x2 = x[:-1]
	dx = x2 - x1
	h = np.zeros((len(x),len(x)))
	
	if k == 'a':
		
		h[0][0]   = 1
		
		h[-1][-1] = 1
	
	elif k == 'b':
	
		h[0][0]   = 2 * dx[0]
		h[0][1]   = dx[0]
		
		h[-1][-1] = 2 * dx[-1]
		h[-1][-2] = dx[-1]

	elif k == 'c':

		h[0][0]   = -dx[1]
		h[0][1]   =  dx[0] + dx[1]
		h[0][2]   = -dx[0]
		
		h[-1][-1] = -dx[-2]
		h[-1][-2] =  dx[-1] + dx[-2]
		h[-1][-2] = -dx[-1]
		
	for i in range(len(x1)-1):
		h[i + 1][i]     = dx[i]
		h[i + 1][i + 1] = 2 * (dx[i + 1] + dx[i])
		h[i + 1][i + 2] = dx[i + 1]
	
	z = np.zeros((len(x)))
	for i in range(1, len(x) - 1):
		z = (y[i + 1] - y[i]) / dx[i] - (y[i]- y[i - 1]) / dx[i - 1]
		z = 6 * z
		
	h = list(h)
	h.append(z)
	h = np.array(h)
	
	return h
	
def coe(x, y):
	
	

fn, nm = argv

rea = open(nm).readlines()
nux = np.zeros((len(rea)))
nuy = np.zeros((len(rea)))

for i in range(len(rea)):
	rec = re.split('\t|\n', rea[i])
	
	while '' in rea:
		rea.remove('')
		
	nux[i] = rea[0]
	nuy[i] = rea[-1]
	
nuh = spl(nux, nuy, 'c')
uth = ut(nuh, 1)
soh = sol(uth)