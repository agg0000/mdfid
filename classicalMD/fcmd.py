#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 00:37:29 2017

@author: zongzi
"""
import math
import numpy as np
import matplotlib.pyplot as plt

#输入基本参数
ipt = open('inp').readline().split(",")
xin = float(ipt[0]) #设置初始位置
pin = float(ipt[1]) #设置初始动量
mas = float(ipt[2]) #设置基本质量
ssd = float(ipt[3]) #设置能垒高度
dt  = float(ipt[4]) #设置时间步长

stp = int(ipt[5]) #设置时间步数


#设置能垒函数
def Fu(a):
    tmp = math.exp(-(a - 5) ** 2)
    Eu  = ssd * tmp
    return Eu

def Fup(a):
    tmp = math.exp(-(a - 5) ** 2)
    Eup = -2 * ssd * (a - 5) * tmp
    return Eup

#初始化
x = xin
p = pin

#用于作图的列表
xa = []
pa = []
Ua = []
Ta = []
Ea = []

#步长迭代
for i in range(stp):
    x0 = x
    p0 = p
    
    U  = Fu(x0)
    T  = p0 ** 2 / (2 * mas)
    E  = T + U
    
    x  = x0 + p0 * dt / mas - Fup(x0) * dt ** 2 / (mas * 2)
    
    U1 = Fup(x0)
    U2 = Fup(x)
    
    p  = p0 - dt * (U1 + U2) / 2 
    
    Ua.append(U)
    Ta.append(T)
    Ea.append(E)
    pa.append(p)
    xa.append(x)
    
#用作作图
Uan = np.array(Ua)
Tan = np.array(Ta)
Ean = np.array(Ea)
pan = np.array(pa)
xan = np.array(xa)

#用作势能函数作图
qa = stp * dt
pq = np.arange(0, qa, dt)

qb  = []

for i in pq:
    j = Fu(i)
    qb.append(j)
    
qc  = np.array(qb)

plt.subplot(234)
plt.xlabel('function')
plt.plot(pq, qc)

plt.subplot(231)
plt.title('totalU')
plt.plot(range(stp), Uan)

plt.subplot(232)
plt.title('totalT')
plt.plot(range(stp), Tan)

plt.subplot(233)
plt.title('totalE')
plt.plot(range(stp), Ean)

plt.subplot(235)
plt.xlabel('momentum')
plt.plot(range(stp), pan)

plt.subplot(236)
plt.xlabel('position')
plt.plot(range(stp), xan)

plt.show()
