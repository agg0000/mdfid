# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 20:33:16 2017

@author: zongzi
"""
import numpy as np
import matplotlib.pyplot as plt

rdate = open('inpt').readlines()

rq = np.array(rdate[0].split(','), dtype='float16')
rp = np.array(rdate[1].split(','), dtype='float16')
ma = float(rdate[2])
bh = float(rdate[3])
dt = float(rdate[4])

stp = int(rdate[5])

def plf(mx, ny):
    return bh * np.exp(-mx ** 2 - ny ** 2)

def plp(mx, ny):
    return -2 * mx * bh * np.exp(-mx ** 2 - ny ** 2)

x1 = rq[0]
y1 = rq[1]
px = rp[0]
py = rp[1]

x1a = []
y1a = []
pxa = []
pya = []

for i in range(stp):
    
    x0 = x1
    y0 = y1
    
    px0 = px
    py0 = py
    
    x1 = x0 + px0 * dt / ma - plp(x0, y0)  * dt ** 2 / (ma * 2)
    y1 = y0 + py0 * dt / ma - plp(x0, y0)  * dt ** 2 / (ma * 2)
    
    U1 = plp(x0, y0)
    U2 = plp(x1, y0)
    U3 = plp(x0, y1)
    
    px = px0 - dt * (U1 + U2) / 2
    py = py0 - dt * (U1 + U3) / 2
    
    x1a.append(x1)
    y1a.append(y1)
    pxa.append(px)
    pya.append(py)
    
rxa = np.array(x1a)
rya = np.array(y1a)

tra = plt.figure('tra')
plt.plot(rxa, rya)

a = np.arange(-1.5, 1.5, 0.001)
b = np.arange(-1.5, 1.5, 0.001)
X, Y = np.meshgrid(a, b)

C = plt.contour(X, Y, plf(X, Y), 6, colors = 'black')
plt.clabel(C, linewidth = 0.1)
tra.show()  

xpan = np.array(pxa)
ypan = np.array(pya)

Tn = (xpan ** 2 + ypan ** 2) / (2 * ma)
Un = plf(rxa, rya)
En = Tn + Un

ene = plt.figure('ene')
sup = np.arange(0, dt * stp, dt)

plt.subplot(221)
plt.plot(sup, Tn, 'red')

plt.subplot(222)
plt.plot(sup, Un, 'black')

plt.subplot(223)
plt.plot(sup, En, 'green')

ene.show()
