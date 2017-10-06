# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 11:36:12 2017

@author: zongzi
"""
import math
import numpy as np
import matplotlib.pyplot as plt

rdate = open('inpt').readlines()

rx = np.array(rdate[0].split(','), dtype='float16')
rp = np.array(rdate[1].split(','), dtype='float16')

ms = float(rdate[2])
bi = float(rdate[3])
de = float(rdate[4])

st = int(rdate[5])

class Pot():
    
    def __inti__(self, bg):
        self.bg = bg
        
    def plf(self, mx, ny):
        return self.bg * math.exp(-mx ** 2 - ny ** 2)

    def plp(self, mx, ny):
        return -2 * mx * self.bg * math.exp(-mx ** 2 - ny ** 2)         
        
class Loop():

    def __init__(self, q1, q2, p1, ma, bh, dt, sp):
        self.q1 = q1
        self.p1 = p1
        self.q2 = q2
        self.sp = sp
        self.dt = dt
        self.ma = ma
        self.bh = bh
        
    def vv(self):
        
        bg  = self.bh
        puf = Pot(bg)
        
        x1 = self.q1
        p1 = self.p1
        y0 = self.q2
        
        xa = [x1]
        pa = [p1]
        
        for i in range(self.sp):
            x0 = x1
            p0 = p1
            
            x1 = x0 + p0 * self.dt / self.ma - puf.plp(x0, y0)  * self.dt ** 2 / (self.ma * 2)
            
            U1 = puf.plp(x0, y0)
            U2 = puf.plp(x1, y0)

            p1 = p0 - self.dt * (U1 + U2) / 2
            
            xa.append(x1)
            pa.append(p1)
        
        return xa, pa

lx = Loop(q1 = rx[0], q2 = rx[1], p1 = rp[0], ma = ms, bh = bi, dt = de, sp = st)
rxa, xpa = lx.vv()

ly = Loop(q1 = rx[1], q2 = rx[0], p1 = rp[1], ma = ms, bh = bi, dt = de, sp = st)
rya, ypa = ly.vv()

plt.plot(rxa, rya)
'''
a = np.arange(-3, 3, 0.001)
b = np.arange(-3, 3, 0.001)

X, Y = np.meshgrid(a, b)

#plt.contourf(X, Y, plf(X, Y), np.linspace(0.025, 3, 5), alpha = 0.1, cmap = plt.cm.hot)

C = plt.contour(X, Y, plf(X, Y), np.linspace(0.005, 2.95, 6), 
                colors = 'b')

plt.clabel(C, linewidth = 0.1)


plt.show()
'''