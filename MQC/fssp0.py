# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 19:15:18 2017

@author: zongzi
"""
import numpy as np
import matplotlib.pyplot as plt

class pot():

    def v11g(a):
        b = 0.01 * (1 - np.exp(- 1.6 * a))
        return b

    def v11l(a):
        b = -0.01 * (1 - np.exp(1.6 * a))
        return b

    def v22g(a):
        b = -pot.v11g(a)
        return b

    def v22l(a):
        b = -pot.v11l(a)
        return b

class sh():
    
    def a12():
        
    def a21():        
    

x1 = np.arange(-10, 0, 0.01)
x2 = np.arange(0, 10, 0.01)

y2 = pot.v11g(x2)
z2 = pot.v22g(x2)

y1 = pot.v11l(x1)
z1 = pot.v22l(x1)

plt.plot(x2, y2, 'black')
plt.plot(x2, z2, 'black')
plt.plot(x1, y1, 'black')
plt.plot(x1, z1, 'black')
plt.show()