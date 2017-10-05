# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 04:13:36 2017

@author: zongzi
"""
import math
import numpy as np
import matplotlib.pyplot as plt

fm = -m * k * x ** 2
rm = -2 * m * math.log(rou) / beta


cop = math.exp(- beta * k * x ** 2)

if rou < cop:
    print("T is too small or k is too big")
    
def pv(a):
    vu = k / 2 * a ** 2
    return vu

def pvp(a):
    vup = k * a
    return vup
    

k = 8
x = np.arange(-5, 5, 0.001)
V = k / 2 * x * x

plt.plot(x, V)
plt.show() 
