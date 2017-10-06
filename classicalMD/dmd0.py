# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 11:36:12 2017

@author: zongzi
"""

import numpy as np
import matplotlib.pyplot as plt

def plf(x, y):
    return 3 * np.exp(-x ** 2 - y ** 2)

a = np.arange(-3, 3, 0.001)
b = np.arange(-3, 3, 0.001)

X, Y = np.meshgrid(a, b)

plt.contourf(X, Y, plf(X, Y), np.linspace(0.025, 3, 5), alpha = 0.1, cmap = plt.cm.hot)

C = plt.contour(X, Y, plf(X, Y), np.linspace(0.025, 2.999999, 5), 
                colors = 'black', linewidth = 0.1)

plt.clabel(C, linewidth = 0.1)

plt.show()