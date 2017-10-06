# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 20:33:16 2017

@author: zongzi
"""
import numpy as np
import matplotlib.pyplot as plt

'''
从inpt中读取数据，
第一行为初始坐标x，y
第二行为初始动量px，py
第三行为质量
第四行为势垒高度
第五行为步长
第六行为步数
'''
rdate = open('inpt').readlines()
#读取坐标值
rq = np.array(rdate[0].split(','), dtype='float16')
#读取动量值
rp = np.array(rdate[1].split(','), dtype='float16')

ma = float(rdate[2])  #读取质量
bh = float(rdate[3])  #读取势垒高度
dt = float(rdate[4])  #读取步长

stp = int(rdate[5])   #读取步数

#定义势能函数
def plfu(mx, ny):
    return bh * np.exp(-mx ** 2 - ny ** 2)

#定义势能函数的x偏导数
def plpx(mx, ny):
    return -2 * mx * bh * np.exp(-mx ** 2 - ny ** 2)

#定义势能函数的y偏导数
def plpy(mx, ny):
    return -2 * ny * bh * np.exp(-mx ** 2 - ny ** 2)

#读取坐标参数的x，y
x1 = rq[0]
y1 = rq[1]
#读取动量参数的px，py
px = rp[0]
py = rp[1]

#用于存取步长下的x，y，px，py的数组
x1a = []
y1a = []
pxa = []
pya = []

#步长迭代
for i in range(stp):
    
    #储存x和y的初值    
    x0 = x1
    y0 = y1
    
    #储存px和py的初值
    px0 = px
    py0 = py
    
    #利用泰勒展开，计算x+Δx
    x1 = x0 + px0 * dt / ma - plpx(x0, y0)  * dt ** 2 / (ma * 2)
    y1 = y0 + py0 * dt / ma - plpy(x0, y0)  * dt ** 2 / (ma * 2)
    
    #计算x方向的势能的偏导数的平均值
    Ux1 = plpx(x0, y0)
    Ux2 = plpx(x1, y0)
    
    #计算y方向的势能函数的平均值
    Uy1 = plpy(x0, y0)
    Uy2 = plpy(x0, y1)
    
    #利用不同方向下的势能的平均值计算p+Δp
    px = px0 - dt * (Ux1 + Ux2) / 2
    py = py0 - dt * (Uy1 + Uy2) / 2
    
    #储存每步的x，y，px，py
    x1a.append(x1)
    y1a.append(y1)
    pxa.append(px)
    pya.append(py)
    
#做出势能函数图与x、y的轨迹图    
rxa = np.array(x1a)
rya = np.array(y1a)

tra = plt.figure('tra')
plt.plot(rxa, rya)

a = np.arange(-1.5, 1.5, 0.001)
b = np.arange(-1.5, 1.5, 0.001)
X, Y = np.meshgrid(a, b)

C = plt.contour(X, Y, plfu(X, Y), 6, colors = 'black')
plt.clabel(C, linewidth = 0.1)
tra.show()  

#作出势能、动能与总能随时间（步长）变化图
xpan = np.array(pxa)
ypan = np.array(pya)

Tn = (xpan ** 2 + ypan ** 2) / (2 * ma)
Un = plfu(rxa, rya)
En = Tn + Un

ene = plt.figure('ene')
sup = np.arange(0, dt * stp, dt)

plt.plot(sup, Tn, 'red')
plt.plot(sup, Un, 'black')
plt.plot(sup, En, 'green')

ene.show()
