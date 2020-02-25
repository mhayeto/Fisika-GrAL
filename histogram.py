# -*- coding: utf-8 -*-
"""
Date: 2020/02/18
Description: This program manages the data from the histogram of Q in order to 
             plot it more correctly.
"""

import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import numpy as np
import sys

# -----------------------------------------------------------------------------
#                                LaTex
# -----------------------------------------------------------------------------
matplotlib.rcParams['text.usetex'] = True
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# -----------------------------------------------------------------------------
#                            Open histogram
# -----------------------------------------------------------------------------
#filename = "test"
#with open(filename+".txt") as f:
#    lines = f.readlines()
#x = [float(line.split()[0]) for line in lines]
#y = [float(line.split()[1]) for line in lines]

# -----------------------------------------------------------------------------
#                            Create histogram
# -----------------------------------------------------------------------------
filename = "histogram_T3_p1_d06"
with open(filename+".txt") as f:
    lines = f.readlines()
x = [float(line.split()[0]) for line in lines]

def count_val(x, val):
    c = 0
    for i in range(len(x)):
        if x[i]==val:
            c += 1
    return c

def histogram(x, r):
#   It creates the histogram of the given array.
    x = np.around(x, r)
    x = x.tolist()
    xh = []
    yh = []
    l = len(x)
    for i in range(l):
        xh.append(x[i])
        yh.append(count_val(x, x[i])/l)
    d = {k: v for k, v in zip(xh, yh)}
    xh, yh = zip(*d.items())
    return xh, yh

x, y = histogram(x, 3)

# -----------------------------------------------------------------------------
#                                Filter
# -----------------------------------------------------------------------------
def select_single_sign(x, y, sign):
#   Selects Q values corresponding to a chosen sign (positive or negative)
    xf = []
    yf = []
    if sign=="pos":
        for i in range(len(x)):
            if x[i]>=0.0:
                xf.append(x[i])
                yf.append(y[i])
    elif sign=="neg":
        for i in range(len(x)):
            if x[i]<=0.0:
                xf.append(x[i])
                yf.append(y[i])
    return xf, yf

def take_abs_x(x,y):
#   It takes the abs. value of all the Q values. 
    xf = []
    yf = []
    for i in range(len(x)):
        xf.append(abs(x[i]))
        yf.append(y[i])
    return xf, yf

def select_with_step(x, y, div):
#   Selects Q values dividing its amount by 'div'.
    xf = []
    yf=[]
    for i in range(0, len(x), div):
        xf.append(x[i])
        yf.append(y[i])
    return xf, yf

def remove_any(x, y, xremove):
#   Removes a selected point from the data set.
    ind = x.index(xremove)
    x.remove(xremove)
    y.remove(y[ind])
    
xf, yf = select_single_sign(x, y, 'pos')
#xf, yf = take_abs_x(x, y)
#xff, yff = select_with_step(xf, yf,  1)

# -----------------------------------------------------------------------------
#                          Mirror
# -----------------------------------------------------------------------------
def mirror_graph(x, y):
#   Doubles the data in order to obtain a mirror image.
    x2, y2 = [], []
    xm, ym = [], []
    l = len(x)
    for i in range(l-1):
        x2.append(-x[l-i-1])
        y2.append(y[l-i-1])
    xm = x2 + x
    ym = y2 + y
    return xm, ym

xm, ym = mirror_graph(xf, yf)

# -----------------------------------------------------------------------------
#                       Interpol / Smooth
# -----------------------------------------------------------------------------
def interpol_smooth(x, y, points):
#   Interpolates and smoothes the given data.
    xx = np.linspace(min(x), max(x), points)
    itp = interp1d(x,y, kind='linear')
    window_size, poly_order = 101, 4
    yy_sg = savgol_filter(itp(xx), window_size, poly_order)
    return xx, yy_sg

xx, yy_sg = interpol_smooth(x, y, 325)
xxm, yy_sgm = interpol_smooth(xm, ym, 475) #475

# -----------------------------------------------------------------------------
#                                Plot
# -----------------------------------------------------------------------------
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, axisbelow=True)
ax.plot(x, y, 'c', alpha=0.5, lw=0.0, label="All", marker='o', markersize=4.0)
ax.plot(xx, yy_sg, 'tomato', alpha=0.5, lw=0.9, label="Interp. + Smooth")
ax.plot(xxm, yy_sgm, 'r', alpha=0.7, lw=0.9, label="Mirrored")
ax.set_xlabel('Order Parameter')
ax.set_ylabel('Probability Distribution')
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.set_title('$C/E_{0}=1$, $E_{0}=1$, $k_{B}T = 3E_0$, $d=0.6$')  
#ax.set_xlim(-1,1)
ax.grid(b=True, which='major', c='silver', lw=0.5, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)

for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)

savefile = filename+".png"
plt.savefig(savefile, dpi=500)
plt.show()

# -----------------------------------------------------------------------------
#                       Landau Free Energy
# -----------------------------------------------------------------------------
def landau_free_energy(yy_sg, beta, N, E0):
#   Based on the previously calculated Probability Distribution (histogram), it
#   calculates the Landau Free Energy corresponding to each Q value.
    y_landau = []
    for P in yy_sg:
        y_landau.append(-np.log(P)/(beta*N*E0))
    return y_landau

#y_landau = landau_free_energy(yy_sg, 1.0/4.0, 1000.0, 10.0)    
#
#fig = plt.figure(facecolor='w')
#ax = fig.add_subplot(111, axisbelow=True)
#ax.plot(xx, y_landau, 'orange', alpha=0.7, lw=0.9, label="$f_{L}(Q)$")
#ax.set_xlabel('Order Parameter')
#ax.set_ylabel('Landau Free Energy')
#ax.yaxis.set_tick_params(length=0)
#ax.xaxis.set_tick_params(length=0)
#ax.set_title('$C/E_{0}=1$, $E_{0}=10$, $k_{B}T = 4E_0$')  
#
#ax.grid(b=True, which='major', c='silver', lw=0.5, ls='-')
#legend = ax.legend()
#legend.get_frame().set_alpha(0.5)
#
#for spine in ('top', 'right', 'bottom', 'left'):
#    ax.spines[spine].set_visible(False)
#
#savefile = filename+".png"
#plt.savefig(savefile, dpi=500)
#plt.show()