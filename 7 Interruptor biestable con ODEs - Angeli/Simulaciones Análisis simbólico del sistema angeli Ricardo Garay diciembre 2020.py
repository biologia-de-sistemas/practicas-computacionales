#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 16:34:48 2020

@author: richie
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
alpha1 = 1
alpha2 = 1
beta1 = 200
beta2 = 10
k1 = 30
k2 = 1
v = 1
gamma1 = 4
gamma2 = 4
######################
a1 = -((alpha1**4) * (alpha2 + beta2)*v**(16)) - alpha2*k2*(alpha1 + beta1)**4*v**6
a2 = 1.632240802e9
a3 = -974.47344e6
a4 = 974.4724e6
a5 = -218.2166e6
a6 = 218.1626e6
a7 = -22.896e6
a8 = 21.816e6
a9 = -9.72e6
a10 = 1.62e6
########################################
c = [a1,a2,0,0,a3,a4,0,0,a5,a6,0,0,a7,a8,0,0,a9,a10]
########################################
def re(x):
    t = []
    for i in range(len(x)):
        if x[i].imag == 0:
            r = x[i].real
            t.append(r)
    return t
########################################
p = re(np.roots(c))
l=[]
for i in range(len(p)):
    w = (((p[i] - 1)*alpha2*k2)/(alpha2 - (alpha2 + beta2)*p[i]))**(1/4)
    l.append(w)
y0 = (alpha2*(k2 + 1))/(beta2 + (k2 + 1)*alpha2)
x0 = (alpha1*(k1 + (v**(gamma1))))/(beta1*(v**(gamma1)) + alpha1*(k1 + (v**(gamma1))))
E = [[l[0],p[0]],[l[1],p[1]],[l[2],p[2]]]
#####################
i =2################Selecciona el punto de equilibrio i= 0,1,2,3,4
#####################
x = E[i][0]
y = E[i][1]
print("punto de equilibrio E_%1i es (%3.5f,%3.5f)"%(i,x,y))
print("\n")
A1 = alpha1 + alpha2 + (beta2*x**(gamma2))/(k2 + x**(gamma2)) + (beta1*(v*y)**(gamma1))/(k1 + (v*y)**(gamma1))
A2 = ((beta1*beta2*gamma1*gamma2*k2*(v**(2*gamma1-1))*(x**(gamma2))*(y**(2*gamma1)))/(((k1 + (v*y)**(gamma1))**2)*((k2 + x**(gamma2))**2)))*((v-1) - k1/((v*y)**(gamma1))) + ((alpha1 + ((beta1 * (v*y)**(gamma1)) / (k1 + (v*y)**(gamma1)))) * (alpha2 + ((beta2 * x**(gamma2))/(k2 + x**(gamma2)))))                                                                                                                                                         
s = [1,A1,A2]
print("Polinomio característico es l^2 +%3.3f l + %3.3f"%(A1,A2))
roots = np.roots(s)
if roots[0]*roots[1]>0:
    print("\n")
    print("Estable")
    print("Las raíces del polinomio caracteristico son: ",roots)
else:
    print("No estable")
    print(roots)
##############################################
###                 PLOT                   ###
##############################################
def f(state,t):
    x, y = state
    xdot = alpha1*(1 - x) - (beta1 * x * (v*y)**(gamma1))/(k1 + (v*y)**(gamma1))
    ydot = alpha2*(1 - y) - (beta2 * y * (x**(gamma2)))/(k2 + (x**(gamma2)))
    return xdot,ydot
t = np.arange(0.0, 1000, 0.001)
fig2 = plt.figure(figsize=(8,6))
ax4 = fig2.add_subplot(1,1,1)
    # quiverplot
    # define a grid and compute direction at each point
x = np.linspace(0, 2, 35)
y = np.linspace(0, 2, 35)
X1 , Y1  = np.meshgrid(x, y)                    # create a grid
DX1, DY1 = f([X1, Y1],t)                        # compute growth rate on the grid
M = (np.hypot(DX1, DY1))                        # norm growth rate 
M[ M == 0] = 1.                                 # avoid zero division errors 
DX1 /= M                                        # normalize each arrows
DY1 /= M
ax4.quiver(X1, Y1, DX1, DY1, M, pivot='mid')
############################
y0 = [2,2]#Initial Condition
EP2 = E[2]
ax4.scatter(y0[0], y0[1], color='red', linewidth = 5, zorder=20, label = 'Condición inicial $(%8.3f,%8.3f)$'%(y0[0], y0[1]))
y = odeint(f, y0, t)
ax4.plot(y[:,0], y[:,1], color='dodgerblue', linewidth=4, zorder=10, label = 'Trayectoria')
ax4.scatter(EP2[0], EP2[1], color='green', linewidth = 15, zorder=20, label = 'Punto de Equilibrio $E^*=(%8.4f,%8.0f)$'%(EP2[0], EP2[1]))
#############################
x0 = [2,1.25]
EP1 = E[0]
ax4.scatter(EP1[0], EP1[1], color='green', linewidth = 15, zorder=20, label = 'Punto de Equilibrio $E^*=(%8.0f,%8.4f)$'%(EP1[0], EP1[1]))
x = odeint(f, x0, t)
ax4.scatter(x0[0], x0[1], color='red', linewidth = 5, zorder=20, label = 'Condición inicial $(%8.3f,%8.3f)$'%(x0[0], x0[1]))
ax4.plot(x[:,0], x[:,1], color='dodgerblue', linewidth=4, zorder=10, label = 'Trayectoria')
#############################
j = 1
ax4.scatter(l[j], p[j], color='black', linewidth = 35, zorder=20, label = 'Punto de Equilibrio $E^*=(%8.0f,%8.4f)$'%(l[j], p[j]))
#############################
print(l[j], p[j])