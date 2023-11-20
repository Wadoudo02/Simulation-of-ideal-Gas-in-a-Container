#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 09:57:41 2022

@author: wadoudcharbak
"""
from Ball import Ball
from Simulation import Simulation
import numpy as np
import pylab as pl
from scipy.optimize import curve_fit

'''

The following code was used to create Figure 7.

''' 


#Velocity Test

Velocity_Sim = Simulation(100, ball_radius=0.01, container = Ball(100000,20,[0,0],[0,0]))
Parameters = Velocity_Sim.run(10000,animate = False, plots = False, plotcolour = "r")

vel = []
for vs in Velocity_Sim.Velocities():
    vel.append(np.sqrt(vs[0]**2 + vs[1]**2))
    
#%%
T = np.mean(Parameters[3])

v_x = []
v_y = []

for v in Velocity_Sim.Velocities():
    v_x.append(v[0])
    v_y.append(v[1])
    
x_points = np.linspace(0, 3.5,200)
maxwell = []

proportion_factor = 35

for v in x_points:
    maxwell.append(proportion_factor*v*np.exp(-0.5*v**2 / (1.38e-23*T)))

pl.figure(2)
pl.title("Speed Distribution vs Maxwell-Boltzmann Distribution")
pl.hist(vel,10,color="c",label = "Speed Distribution")

data       = vel
y,binEdges = np.histogram(data,bins=10)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
menStd     = np.sqrt(y)
width      = 0.05
pl.errorbar(bincenters, y, yerr=menStd, fmt="o", mew=2, ms=3, capsize=4,color='orange',)

pl.plot(x_points,maxwell,label = "MB Distribution",color="r")
pl.xlabel("Velocity m/s")
pl.ylabel("Relative Intensity (Arb units)")
pl.legend()
pl.grid()
#pl.savefig("Fig.png", dpi=300)



pl.figure(3)
pl.title("Velocity_X Distribution")
pl.hist(v_x,10,color="r")
pl.xlabel("Velocity m/s")
pl.ylabel("Relative Intensity (Arb units)")
pl.grid()

pl.figure(4)
pl.title("Velocity_Y Distribution")
pl.hist(v_y,10,color="r")
pl.xlabel("Velocity m/s")
pl.ylabel("Relative Intensity (Arb units)")
pl.grid()

