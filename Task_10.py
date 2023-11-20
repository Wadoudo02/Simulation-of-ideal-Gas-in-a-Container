#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 23:05:38 2022

@author: wadoudcharbak
"""

'''

The following code below was used for the investigation that produced figures 1,2,3 & 4.

'''

from Ball import Ball
from Simulation import Simulation
import numpy as np
import pylab as pl

Sim = Simulation(100, ball_radius=0.01, container = Ball(100000,20,[0,0],[0,0]))
Parameters = Sim.run(100,animate = False, plots = True, plotcolour = "r")

#%%

'''
vel = []
for vs in Sim.Velocities():
    vel.append(np.sqrt(vs[0]**2 + vs[1]**2))

'''

plotcolour = 'm'

pl.figure(1)
pl.ylim(0,100)
pl.xlabel("Time / s")
pl.ylabel("Kinetic Energy / J")
pl.plot(Parameters[0],Parameters[1],color="g")
pl.grid()
pl.title("Total KE verses time")
#pl.savefig("Fig 3.png", dpi=300)


pl.figure(2)
pl.ylim(0,150)
pl.xlabel("Time / s")
pl.ylabel("Momentum / Kgms^-1")
pl.plot(Parameters[0],Parameters[4],color=plotcolour)
pl.grid()
pl.title("Total Momentum verses time")
#pl.savefig("Fig 4.pdf", dpi=300)


pl.show()

# [t,KE,P,T,Mom]