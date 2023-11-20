#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 22:20:51 2022

@author: wadoudcharbak
"""

from Ball import Ball
from Simulation import Simulation
import numpy as np
import pylab as pl
from scipy.optimize import curve_fit

'''

The following code was used to create Figure 8. 


'''
simulation_of_balls_number = []
radius_10_times = []
radius_100_times = []

parameters_number = []
parameters_10 = []
parameters_100 = []

start = 1
end = 201

frames = 50

TC = 10

for i in range(start,end,5):
    #print(i)
    simulation_of_balls_number.append(Simulation(i, temp_coef=TC, ball_radius=0.01, container = Ball(100000,20,[0,0],[0,0])))
    radius_10_times.append(Simulation(i, temp_coef=TC, ball_radius=0.1, container = Ball(100000,20,[0,0],[0,0])))
    radius_100_times.append(Simulation(i, temp_coef=TC, ball_radius=1, container = Ball(100000,20,[0,0],[0,0])))
    
for i in range(len(simulation_of_balls_number)):
    parameters_number.append(simulation_of_balls_number[i].run(frames,animate = False, plots = False, plotcolour = "r"))
    parameters_10.append(radius_10_times[i].run(frames,animate = False, plots = False, plotcolour = "r"))
    parameters_100.append(radius_100_times[i].run(frames,animate = False, plots = False, plotcolour = "r"))
    


#ORDER [t,KE,P,T,Mom]
#%%
def F_ideal_gas(N,T):
    P = N*1.38e-23*T/(np.pi*20**2)
    return P

def Van_der_walls(N,T,a,b):
    P = N*1.38e-23*T/(np.pi*20**2 - N*b) - a*((N/np.pi*20**2)**2)
    return P

N = list(range(start,end,5))
P_theory = []

P_sim = [] 
T_sim = []

P_10_sim = []
T_10_sim = []

P_100_sim = []
T_100_sim = []

for g in parameters_number:
    P_sim.append(g[2])
    T_sim.append(np.mean(g[3]))


for f in parameters_10:
    P_10_sim.append(f[2])
    T_10_sim.append(np.mean(f[3]))
    
for h in parameters_100:
    P_100_sim.append(h[2])
    T_100_sim.append(np.mean(h[3]))

#for n in N:
#    P_theory.append(F_ideal_gas(n, 6e24))

param, param_cov = curve_fit(Van_der_walls, N, P_100_sim,[3e26,0,1],absolute_sigma = True)

Fit = []

for n in N:
    Fit.append(Van_der_walls(n, param[0],param[1],param[2]))

perr = np.sqrt(np.diag(param_cov))

P_error = np.array(Fit)*perr[0]*1.38e-23/(np.pi*20**2)

pl.figure(1)
pl.title("Pressure vs Number of Balls")
pl.plot(N,P_sim,"o",label = "Simulation",color = "m")
pl.plot(N,P_10_sim,"o",label = "Simulation 10r",color = "g")
pl.plot(N,P_100_sim,"o",label = "Simulation 100r",color = "r")
#pl.errorbar(N,Fit,yerr=P_error,label = "Theory",mew=2, ms=3, capsize=4,color="c")
pl.plot(N,Fit,label = "Van der Waal's Law",color="c")
pl.xlabel("Number of Balls")
pl.ylabel("Pressure / Pa")
pl.legend()
pl.grid()
pl.savefig("Fig 8.png", dpi=300)
pl.show()


print("The value for temperature here is", param[0],"±", perr[0],"K.")
print("The value for a here is", param[1],"±", perr[1])
print("The value for b here is", param[2],"±", perr[2])
