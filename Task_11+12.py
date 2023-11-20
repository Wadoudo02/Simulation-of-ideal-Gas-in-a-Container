#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 10:13:23 2022

@author: wadoudcharbak
"""

from Ball import Ball
from Simulation import Simulation
import numpy as np
import pylab as pl
from scipy.optimize import curve_fit

'''

The following code was used to create Figures 5 & 6. 

If you want to reproduce figure 5, remove the arrays for radius_10 and radius_100.

'''

#r_20 = Simulation(5,container = Ball(100000,20,[0,0],[0,0]))
#r_30 = Simulation(5,container = Ball(100000,30,[0,0],[0,0]))
'''
n_5 = Simulation(5, ball_radius=0.1, container = Ball(100000,20,[0,0],[0,0]))
n_10 = Simulation(10, ball_radius=0.1, container = Ball(100000,20,[0,0],[0,0]))
n_15 = Simulation(15, ball_radius=0.1, container = Ball(100000,20,[0,0],[0,0]))
n_20 = Simulation(20, ball_radius=0.1, container = Ball(100000,20,[0,0],[0,0]))
n_25 = Simulation(25, ball_radius=0.1, container = Ball(100000,20,[0,0],[0,0]))
n_30 = Simulation(30, ball_radius=0.1, container = Ball(100000,20,[0,0],[0,0]))

Parameter_5 = n_5.run(100,animate = False, plots = False, plotcolour = "b")
Parameter_10 = n_10.run(100,animate = False, plots = False, plotcolour = "r")
Parameter_15 = n_15.run(100,animate = False, plots = False, plotcolour = "r")
Parameter_20 = n_20.run(100,animate = False, plots = False, plotcolour = "r")
Parameter_25 = n_25.run(100,animate = False, plots = False, plotcolour = "r")
Parameter_30 = n_30.run(100,animate = False, plots = False, plotcolour = "r")

This was my origionl way of doing things but it was inefficient and would take way toom long to scale up,
for this reason I transfered over to loops.
'''


simulation_of_balls_number = []
radius_10_times = []
radius_100_times = []

parameters_number = []
parameters_10 = []
parameters_100 = []

start = 1
end = 11

frames = 30

ball_number = 50

for i in range(start,end,1):
    #print(i)
    simulation_of_balls_number.append(Simulation(ball_number,temp_coef=i, ball_radius=0.02, container = Ball(100000,20,[0,0],[0,0])))
    radius_10_times.append(Simulation(ball_number,i, ball_radius=0.2, container = Ball(100000,20,[0,0],[0,0])))
    radius_100_times.append(Simulation(ball_number,i, ball_radius=2, container = Ball(100000,20,[0,0],[0,0])))

'''
for i in range(start,end,5):
    #print(i)
    simulation_of_balls_number.append(Simulation(i, ball_radius=0.01, container = Ball(100000,20,[0,0],[0,0])))
    radius_10_times.append(Simulation(i, ball_radius=0.1, container = Ball(100000,20,[0,0],[0,0])))
    radius_100_times.append(Simulation(i, ball_radius=1, container = Ball(100000,20,[0,0],[0,0])))
'''    

for i in range(len(simulation_of_balls_number)):
    parameters_number.append(simulation_of_balls_number[i].run(frames,animate = False, plots = False, plotcolour = "r"))
    parameters_10.append(radius_10_times[i].run(frames,animate = False, plots = False, plotcolour = "r"))
    parameters_100.append(radius_100_times[i].run(frames,animate = False, plots = False, plotcolour = "r"))
    


#ORDER [t,KE,P,T,Mom]
#%%
V = (np.pi*20**2)

def F_ideal_gas(T,V):
    P = ball_number*1.38e-23*T/V
    return P

N = list(range(start,end,1))
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

param1, param_cov1 = curve_fit(F_ideal_gas, T_sim, P_sim,(np.pi*20**2),absolute_sigma = True)

Fit1 = []

for t1 in T_sim:
    Fit1.append(F_ideal_gas(t1, param1[0]))

param10, param_cov10 = curve_fit(F_ideal_gas, T_10_sim, P_10_sim,(np.pi*20**2),absolute_sigma = True)

Fit10 = []

for t10 in T_10_sim:
    Fit10.append(F_ideal_gas(t10, param10[0]))

param100, param_cov100 = curve_fit(F_ideal_gas, T_100_sim, P_100_sim,(np.pi*20**2),absolute_sigma = True)

Fit100 = []

for t100 in T_100_sim:
    Fit100.append(F_ideal_gas(t100, param100[0]))

pl.figure(1)
pl.title("Pressure vs Temperature")
pl.plot(T_sim,P_sim,"o",label = "Simulation",color = "c")
pl.plot(T_sim,Fit1,label = "Fit 1",color="c")

pl.plot(T_10_sim,P_10_sim,"o",label = "Simulation 10r",color = "g")
pl.plot(T_10_sim,Fit10,label = "Fit 2",color="g")

pl.plot(T_100_sim,P_100_sim,"o",label = "Simulation 100r",color = "r")
pl.plot(T_100_sim,Fit100,label = "Fit 3",color="r")

#pl.errorbar(N,Fit,yerr=P_error,label = "Theory",mew=2, ms=3, capsize=4,color="c")
pl.xlabel("Temperature / K")
pl.ylabel("Pressure / Pa")
pl.legend()
pl.grid()
pl.savefig("Fig.png", dpi=300)

pl.show()


'''

param, param_cov = curve_fit(F_ideal_gas, N, P_100_sim,3e24,absolute_sigma = True)

Fit = []

for n in N:
    Fit.append(F_ideal_gas(n, param[0]))

perr = np.sqrt(np.diag(param_cov))

P_error = np.array(Fit)*perr[0]*1.38e-23/(np.pi*20**2)


pl.figure(1)
pl.title("Pressure vs Number of Balls")
pl.plot(N,P_sim,"o",label = "Simulation",color = "m")
pl.plot(N,P_10_sim,"o",label = "Simulation 10r",color = "g")
pl.plot(N,P_100_sim,"o",label = "Simulation 100r",color = "r")
#pl.errorbar(N,Fit,yerr=P_error,label = "Theory",mew=2, ms=3, capsize=4,color="c")
pl.plot(N,Fit,label = "Theory",color="c")
pl.xlabel("Number of Balls")
pl.ylabel("Pressure / Pa")
pl.legend()
pl.grid()
pl.show()

print("The value for temperature here is", param[0],"±", perr[0],"K.")
'''  

#%%


