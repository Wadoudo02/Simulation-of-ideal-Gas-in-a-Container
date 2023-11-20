#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 10:11:38 2022

@author: wadoudcharbak
"""


'''
%matplotlib inline – for rendering in the plot browser.
%matplotlib auto – for rendering in a window. This mode enables animations to play.
'''

import numpy as np
import pylab as pl
import random as rd
from Ball import Ball


        
        
class Simulation(Ball):
    
    #__doc__ = 
    '''
    This is the class contains and runs a simulation of balls in a container,
    it has attributes of Number of Balls, Temperature coef, Ball radius and container.
    '''
    
    container_collisions = 0
    container_impulse = 0
    time_running = 0
    
    def __init__(self, num_of_balls, temp_coef = 1, ball_radius = 1, container = Ball(mass = 1000000, radius = 20, v = [0,0], r = [0,0])):
        '''
        Initialises the similation.

        Parameters
        ----------
        num_of_balls : Integer 
            This is the number of balls for the simulation to run, these are arranged systematically depending on the number of balls.
        temp_coef : Float, optional
            A coefficient that changes the width of the gaussian of the normal distribution of velocity. It changes the velocity of the balls whilst keeping the average
            of thier speeds 0, changing the Temperature of the whole system. The default is 1.
        ball_radius : Float, optional
            Controls the radius of the balls in the simulation.. The default is 1.
        container : Ball, optional
            The container you want to use for your experiment. The default is Ball(mass = 1000000, radius = 20, v = [0,0], r = [0,0]).

        Returns
        -------
        None.

        '''
        
        #num_of_balls = len(total_balls)
        angle_basis = 2*np.pi / num_of_balls
        
        self.__temp_coef = temp_coef
        self._ball =  []
        
        random_v = np.random.normal(0, temp_coef,size=(2,num_of_balls)) #Random Velocities with a mean of 0.
      
        
        
        for i in range(num_of_balls):
            self._ball.append(Ball(mass = 1, radius = ball_radius, v = [random_v[0][i],random_v[1][i]], r = 10*np.array([np.cos(i*angle_basis),np.sin(i*angle_basis)])))
        
        self._container = container

    def __repr__(self):
        return "Simulation of %s balls with temperture coefficient = %s." % (len(self._ball), self.__temp_coef)
    
    def Velocities(self):
        '''
        Returns an array of all the velocities of all the balls in the simulaton. Used for histogram analysis.

        Returns
        -------
        v : Array
            An array of all the velocities of all the balls in the simulaton.

        '''
        v=[]
        for ball in self._ball:
            v.append(ball.vel())
        return v

    def next_collision(self):

        '''
        This calculates the time to the next collision of the balls, moves them to that position and collides them.
        
        Changes
        -------
        
        This changes the velocities, positions of all the balls and container. 

        Returns
        -------
        None.

        '''
        time = [[] for _ in range(len(self._ball))] #2D Array with each position being collision times of each ball to every other ball.
        for i in range(len(self._ball)):
            for j in range(len(self._ball)):
                if i == j: # This is to avoid the balls attempting to collide with themselves  
                    time[i].append(10000) 
                else: 
                    time[i].append(self._ball[i].time_to_collision(self._ball[j]))
            time[i].append(self._ball[i].time_to_collision(self._container, True))
            #print("---")
        time = np.array(time)
        time[time<0] = 10000
        time_min = np.nanmin(time)
        #print(time_min)
        minpos = [np.where(time == time_min)[0][0], np.where(time == time_min)[1][0]] # this find the minimum time position  
        '''
        time_min = 999.9
        minpos = [0,0]
        for i in range(len(time)):
            if np.nanmin(time[i]) < time_min  and np.nanmin(time[i]) > 0:
                minpos = [i,time[i].index(min(time[i]))] #Find minimum position
                time_min = np.nanmin(time[i])
        '''
        
        
        #print(time)
        #print(time_min)
        #print(minpos)
        
      
        for ball in self._ball:
            ball.move(time_min)
        self._container.move(time_min)
        
        Simulation.time_running += time_min
        
        if minpos[1] == (len(self._ball)):
            c_mom_before = self._container.momentum()
            self._ball[minpos[0]].collision(self._container)
            c_mom_after= self._container.momentum()
            #Calculating the momentum before and after to calculate the change in momentum
            delta_m = c_mom_before - c_mom_after
            Simulation.container_impulse += (np.sqrt(delta_m[0]**2 + delta_m[1]**2))
            Simulation.container_collisions += 1
            
        else:
            self._ball[minpos[0]].collision(self._ball[minpos[1]])
        
    
    def run(self, num_frames, animate = False, plots = True, plotcolour = "b"):

    
        '''
        This actually runs and animates the simulation.

        Parameters
        ----------
        num_frames : int
            This controls the number of frames which corresponds to the number of collisions.
            This also controls number of frames animated if desired.
        animate : Boolean, optional
            This controls if the balls are animated on a plot. The default is False.
        plots : Boolean, optional
            This controls if plots are animated or not. The default is True.
        plotcolour : String, optional
            This controls the plot colour. The default is "b".

        Returns
        -------
        An animation of the balls in the container if animate is set to True. 
        
        Returns plots for:
            - Distance between each pair of balls
            - Distance of Balls from centre of Container
            - KE verses time
            - Momentum verses time
            
        if plots is set to True.
        
        Returns a list of parameters of the frame number, KE, Pressure, Temperature and Momentum of the system in the Simulation.

        '''
        t = []
        KE = []
        #P = []
        Mom = []
        h = []  
        
        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-self._container.radius-1, 
            self._container.radius+1), ylim=(-self._container.radius-1, self._container.radius+1))
            ax.set_aspect('equal', adjustable='box')
            ax.add_artist(self._container.get_patch(True))
            colours = ['r','y','g','w','b'] 
            for i in range(len(self._ball)):
                ax.add_patch(self._ball[i].get_patch())
        
        for frame in range(num_frames):
              
            
            total_KE = 0
            total_Mom = 0
            for ball in self._ball:
                total_KE += ball.KE()
                total_Mom += (np.sqrt(ball.v[0]**2 + ball.v[1]**2) * ball.mass)
                h.append(ball.distance_from_centre())
            KE.append(total_KE)
            t.append(frame)
            Mom.append(total_Mom)
            #P.append((np.sqrt(self._container.v[0]**2 + self._container.v[1]**2) * self._container.mass) / Simulation.time_running)
            
            if animate:
                for i in range(len(self._ball)):
                    ax.add_patch(self._ball[i].get_patch())
                pl.pause(0.1)
            
            self.next_collision() 
            
        T=[]
        for e in KE:
            T.append(e/(len(self._ball)*1.38e-23)) # From Ideal Gas equation with 2 degrees of freedom
        #print("The Temperature is", T)
        
          
        if animate:
            
            pl.show()
            
        '''
        for i in range(len(self._ball)):
            h.append(self._ball[i].distance_from_centre())
        '''
        d = []
        for i in range(int(len(self._ball))):
            for j in range(len(self._ball)):
                if i == j:
                    d.append(-5) 
                else: 
                    d.append(self._ball[i].distance_from_ball(self._ball[j]))   
              
        if plots:
            pl.figure(2)
            pl.title("Distance of Balls from Centre of Container")
            pl.hist(h,color=plotcolour)
            pl.ylabel("Number of Balls")
            pl.xlabel("Distance (m)")
            pl.grid()
            pl.savefig("Fig 1.png", dpi=300)
            
            
            for i in d:
                if i<0:
                    d.remove(i)
                    
                    
            pl.figure(3)
            pl.hist(d,20,color=plotcolour)
            pl.ylabel("Number of Balls")
            pl.xlabel("Distance (m)")
            pl.title("Distance between each Pair of Balls")
            pl.grid()
           # pl.savefig("Fig 2.png", dpi=300)
            
            
            '''
            
            pl.figure(4)
            pl.ylim(0,0.5)
            pl.xlabel("Time / s")
            pl.ylabel("Kinetic Energy / J")
            pl.plot(t,KE,color=plotcolour)
            pl.grid()
            pl.title("KE verses time")
            
         
            pl.figure(6)
            #pl.ylim(0,0.5)
            pl.xlabel("Time / s")
            pl.ylabel("Momentum / Kgms^-1")
            pl.plot(t,Mom,color=plotcolour)
            pl.grid()
            pl.title("Momentum verses time")
            '''
            
            pl.show()
            
       # print(Simulation.time_running)
        
        P = Simulation.container_impulse / Simulation.time_running
        #print((np.sqrt(self._container.v[0]**2 + self._container.v[1]**2) * self._container.mass) / Simulation.time_running)
        
        #print(P)
        
        return [t,KE,P,T,Mom]
   
#%%

#Tests for code

'''
test = Simulation(5)

test.run(20,animate = False,plotcolour="r")

        
#%%
m = Ball(1,1,[1,0],[-5,0])
n = Ball(1,1,[-1,0],[5,0])

m.time_to_collision(n)
  '''  
#%%
'''
r_2r = Ball(1,2,[1,6],[-10,-10])
r = Ball(1,1,[-1,-1],[5,5])

a = Ball(1,1,[1,0],[-3,0])
b = Ball(1,1,[-1,-1],[3,3])

a.time_to_collision(b,True)

ball = Ball(1,1,[-1,-2.5],[0,-10])

q = Ball(1,1,[1,0],[-5,0])
r = Ball(1,1,[-1,0],[5,0])
s = Ball(1,1,[0,1],[0,-5])
t = Ball(1,1,[0,-1],[0,5])

fourballtest = [q,r,s,t]

m = Ball(1,1,[1,0.1],[-5,0])
n = Ball(1,1,[-1,0],[5,0])

twoballtest = [m,n]

sus = [a,b]
sus = [r,r_2r,a,b]

container = Ball(100000,20,[0,0],[0,0])

test = Simulation(3,container)

test.run(15,animate = True)


'''
'''
Numpy precision, if 2 happen AT THE SAME TIME then NEED TO MAKE COLLISIONS SIMULTANIUS
'''


test = Simulation(5)

test.run(20,animate = True,plots=False)


    

