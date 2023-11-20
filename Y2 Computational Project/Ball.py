#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 10:11:38 2022

@author: wadoudcharbak
"""

import numpy as np
import pylab as pl
import random as rd

class Ball:
    #__doc__ = 
    '''
    This is the class that hold infomation for a ball,
    it has attributs of Mass, Radius, Velocity and Position.
    '''
    
    
    def __init__(self, mass = 1, radius = 1, v = [0,0], r = [0,0]):
        '''
        Initialises a ball.

        Parameters
        ----------
        mass : Float, optional
            The Mass of the Ball. The default is 1.
        radius : Float, optional
            The Radius of the Ball. The default is 1.
        v : Array, optional
            An array of velocity, the first value is the velocity in the x direction. 
            The second is the veloity in the y direction. The default is [0,0].
        r : Array, optional
            An array of position, the first value is the position in x. 
            The second is the position in  y. The default is [0,0].

        Raises
        ------
        Exception
            If the r vector or v vector has an array larger 1 x 2, i.e it is in a Higher dimention than 2D.

        Returns
        -------
        None.

        '''
        
        if len(r) > 2:
            raise Exception("Parameter r has incorrect size, please go back to 2D.")
        if len(v) > 2:
            raise Exception("Parameter v has incorrect size, please go back to 2D.")
            
        self.mass = mass
        self.radius = radius
        self.r = np.array(r,dtype = float)
        self.v = np.array(v,dtype = float)
        
    def __repr__(self):
        return "M = %g, R = %g, v = %s, r = %s" % (self.mass, self.radius,self.v,self.r)
        
    def __str__(self):
        return "Mass = %g, Radius = %g, Velocity = array(%s), Position = array(%s)" % (self.mass, self.radius,self.v,self.r)
    
    
    def pos(self):
        '''
        Returns the position of the ball.

        Returns
        -------
        Array
            Returns a 2D array of the position of the ball.
            Position 0 of the array is the ball's position in the x axis.
            Position 1 of the array is the ball's position in the y axis.

        '''
        return self.r
    
    def vel(self):
        '''
        Returns the velocity of the ball.

        Returns
        -------
        Array
            Returns a 2D array of the velocity of the ball.
            Position 0 of the array is the ball's velocity in the x direction.
            Position 1 of the array is the ball's velocity in the y direction.

        '''
        return self.v
    
    def move(self, dt):
        '''
        Moves the ball by a distance, dt * v.

        Parameters
        ----------
        dt : Float
            The time dt you want to move the ball by.

        Changes
        -------
        Changes the position of the ball according to the time dt and the velocity of the ball.

        '''
        self.r += self.v*dt
        
    def time_to_collision(self, other, other_is_container = False):
        
        '''
        Calculates the time taken for 2 balls (or container) to collide.
        
        Parameters
        ----------
        other : Ball
            This is the other ball (or container) you wish to calculate the time taken to collide with the initial ball.
            
        other_is_container : Boolean, optional
            This asks you if the other ball you are colliding with is the container. 
            This is important as the mathmatical calculations differ depending on if the ball is colliding within a container or to another ball.
            The default is False.

        Returns
        -------
        Float Variable
            The returned value will be the time taken for the balls to colide.
            
        Time to collide should return a positive real value, as a negitive time cannot exist, 
        to solve this it must take the maximum of the +/- values of the square root.

        '''
        
        if other_is_container == False:
            Radius_col = self.radius + other.radius
        else:
            Radius_col = self.radius - other.radius
            
        r_col = self.r - other.r
        v_col = self.v - other.v
        
        v_sqrd = np.dot(v_col,v_col)
        r_sqrd = np.dot(r_col,r_col)
        
        if round(r_sqrd,3) == (self.radius + other.radius)**2:
            return 100000
        
        #print(r_col,v_col,r_sqrd,v_sqrd)
        
        
        #print((-np.dot(r_col,v_col)+np.sqrt(np.dot(r_col,v_col)**2 \
        #        - v_sqrd * (r_sqrd - Radius_col**2)))/v_sqrd)
        #print("===")
        
        
    
        positive_root = (-np.dot(r_col,v_col)+np.sqrt(np.dot(r_col,v_col)**2 \
                                                 - v_sqrd * (r_sqrd - Radius_col**2)))/v_sqrd
            
        negative_root = (-np.dot(r_col,v_col)-np.sqrt(np.dot(r_col,v_col)**2 \
                                                 - v_sqrd * (r_sqrd - Radius_col**2)))/v_sqrd
            
        lst = [positive_root,negative_root]
        
        if abs(positive_root) < 10**-7 or abs(negative_root) < 10**-7:
            return 10000
        
        if positive_root >= 0 and negative_root >= 0:  
            return min(n for n in lst  if n>0)
        else:
            return positive_root
        
    
    def collision(self,other):
        '''
        This collides one ball with another ball. It changes the velocity of the balls from beofre their collision to after their collision.

        Parameters
        ----------
        other : Ball
            This function collides the ball initialised with the other ball. 
            It takes into account thier mass, position and current velocites and calculates new velocities.
            
        Changes
        -------
        
        self and other: Ball
            The function changes the velocites of the called ball and the other ball. 

        Returns
        -------
        None.

        '''
        
        normal = other.r - self.r
        n = np.sqrt(normal[0]**2 + normal[1]**2)
        
        r1 = self.r
        r2 = other.r
        
        u1 = self.v
        u2 = other.v
        
        m1 = self.mass
        m2 = other.mass
            
        v1 = u1 - 2*m2 * np.dot(u1 - u2, r1 - r2) * (r1 - r2) / (n**2 * (m1 + m2)) 
        v2 = u2 - 2*m1 * np.dot(u2 - u1, r2 - r1) * (r2 - r1) / (n**2 * (m1 + m2)) 
        
        self.v = v1
        other.v = v2
        
    def get_patch(self,container = False,colour = 'r'):
        '''
        This returns the patch of the ball which allows the ball to be drawn.
        
        Parameters
        ----------
        container : Ball, optional
            If the container is called, different drawing parameters are required.
            The default is False.

        Returns
        -------
        patch : patch
            This reutrns a value for the patch of the balls. 
            This allows them to be drawn and animated in Matplotlib.

        '''
        if container == True:
            patch =  pl.Circle(self.pos(), self.radius, ec='b', fill=False, ls='solid')
        else:
            #colours = ['r','y','g','w','b']
            #colours[rd.randint(0, 4)]
            patch = pl.Circle(self.pos(), self.radius, fc=colour, ec='r') 
        return patch
              
    def KE(self):
        '''
        

        Returns
        -------
        Float
            The Kinetic energy of the ball.

        '''
        #velo = np.sqrt(other.v[0]**2 + other.v[1]**2)
        #Eo = 0.5*other.mass*(velo**2)
        #print("KE  =", Es , "for ball &", Eo , "for container.")
        #print("Total KE =", Eo+Es) 
        
        vels = np.sqrt(self.v[0]**2 + self.v[1]**2)
        return 0.5*self.mass*(vels**2)
    
    def momentum(self):
        '''
        Returns the momentum of the ball in a vector format.

        Returns
        -------
        Float
            The momentum of the ball.

        '''
        return self.mass * self.v
    
    def distance_from_centre(self):
        '''
        Returns the distance of the ball from the centre of the container in a 1D format.

        Returns
        -------
        Float
            The distance of the ball to the centre of the container.

        '''
        return np.sqrt(self.r[0]**2 + self.r[1]**2)
    
    def distance_from_ball(self,other):
        '''
        Calculates the distance between two balls.

        Parameters
        ----------
        other : Ball
            The other ball to make the distance comparison to.

        Returns
        -------
        Float
            The distance from the "self" ball to the other ball.

        '''
        r_col = self.r - other.r
        return np.sqrt(r_col[0]**2 + r_col[1]**2)
        
    
    
        