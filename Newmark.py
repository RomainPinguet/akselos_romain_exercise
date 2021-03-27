#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
author : Romain Pinguet
"""


import numpy as np
from math import pi
from numpy.linalg import inv

def Newmark(M,K,C,F,deltat,duration,beta,gamma):
    
    """
    This function aims to solve dynamically the system of equations :

            Mu''(k+1) + Cu'(k+1) + Ku(k+1) = F
            
            u(t0)=u0
            u'(t0)=u'0
            
    with a Newmark algorithm for a 2 beam element problem, where a constant 
    vertical force is applied at the extremity of the second element
    
    Inputs :    M = Mass matrix 
                K = Stiffness
                C = Damping matrix
                F = Vertical force at the extremity (N)
                deltat = time step (s)
                duration = duration of the simulation (s)
                beta = parameter of the time integration algo (displacement)
                gamma = parameter of the time integration algo (velocity)
                
    Output :    time = simulation time 
                wc = vertical displacement at at the extremity
    
    """
        
    # initialize the applied force
    
    f0=np.array([[0,0,0,0,F,0]]).transpose()
    
    # initialize displacement, velocity and acceleration
    
    u0=np.zeros((6,1))    
    udot0=np.zeros((6,1))    
    udotdot0=inv(M).dot(f0)

    
    # initialize loop 
    
    u=u0
    udot=udot0
    udotdot=udotdot0
    t=0
    
    #output variables 
    
    # vertical position of node C
    wc=[u[4]]  
    # simulation time
    time=[t]
        
    Niter=int(duration/deltat)
    
    for i in range(Niter):
        
    # step 1  Calculation of the predictor
        u_tild=u+udot*deltat+udotdot*(0.5-beta)*deltat**2
        udot_tild=udot+udotdot*(1-gamma)*deltat
    
    # step 2 Solution of the linear problem 
        udotdot=inv(M+C*gamma*deltat+K*beta*deltat**2).dot(f0-C.dot(udot_tild)-K.dot(u_tild))
    
    #step 3 Calculation of the corrector
        udot=udot_tild+udotdot*gamma*deltat 
        u=u_tild+udotdot*beta*deltat**2
        
    #store time and wc
        t=t+deltat    
        time.append(t)
        wc.append(u[4])
                
    wc=np.array(wc)
    return time,wc

