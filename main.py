#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
author : Romain Pinguet    
    
This code aims to solve the dynamic structural problem of two orthogonal beams
of same length where the extremity of the first beam (A) is fully constrained 
and a vertical force F is applied to the extremity of the second beam (C).

     _________C
   B|         | Force F
    |         v
    |
    |
    | A (fully constrained)
 -------     
 //////    
 
This program calls :
    Newmark.py
    matrix_assembly.py
    
The input parameters section below can be modified by the user 

"""

import matplotlib.pyplot as plt
import numpy as np
from math import pi
from matrix_assembly import *
from Newmark import *




def getdisplacement(F, deltat=0.01):
    
    """
    Get the displacement of point C and the simulation time
    
    Inputs :    F = Vertical force applied in C (N)
                deltat = tstep of the Newmark method (s)  (default 0.01 s)
                
    Output :    time = simulation time (s)
                wC = vertical displacement (m)
                    
    """
    
    
    ###################          Input parameters       #######################
    
    ###########################        Mechanical inputs
    
    # Young modulus (Pa)
    E=200e9
    #beam Length (m)
    L=6
    # mass density (kg/m3)
    rho=7850
    #cross section area (m2)
    A=2827.43*1e-6
    #cross section moment (m4)
    I=2.89812e6*1e-12
    #coeff lambda for the damping matrix C=lambda*K (s)
    lamb=0.03
    
    
    ############################       Numerical Inputs for Newmark method
    
    # duration of the simulation (s)
    duration=10   
    
    # gamma and beta for Newmark method 
    gamma=0.5
    beta=0.25 # these parameters are equivalent to the trapezoidal rule
    
    
    ##############################################################################
    
    # step 1 : Construct the matrices of mass and stiffness
    
    #Element 1
    
    #Mass Matrix
    M1=MassMat(rho,A,L)
    #Stiffness Matrix
    K1=StiffnessMat(E,A,L,I)
    
    
    #Element 2
    
    #Mass Matrix
    M2=MassMat(rho,A,L)
    #Stiffness Matrix
    K2=StiffnessMat(E,A,L,I)
    
    
    # step 2 : Transform the matrices from local to global coordinate system
    
    #Element 1
    theta1=pi/2
    #Mass Matrix
    M1g=TransformMat(M1,theta1)
    #Stiffness Matrix
    K1g=TransformMat(K1,theta1)
    
    
    #Element 2
    theta2=0
    #Mass Matrix
    M2g=TransformMat(M2,theta2)
    #Stiffness Matrix
    K2g=TransformMat(K2,theta2)
    
    
    # step 3 : Assemble the global matrices
    
    #Global Mass matrix
    Ma=assemblyMat(M1g,M2g)
    
    #Global Stiffness Matrix
    Ka=assemblyMat(K1g,K2g)
    
    
    # step 4 : Reduce the order of the matrices
    
    #list of constrained DoFs. DoFs are defined as A [0,1,2] B [3,4,5] and C [6,7,8]
    constrained_DoF=[0,1,2]   
    
    #Reduced Global Mass matrix
    M=reduceMorder(Ma,constrained_DoF)
    
    #Reduced Global Stiffness matrix
    K=reduceMorder(Ka,constrained_DoF)
    
    #Reduced Global Damping matrix
    C=lamb*K
    
    
    # step 5 : Solve the dynamic problem with Newmark method 
    
    
    #Equation solving
    time,wC=Newmark(M,K,C,F,deltat,duration,beta,gamma)

    return time,wC


def plotdisplacement(time,wC,output_path,deltat=0.01):
    
    """
   plot the displacement of point C as a function of time
    
    Inputs :    time = simulation time (s)
                wC = vertical displacement (m)
                output_path = path of the output figure 
                deltat = tstep of the Newmark method (s)  (default 0.01 s)
                    
    """

    ######     Plot the time series of the vertical position of node C 
    plt.figure(1)
    plt.plot(time,wC*1e3,label="time step = %s s" %deltat)
    plt.xlabel("time [s]")
    plt.ylabel("displacement [mm]")
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.grid()
    plt.title("Vertical displacement of node C")
    plt.legend(loc=0)
    plt.savefig(output_path,bbox_inches='tight')

        
#
######     Time step convergence analysis  (plot)

def plotTstepconv(F,output_path):
    
    """
   plot time series of displacement of node C for different time steps
    
    Inputs :    F = input force (N)
                output_path = path of the output figure 
                
    (This function could be optimized by avoiding iterations on step 1 to 4 from
    the function getdisplacement)
    """

    
    tsteps=[0.2,0.1,0.05,0.01,0.005,0.001]
    plt.figure(2)
    for deltat in tsteps: 
        
        #Equation solving
        time,wC=getdisplacement(F,deltat)
       
        ######     Plot the time series of the vertical position of node C 
   
        plt.plot(time,wC*1e3,label="time step = %s s" %deltat)
    plt.xlabel("time [s]")
    plt.xlim(0,2)
    plt.ylabel("displacement [mm]")
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.title("Time step convergence analysis")
    plt.legend(loc=0)
    plt.grid()
    plt.savefig(output_path,bbox_inches='tight')


if __name__ == "__main__":
    F=-1
    time,wC=getdisplacement(F)
    plotdisplacement(time,wC,"Vertical_displacement_C.png")
    plotTstepconv(F,"convergence_tstep.png")