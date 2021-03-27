#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
author: Romain Pinguet
"""

import matplotlib.pyplot as plt
import numpy as np
from math import pi
from numpy.linalg import inv


def StiffnessMat(E,A,L,I):
    
    """
    Generate the stiffness matrix 
    
    Inputs :    E = Young Modulus (Pa)
                A = Cross section area (m^2)
                L = Length of the element (m)
                I = Cross section second moment of area (m^4)
                
    Output :    K = Stiffness Matrix of Euler-Bernoulli Straight Beam Element
    
    """
    #init matrix
    K=np.zeros((6,6))
    
    #value of the matrix elements
    k11=E*A/L
    k22=12*E*I/(L**3)
    k23=6*E*I/(L**2)
    k33=4*E*I/L
    k36=2*E*I/L
    
    #Fill the stiffness matrix
    K=[
       [k11,0,0,-k11,0,0],
       [0,k22,k23,0,-k22,k23],
       [0,k23,k33,0,-k23,k36],
       [-k11,0,0,k11,0,0],
       [0,-k22,-k23,0,k22,-k23],
       [0,k23,k36,0,-k23,k33]
       ]
    #Convert to array
    K=np.array(K)
    
    return K

def MassMat(rho,A,L):
    
    """
    Generate the mass matrix 
    
    Inputs :    rho = Mass density (kg/m^3)
                A = Cross section area (m^2)
                L = Length of the element (m)
                
    Output :    M = Mass Matrix of Euler-Bernoulli Straight Beam Element
                    associated with translational inertia
                    
    """
    
  
    #initialize the matrix
    M=np.zeros((6,6))
    #fill the matrix
    M=[
       [140,0,0,70,0,0],
       [0,156,22*L,0,54,-13*L],
       [0,22*L,4*L**2,0,13*L,-3*L**2],
       [70,0,0,140,0,0],
       [0,54,13*L,0,156,-22*L],
       [0,-13*L,-3*L**2,0,-22*L,4*L**2],
      ]
    #convert to array
    M=(rho*A*L*1/420)*np.array(M)
    
    return M


def TransformMat(Matrix,theta):

    """  
    
    Transform the Matrix M in the global coordinate system
    
    Inputs :    Matrix = 6x6 matrix of Euler-Bernoulli Straight Beam Element
                theta = angle between the element and the horizontal axis
                        of the global coordinate system (rad)
                
                
    Output :    Matrix = Matrix in the global coordinate system
    
    """    
    # Initialize the rot matrix
    rot=np.zeros((6,6))
    
    # Fill the rot matrix
    rot=[
         [np.cos(theta),np.sin(theta),0,0,0,0],
         [-np.sin(theta),np.cos(theta),0,0,0,0],
         [0,0,1,0,0,0],
         [0,0,0,np.cos(theta),np.sin(theta),0],
         [0,0,0,-np.sin(theta),np.cos(theta),0],
         [0,0,0,0,0,1]
         ]
    
    # Convert rot matrix to array
    rot=np.array(rot)
    rot_t=np.transpose(rot)
    
    # Apply the rotation to the input matrix
    Matrix=np.dot(rot_t,np.dot(Matrix,rot))
    
    return Matrix


def assemblyMat(Melem1,Melem2):
    
    """  
    
    Assemble the matrices of two elements (the second node of element 1 is 
    the first node of element 2)
    
    Inputs :    Melem1 = Matrix of the first Element 
                Melem2 = Matrix of the second Element
                
                
    Output :    Matrix = Global matrix of the system 
    
    """        

# Stiffness matrix 
    c0=np.zeros((6,1))
    l0=np.zeros((1,9))
       
    Melem1=np.hstack((Melem1,c0,c0,c0))
    Melem1=np.vstack((Melem1,l0,l0,l0))

    Melem2=np.hstack((c0,c0,c0,Melem2))
    Melem2=np.vstack((l0,l0,l0,Melem2))
    
    Matrix=Melem1+Melem2

    return Matrix
    

def reduceMorder(M,constrained_DoF):
    
    """
    
    Reduce the order of the matrix by removing the line and column of the 
    constrained DoF
    
    Inputs :    M = Global matrix of the system
                constrained_DoF = constrained DoF position in the global
                                  displacement vector 
                
    Output :    M = Reduced matrix 

    """         
    
    M = np.delete(M, constrained_DoF, 0)
    M = np.delete(M, constrained_DoF, 1)
    
    return M
    
