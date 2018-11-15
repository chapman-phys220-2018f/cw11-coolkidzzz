#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###
# Name:Gabriella Nutt and Raha
#Student ID: 2307512
#Email: nutt@chapman.edu
#Course: PHYS220/MATH220/CPSC220 Fall 2018
#Assignment: HW11
###


import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


J = np.array([[0, 1],[ -1, 0]])

def drdt(r):
    """Returns the derivative of r by multiplying it by the matrix of J"""

    return J@ r



def euler(r0,N):
    """The function euler should return an array of points of the  length ...."""
    dt = 2*np.pi/N

    t = np.arange(0,5*(2*np.pi),dt)

    e = np.zeros((2,len(t))) #makes it a 2 by n matrix
    e[:,0]=r0
    for i in range(1,len(t)):
        e[:,i]= e[:,i-1]+ dt*drdt(e[:,i-1]) #redefining r in drdt
    return e

def eplot(e,t):
    """"this plots the euler function!"""
    sin = -1*(np.sin(t))
    cos = np.cos(t)

    plt.plot(t,e[0,:],color="k",label='e vs x') # x is the 0 row of e
    plt.plot(t,e[1,:],color="green",label='e vs y') # y is the 1st row of e
    plt.plot(t,sin,'r--',color="orange",label='-sinx')
    plt.plot(t,cos,'r--',color="yellow",label='cosx')
    plt.legend(loc='upper left') #puts legend in upper left corner
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Plot of euler')
    plt.ylim((-5,5))
    plt.show()

def heun(r0,N):
    """heun ....."""
    dt = 2*np.pi/N
    t = np.arange(0,5*(2*np.pi),dt)
    h = np.zeros((2,len(t)))
    h[:,0]=r0
    for i in range(1,len(t)):
        u = h[:,i-1]+ dt*drdt(h[:,i-1])
        h[:,i]= h[:,i-1]+ (dt/2)*drdt(h[:,i-1])+(dt/2)*drdt(u)
    return h

def hplot(h,t):
    """"this plots the heun function!"""
    sin = -1*(np.sin(t))
    cos = np.cos(t)

    plt.plot(t,h[0,:],color="k",label='h vs x') # x is the 0 row of h
    plt.plot(t,h[1,:],color="pink",label='h vs y') # y is the 1st row of h
    plt.plot(t,sin,'r--',color="orange",label='-sinx')
    plt.plot(t,cos,'r--',color="yellow",label='cosx')
    plt.legend(loc='upper left') #puts legend in upper left corner
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Plot of heun')
    plt.ylim((-5,5))
    plt.show()


def r2(r0,N):
    """function for 2nd-order Runge-Kutta Method"""
    dt = 2*np.pi/N
    t = np.arange(0,5*(2*np.pi),dt)
    r2 = np.zeros((2,len(t)))
    r2[:,0]=r0
    for i in range(1,len(t)):
        K1= dt*drdt(r2[:,i-1])
        K2 = dt*drdt(r2[:,i-1]+K1/2)
        r2[:,i]= r2[:,i-1]+K2
    return r2

def r2plot(r2,t):
    """"this plots the 2nd-order Runge-Kutta Method!"""
    sin = -1*(np.sin(t))
    cos = np.cos(t)

    plt.plot(t,r2[0,:],color="k",label='r2 vs x') # x is the 0 row of r2
    plt.plot(t,r2[1,:],color="purple",label='r2 vs y') # y is the 1st row of r2
    plt.plot(t,sin,'r--',color="orange",label='-sinx')
    plt.plot(t,cos,'r--',color="yellow",label='cosx')
    plt.legend(loc='upper left') #puts legend in upper left corner
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Plot of 2nd-Runge Kutta')
    plt.ylim((-5,5))
    plt.show()

def r4(r0,N):
    """function for 4th-order Runge-Kutta Method"""
    dt = 2*np.pi/N
    t = np.arange(0,5*(2*np.pi),dt)
    r4 = np.zeros((2,len(t)))
    r4[:,0]=r0
    for i in range(1,len(t)):
        K1= dt*drdt(r4[:,i-1])
        K2 = dt*drdt(r4[:,i-1]+K1/2)
        K3 = dt*drdt(r4[:,i-1]+K2/2)
        K4 = dt*drdt(r4[:,i-1]+K3)
        r4[:,i]= r4[:,i-1]+(K1+2*K2+2*K3+K4)/6
    return r4

def r4plot(r4,t):
    """"this plots the 4th-order Runge-Kutta Method!"""
    sin = -1*(np.sin(t))
    cos = np.cos(t)

    plt.plot(t,r4[0,:],color="k",label='r4 vs x') # x is the 0 row of r4
    plt.plot(t,r4[1,:],color="blue",label='r4 vs y') # y is the 1st row of r4
    plt.plot(t,sin,'r--',color="orange",label='-sinx')
    plt.plot(t,cos,'r--',color="yellow",label='cosx')
    plt.legend(loc='upper left') #puts legend in upper left corner
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Plot of 4th-Runge Kutta')
    plt.ylim((-5,5))
    plt.show()