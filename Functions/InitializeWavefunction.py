# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 10:47:11 2014

@author: thomas
"""
import numpy as np
from math import *

def InitialWavefunction1D(x,sigma):
    Psi = ((abs(x)<5e-6).astype(int))
    Psi = Psi*np.exp(-(x**2)/(2*sigma**2))
    return Psi
    
def SimpleFunction1D(x):
    Psi = ((abs(x)<5e-6).astype(int))
    return Psi
    
def InitialWavefunction2D(x,y,sigma):
    xyGrid = (x[:,np.newaxis]**2+y[np.newaxis,:]**2)**(1./2.)
    Psi = ((xyGrid<5e-6).astype(int))
    Psi = Psi*np.exp(-(xyGrid**2)/(2*sigma**2))
    return Psi
    
def SimpleFunction2D(x,y):
    xyGrid = (x[:,np.newaxis]**2+y[np.newaxis,:]**2)**(1./2.)
    Psi = ((xyGrid<5e-6).astype(int))
    return Psi
    
def InitialWavefunction3D(x,y,z,sigma):
    xyzGrid = (x[:,np.newaxis,np.newaxis]**2+y[np.newaxis,:,np.newaxis]**2+z[np.newaxis,np.newaxis,:]**2)**(1./2.)
    Psi = ((xyzGrid<5e-6).astype(int))
    Psi = Psi*np.exp(-(xyzGrid**2)/(2*sigma**2))
    return Psi
    
def SimpleFunction3D(x,y,z):
    xyzGrid = (x[:,np.newaxis,np.newaxis]**2+y[np.newaxis,:,np.newaxis]**2+z[np.newaxis,np.newaxis,:]**2)**(1./2.)
    Psi = ((xyzGrid<5e-6).astype(int))
    return Psi   
    
def NarrowInitialMomentumWavefunction(k,m):
    dk = k[2]-k[1]
    phi = np.exp(-k**2/(2*(m*dk)**2))
    psi = np.fft.fftshift(np.fft.fft(phi))
    return psi