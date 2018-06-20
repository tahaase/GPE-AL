"""
Created on Thu Sep 25 10:58:05 2014

@author: thomas
"""

import numpy as np
from math import *

class InitializeParameters():
    def __init__(self,Ox,Oy,Oz,d,l,gp,aN):
        self.Omega_x = Ox
        self.Omega_y = Oy
        self.Omega_z = Oz
        self.Dimensions = d
        self.Length = l
        self.GridPoints = gp
        self.Initialized = False
        self.AtomNumber = aN
        self.hBar = 1.054e-34
        
   
    def SetUpArray(self):
        L = self.Length
        gridPoints = self.GridPoints
#        X = np.linspace(-L/2,L/2,gridPoints)
        X = np.arange(-gridPoints/2.,gridPoints/2.)/gridPoints*L
        return X
    
    def SetUpPhaseSpaceArray(self):
        L = self.Length
        gridPoints = self.GridPoints
#        K = np.linspace(-gridPoints/2,gridPoints/2,gridPoints)*2*pi/L
        K = np.arange(-gridPoints/2.,gridPoints/2.)*2.*pi/L
        return K