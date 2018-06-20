# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 09:55:34 2015

@author: thomas
"""
import numpy as np
import pickle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, FuncFormatter
import matplotlib.mlab as ml
import matplotlib.pyplot as plt

from Functions import Definitions as Def 
from Functions import TimeEvolveFunction1D as TE
from Functions import Constants
c = Constants.Constants()
from Functions import FileReader as FR

""" *********************************************************************** """
gridPoints = 4096
Pot = 10
""" *********************************************************************** """
SFreq = 200
CutFreq = 10
Alpha = 0.7
""" *********************************************************************** """
deltaCor = 0.02
sigmaCor = 1.5
""" *********************************************************************** """
VPink = Def.PinkPotential(SFreq,CutFreq,Alpha,gridPoints)
VPinkAv = np.average(VPink)
VPinkVar = np.average(VPink**2)


VCor = Def.CorelatedRandomSave(deltaCor,sigmaCor,gridPoints,None)
VCorAv = np.average(VCor)
VCorVar = np.average(VCor**2)

facPink = (Pot/VPinkVar)**(1./2.)
facCor = (Pot/VCorVar)**(1./2.)

VPinkMod = VPink*facPink
VPinkModAv = np.average(VPinkMod)
VPinkModVar = np.average(VPinkMod**2)

VCorMod = VCor*facCor
VCorModAv = np.average(VCorMod)
VCorModVar = np.average(VCorMod**2)
