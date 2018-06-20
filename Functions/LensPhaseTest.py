# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 16:31:34 2015

@author: thomas
"""
import numpy as np
from math import pi
grid = 256
rad = grid/2
L = 150e-6

fac = 4
sig = 5e-6
sep = sig*fac

lens = np.zeros(grid)
lens[grid/2.-rad/2.:grid/2.+rad/2.] = 1

angles = np.random.normal(0,1,grid)*pi*0.8
phase = np.exp(1j*angles)

dLens = lens*phase

x = np.arange(-grid/2.,grid/2.)/grid*L
Gaus = np.exp(-(x-sep/2)**(2.)/(2*sig)**(2.))+np.exp(-(x+sep/2)**(2.)/(2*sig)**(2.))

Id = dLens*Gaus
IdAbs = np.abs(Id)
If = np.abs(np.fft.fftshift(np.fft.fft(Id)))
#If = (If - np.average(If))/(np.max(If)-np.min(If))
If = (If/np.average(If))-1.
