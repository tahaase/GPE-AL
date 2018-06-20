# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 13:05:37 2015

@author: thomas
"""
from math import pi
import numpy as np

Vo = 100
sigma = 1e-6
gridPoints = 512
L = 600e-6
I = np.random.exponential(Vo,gridPoints)
k = np.arange(-gridPoints/2.,gridPoints/2.)*2.*pi/L
Fil = np.zeros(gridPoints) 
Fil[(-1/sigma<k) & (k<1/sigma)] = 1

def autocorrelate(x):
    result = np.correlate(x,x,mode='full')
    return result[result.size/2:]

If = np.fft.fftshift(np.fft.fft(I))
If = If*Fil
I = np.fft.ifft(np.fft.ifftshift(If))

ACor = autocorrelate(I)
I = np.real(I)
I = I - np.average(I)
Avg = np.average(I)
AvgSq = np.average(I**2)