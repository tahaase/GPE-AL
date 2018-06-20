# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 11:18:01 2014

@author: thomas
"""
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.mlab as ml
import matplotlib.pyplot as plt
import pickle
import os

plt.close('all')

'''Number of peaks and peak width control'''
pNumb = 2000
var = 9e-8

'''Dissorder control'''
sepF = 400
noiseStr = 0.1
'''Simulation Parameters'''
gP = 4096
L = 800e-6

'''Generation'''

'''1. Definitions'''
def gaussian(x,a,c):
    f = np.exp(-(x-a)**2/(2*c**2))
    return f
    
'''Code'''

x = np.arange(-gP/2.,gP/2.)/gP*L
xi = []
Vbase = np.zeros(gP)
Vrand = np.zeros(gP)

if (pNumb%2 == 0):
    fac = 1./(2.*pNumb)
    for i in range(1,pNumb+1):
        a = fac*(2*i-1)*gP
        xi.append(x[int(a)])
else:
    fac = 1./(pNumb+1.)
    for i in range(1,pNumb+1):
        a = fac*i*gP
        xi.append(x[int(a)])  

if(sepF >0):
    xir = xi + (L/sepF)*(np.random.rand(len(xi))-1./2.)
else: 
    xir = xi

for i in range(len(xi)):
    Vbase = Vbase+(gaussian(x,xi[i],var))
    Vrand = Vrand+(gaussian(x,xir[i],var))

Vbase = Vbase - Vbase.min()
Vbase = Vbase/Vbase.max()
Vrand = Vrand + (1+noiseStr*np.random.rand(gP))
Vrand = Vrand - Vrand.min()
Vrand = Vrand/Vrand.max()

fig,axes = plt.subplots()
axes.plot(x, Vbase,'r',x,Vrand)
axes.set_xlabel('x')
axes.set_xlim(x.min(),x.max())
axes.set_ylabel('Random Potential')
axes.set_ylim(Vbase.min(),Vbase.max())
axes.set_title('RandomPotential')
plt.figtext(0.05,0.9,'Peak Number = '+str(pNumb)+'\n'+'Peak variance = '+str(var)+'\n'+'Peak Separation factor = '+str(sepF)+'\n'+'Noise strength = '+str(noiseStr))
plt.savefig('RandomPotential.png',dpi = 250)
plt.show(block = False)

Vrandf = np.real(np.fft.fftshift(np.fft.fft(Vrand)))
Vbasef = np.real(np.fft.fftshift(np.fft.fft(Vbase)))
fig,axes = plt.subplots()
axes.plot(x,Vbasef,'r',x,Vrandf)


random_Potential_1D = open('Random_Potential_1D_creator.pkl','wb')
pickle.dump(Vrand,random_Potential_1D)
random_Potential_1D.close()
#random_Potential_1D = open('1D/Random_Potential_1D_creator.pkl','wb')
#pickle.dump(Vrand,random_Potential_1D)
#random_Potential_1D.close()
print 'Potential saved to file.'
        
cor = np.correlate(Vbase,Vrand,"same")

fig, axes = plt.subplots()
axes.plot(x, cor)
axes.set_xlabel('x')
axes.set_xlim(x.min(),x.max())
axes.set_ylabel('Cor')
axes.set_ylim(cor.min(),cor.max())
axes.set_title('Correlation')   
plt.show(block = False)