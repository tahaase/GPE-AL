"""
This script is set-up to create the harmonic potential of the trap. It requires the simulation parameters to be initialized.
The function returns the value of the harmonic potential at a certain point in the array. 
@author: Thomas Haase
"""
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
#import matplotlib.pyplot as plt
#import numpy as np
import numpy as np
from math import *

hBar = 1.054e-34
m = 87*1.66e-27
a = 5.186e-9

def GetHarmonicPotential1D(x,Ox):
    return 0.5*m*(Ox**2*x**2)

def GetHarmonicPotential2D(x,y,Ox,Oy):   
    return 0.5*m*(Ox**2*x[:,np.newaxis]**2+Oy**2*y[np.newaxis,:]**2)
    
def GetHarmonicPotential3D(x,y,z,Ox,Oy,Oz):
    return 0.5*m*(Ox**2*x[:,np.newaxis,np.newaxis]**2+Oy**2*y[np.newaxis,:,np.newaxis]**2+Oz**2*z[np.newaxis,np.newaxis,:]**2)

def GetKineticEnergyGrid1D(x):
    return 0.5*(hBar**2/m)*x**2
   
def GetKineticEnergyGrid2D(x,y):    
    return 0.5*(hBar**2/m)*(x[:,np.newaxis]**2+y[np.newaxis,:]**2)

def GetKineticEnergyGrid3D(x,y,z):
    return 0.5*(hBar**2/m)*(x[:,np.newaxis,np.newaxis]**2+y[np.newaxis,:,np.newaxis]**2+z[np.newaxis,np.newaxis,:]**2)
    
def GetInteractionPotential1D(omega_x):
    intPot = 4*pi*a*(hBar**2)/m
    intPot = intPot/((hBar/m*omega_x)**(1./2.))
    return intPot
    
def GetInteractionPotential2D(omega_x):
    intPot = 4*pi*a*(hBar**2)/(m)
    intPot = intPot/((hBar/(m*omega_x))**(1./2.))
    return intPot
    
def GetInteractionPotential3D():
    intPot = 4*pi*a*(hBar**2)/m
    return intPot
    
def GetDoubleWellPot2D(x,y,Ox,Oy,L):
    VHarmy = 0.5*m*Oy**2*y**2    
    VDubx = 1e-6*m*(0.045*(x*12/L)**4-1.5*(x*12/L)**2)
    VDubx = VDubx - min(VDubx)
    
    DoublW = (VDubx[np.newaxis,:]+VHarmy[:,np.newaxis])*100.
    return DoublW


#V=[]
#for i in range(0,gridPoints):
#    V.append([])
#    for j in range(0,gridPoints):
#        V[i].append(GetHarmonicPotential(i,j,64))
#
#fig = plt.figure()
#fig.suptitle('Harmonic Potential')
#ax = fig.gca(projection='3d')
#X,Y = np.meshgrid(x,y)
#Va = np.asarray(V)
#surf = ax.plot_surface(X,Y,Va,rstride=1, cstride=1, cmap=cm.coolwarm,
#        linewidth=0, antialiased=False)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#ax.set_zlabel('Potential Energy')
#ax.xaxis.labelpad=20
#ax.zaxis.labelpad=100
#fig.colorbar(surf, shrink=0.5, aspect=5)
#plt.show()