# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 09:26:30 2015

@author: thomas
"""
import pickle
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, FuncFormatter
import matplotlib.mlab as ml
import matplotlib.pyplot as plt


minCut = 8./1000.
maxCut = 2./10.

gridPoints = 512
L = 600e-6
xAr = np.arange(-gridPoints/2.,gridPoints/2.)/gridPoints*L

def ExponentialFit(Psi,x,nu1,nu2):
    logD = np.log(Psi[nu1:nu2])
    X = x[nu1:nu2]
    m = np.polyfit(X,logD,1)
    return m

#__path__  = '/home/thomas/Documents/Simulations/GPE/Results/2015-02-05_21-03-57/Averages'
#__path__ = '/home/thomas/Docume!nts/Simulations/GPE/Results/2015-02-09_10-20-53/Averages'
#__path__ = '/home/thomas/Documents/Simulations/GPE/Results/2015-02-26_11-41-02/Averages'
#__path__ = '/home/thomas/Documents/Simulations/GPE/Results/2015-02-26_16-06-07/Averages'
#__path__ = '/home/thomas/Documents/Simulations/GPE/Results/2015-02-27_09-37-45/Averages'
#__path__ = '/home/thomas/Documents/Simulations/GPE/Results/2015-02-27_16-14-10/Averages'
__path__ = '/home/thomas/Documents/Simulations/GPE/Results/2015-03-01_19-25-24/Averages'
#__path__ = '/home/thomas/Documents/Simulations/GPE/Results/2015-03-01_20-14-59/Averages'
#__path__ = '/home/thomas/Documents/Simulations/GPE/Results/2015-03-03_09-48-56/Averages'
#__path__ = '/home/thomas/Documents/Simulations/CNFD/Results'
#__path__ ='/home/thomas/Documents/Simulations/GPE/Results/2015-04-09_09-29-37/Averages'
with open(__path__ +'/PsiFinal.pkl','rb') as xfile:
    Psi = pickle.load(xfile)
    
with open ("/home/thomas/Documents/Simulations/GPE/Results/2015-02-26_08-37-58/Results_Run_1/Psi_Step_250000.pkl") as xfile:
    Psi = pickle.load(xfile)
APsi = np.abs(Psi)
logD  = np.log(APsi)
gridPoints = np.size(Psi)
numbsa = ((APsi > (APsi.max()*maxCut)).astype(int)).nonzero()
numbsb = ((APsi > (APsi.max()*minCut)).astype(int)).nonzero()
n1 = np.min(numbsa)
n2 = np.max(numbsa)
n3 = np.min(numbsb)
n4 = np.max(numbsb)

m1,b1 = ExponentialFit(APsi,xAr,n3,n1)
m2,b2 = ExponentialFit(APsi,xAr,n2,n4)

Exp1 = np.exp(m1*xAr+b1)
Exp2 = np.exp(m2*xAr+b2)
        
fig,ax1 = plt.subplots()
fig.suptitle('Final Density Profile - Log scale',fontsize = 20)
ax1.plot(xAr,APsi,'b',linewidth=1)
#ax1.plot(xAr[n3:n1],Exp1[n3:n1],'r',xAr[n2:n4],Exp2[n2:n4],'r',linewidth=3)
#ax1.plot((xAr[n1],xAr[n1]),(logD.min(),logD.max()),'g--')
#ax1.plot((xAr[n2],xAr[n2]),(logD.min(),logD.max()),'g--')
#ax1.plot((xAr[n3],xAr[n3]),(logD.min(),logD.max()),'m--')
#ax1.plot((xAr[n4],xAr[n4]),(logD.min(),logD.max()),'m--')
ax1.set_xlabel('x [$\mu$m]',fontsize=18)
#ax1.set_xlim(np.min(xAr),np.max(xAr))
ax1.set_xlim(xAr[0],xAr[gridPoints-1])
ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, pos:('%1.f')%(x*1e6)))
ax1.set_yscale('log')
ax1.set_ylabel(' Density $|\psi|^2$',fontsize=18)
ax1.set_ylim(APsi.min()*0.1,APsi.max()*100)

fig,ax2 = plt.subplots()
fig.suptitle('Final Density Profile',fontsize = 20)
ax2.plot(xAr,APsi,'b',linewidth=1)
#ax2.plot(xAr[n3:n1], Exp1[n3:n1],'r',linewidth=3)
#ax2.plot(xAr[n2:n4], Exp2[n2:n4],'r',linewidth=3)
#ax2.plot((xAr[n1],xAr[n1]),(0,APsi.max()),'g--')
#ax2.plot((xAr[n2],xAr[n2]),(0,APsi.max()),'g--')
#ax2.plot((xAr[n3],xAr[n3]),(0,APsi.max()),'m--')
#ax2.plot((xAr[n4],xAr[n4]),(0,APsi.max()),'m--')
ax2.set_ylim(0,APsi.max()*1.05)
ax2.xaxis.set_major_formatter(FuncFormatter(lambda x, pos:('%1.f')%(x*1e6)))
ax2.set_xlim(np.min(xAr)*0.5,np.max(xAr)*0.5)
ax2.set_xlabel('x [$\mu$m]',fontsize=18)
ax2.set_ylabel('Density $|\psi|^2$',fontsize=18)