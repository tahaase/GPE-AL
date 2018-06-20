# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 14:23:03 2014

@author: thomas
"""
from math import *
import numpy as np
import os
import pickle    
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, FuncFormatter
import matplotlib.mlab as ml
import matplotlib.pyplot as plt
    

def MakeCoordinateArray(gridPoints,L):
    return np.arange(-gridPoints/2.,gridPoints/2.)/gridPoints*L
    
def MakeCoordinateArray2(gridPoints,L):
    return np.linspace(-L/2,L/2.,gridPoints)

def MakeMomentumArray(gridPoints,L):
    return np.arange(-gridPoints/2.,gridPoints/2.)*2.*pi/L
    
def MakeMomentumArray2(gridPoints,L):
    return np.linspace(-gridPoints/2.,gridPoints/2,gridPoints)*2.*pi/L

def HarmonicPotential(m,Ox,x):
    return 0.5*m*(Ox**2*x**2)
    
def KineticEnergy(hBar,m,x):
    return 0.5*(hBar**2/m)*x**2
    
def InteractionPotential(a,hBar,m,omega_x):
    intPot = 4*pi*a*(hBar**2)/m
    intPot = intPot/((hBar/m*omega_x)**(1./2.))
    return intPot
    
def RandomPotential(pNumb,var,sepF,noiseStr,gP,x):

    return RP
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
    
def PinkPotential(sfreq,cut,alpha,gp):
    wn=np.random.normal(0.,np.sqrt(sfreq),gp)
    s = np.fft.rfft(wn)
    f = np.fft.fftfreq(gp, d=1./sfreq)[:len(s)]
    f[-1] = np.abs(f[-1])
    sfft = s*one_over_f(f, cut, alpha)
    pn = np.fft.irfft(sfft)
#    pn = (pn-pn.min())/pn.max()
    pn = pn/pn.max()
    return pn    

def CorelatedRandom(delta,sigma,gp,subpath):
    if (delta != 0):   
        a = int((1/delta)*4)
    else:
        a = 1
    cr = np.zeros(gp+a)
    cr[0] = np.random.rand()
    for i in range (1,len(cr)):
        cr[i] = (1-delta)*cr[i-1]+delta*np.random.normal(0,sigma)
#    cr = cr - cr.min()
    return (cr[a-1:gp+a-1]-np.average(cr[a-1:gp+a-1]))
    
def DiffusivePlateIntensity(sepFac,grid,x):
    sig = 100e-6
    sep = sig*sepFac
    rad = grid/2
    lens = np.zeros(grid)
    lens[grid/2.-rad/2.:grid/2.+rad/2.] = 1
    Gauss = np.exp(-(x-sep/2)**(2.)/(2*sig)**(2.))+np.exp(-(x+sep/2)**(2.)/(2*sig)**(2.))
    angles = np.random.normal(0,1,grid)*pi*0.8
    phase = np.exp(1j*angles)
    dPlate = Gauss*phase
    Id = lens*dPlate
    If = np.abs(np.fft.fftshift(np.fft.fft(Id)))
    If = (If/np.average(If))-1.    
    return If
    
def SpeckleIntensity(kvec,corLen,gP):
    rand = np.random.exponential(1,gP)
    fil = np.zeros(gP)
    fil[abs(kvec)< 0.5/corLen] = 1    
    randf = np.fft.fftshift(np.fft.fft(rand))
    randf = randf*fil
    rand = np.fft.ifft(np.fft.ifftshift(randf))

    return np.abs(rand)
    
def one_over_f(f, cut, alpha):
    Frq = np.ones_like(f)
    Frq[f<cut] = np.abs((f[f<cut]/cut)**(-alpha))
    Frq[0] = 1
    return Frq
    
def PickleSave(path,filename,array):
    if not os.path.exists(path):os.makedirs(path)
    fileName = path+'/'+filename+'.pkl'
    state = open(fileName,'wb')
    pickle.dump(array,state)
    state.close()
    
def PickleOpen(filepath):
    state = open(filepath,'rb')
    array = pickle.load(state)
    state.close()
    return array
#def RescalePotential(Vin,nKvalue,kb):
#    B = kb*(nKvalue*1e-9)
#    A = np.average(Vin)
#    return Vin*(B/A)
    
def RescalePotential(Vin,nK,kb):
    B = kb*(nK*1e-9)
    A = np.average(Vin**2)
    Vout = Vin*(B/(A**(0.5)))
    return Vout
    
def AddMomentum(psi,beta,hb,x,pR):
    return psi*np.exp(-1j*beta*pR*x/hb)
    
def CreateNarrowMomentumGroundState(n,k,m,aN,dx,kRecoil):
    phi = (1/((2*pi)**(0.5)*(n*kRecoil)))*np.exp(-(k-m*kRecoil)**2/(2*(n*kRecoil)**2))
    phi = np.fft.ifftshift(phi)    
    psi = np.fft.ifft((phi))
    psi = np.fft.fftshift(psi)
    psi = psi*(1/(sum(abs(psi)**2*dx)))**(1./2.) 
    return psi

#    
def SaveFinal(psi,subpath,totalTime):
    final_state_function = open(subpath+'/Psi_Step_'+str(totalTime)+'.pkl','wb')
    pickle.dump(psi,final_state_function)
    final_state_function.close()
    
def GetMaxDistance(psi,cut,xAr):
    numbs = ((abs(psi**2)>(abs(psi**2).max())*cut).astype(int)).nonzero()
    return xAr[np.max(numbs)]

def Gaussian(x,a,c):
    return np.exp(-(x-a)**2/(2*c**2))

def ReadGroundState(path,subpath,xAr):
    with open(path+'/GroundState.pkl','rb') as ground:
        Psi = pickle.load(ground)
    PsiA = abs(Psi**2)
    fig, ax1 = plt.subplots()
    ax1.plot(xAr,PsiA)
    ax1.text(xAr[10],PsiA.max()*(1.1),'Elapsed time = 0s')
    ax1.set_ylabel('$|\psi|^2$')
    ax1.set_ylim(0,PsiA.max()*1.2)
    ax1.set_xlabel('x [$\mu$m]')
    ax1.set_xlim(xAr.min(),xAr.max())
    ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, pos:('%.1f')%(x*1e6)))
    plt.subplots_adjust(top=0.85)
    plt.tight_layout()
    plt.savefig(subpath+'/Psi_Step_'+str(0)+'_.png',dpi = 250)
    plt.close('all')
    return Psi

def Autocorr(x):
    result = np.correlate(x,x,mode='full')
    return result[result.size/2:]
    
def NotifyEnd():
    import smtplib
    server = smtplib.SMTP('mailhost.auckland.ac.nz', port = 25)
    server.ehlo()
    server.starttls()
    server.login("thaa191","tra50546")
    msg = "\r\n".join([
      "From: thaa191@aucklanduni.ac.nz",
      "To: thaa191@aucklanduni.ac.nz",
      "Subject: Simulation",
      "",
      "Simulation is complete."
      ])
    server.sendmail("thaa191@aucklanduni.ac.nz","thaa191@aucklanduni.ac.nz",msg)
    server.quit()
