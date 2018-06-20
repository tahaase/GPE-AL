# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 15:42:46 2014

@author: thomas
"""
import numpy as np

def GetCorrelation(x,x0,beta,charL):
    return 0.5*(1+((x-x0)**2)/(charL**2))**(-beta/2)*1e-20
    
def GetRandomPotential(corG,dx,Vo,gp):
        a = np.random.rand(gp)/dx
        Vk = 1e-30*np.multiply(a,corG**(1./2.))
        Vtemp = np.real(np.fft.ifft(np.fft.ifftshift(Vk)))
        V = np.fft.fftshift( Vtemp - min(Vtemp) )
        return V

def GetSuperRandomPotential(gp,Vo,mu,sigma):
    a = np.random.normal(mu,sigma,gp)+1j*np.random.normal(mu,sigma,gp)
    b = abs(a)**2/((abs(np.mean(a))**2)-1)
    b = b/b.max()
    Vr = Vo*b
    Vr = Vr+abs(Vr.min())
    return Vr