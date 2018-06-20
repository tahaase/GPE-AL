# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 16:50:51 2014

@author: thomas
"""
import numpy as np
import pickle

def RealTimeEvolve(psi,timeSteps,freetime,dt,intPot,KEn,aN,dx,dW,sampStep,ramp,path):
    hBar = 1.054e-34 
    factor = np.linspace(0,1,ramp+1)
    psiavg = 0
    phi = np.fft.ifftshift(np.fft.ifft(psi))
    phi = phi*np.exp(-1j*(KEn/hBar)*dt*freetime)
    psi = np.fft.fft(np.fft.fftshift(phi))
    psi = psi*(1/(sum(abs(psi)**2*dx)))**(1./2.)
    for i in xrange(1,timeSteps+1):    
        #Real time split step evolution
        if (i < ramp):
            fac = factor[i+1]
        else:
            fac = 1   
            
        phi = np.fft.ifft(psi)                            #Momentum space wave function
        phi = np.fft.ifftshift(phi)  
        phi = phi*np.exp(-1j*(KEn/hBar)*dt*0.5)                #Kinetic energy evolution
        phi = np.fft.fftshift(phi)
        psi = (np.fft.fft(phi))                             #Coordinate space wave function
        psi = psi*np.exp(-1j*(((aN*intPot)*(abs(psi)**2)+fac*dW)/hBar)*dt) #Interaction potential evolution        
        phi = np.fft.ifft(psi)                            #Momentum space wave function
        phi = np.fft.ifftshift(phi)  
        phi = phi*np.exp(-1j*(KEn/hBar)*dt*0.5)                #Kinetic energy evolution
        phi = np.fft.fftshift(phi)
        psi = (np.fft.fft(phi))  
        
        psi = psi*(1/(sum(abs(psi)**2*dx)))**(1./2.)       #Renormalizing step
        
        if((i)%10000 == 0):
            print "Real time step "+ str(i)+" out of "+str(timeSteps)
#            print np.max(fac*dW)
           
        if((i)%sampStep == 0):
            fileName = path+'/Psi_Step_'+str(i)+'.pkl'            
            densityPsiFile = open(fileName,'wb')
            pickle.dump(psi,densityPsiFile)
            densityPsiFile.close()
            print('Wave function saved on step '+str(i))
        
        if (i == timeSteps+1-1000):
            psiavg = psi/1000
        if (i> timeSteps+1-1000):
            psiavg = psiavg+(psi/1000)
            psiavg = psiavg*(1/(sum(abs(psiavg)**2*dx)))**(1./2.)
    return psi, psiavg
  
def ImaginaryTimeEvolve(xAr,timeSteps,dt,HPot,intPot,KEn,aN,dx,hBar,path):
    psi = ((abs(xAr)<5e-6).astype(int))
    dt = dt*(-1j)    
    for i in xrange(0,timeSteps):
        phi = np.fft.fft(psi)
        phi = np.fft.fftshift(phi)
        phi = phi*np.exp(-1j*(KEn/hBar)*dt)
        phi = np.fft.ifftshift(phi)
        psi = np.fft.ifft(phi)
        psi = psi*np.exp((-1j*(HPot+intPot*(abs(psi)**2))/hBar)*dt)
      
        psi = psi*(aN/(sum(abs(psi)**2*dx)))**(1./2.)    
        
        if(i%10000 == 0):
            print "Imaginary time step "+ str(i)+" out of "+str(timeSteps)
            
    fileName = path+'/GroundState.pkl'
    groundState = open(fileName,'wb')
    pickle.dump(psi,groundState)
    groundState.close()
    return psi