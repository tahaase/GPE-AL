# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 13:26:51 2015

@author: thomas
"""
import numpy as np 
import pickle
class TimeEvolve():
    def __init__(self,timeSteps,sampTime,imaginaryTime,freeTime,rampTime,dt,dx,aNumb):
        self.timeSteps = timeSteps
        self.imaginaryTimeSteps = imaginaryTime
        self.sampleTime = sampTime
        self.freeExpansionTime = freeTime
        self.potentialRampTime = rampTime
        self.dt = dt
        self.dx = dx
        self.atomNumber = aNumb
        
    def RealTimeEvolve(self,psi,intPot,KEn,Vrand,__subpath__):
        hBar = 1.054e-34 
        factor = np.linspace(0,1,self.potentialRampTime+1)
        psiavg = 0
        phi = np.fft.ifftshift(np.fft.ifft(psi))
        phi = phi*np.exp(-1j*(KEn/hBar)*self.dt*self.freeExpansionTime)
        psi = np.fft.fft(np.fft.fftshift(phi))
        psi = psi*(1/(sum(abs(psi)**2*self.dx)))**(1./2.)
        for i in xrange(1,self.timeSteps+1):    
            #Real time split step evolution
            if (i < self.potentialRampTime):
                fac = factor[i+1]
            else:
                fac = 1   
                
            phi = np.fft.ifft(psi)                            #Momentum space wave function
            phi = np.fft.ifftshift(phi)  
            phi = phi*np.exp(-1j*(KEn/hBar)*self.dt*0.5)                #Kinetic energy evolution
            phi = np.fft.fftshift(phi)
            psi = (np.fft.fft(phi))                             #Coordinate space wave function
            psi = psi*np.exp(-1j*(((self.atomNumber*intPot)*(abs(psi)**2)+fac*Vrand)/hBar)*self.dt) #Interaction potential evolution        
            phi = np.fft.ifft(psi)                            #Momentum space wave function
            phi = np.fft.ifftshift(phi)  
            phi = phi*np.exp(-1j*(KEn/hBar)*self.dt*0.5)                #Kinetic energy evolution
            phi = np.fft.fftshift(phi)
            psi = (np.fft.fft(phi))  
            
            psi = psi*(1/(sum(abs(psi)**2*self.dx)))**(1./2.)       #Renormalizing step
            
            if((i)%10000 == 0):
                print "Real time step "+ str(i)+" out of "+str(self.timeSteps)
    #            print np.max(fac*dW)
               
            if((i)%self.sampleTime == 0):
                fileName = __subpath__+'/Psi_Step_'+str(i)+'.pkl'            
                densityPsiFile = open(fileName,'wb')
                pickle.dump(psi,densityPsiFile)
                densityPsiFile.close()
                print('Wave function saved on step '+str(i))
            
            if (i == self.timeSteps+1-1000):
                psiavg = psi/1000
            if (i> self.timeSteps+1-1000):
                psiavg = psiavg+(psi/1000)
                psiavg = psiavg*(1/(sum(abs(psiavg)**2*self.dx)))**(1./2.)
        return psi, psiavg
        
    def ImaginaryTimeEvolve(self,xAr,HPot,intPot,KEn,Vrand,path):
        psi = ((abs(xAr)<5e-6).astype(int))
        dt = self.dt*(-1j)
        hBar = 1.054e-34 
        for i in xrange(0,self.imaginaryTimeSteps):
            phi = np.fft.fft(psi)
            phi = np.fft.fftshift(phi)
            phi = phi*np.exp(-1j*(KEn/hBar)*dt)
            phi = np.fft.ifftshift(phi)
            psi = np.fft.ifft(phi)
            psi = psi*np.exp((-1j*(HPot+Vrand+intPot*(abs(psi)**2))/hBar)*dt)
      
            psi = psi*(self.atomNumber/(sum(abs(psi)**2*self.dx)))**(1./2.)    
            
            if(i%10000 == 0):
                print "Imaginary time step "+ str(i)+" out of "+str(self.imaginaryTimeSteps)
            
        fileName = path+'/GroundState.pkl'
        groundState = open(fileName,'wb')
        pickle.dump(psi,groundState)
        groundState.close()
        return psi
                    