# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 16:55:36 2014

@author: thomas
"""
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, FuncFormatter
import matplotlib.mlab as ml
import matplotlib.pyplot as plt
import numpy as np
import pickle
from PIL import Image
from math import *
import natsort

from Im2gif import writeGif
import os

plt.close('all')
def CreatePlots(xAr,dx,gP,L,samptime,dt,subpath):
    xAv = []
    xVar = []
    plotNumb = 0
    for j in os.listdir(subpath):
        if(j.endswith('.pkl')):
            plotNumb = plotNumb+1
    tAr = np.linspace(1,plotNumb,plotNumb-1)*dt*samptime
    plt.ioff()
    for i in range(1,plotNumb):
        f = i*samptime
        fileName = subpath+'/Psi_Step_'+str(f)+'.pkl'
        State = open(fileName,'rb')
        psi = pickle.load(State)
        State.close()
        print 'Data file for time step '+str(f)+' loaded.'
        
        Psi = abs(psi**2)
                
        xAv.append(sum(np.conjugate(psi)*xAr*psi*dx)/sum(Psi*dx))
        xVar.append(sum(np.conjugate(psi)*xAr**2*psi*dx)/sum(Psi*dx))
        
        fig, ax1 = plt.subplots()
        ax1.plot(xAr,Psi)
        ax1.text(xAr[10],Psi.max()*(1.1),'Elapsed time = '+str(f*dt)+'s')
        ax1.set_ylabel('$|\psi|^2$')
        ax1.set_ylim(0,Psi.max()*1.2)
        ax1.set_xlabel('x [$\mu$m]')
        ax1.set_xlim(xAr.min(),xAr.max())
        ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, pos:('%.1f')%(x*1e6)))
        
        plt.subplots_adjust(top=0.85)
        plt.tight_layout()
        plt.savefig(subpath+'/Psi_Step_'+str(i)+'_.png',dpi = 250)
        plt.close('all')
        
    fig, (ax1,ax2) = plt.subplots(1,2)
    ax1.plot(tAr,np.real(xAv))
    ax1.set_ylim(-L/6,L/6)
    ax2.plot(tAr,np.real(xVar))
    ax2.set_ylim(0,np.max(np.real(xVar))*1.1)
    plt.savefig(subpath+'/Variance.png',dpi=250)
    plt.close('all')

def CreateLogPlots(xAr,dx,gP,L,maxCut,minCut,samptime,dt,subpath):
    xAv = []
    xVar = []
    plotNumb = 0
    for j in os.listdir(subpath):
        if(j.endswith('.pkl')):
            plotNumb = plotNumb+1
    tAr = np.linspace(1,plotNumb,plotNumb-1)*dt*samptime
    plt.ioff()
    for i in range(1,plotNumb):
        f = i*samptime
        fileName = subpath+'/Psi_Step_'+str(f)+'.pkl'
        
        State = open(fileName,'rb')
        psi = pickle.load(State)
        State.close()
        print 'Data file for time step '+str(f)+' loaded.'
        
        Psi = abs(psi**2)
        logD = np.log(Psi)        
        xAv.append(sum(np.conjugate(psi)*xAr*psi*dx)/sum(Psi*dx))
        xVar.append(sum(np.conjugate(psi)*xAr**2*psi*dx)/sum(Psi*dx))
      
        numbsa = ((Psi > (Psi.max()*maxCut)).astype(int)).nonzero()
        numbsb = ((Psi > (Psi.max()*minCut)).astype(int)).nonzero()
        n1 = np.min(numbsa)
        n2 = np.max(numbsa)
        n3 = np.min(numbsb)
        n4 = np.max(numbsb)
            
        m1,b1 = ExponentialFit(Psi,xAr,n3,n1)
        m2,b2 = ExponentialFit(Psi,xAr,n2,n4)
            
        Exp1 = m1*xAr+b1
        Exp2 = m2*xAr+b2
                
        fig, (ax1,ax2) = plt.subplots(1,2)
        ax1.plot(xAr,logD,'b',xAr,Exp1,'r',xAr,Exp2,'r')
        ax1.plot(xAr,logD)
        ax1.set_xlabel('x [$\mu$m]')
        ax1.set_xlim(np.min(xAr),np.max(xAr))
        ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, pos:('%1.f')%(x*1e6)))
        ax1.set_ylabel('Log $|\psi|^2$')
        ax1.set_ylim(10,logD.max()+1)
        
        ax2.plot(xAr,Psi)
        ax2.text(xAr[10],Psi.max()*(1.1),'Elapsed time = '+str(f*dt)+'s')
        ax2.set_ylabel('$|\psi|^2$')
        ax2.set_ylim(0,Psi.max()*(1.2))
        ax2.set_xlabel('x [$\mu$m]')
        ax2.set_xlim(xAr.min(),xAr.max())
        ax2.xaxis.set_major_formatter(FuncFormatter(lambda x, pos:('%.1f')%(x*1e6)))
        
        plt.subplots_adjust(top=0.85)
        plt.tight_layout()
        plt.savefig(subpath+'/Psi_Step_'+str(i)+'_.png',dpi = 250)
        plt.close('all')
    
    fig, (ax1,ax2) = plt.subplots(1,2)
    ax1.plot(tAr,np.abs(xAv))
    ax1.set_ylim(-L/4,L/4)
    ax2.plot(tAr,np.abs(xVar))
    ax2.set_ylim(0,np.max(np.abs(xVar))*1.1)
    plt.savefig(subpath+'/Variance.png',dpi=250)
    plt.close('all')
        
def CreateGif(subpath):
    file_names = natsort.natsorted((fn for fn in os.listdir(subpath) if fn.endswith('.png')),key = lambda y:y.lower())
    plots = [Image.open(subpath+'/'+fn) for fn in file_names]
    
    size = (1024,1024)
    
    for im in plots:
        im.thumbnail(size, Image.ANTIALIAS)
        
    filename = subpath+'/Plot.gif'
    writeGif(filename,plots, duration = 0.8)

def LogPlot(Psi):
    loc = 0
    with open(path +'/PsiAverage.pkl','wb') as xfile:
        Psi = pickle.load(xfile)
    return loc
    
def ExponentialFit(Psi,x,nu1,nu2):
    logD = np.log(Psi[nu1:nu2])
    X = x[nu1:nu2]
    m = np.polyfit(X,logD,1)
    return m
    
def ReadAverages(path):
    plotNumb = 0
    psiavg = 0
    for j in os.listdir(path):
        if(j.endswith('.pkl')):
            plotNumb = plotNumb+1
    for i in range (1,plotNumb+1):
        fileName = path+'/Psi'+str(i)+'.pkl'
        State = open(fileName,'rb')
        psi = pickle.load(State)
        State.close()
        if (i==1):
            psiavg = psi
        else:
            psiavg = psiavg + psi
        psiavg = psiavg / plotNumb
    print "Averaged over runs" 
    return psiavg