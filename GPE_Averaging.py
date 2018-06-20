import numpy as np
import os,datetime
import pickle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, FuncFormatter
import matplotlib.mlab as ml
import matplotlib.pyplot as plt

from Functions import Definitions as Def 
from Functions import TimeEvolveFunction1D as TE
from Functions import Constants
c = Constants.Constants()
from Functions import FileReader as FR
""" *********************************************************************** """
rNumb = 100
vFac = 10
dfac = 0.25
sepFac = 3
kat = 0.2*c.pRecoil
""" *********************************************************************** """
timeSteps = 5000000
sampTime = (timeSteps/20)
freeExpansionTimeSteps = 0
imaginaryTimeSteps = 500000
potRampTime = 0
dt = 1e-8 
""" *********************************************************************** """
gridPoints = 256
L = 600e-6
atomNumber = 30000
""" *********************************************************************** """
__path__ = os.getcwd()+'/Results/'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
#__path__ = os.getcwd()+'/Results/'+'test'
if not os.path.exists(__path__):os.makedirs(__path__)
__support__ = __path__+'/SupportFiles'
if not os.path.exists(__support__):os.makedirs(__support__)
__averages__ = __path__+'/Averages'
if not os.path.exists(__averages__):os.makedirs(__averages__)
""" *********************************************************************** """    
xAr = Def.MakeCoordinateArray(gridPoints,L)
dx = xAr[2]-xAr[1]
kxAr = Def.MakeMomentumArray(gridPoints,L)
dk = kxAr[2]-kxAr[1]
HPot = Def.HarmonicPotential(c.mRb,c.X_Trap_Freq,xAr)
KEn = Def.KineticEnergy(c.hb,c.mRb,kxAr)
IPot = Def.InteractionPotential(c.a,c.hb,c.mRb,c.X_Trap_Freq)
""" *********************************************************************** """
Psi = Def.CreateNarrowMomentumGroundState(dk,6,kxAr,kat,atomNumber,dx,None)
Def.PickleSave(__support__,'GroundState',Psi)
""" *********************************************************************** """
#Vrand = Def.CorelatedRandom(dfac,1,gridPoints,None)
Vrand = Def.DiffusivePlateIntensity(sepFac,gridPoints,xAr) 
Def.PickleSave(__support__,'RandomPotential',Vrand)
Vrand = Def.RescalePotential(Vrand,vFac,c.Kb)
""" *********************************************************************** """
for i in xrange(0,rNumb):
    print "Starting run "+str(i+1)+" out of "+str(rNumb)+"." 
#    Vrand = Def.CorelatedRandom(dfac,1,gridPoints,None) 
#    Def.PickleSave(__support__,'RandomPotential'+str(i+1),Vrand)
#    Vrand = Def.RescalePotential(Vrand,vFac,c.Kb)
    __subpath__ = __path__ + '/Results_Run_'+str(i+1)
    if not os.path.exists(__subpath__):os.makedirs(__subpath__)
    Psi = Def.PickleOpen(__support__+'/GroundState.pkl')    
    Psi,PsiAvg = TE.RealTimeEvolve(Psi,timeSteps,freeExpansionTimeSteps,dt,IPot*0,KEn,atomNumber,dx,Vrand,sampTime,potRampTime,__subpath__)
    FR.CreatePlots(xAr,dx,gridPoints,L,sampTime,dt,__subpath__)
    FR.CreateGif(__subpath__)
    Def.PickleSave(__averages__,'Psi'+str(i+1),PsiAvg)
    print "Average wafefunction saved to file"
    
    with open(__subpath__+'/SimulationInfo.txt','wb') as infofile:
        infofile.write('Simulation Conducted with the following parameters:\n'
        +'Number of Atoms = '+str(atomNumber)+'.\n'
        +'Grid Length is = '+str(L)+'m.\n'
        +'Time steps = '+str(timeSteps)+'.\n'
        +'Time step value = '+str(dt)+'s.\n'
        +'Potential ramp steps = '+str(potRampTime)+'.\n'
        +'Grid Points = '+str(gridPoints)+'.\n'
        +'An initial momentum is added to the wavefunction of '+str(kat/c.pRecoil)+' photon recoils.\n'
        +'\n'
        +'The potential used was corelated noise with:\n'
        +'deltaCorFactor =' +str(dfac)+'\n'
        )

Psi = FR.ReadAverages(__averages__)
Def.PickleSave(__averages__,'PsiFinal',Psi)  
    