"""
The Following script GPE_Averaging_Par is written bassed on the multiprocessing package do allow each process to solve for one 
time evolution run. The scripts creates a pool of 4 processes and uses the apply call on the pool to manage the primary for loop. 

After concluded the script reads the results of each average and saves the result into the numpy arrays Psi, APsi and LogPsi. 
A pickle file with the information is also created. 

Author: Thomas Haase
"""
import numpy as np
import os,datetime
#import pickle
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter, FuncFormatter
#import matplotlib.mlab as ml
#import matplotlib.pyplot as plt
import multiprocessing as mp

from Functions import Definitions as Def 
from Functions import TimeEvolveFunction
from Functions import Constants
c = Constants.Constants()
from Functions import FileReader as FR
""" *********************************************************************** """
"""
Definitions and path creation block. Initializes all required values. 
"""
""" *********************************************************************** """
rNumb = 4
vFac = 0.01
dfac = 0.25
sepFac = 2
kat = 0
corLen = 0.15e-6
""" *********************************************************************** """
time = 0.1
dt = 1e-8
timeSteps = int(time/dt)
sampTime = int((timeSteps/20))
freeExpansionTimeSteps = 0
imaginaryTimeSteps = 800000
potRampTime = 0

""" *********************************************************************** """
gridPoints = 4096
#L = 1000e-6
L = 400e-6
#L = 300e-6
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
#Vrand = Def.CorelatedRandom(dfac,1,gridPoints,None)
#Vrand = Def.DiffusivePlateIntensity(sepFac,gridPoints,xAr) 
Vrand = Def.SpeckleIntensity(kxAr,corLen,gridPoints)
Def.PickleSave(__support__,'RandomPotential',Vrand)
Vrand = Def.RescalePotential(Vrand,vFac,c.Kb)
""" *********************************************************************** """
"""
Second block initializes the time evolution class. 
"""
""" *********************************************************************** """
TE = TimeEvolveFunction.TimeEvolve(timeSteps,sampTime,imaginaryTimeSteps,freeExpansionTimeSteps,potRampTime,dt,dx,atomNumber)
""" *********************************************************************** """
#Psi = TE.ImaginaryTimeEvolve(xAr,HPot,IPot,KEn,Vrand*0,__support__)
Psi = Def.CreateNarrowMomentumGroundState(0.01,kxAr,kat,atomNumber,dx,c.kRecoil)
Def.PickleSave(__support__,'GroundState',Psi)
""" *********************************************************************** """
"""
Defines the main loop from which the time evolution is conducted. 
"""
""" *********************************************************************** """
def Main_Loop(i,rNumb,__path__,__support__):
    print "Starting run "+str(i+1)+" out of "+str(rNumb)+"." 
#    Vrand = Def.CorelatedRandom(dfac,1,gridPoints,None) 
#    Vrand = Def.DiffusivePlateIntensity(sepFac,gridPoints,Def.MakeCoordinateArray(gridPoints,1e-4)) 
#    Def.PickleSave(__support__,'RandomPotential'+str(i+1),Vrand)
#    Vrand = Def.RescalePotential(Vrand,vFac,c.Kb)
    __subpath__ = __path__ + '/Results_Run_'+str(i+1)
    if not os.path.exists(__subpath__):os.makedirs(__subpath__)
    Psi = Def.PickleOpen(__support__+'/GroundState.pkl')    
    Psi, PsiAvg = TE.RealTimeEvolve(Psi,IPot*0,KEn,Vrand,__subpath__)
    Def.PickleSave(__averages__,'Psi'+str(i+1),PsiAvg)
#    FR.CreatePlots(xAr,dx,gridPoints,L,sampTime,dt,__subpath__)
#    FR.CreateGif(__subpath__)
#    Def.PickleSave(__averages__,'Psi'+str(i+1),PsiAvg)
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
""" *********************************************************************** """        
"""
Defines the pool for multiple processes. Uses the default number of cores of the computer. Each core does one run of the rNumb simulations. 
""" 
""" *********************************************************************** """
if __name__ == '__main__':
    pool = mp.Pool(processes=4)    
    [pool.apply_async(Main_Loop,args = (x,rNumb,__path__,__support__,)) for x in xrange(rNumb)]
    pool.close()
    pool.join()


""" *********************************************************************** """ 
"""
Final Block. Averaging and Analysis of simulations. 
""" 
""" *********************************************************************** """
#Def.NotifyEnd()
Psi = FR.ReadAverages(__averages__)
Def.PickleSave(__averages__,'PsiFinal',Psi)
APsi = np.abs(Psi)**2
LogPsi = np.log(APsi)
""" *********************************************************************** """  