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
#factors = np.arange(10,50,10)
factors = [20]
#factors = [0]
doLogPlot = False
""" *********************************************************************** """
#beta = 0.1
beta = 0
BreakingDistancePlot = False
""" *********************************************************************** """
timeSteps = 10000000
sampTime = (timeSteps/20)
freeExpansionTimeSteps = 0
""" *********************************************************************** """
imaginaryTimeSteps = 500000
potRampTime = 1
dt = 4e-8
""" *********************************************************************** """
maxCut = 70./100.
lowCut = 20./100.
""" *********************************************************************** """
SFreq = 200
CutFreq = [150]
Alpha = 0.7
""" *********************************************************************** """
#deltaCorFactors = [0.010,0.05,0.1,0.2]
deltaCorFactors = [0.1]
sigmaCor = 1
""" *********************************************************************** """
gridPoints = 256
#L = 1500e-6
#L = 150e-6
L = 600e-6
atomNumber = 30000
""" *********************************************************************** """
#sepFactors = [2]
""" *********************************************************************** """
__path__ = os.getcwd()+'/Results/'+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
#__path__ = os.getcwd()+'/Results/'+'test'
if not os.path.exists(__path__):os.makedirs(__path__)
""" *********************************************************************** """
xAr = Def.MakeCoordinateArray(gridPoints,L)
dx = xAr[2]-xAr[1]
kxAr = Def.MakeMomentumArray(gridPoints,L)
dk = kxAr[2]-kxAr[1]
xMax = np.zeros(len(factors))
VAvg = np.zeros(len(factors))
LocLen = np.zeros(len(factors))
""" *********************************************************************** """
HPot = Def.HarmonicPotential(c.mRb,c.X_Trap_Freq,xAr)
KEn = Def.KineticEnergy(c.hb,c.mRb,kxAr)
IPot = Def.InteractionPotential(c.a,c.hb,c.mRb,c.X_Trap_Freq)
""" *********************************************************************** """
#random_Potential_1D = open('Random_Potential_1D_creator.pkl','rb')
#VRand = pickle.load(random_Potential_1D)
#random_Potential_1D.close()
#VRand = Def.PinkPotential(SFreq,CutFreq,Alpha,gridPoints)
#VRand = Def.CorelatedRandomSave(deltaCor,sigmaCor,gridPoints,None)
#with open (__path__+'/Random_Potential.pkl','wb') as vrfile:
#    pickle.dump(VRand,vrfile)
#print 'Random Potential saved to file.'
""" *********************************************************************** """
#TE.ImaginaryTimeEvolve(xAr,imaginaryTimeSteps,dt,HPot,IPot,KEn,atomNumber,dx,c.hb,__path__)
Def.PickleSave(__path__,'GroundState',Def.CreateNarrowMomentumGroundState(dk,6,2*kxAr,dk,atomNumber,dx,__path__))
""" *********************************************************************** """
#for i in range(0,len(factors)):
#    print 'Starting run '+str(i+1)+' out of '+str(len(factors))+'.'    
#    fac = factors[i]
#    __subpath__ = __path__+'/Potential_'+str(fac)+'nK'
#    if not os.path.exists(__subpath__):os.makedirs(__subpath__)
##    with open(__path__+'/GroundState.pkl','rb') as ground:
##        Psi = pickle.load(ground)
#    Psi = Def.ReadGroundState(__path__,__subpath__,xAr)
#    KenN = np.max(((abs(Psi**2)>0.001).astype(int)).nonzero())
#    VRandS = Def.RescalePotential(VRand,fac,c.Kb)
##    Psi = Def.AddMomentum(Psi, beta, c.hb, xAr,c.pRecoil)
#    Psi,PsiAvg = TE.RealTimeEvolve(Psi,timeSteps,freeExpansionTimeSteps,dt,IPot*0,KEn,atomNumber,dx,VRandS*(-1),sampTime,potRampTime,__subpath__)
#    FR.CreatePlots(xAr,dx,gridPoints,L,sampTime,dt,__subpath__)
##    FR.CreateLogPlots(xAr,dx,gridPoints,L,maxCut,lowCut,sampTime,dt,__subpath__)    
#    FR.CreateGif(__subpath__)
#    if (BreakingDistancePlot == True):    
#        xMax[i]=(Def.GetMaxDistance(Psi,lowCut,xAr))
#    if (doLogPlot == True):
#        LocLen[i] = FR.LogPlot(PsiAvg,__subpath__)
""" *********************************************************************** """
#for j in range(0,len(CutFreq)):
for j in range(0,len(deltaCorFactors)):
#for j in range(0,len(sepFactors)):
    dfac = deltaCorFactors[j]
    __subpath1__ = __path__ +'/RandomCorelation_'+str(dfac)
#    cfreqfac = CutFreq[j]
#    __subpath1__ = __path__ +'/RandomCorelation_'+str(cfreqfac)
#    sepFac = sepFactors[j]
#    __subpath1__ = __path__+'/Separation_Factor_'+str(sepFac)
    if not os.path.exists(__subpath1__): os.makedirs(__subpath1__)
    VRand = Def.CorelatedRandom(dfac,sigmaCor,gridPoints,None)
#    VRand = Def.PinkPotential(SFreq,cfreqfac,Alpha,gridPoints)
#    VRand = Def.DiffusivePlateIntensity(sepFac,gridPoints,xAr)
#    VRand = Def.PickleOpen('/home/thomas/Documents/Simulations/GPE/Results/2015-02-04_11-23-53/RandomCorelation_0.2/Random_Potential.pkl')
    with open (__subpath1__+'/Random_Potential.pkl','wb') as vrfile:
        pickle.dump(VRand,vrfile)
    print 'Random Potential saved to file. Starting set '+str(j+1)+' of ' + str(len(deltaCorFactors))+'.'
#    print 'Random Potential saved to file. Starting set '+str(j+1)+' of ' + str(len(sepFactors))+'.'
    
    for i in range(0,len(factors)):
        print 'Starting run '+str(i+1)+' out of '+str(len(factors))+'.'    
        fac = factors[i]
        __subpath2__ = __subpath1__+'/Potential_'+str(fac)+'nK'
        if not os.path.exists(__subpath2__):os.makedirs(__subpath2__)
        Psi = Def.ReadGroundState(__path__,__subpath2__,xAr)
        KenN = np.max(((abs(Psi**2)>0.001).astype(int)).nonzero())
        VRandS = Def.RescalePotential(VRand,fac,c.Kb)
        Psi = Def.AddMomentum(Psi, beta, c.hb, xAr,c.pRecoil)
        Psi,PsiAvg = TE.RealTimeEvolve(Psi,timeSteps,freeExpansionTimeSteps,dt,IPot*0,KEn,atomNumber,dx,VRandS*(-1),sampTime,potRampTime,__subpath2__)
        FR.CreatePlots(xAr,dx,gridPoints,L,sampTime,dt,__subpath2__)
        #    FR.CreateLogPlots(xAr,dx,gridPoints,L,maxCut,lowCut,sampTime,dt,__subpath__)    
        FR.CreateGif(__subpath2__)
        if (BreakingDistancePlot == True):    
            xMax[i]=(Def.GetMaxDistance(Psi,lowCut,xAr))
        if (doLogPlot == True):
            LocLen[i] = FR.LogPlot(PsiAvg,__subpath2__)
        
        with open(__subpath2__+'/SimulationInfo.txt','wb') as infofile:
            infofile.write('Simulation Conducted with the following parameters:\n'
            +'Number of Atoms = '+str(atomNumber)+'.\n'
            +'Grid Length is = '+str(L)+'m.\n'
            +'Time steps = '+str(timeSteps)+'.\n'
            +'Time step value = '+str(dt)+'s.\n'
            +'Potential ramp steps = '+str(potRampTime)+'.\n'
            +'Grid Points = '+str(gridPoints)+'.\n'
            +'An initial momentum is added to the wavefunction of '+str(beta)+' photon recoils.\n'
            +'\n'
            +'The potential used was corelated noise with:\n'
            +'deltaCorFactor =' +str(dfac)+'\n'
            +'sigmaCor =' +str(sigmaCor)+'\n'
            )
""" *********************************************************************** """

#if (BreakingDistancePlot == True):
#    VAvg = factors/c.Kb
#    plt.ion()
#    fig, ax = plt.subplots()
#    ax.plot(xMax,VAvg)
#    ax.set_xlabel('x')
#    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos:('%.01f')%(x*1e6)))
#    ax.set_ylabel('Potential [nK]')
#    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos:('%0.1f')%(x*1e9)))
#    plt.savefig(__path__+'/Distance.png',dpi = 250)
#    plt.show()
#    
#    with open(__path__ +'/xMax.pkl','wb') as xfile:
#        pickle.dump(xMax, xfile)
#    with open(__path__ + '/VAvg.pkl','wb') as vfile:
#        pickle.dump(VAvg, vfile)
    
#with open(__path__+'/SimulationInfo.txt','wb') as infofile:
#    infofile.write('Simulation Conducted with the following parameters:\n'
#    +'Number of Atoms = '+str(atomNumber)+'.\n'
#    +'Grid Length is = '+str(L)+'m.\n'
#    +'Time steps = '+str(timeSteps)+'.\n'
#    +'Time step value = '+str(dt)+'s.\n'
#    +'Potential ramp steps = '+str(5000)+'.\n'
#    +'Grid Points = '+str(4096)+'.\n'
#    +'An initial momentum is added to the wavefunction of '+str(beta)+' photon recoils.\n'
#    +'\n'
#    +'The potential used was pink noise with:\n'
#    +'Frequency Range = ' +str(SFreq)+'.\n'
#    +'Frequency cut off = ' +str(CutFreq)+'.\n'
#    +'Alpha = '+str(Alpha)+'.\n'
#    +'\n'
#    +'deltaCor =' +str(deltaCor)+'\n'
#    +'sigmaCor =' +str(sigmaCor)+'\n'
#    )
#Def.NotifyEnd()