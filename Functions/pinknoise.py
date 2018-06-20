# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 16:06:33 2014

@author: thomas
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
from math import pi

def one_over_f(f, knee, alpha):
    desc = np.ones_like(f)
    desc[f<KNEE] = np.abs((f[f<KNEE]/KNEE)**(-alpha))
    desc[0] = 1
    return desc
plt.close('all')
white_noise_sigma =  1 #mK * sqrt(s)
 
SFREQ = 200 #Hz
#KNEE = 1000 / 1e3 #Hz
KNEE = 10
ALPHA = .7
#N = SFREQ * 3600 * 2 # 4 hours
N = 4096 
#generate white noise in time domain
wn=np.random.normal(0.,white_noise_sigma*np.sqrt(SFREQ),N)
 
#shaping in freq domain
s = np.fft.rfft(wn)
f = np.fft.fftfreq(N, d=1./SFREQ)[:len(s)]
f[-1]=np.abs(f[-1])
fft_sim = s * one_over_f(f, KNEE, ALPHA)
T_sim = np.fft.irfft(fft_sim)
T_sim = (T_sim-T_sim.min())/T_sim.max() 
#PSD - 1 hour window
NFFT = int(SFREQ*60*60*1)
s_sim, f_sim  = mlab.psd(T_sim, NFFT=NFFT, Fs=SFREQ, scale_by_freq=True)
x = np.linspace(0,len(wn),len(wn))

#plot
plt.figure()
plt.plot(f_sim, np.sqrt(s_sim), label='sim')
plt.loglog(f_sim, one_over_f(f_sim, KNEE, ALPHA) * white_noise_sigma*1e3*np.sqrt(2), 'r',label='noise model')
plt.vlines(KNEE,*plt.ylim())
plt.grid(); plt.xlabel('Freq'); plt.title('Amplitude spectrum'); plt.legend()

plt.figure()
plt.plot(x,T_sim)