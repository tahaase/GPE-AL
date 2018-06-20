
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 14:37:54 2014

@author: thomas
"""
from math import pi

class Constants():
    def __init__(self):
        self.hb = 1.054e-34
        self.Kb = 1.3806e-23
        self.h = 6.626e-34
        self.mRb = 87*1.66e-27
        self.a = 5.186e-9
        self.X_Trap_Freq = 150
        self.pRecoil = self.h/780.2e-9
        self.kRecoil = 2*pi/780.2e-9