# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 10:08:54 2014

@author: thomas
"""

def parabolicEdgeFilter(xAr,L):
    fil = -(xAr-L/2.)*(xAr+L/2.)
    fil = fil/(fil.max()/10)
    m = (fil<1).astype(int)
    n = (fil>1).astype(int)
    fil = fil*m+n
    return fil