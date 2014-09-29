# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 12:33:40 2014

@author: Mark
"""

#import scipy
import numpy as np
from scipy.special import fresnel
import matplotlib
import matplotlib.pyplot as plt

def LocBeamInt(x, z): # function works out beam intensity as a function of distance from slit (z) and position along axis of slit (x)
    alpha_1 = (x+d/2)*(np.sqrt(2/(lam*z)))
    alpha_2 = (x-d/2)*(np.sqrt(2/(lam*z)))
    s1, c1 = fresnel(alpha_1) #Fresnel integrals
    s2, c2 = fresnel(alpha_2) #Fresnel integrals
    I = 0.5*((c2-c1)**2+(s2-s1)**2)    
    return I;

d = 1.0             #slitwidth mm
lam = 2.66e-4       #wavelength mm

xposlim = 10 #number of x points
xpos = np.linspace(-d,d,xposlim)

zposlim = 10 #number of z points
zpos = np.linspace(0,10,zposlim)

LocalIntensity = LocBeamInt(xpos,zpos)

# plot limits
x_min = np.min(xpos)
x_max = np.max(xpos)
z_min = np.min(zpos)
z_max = np.max(zpos)

fig3 = plt.figure(1)
plt.imshow(LocalIntensity, cmap=matplotlib.cm.gray, interpolation='nearest', aspect='auto', origin='lower', extent = [ (z_min), (z_max), x_min, x_max,]) 
plt.title("Beam intensity")
plt.xlabel('Distance along slit (mm)')
plt.ylabel('Distance from slit (mm)')
plt.show()
