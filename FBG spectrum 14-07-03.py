# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 12:33:40 2014

@author: Mark
"""

#import scipy
import numpy as np
from scipy.special import fresnel
#import matplotlib
import matplotlib.pyplot as plt

def LocBeamInt(x, z):
    alpha_1 = (x+d/2)*(np.sqrt(2/(lam*10**z)))
    alpha_2 = (x-d/2)*(np.sqrt(2/(lam*10**z)))
    s1, c1 = fresnel(alpha_1) #Fresnel integrals
    s2, c2 = fresnel(alpha_2) #Fresnel integrals
    I = 0.5*((c2-c1)**2+(s2-s1)**2)    
    return I;

def LocBeamInt2(x, z):
    s3, c3 = fresnel(x*(np.sqrt(2/(lam*10**z))))
    I = 0.5*((c3[0:steplen]-c3[steplen:])**2+(s3[0:steplen]-s3[steplen:])**2)    
    return I;
    
def BeamShape(x):
    A = 1.0                             #Beam intensity normalised to maximum of 1
    B = GratLenOutside/2                #Maximum intensity assumed to be at the centre of the fibre
    C = 500.0                           #sigma of the beam
    I = A*np.exp(-((x-B)**2/(2*C**2)))  #Intensity of beam
    return I;

def CompRef(walam, xpoz, zpoz):
    Inn = 1.5               #base refractive index 
    xp = Period/2           #length of half period
    xn = xp                 #xn and xp same length

    GratLen = GratLenOutside
    Intermin = LocBeamInt(xpoz, zpoz).real
    
    refp = 0.0
    for x in range(0,GratLen-1,1):      
        
        Inp = Inn + 0.0001*Intermin[x]  #high refractive index is small change on base refractive index
        
        rn = (Inp - Inn)/(Inp+Inn)
        rp = (Inn - Inp)/(Inp+Inn)
        
        thetap = (2*np.pi*Inp*xp)/walam
        thetan = (2*np.pi*Inn*xn)/walam
        
        refn = (rn + refp*np.exp(-2.0j * thetan)) / (1+ rn*refp*np.exp(-2.0j * thetan))
        refp = (rp + refn*np.exp(-2.0j * thetap)) / (1+ rp*refn*np.exp(-2.0j * thetap))
        #print x
    return refp*np.conj(refp);  #returns reflection value from reflectivity  


Period = 0.5 # period length microns
GratLenOutside = 2000 #2000 #periods
GratLimP = Period*GratLenOutside*10**-3 ##length of grating in mm
pointval = 500 #number of points
d = 1.0#0.75        #slitwidth mm
lam = 2.66e-4       # wavelength mm
zposlim = 500
RefSpace = np.empty([pointval, zposlim]) #Array of Reflection values
IntSpace = np.empty([pointval, zposlim])
xpos1 = np.linspace(-GratLimP/2,GratLimP/2,GratLenOutside)# distance in x direction mm
xpos2 = np.linspace(-GratLimP/2,GratLimP/2,pointval)
zpos = np.linspace(0.5,3.5,zposlim) # log10 of the distance to the screen mm
lampos = np.linspace(1.498,1.502,pointval) #wavelength range microns



dpos = np.linspace(0.5,1.5,zposlim)
for y in range(0,zposlim,1):
    d = dpos[y]
    #steplen = (GratLenOutside+(d/Period)*10**3)/2
    #xpos3 = np.linspace(-(GratLimP+d)/2, (GratLimP+d)/2, steplen*2)
    for x in range(0,pointval,1):
        RefSpace[x,y] = CompRef(lampos[x], xpos1, 2.5).real
        #IntSpace[x,y] = LocBeamInt(xpos2[x], 2.5)

#for tep in range(1,9,1):
#    y = 10*tep
##for y in range(0,zposlim,1):
#    for x in range(0,pointval,1):
#        RefSpace[x,y] = CompRef(lampos[x], xpos1, zpos[y]).real
#        IntSpace[x,y] = LocBeamInt(xpos2[x], zpos[y])
#    
#    
#    yshow = y
##Wavelength profile of grating
#    fig1 = plt.figure(num=1, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
#    plt.plot(xpos2, IntSpace[:,yshow]/max(IntSpace[:,yshow]))
#    plt.title(u'Incident beam spatial profile', fontsize = 20.0)
#    plt.ylabel(u'Intensity (a.u.)', fontsize = 20.0)
#    plt.xlabel(u'Position (mm)', fontsize = 20.0)
#
##Wavelength profile of grating
#    fig2 = plt.figure(num=2, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
#    plt.plot(lampos*10**3, RefSpace[:,yshow]*100)
#    plt.title(u'Grating wavelength profile', fontsize = 20.0)
#    plt.ylabel(u'Reflectivity (%)', fontsize = 20.0)
#    plt.xlabel(u'wavelength (nm)', fontsize = 20.0)

# plot limits
x_min = np.min(lampos)
x_max = np.max(lampos)
z_min = np.min(dpos)
z_max = np.max(dpos)

zplot = 3#310 # distance at which you would  like to plot the intesity pattern mm

fig3 = plt.figure(3)
plt.imshow(RefSpace, cmap=matplotlib.cm.gray, interpolation='nearest', aspect='auto', origin='lower', extent = [ (z_min), (z_max), x_min, x_max,] )# [ 10**(z_min), 10**(z_max), x_min, x_max,] )
#plt.xscale('log')
plt.title("Wavelength reflectivity intensity " + str(zposlim) + "pts")
#plt.xlabel('z (mm)')
plt.xlabel('Slitwidth (micron)')
plt.ylabel('Wavelength (mm)')
plt.show()