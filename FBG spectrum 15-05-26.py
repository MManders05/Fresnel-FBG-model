# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 12:33:40 2014

@author: Mark
"""


import numpy as np
import matplotlib.pyplot as plt
import time as time
from numpy.lib.scimath import sqrt as csqrt
from scipy.special import fresnel

plt.close()
timepos = time.time()

#def LocBeamInt(x, z, lam, d):
#    alpha_1 = (x+d/2)*(np.sqrt(2/(lam*z)))
#    alpha_2 = (x-d/2)*(np.sqrt(2/(lam*z)))
#    s1, c1 = fresnel(alpha_1) #Fresnel integrals
#    s2, c2 = fresnel(alpha_2) #Fresnel integrals
#    I = 0.5*((c2-c1)**2+(s2-s1)**2)    
#    return I;

def SlitPhaseP(SlitPos, x, z, k):
    rpos = z + ((x - SlitPos)**2)/(2*z)
    phase = np.exp(1.0j*k*rpos)/rpos
    return phase

def SlitPhase(x, z, lam, d): #code to work out phase at grating
    IntMax = MaskNumGlobal #10**4
    k = 2*np.pi/lam
    SlitPoint = np.linspace(-d/2,d/2,IntMax)
    ArrayPhase = SlitPhaseP(SlitPoint, x, z, k)
    totPhase = np.sum(ArrayPhase)
    totPhase = (totPhase/IntMax)
    return totPhase;

def PhasBeamAmp(MaskPos, MaskPeriod, FibPos, Knum):
    #Additional phase
    BonPhase = SlitPhase(MaskPos, SlitDist, ConLam, Slitwidth)
    
    #MaskLow
    LowPos = MaskPos - MaskPeriod/4.0
    LowR = np.sqrt((LowPos - FibPos)**2 + FibrePos**2)
    LowE = BonPhase*(np.exp(1.0j*(Knum*LowR))/(LowR))#**2)
    #MaskHi 
    HiPos = MaskPos + MaskPeriod/4.0
    HiR = np.sqrt((HiPos - FibPos)**2 + FibrePos**2)
    HiE = BonPhase*(np.exp(1.0j*(Knum*HiR+np.pi))/(HiR))#**2)
                
    LocE = (LowE + HiE)/2
    return LocE;

def PhasBeamInt(FibrePoint, MaskPeriod, MaskLen, ConLam, FibrePos):
    #Code for working out local intensity at fibre
    zrange = len(FibrePoint)
    I = np.zeros(zrange)
    
    MaskNum = int(MaskLen/MaskPeriod)
    MaskPoint = np.linspace(-MaskLen/2, MaskLen/2, MaskNum)
    Knum = 2*np.pi/ConLam
    
    for y in range(0,zrange,1):
        FibPos = FibrePoint[y]
        LocE = PhasBeamAmp(MaskPoint, MaskPeriod, FibPos, Knum)      
        FibE = np.sum(LocE)
        LocI = FibE*np.conj(FibE)
        I[y] = LocI.real
    return I;

def TMM(period, ref_mod, delta, alpha, dz, kapaac):

    T11 = np.cosh(alpha*dz)-1j*(delta/alpha)*np.sinh(alpha*dz)
    T22=np.cosh(alpha*dz)+1j*(delta/alpha)*np.sinh(alpha*dz)
    T12=-1j*(kapaac/alpha)*np.sinh(alpha*dz)
    T21=1j*(kapaac/alpha)*np.sinh(alpha*dz)

    return np.array([[T11,T12],[T21,T22]])

# constants - all dimensional units in m

r_core = 5e-6 # radius of the core
lambda_cut = 700e-9 # wavelength in mm
n_core = 1.445 # refractive idnex of the core
n_clad = np.sqrt(-(2.405**2*lambda_cut**2)/(4*np.pi**2*r_core**2)+n_core**2) # calculate cladding index given cut off wavelength and core radius and refractive index
d_n_max = 1e-4 # magnitude of index modulation


Bragg_wavelength = 1554e-9#750e-9 #this must be in units of meters
period = Bragg_wavelength/(2*n_core)
GratL = 2e-3 #m #length of the grating
no_sections = 10**4#10**4 #The number of points along the grating looked at
section_length = GratL/no_sections

lambda_min = 1548e-9#748.3e-9 #grating miniumum length
lambda_points = 500 #The number of wavelength points looked at
lambda_range = 10e-9#3e-9  #The range of wavelengths looked at
lambda_max = lambda_min + lambda_range

wavelength = np.linspace(lambda_min,lambda_max,lambda_points)

ko = 2*np.pi/wavelength #wavevector of that wavelength

#calculating effective index of fibre mode
vf=((2*np.pi/wavelength)*r_core*np.sqrt(n_core**2-n_clad**2))

u = (1+np.sqrt(2))*vf/(1+(4+vf**4)**0.25)
bfib = 1-u**2/vf**2

n_eff = np.sqrt(bfib*(n_core**2-n_clad**2)+n_clad**2) #effective refractive index of the mode
R=np.zeros((lambda_points))
dz_or = section_length #np.int(section_length/period)*period #this value is currently zero which doesn't work 
dz = dz_or
Length = dz_or*no_sections

## Mark's addition
# Adding in the modification from changing beam intensity
Slitwidth = 700e-6#1.0e-3 #slitwideth is 1 mm
ConLam = 2.48e-9#2.66e-9 #construction wavelength is 2.48nm
SlitDist = 0.2#1 #distance from slit to fibre is 20cm
z = np.linspace(-Length/2,Length/2,no_sections) #points along fibre

#Input values
MaskPeriod = 550e-9 #m #grating period of the mask
MaskLen = 1e-2 #m #length of the phase mask
MaskNumGlobal = int(MaskLen/MaskPeriod) #Global mask point number
#print MaskNumGlobal
FibrePos = 1e-2 #m #position of fibre from mask
zhalf = np.linspace(-Length/2,0,no_sections/2) # creates half the fibre length

FibIntp = PhasBeamInt(z, MaskPeriod, MaskLen, ConLam, FibrePos) #The intensity of light after the Phase Mask incident on the fibre
FibInt = FibIntp/max(FibIntp) #Normalising FibIntp
delta_n = FibInt #makes sure delta_n is the same length as FibInt, it is the refractive index modulation along the fibre

print time.time() - timepos

#Threshold code
for i in range(0,len(FibInt)):
    if FibInt[i] > 0.82: #threshold set to 80%
        delta_n[i] = d_n_max*FibInt[i] #n_core + d_n_max*FibInt[i] #FibInt[i]
    else:
        delta_n[i] = 0.0 #n_core #0
        
print time.time() - timepos

##Functioning bit of code
for i in range(0,lambda_points):
    MAT = np.eye(2)

    for m in range(no_sections-1,0,-1):
        kapadc = 4*np.pi*delta_n[m]/wavelength[i]
        kapaac = kapadc/2
        delta = kapadc+0.5*(2*ko[i]*n_eff[i]-2*np.pi/period)
        alpha = csqrt(kapaac**2-delta**2)
                
        T = TMM(period, delta_n[m], delta, alpha, dz, kapaac)
        
        MAT = np.dot(MAT,T)
        
        R[i] = np.abs(MAT[0,1]/MAT[0,0])

plt.figure(1)
plt.plot(z*10**3,delta_n/max(delta_n))
plt.title('Refractive index modulation of the fibre')
plt.ylabel('Normalised refractive index modulation (a.u.)', fontsize = 20.0)
plt.xlabel('Position along fibre (mm)', fontsize = 20.0)
plt.savefig('Refractive index modulation 150526.png')

plt.figure(2)
plt.plot(wavelength*10**6,R*10**2)
plt.title('Grating reflection intensity as a function of wavelength')
plt.ylabel('Reflection (%)', fontsize = 20.0)
plt.xlabel(u'Wavelength (\u00B5m)', fontsize = 20.0)
plt.savefig('Wavelength reflection 150526.png')

plt.show()

print time.time() - timepos