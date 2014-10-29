# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 11:47:19 2014

@author: Matthew
"""

#import things that I need
import numpy as np
from scipy.special import fresnel
from scipy import interpolate
import matplotlib
import matplotlib.pyplot as plt
import Tkinter
import Tkconstants
import tkFileDialog
import csv
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg


#variables
#slit_width = 0.1 #slitwidth mm
#lam = 0.00024 #wavelength mm 
#xposlim = 10000 #number of x points
#xpos = np.linspace(-0.6*slit_width,0.6*slit_width,xposlim) #Moved down
#zposlim = 1000 #number of z points
#zpos = np.linspace(0,0.5,zposlim)
#LocBeamInt(xpos,zpos)


#frame
class TkDialog(Tkinter.Frame):
    def __init__(self, root):
        Tkinter.Frame.__init__(self, root)
        root.title("Spectral Processing")

        self.plotnum = Tkinter.IntVar()
        self.plotnum.set(0)
        self.slide_len = int(0)
                
        # options for a button
        button_opt = {'fill': Tkconstants.BOTH, 'padx': 5, 'pady': 1, 'side': Tkconstants.TOP}
        radio_opt = {'fill': Tkconstants.X, 'padx': 5, 'pady': 1, 'side': Tkconstants.TOP, 'anchor': Tkconstants.W}

        #Variables 
        self.Varaibles = Tkinter.LabelFrame(self, text='Variables')
        self.Varaibles.pack(**radio_opt)      

        self.slit_width_len = Tkinter.IntVar()
        self.slit_width_len.set(0.1)
        Tkinter.Label(self.Varaibles, text='Slit width (mm)').pack()
        self.slit_width = Tkinter.Entry(self.Varaibles, textvariable=self.slit_width_len)
        self.slit_width.pack(**radio_opt)

        self.lam_len = Tkinter.IntVar()
        self.lam_len.set(0.000266)
        Tkinter.Label(self.Varaibles, text='Wavelength (mm)').pack()
        self.lam = Tkinter.Entry(self.Varaibles, textvariable=self.lam_len)
        self.lam.pack(**radio_opt)

        #Parmeters 
        self.Parameters = Tkinter.LabelFrame(self, text='Model parameters')
        self.Parameters.pack(**radio_opt)      

        self.xposlim_len = Tkinter.IntVar()
        self.xposlim_len.set(10000)
        Tkinter.Label(self.Parameters, text='Slit width bins').pack()
        self.xposlim = Tkinter.Entry(self.Parameters, textvariable=self.xposlim_len)
        self.xposlim.pack(**radio_opt)

        self.zposlim_len = Tkinter.IntVar()
        self.zposlim_len.set(1000)
        Tkinter.Label(self.Parameters, text='Screen distance bins').pack()
        self.zposlim = Tkinter.Entry(self.Parameters, textvariable=self.zposlim_len)
        self.zposlim.pack(**radio_opt)

        self.zpos_max_len = Tkinter.IntVar()
        self.zpos_max_len.set(0.5)
        Tkinter.Label(self.Parameters, text='Screen distance max (mm)').pack()
        self.zpos_max = Tkinter.Entry(self.Parameters, textvariable=self.zpos_max_len)
        self.zpos_max.pack(**radio_opt)

        self.zpos_min_len = Tkinter.IntVar()
        self.zpos_min_len.set(0.01)
        Tkinter.Label(self.Parameters, text='Screen distance min (mm)').pack()
        self.zpos_min = Tkinter.Entry(self.Parameters, textvariable=self.zpos_min_len)
        self.zpos_min.pack(**radio_opt)

        #Grating length 
        self.Gratlen_len = Tkinter.IntVar()
        self.Gratlen_len.set(0.1)
        Tkinter.Label(self.Parameters, text='Grating length (mm)').pack()
        self.Gratlen = Tkinter.Entry(self.Parameters, textvariable=self.Gratlen_len)
        self.Gratlen.pack(**radio_opt)


        # creating new frame for containing plot and toolbar
        self.frame = Tkinter.Frame(root)
        self.frame.pack(side=Tkconstants.RIGHT)
        
        # creating canvas for plot
        self.f = plt.Figure()
        self.canvas = FigureCanvasTkAgg(self.f, master=self.frame)
        self.canvas.get_tk_widget().pack(side=Tkconstants.BOTTOM, fill=Tkconstants.BOTH, expand=1)
        self.canvas.show()
       
        # creating toolbar for plot
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.frame)
        self.toolbar.pack(side=Tkconstants.BOTTOM, fill=Tkconstants.BOTH, expand=1)
        self.toolbar.update()        

        self.run = Tkinter.Button(self, text='Run', command=self.run)
        self.run.pack(**radio_opt)

        self.run = Tkinter.Button(self, text='Save slice', command=self.save)
        self.run.pack(**button_opt)

    
     # initializing plot
    def initPlot(self):
        
        mutable_object = {} 

        def onclick(event):
            global chosen
            zposlim = int(self.zposlim.get())
            zpos_max = float(self.zpos_max.get())

            if event.xdata != None and event.ydata != None:
    
                #print xpos
                self.fig1.cla() #Clear previous to prevent stacking
                chosen = (event.xdata/zpos_max)*zposlim
                self.fig1.plot(xpos, Intprof_newaxis.T[chosen]) #raw scan
                #FBG spectrum
                self.fig3.cla() #clear previous
                self.fig3.plot(Wav, (FBGint_newaxis.T[chosen]/max(FBGint_newaxis.T[chosen])))
                w_min = np.min(Wav)
                w_max = np.max(Wav)
                self.fig3.axis([w_min, w_max, 0.0, 1.0])
                
                self.canvas.draw()

        self.f.clear()
        self.fig1 = self.f.add_subplot(312)
        self.fig2 = self.f.add_subplot(313)
        self.fig3 = self.f.add_subplot(311) 
        cid = self.canvas.mpl_connect('motion_notify_event', onclick)


    def save(self):
        data = Intprof_newaxis.T[chosen]

        savefilename = tkFileDialog.asksaveasfile(mode='w',defaultextension=".txt")
        for line in data:
            line = line[0]
            savefilename.write("%s\n" % str(line))
        savefilename.close()
        
        

    def run(self):
        
        global xpos
        global Intprof_newaxis
        #added global array Wav & FBGint_newaxis
        global Wav
        global FBGint_newaxis

        slit_width = float(self.slit_width.get())
        lam = float(self.lam.get())
        xposlim = int(self.xposlim.get())
        zposlim = int(self.zposlim.get())
        zpos_max = float(self.zpos_max.get())
        zpos_mim = float(self.zpos_min.get())
        Gratlen = float(self.Gratlen.get())

        zpos = np.linspace(zpos_mim,zpos_max,zposlim)

        xpos = np.linspace(-1.2*slit_width,1.2*slit_width,xposlim)

        def LocBeamInt(x, z, lam): # function works out beam intensity as a function of distance from slit (z) and position along axis of slit (x)
            alpha_1 = (x+slit_width/2)*(np.sqrt(2/(lam*z)))
            alpha_2 = (x-slit_width/2)*(np.sqrt(2/(lam*z)))
            s1, c1 = fresnel(alpha_1) #Fresnel integrals
            s2, c2 = fresnel(alpha_2) #Fresnel integrals
            I = 0.5*((c2-c1)**2+(s2-s1)**2)    
            return I;
            
        def FBGref(d, lam, z, GratLen, BasWav, Wav): #function produces a grating from the fresnel beam intensity
            #This creates the beam intensity profile from the fresnel pattern
            Intpoints = 1000
            Slitpoints = np.linspace(-10,10,Intpoints)
            Intval = LocBeamInt(Slitpoints, z, lam)
            Intval = Intval/max(Intval) #Normalise by maximum
            tck = interpolate.splrep(Slitpoints, Intval, s=0)
            #This creates the grating period and properties
            Inn  = 1.45              #base refractive index
            deltaN = 0.0002#1.6*10**-4            #change in refractive index
            Period = (BasWav)/(2*Inn) #Calculates period from base wavelength
            PerNum = int((GratLen*10**3)/Period) #Calculates the number of periods
            xp = Period/2           #positive index length
            xn = xp                 #negative index length    
            
            #This does a spline interpolation to convert the fresnel intensity into a local intensity for the grating
            xnew = np.linspace(-GratLen/2,GratLen/2,PerNum)
            LocInt = interpolate.splev(xnew, tck, der=0)
            LocInt = LocInt/max(LocInt)    
            
            #print "Base wavelength =", 2*(Inn+deltaN/2)*Period
            #print "Bandwidth =", (2*deltaN*Inn/3.14)*2*(Inn+deltaN/2)*Period
            
            refp = 0.0    
            for x in range(0,PerNum-1,1):      
                
                Inp = Inn + deltaN*LocInt[x] #high refractive index is small change on base refractive index
                
                rn = (Inp - Inn)/(Inp+Inn)
                rp = (Inn - Inp)/(Inp+Inn)
                
                thetap = (2*np.pi*Inp*xp)/Wav
                thetan = (2*np.pi*Inn*xn)/Wav
                
                refn = (rn + refp*np.exp(-2.0j * thetan)) / (1+ rn*refp*np.exp(-2.0j * thetan))
                refp = (rp + refn*np.exp(-2.0j * thetap)) / (1+ rp*refn*np.exp(-2.0j * thetap))
            
            return refp*np.conj(refp);  #returns reflection value from reflectivity 

        # plot limits
        x_min = np.min(xpos)
        x_max = np.max(xpos)
        z_min = np.min(zpos)
        z_max = np.max(zpos)

                
        Wav = np.linspace(1.5,1.6,10**3) #Grating wavelength in microns      
        
        Intprof = np.empty([xposlim, zposlim])
        FBGint = np.empty([len(Wav), zposlim])
        
        for x in range(0,zposlim,1):
            Intprof[:,x] = LocBeamInt(xpos,zpos[x], lam).real
            FBGint[:,x] = FBGref(slit_width,lam,zpos[x],Gratlen,1.55,Wav).real
        
        #d = 0.15#0.75        #slitwidth mm
        #Gratlen = 0.2       #grating length mm
        #lam = 2.66e-6       #Writing wavelength mm 0.00024
        #z = 400.0            #position of screen mm
        #BasWav = 1.55        #base wavelength in microns
        #FBGint = FBGref(d,lam,400.0,Gratlen,BasWav,Wav)
        
        print "Run complete!"
        self.initPlot()
        self.fig1 = self.f.add_subplot(312)
        self.fig2 = self.f.add_subplot(313)
        self.fig3 = self.f.add_subplot(311)
                
        im = self.fig2.imshow(Intprof, cmap=matplotlib.cm.gray, interpolation='nearest', aspect='auto', origin='lower', extent = [ (z_min), (z_max), x_min, x_max,])    
        
        self.fig2.axis([z_min, z_max, x_min, x_max])

        self.canvas.draw()
        Intprof_newaxis = np.array(Intprof)[np.newaxis]
        FBGint_newaxis = np.array(FBGint)[np.newaxis]
        #print len(Wav), len (FBGint_newaxis)

        return Intprof_newaxis, xpos, FBGint_newaxis, Wav;



if __name__ == '__main__':
    root = Tkinter.Tk()
    root.resizable(width=False, height=False)
    stat = Tkinter.BooleanVar(root)
    stat.set(True)
    TkDialog(root).pack()
    root.mainloop()