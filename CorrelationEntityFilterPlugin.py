#!/usr/bin/env python
from cPlugin import FilterPlugin
import numpy as np
import os
from Tkinter import *
import Pmw

class CorrelationEntityPlugin(FilterPlugin):
    """
    Depends:  get_calphas
    Provides: 
    """
    def __init__(self):                 # init must have no arguments
        FilterPlugin.__init__(self)
        # data
        # parameters
        self.max = 1.0              # the matrix file
        self.min = 0.4               # cutoff

# --------------
    def get_filter(self, data, reference, *args, **kwargs):
        """
        Generate a filtering matrix of type=np.bool, to filter
        *data*, a matrix referring to the *reference* object in PyMOL.
        """
        if not self.is_activated:
            return np.ones(data.shape, dtype=np.bool)
        
        # retrieve parameters from gui
        self.min = self.minf.get()
        self.max = self.maxf.get()
        
        return self.filter(data, self.min, self.max)

# --------------
    def filter(self, data, min, max):
        return (abs(data)>=min) & (abs(data)<=max)

# --------------
    def GUI(self, page):
        """
        Add the specific GUI elements for this plugin to the provided Tk widget.
        """
        group = Pmw.Group(page,tag_text='Matrix values')
        group.pack()
        frame = Frame(group.interior())
        frame.pack()
        self.minf = self.makeDoubleVar(frame, "Low cut-off", self.changeCcol, tooltip="Only correlations with absolute value\ngreater than this cutoff will\nbe stored.")
        self.minf.set(self.min)
        self.maxf = self.makeDoubleVar(frame, "High cut-off", self.changeCcoh, tooltip="Only correlations with absolute value\nsmaller than this cutoff will\nbe stored.")
        self.maxf.set(self.max)
        # self.balloon = Pmw.Balloon(frame)      
        # self.balloon.bind(self.cco_mode, """Correlation cutoff mode""")
        return group

        
    def changeCcol(self, a, *args):
        if len(args)!=0: 
            a=args[0]
        val = max(0.0, float(self.minf.get()) + float(a)*0.05)
        self.minf.set(val)
        # self.msg("Correlation Cutoff (low) changed to %f"%val)
   
    def changeCcoh(self, a, *args):
        if len(args)!=0: 
            a=args[0]
        val = max(0.0, float(self.maxf.get()) + float(a)*0.05)
        self.maxf.set(val)
        # self.msg("Correlation Cutoff (high) changed to %f"%val)

#---- makevar methods copied from cPlot: may be implemented as services
    def makeDoubleVar(self, parent, label, callback, tooltip=None):
        var = DoubleVar()
        self.makeVar(var, parent, label, callback, tooltip=tooltip)
        return var
    
    def makeVar(self, var, parent, label, callback, tooltip=None):
        fr = Frame(parent)
        loc = Entry(fr,textvariable=var,bg='black',fg='green')
        scr = Scrollbar(fr,orient="horizontal",command=callback)
        lab = Label(fr, text=label)
        scr.pack(side=LEFT)
        loc.pack(side=LEFT)
        lab.pack(side=LEFT)
        fr.pack(fill='x',padx=4,pady=1) # vertical
        # if tooltip:
        #     self.balloon.bind(fr, tooltip)
    
# --------------
    def check(self):
        """
        Check wether the plugin is ready for a *get_filter* call (e.g. all GUI inputs
        are valid, all information is available, ecc.).

        Inputs are checked on change, and have default values: return True
        """
        try:
            minf = self.minf.get()
        except:
            return (False,"The lower cut-off value was not correctly specified")
        try:    
            maxf = self.maxf.get()
        except:
            return (False,"The higher cut-off value was not correctly specified")
        #if maxf > 1.0 or maxf < 0.0 or minf > 1.0 or minf < 0.0:
        #    return (False,"The cut-off values must be between 0.0 and 1.0")
        if maxf <= minf:
            return (False,"The higher cut-off value must be higher than the lower cut-off value")
        return (True,"")

        

