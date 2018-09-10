#!/usr/bin/env python
#  File: Distance3DFilterPlugin.py 
#
#  Copyright (C) 2011 Marco Pasi <mf.pasi@gmail.com> 
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# "Distance3DFilterPlugin.py v0.1 (C) 2011 Marco Pasi"
#
from cPlugin import FilterPlugin
import numpy as np
import os
from Tkinter import *
import Pmw
from math import floor

class SequenceDistancePlugin(FilterPlugin):
    """
    Depends:  get_calphas, parse_dat
    Provides: 
    """
    def __init__(self):                 # init must have no arguments
        FilterPlugin.__init__(self)
        # data
        # self.dist       = 0          # cutoff distance
        # parameters
        self.seq_co = 12              # the matrix file
        self.mode = 0            # cutoff

# --------------
# Sequence distance filter
    def filter(self, mode, dist, reference):
        """
        Sequence distance filter

        dist: sequence distance cut-off
        mode: 0/1. Switch filtering mode:
            0 filters residues WITHIN distance dist (abs(i-j)<dist), 
            1 filters residues OUTSIDE sequence distance dist (abs(i-j)>dist),
        """
        cas = self.services.provideServiceByName("get_calphas")(reference)
        chains = [ i.chain for i in cas ]
        filter_matrix=np.ones((len(cas),len(cas)), dtype=np.bool)
        if mode == 0: 
            for i in range(len(cas)):
                for j in range(len(cas)):
                    if abs(i-j)>dist and chains[i] == chains[j]:
                        filter_matrix[i][j] = False
        else:
            for i in range(len(cas)):
                for j in range(len(cas)):
                    if abs(i-j)<=dist and chains[i] == chains[j]:
                        filter_matrix[i][j] = False
        return filter_matrix

# --------------
    def get_filter(self, data, reference, *args, **kwargs):
        """
        Generate a filtering matrix of type=np.bool, to filter
        *data*, a matrix referring to the *reference* object in PyMOL.
        """
        if not self.is_activated:
            return np.ones(data.shape, dtype=np.bool)

        # retrieve parameters from gui
        self.seq_co = self.seq_cof.get()
        #a if b else c
        self.mode = 0 if self.co_modef.getvalue() == "within distance" else 1        

        return self.filter(self.mode, self.seq_co, reference)

# --------------
    def GUI(self, page):
        """
        Add the specific GUI elements for this plugin to the provided Tk widget.
        """
        seqdistG = Pmw.Group(page,tag_text='Sequence Distance')
        seqdistG.pack()
        seqdistF = Frame(seqdistG.interior())
        seqdistF.pack()

        # input file
        # distance cut-off
        self.seq_cof = self.makeDoubleVar(seqdistF, "Distance cut-off", self.change_seq_co, tooltip="""\
Only correlations between residues with
sequence distance lower/higher than this cutoff
will be plotted.""")
        self.seq_cof.set(self.seq_co)
        # distance cut-off mode
        self.co_modef = Pmw.OptionMenu(seqdistF,
                                       labelpos = 'w',
                                       label_text = 'Cut-off mode',
                                       items = ["within distance","over distance"],
                                       initialitem = "within distance"
                                       )
        self.co_modef.pack(fill='x',expand=1,padx=4)
        #self.balloon.bind(self.co_modef, """Sequence distance cutoff mode""")
        return seqdistG

#----
    def change_seq_co(self, a,*args):
        if len(args)!=0: 
            a=args[0]
        val = int(self.seq_cof.get()) + int(a)
        print val
        self.seq_cof.set(val)
        #self.msg("Sequence Distance Cutoff changed to %d"%val)
        
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
        #chains=list(set([ i.chain for i in self.services.provideServiceByName("get_calphas")(reference) ]))
        try:
            seq_co = self.seq_cof.get()
        except:
            return (False,"The lower cut-off value was not correctly specified")
        if seq_co != floor(seq_co):
            return (False,"The cut-off value must be an integer")
        if seq_co < 1:
            return (False,"The sequence cut-off must be > 0")
        return (True,'')

        

