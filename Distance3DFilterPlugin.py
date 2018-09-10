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

import tkFileDialog
import numpy as np
import os, sys
from Tkinter import *
import Pmw
import Tkinter


class Distance3DFilterPlugin(FilterPlugin):
    """
    Depends:  get_calphas, parse_dat
    Provides: distance_matrix
    
    The "get_dm" method will be provided using the capability name "distance_matrix".
    """
    def __init__(self):                 # init must have no arguments
        FilterPlugin.__init__(self)
        # data
        self.dm       = None            # the matrix
        # parameters
        self.fname = ""              # the matrix file
        self.cutoff_l = 0.4          # lower cutoff
        self.cutoff_h =1.0          # higher cutoff
        # parameters used to compute current dm, for caching
        self._fname = ""
        self._reference= None

# --------------
    def compute_dm(self, reference):
        """
        Compute distance matrix.  The *calphas* array is obtained using the
        "get_calphas" capability, provided by the cPlot class.
        """
        cas = self.services.provideServiceByName("get_calphas")(reference)
        coord = np.array([ca.coord for ca in cas]) # coordinates in sane format
        n = coord.shape[0]
        d = np.zeros((n,n,3), np.float)
        for i in range(3):
            d[:,:,i] = coord[:,i,None] - coord[:,i]
        return (np.sqrt(np.sum(d**2,2)))/10.0

# --------------
    def parse_dm(self, filen):
        """
        Parse a distance matrix.  The .dat file parsing (parse_dat) capability is
        provided by the cPlot class.
        """
        return np.array(self.services.provideServiceByName("parse_dat")(filen))
        
        
# --------------
    def get_dm(self, reference):
        """
        Get the Distance matrix, either from file, or compute it from reference.
        This method caches the matrix for further calls using the same parameters.
        """
        if self.fname:                       # always use file if available
            if self.fname != self._fname:
                # if cached version is not valid, re-compute
                # after saving the new parameters
                self._fname = self.fname
                self.dm = self.parse_dm(self.fname)
        else:
            if self._reference != reference:
                # if cached version is not valid, re-compute
                # after saving the new parameters
                self._reference = reference
                self.dm = self.compute_dm(reference)
        return self.dm

# --------------
# 3D distance filter
    def filter(self, dm, cutoff_l, cutoff_h):
        """
        Filter on the basis of spatial distance.

        dm: a distance matrix
        cutoff_l: distance lower limit
        cutoff_h: distance high limit
        
        """
        return (dm >= cutoff_l) & (dm < cutoff_h)

# --------------
    def get_filter(self, data, reference, *args, **kwargs):
        """
        Generate a filtering matrix of type=np.bool, to filter
        *data*, a matrix referring to the *reference* object in PyMOL.
        """
        if not self.is_activated:
            return np.ones(data.shape, dtype=np.bool)

        # retrieve parameters from gui
        self.fname = None
        if self.ff.valid():
            self.fname = self.ff.get()
        self.cutoff_l= self.colf.get()
        self.cutoff_h= self.cohf.get()

        self.dm = self.get_dm(reference) # get_dm handles caching

        return self.filter(self.dm, self.cutoff_l, self.cutoff_h)

# --------------
    def GUI(self, page):
        """
        Add the specific GUI elements for this plugin to the provided Tk widget.
        """
        distG = Pmw.Group(page,tag_text='Distance matrix')
        distG.pack()
        distF = Frame(distG.interior())
        distF.pack()

        # input file
        self.ff = Pmw.EntryField(distF,
                                 labelpos='w',
                              #   label_pyclass = FileDialogButtonClassFactory.get(lambda x: self.ff.setentry(x), filter=("*.dat")),
                                 validate = {'validator': self.validate_file,},
                                 value = self.fname,
                                 label_text = 'Choose Distance Matrix file:')
        self.ff.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        def open_distance_matrix():
            self.fname = tkFileDialog.askopenfilename()
            self.ff.setentry(self.fname)
        trial=Button(distF,command=open_distance_matrix,text="Choose distance Matrix file")
        trial.pack(side=LEFT)
        # distance cut-off
        self.colf = self.makeDoubleVar(distF, "Low cut-off (nm)", self.change_col,
                                      tooltip="""\
Only correlations between residues with
distance higher than this cutoff
will be plotted.""")
        self.colf.set(self.cutoff_l)

        self.cohf = self.makeDoubleVar(distF, "High cut-off (nm)", self.change_coh,
                                      tooltip="""\
Only correlations between residues with
distance lower than this cutoff
will be plotted.""")
        self.cohf.set(self.cutoff_h)
        return distG

#----
    def change_col(self, a, *args):
        if len(args)!=0: 
            a=args[0]
        val = float(self.colf.get()) + 0.05*float(a)
        self.colf.set(val)
        
    def change_coh(self, a, *args):
        if len(args)!=0: 
            a=args[0]
        val = float(self.cohf.get()) + 0.05*float(a)
        self.cohf.set(val)
        
    def validate_file(self, s, *args):
        if len(args)!=0: 
            a=args[0]
        if s == '':
            return Pmw.PARTIAL
        elif os.path.isfile(s):
            #self.msg("Distance Matrix set to %s"%s)
            return Pmw.OK
        elif os.path.exists(s):
            return Pmw.PARTIAL
        else:
            return Pmw.PARTIAL

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
            colf = self.colf.get()
            cohf = self.cohf.get()
        except:
            return (False,"The distance cut-off value was not correctly specified")
        if colf > 0.0 and cohf > 0.0:
            if self.ff.valid():
                return (True,"Distance Matrix read from file")
            return (True,"Distance Matrix calculated from pdb")
        return (False, "The distance cut-off must be > 0.0")
        
class FileDialogButtonClassFactory:
    def get(fn,filter='*'):
        """This returns a FileDialogButton class that will
        call the specified function with the resulting file.
        """
        class FileDialogButton(Tkinter.Button):
            # This is just an ordinary button with special colors.
            def __init__(self, master=None, cnf={}, **kw):
                '''when we get a file, we call fn(filename)'''
                self.fn = fn
                self.__toggle = 0
                apply(Tkinter.Button.__init__, (self, master, cnf), kw)
                self.configure(command=self.set)
            def set(self):
                fd = PmwFileDialog(self.master,filter=filter)
                fd.title('Please choose a file')
                n=fd.askfilename()
                if n is not None:
                    self.fn(n)
        return FileDialogButton
    get = staticmethod(get)
        

