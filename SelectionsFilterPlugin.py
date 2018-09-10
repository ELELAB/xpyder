#!/usr/bin/env python
#  File: SelectionsFilterPlugin.py 
#
#  Copyright (C) 2011
#        Matteo Tiberti <matteo.tiberti@gmail.com>
#        Marco Pasi <mf.pasi@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
# "SelectionsFilterPlugin.py v0.2 (C) 2011 Matteo Tiberti, Marco Pasi"

from pymol import cmd
from cPlugin import FilterPlugin
import numpy as np
import os
from Tkinter import *
import Pmw

class SelectionsFilterPlugin(FilterPlugin):
    """
    Depends:  get_calphas, resi2i
    Provides: 
    """
    def __init__(self):                 # init must have no arguments
        FilterPlugin.__init__(self)
        # data
        # parameters
        self.single = True              # selection mode
        self.sel1=""                    # first selection 
        self.sel2=""                    # second selection

# --------------
    def get_filter(self, data, reference, *args, **kwargs):
        """
        Generate a filtering matrix of type=np.bool, to filter
        *data*, a matrix referring to the *reference* object in PyMOL.
        """
        blank = np.ones(data.shape, dtype=np.bool)

        if not self.is_activated:
            return blank
        
        self.sel1 = self.sel1f.getvalue()
        self.sel2 = self.sel2f.getvalue()

        if self.single:
            self.sel2 = self.sel1

        return self.filter(data)

# --------------
    def filter(self, data):
        """
        Filter correlations according to 1 or 2 selections.
        """
        out = np.zeros(data.shape, dtype=np.bool)

        cas1 = self.services.provideServiceByName("get_calphas")(self.sel1)
        if self.sel2 == self.sel1:
            cas2 = cas1
        else:
            cas2 = self.services.provideServiceByName("get_calphas")(self.sel2)
        

        for i in range(len(cas1)):
            for j in range(len(cas2)):
                # self.msg("1: %d: %s%s = %d"%(i,cas1[i].resi,cas1[i].chain,self.resi2i(cas1[i])))
                # self.msg("2: %d: %s%s = %d"%(i,cas2[j].resi,cas2[j].chain,self.resi2i(cas2[j])))
                idx1 = self.services.provideServiceByName("resi2i")(cas1[i])
                idx2 = self.services.provideServiceByName("resi2i")(cas2[j])
                # print "%s  %s" % (idx1,idx2)
                out[idx1][idx2] = True
                out[idx2][idx1] = True
        return out


# --------------
    def selectionValidation(self, s, mode='pmw', second=False):
        objs = self.services.provideServiceByName("get_objects_names")()+["all"]
        if second and self.single:
            if mode == 'pmw':
                return Pmw.PARTIAL
            else:
                return False
        if s in objs:
            if mode == 'pmw':
                return Pmw.OK
            return True
        if mode == 'pmw':
            return Pmw.PARTIAL
        return False

# --------------
    def selection_mode(self,var):
        self.single = var=="single"

# --------------
    def GUI(self, page):
        """
        Add the specific GUI elements for this plugin to the provided Tk widget.
        """
        group = Pmw.Group(page,tag_text='Selections')
        group.pack()
        frame = Frame(group.interior())
        frame.pack()

        self.modesel= Pmw.RadioSelect(frame,
            buttontype = 'radiobutton',
            orient = 'horizontal',
            labelpos = 'n',
            command = self.selection_mode,
            label_text = 'Selection mode',
            hull_borderwidth = 2,
            hull_relief = 'ridge',
        )
        self.modesel.pack(fill=X,side = 'top', expand = 1, padx = 10, pady = 10)
        self.modesel.add('single')
        self.modesel.add('double')
        self.modesel.setvalue(['double','single'][0+self.single]) # ;(
        self.sel1f = Pmw.EntryField(frame,
                                         labelpos='w',
                                         label_text='Selection 1: ',
                                         value=self.sel1,
                                         validate = {'validator':self.selectionValidation},
                     )
        self.sel2f = Pmw.EntryField(frame,
                                         labelpos='w',
                                         label_text='Selection 2: ',
                                         value=self.sel2,
                                         validate = {'validator':self.selectionValidation,
                                                     'second': True},
                     )
    
        self.sel1f.pack(fill='both',padx=4,pady=1,expand=0) # vertical
        self.sel2f.pack(fill='both',padx=4,pady=1,expand=0) # vertical
        # self.balloon.bind(self.sel1f, """Selection to plot correlations to.""")
        # self.balloon.bind(self.sel2f, """Selection to plot correlations to.""")       
        return group

 
# --------------
    def check(self):
        """
        Check wether the plugin is ready for a *get_filter* call (e.g. all GUI inputs
        are valid, all information is available, ecc.).

        Inputs are checked on change, and have default values: return True
        """
        ref    = self.services.provideServiceByName("get_reference")()
        refids = cmd.identify(ref,1)
        if not refids:
            return (False,"The reference set must have at least one atom")
        refobjs= list(set(zip(*refids)[0]))
        refset = set(zip(*refids)[1])

        sel1 = self.sel1f.getvalue()
        if not self.selectionValidation(sel1, mode='pymol'):
            return (False,"Selection 1 is not valid")
        
        sel1ids  = cmd.identify(sel1,1)
        if not sel1ids:
            return (False,"Selection 1 must have at least 1 atom")
        sel1objs = list(set(zip(*sel1ids)[0]))
        sel1set  = set(zip(*sel1ids)[1])
        
        if len(sel1objs) == 1:
            if sel1objs[0] not in refobjs or sel1set > refset:
                return (False, "Selection 1 must be part of the reference <%s>"%ref)
        else:
            return (False, "Selection 1 must be part of a single object")
        
        if not self.single:
            sel2 = self.sel2f.getvalue()
            if not self.selectionValidation(sel2, mode='pymol'):
                return (False,"Selection 2 is not valid")
            
            sel2ids = cmd.identify(sel2,1)
            if not sel2ids:
                return (False,"Selection 1 must have at least 1 atom")
            sel2objs = list(set(zip(*sel2ids)[0]))
            sel2set  = set(zip(*sel2ids)[1])
            
            if len(sel2objs) == 1:
                if sel2objs[0] not in refobjs or sel2set > refset:
                    return (False, "Selection 2 must be part of the reference <%s>"%ref)
            else:
                return (False, "Selection 2 must be part of a single object")
            
            if not ( len(sel1set & sel2set) == len(sel1set | sel2set) or len(sel1set & sel2set) == 0 ) :
                return (False, "Selections must be identical or non-overlapping")
        return (True,"")
