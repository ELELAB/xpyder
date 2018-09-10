# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 16:12:41 2011

@author: teo
"""

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

class PosNegFilterPlugin(FilterPlugin):
    """
    Depends:  get_calphas, parse_dat
    Provides: 
    """
    def __init__(self):                 # init must have no arguments
        FilterPlugin.__init__(self)
        # data
        self.positive       = True          # cutoff distance
        # parameters
# --------------
# Sequence distance filter
    def filter(self, data, ispositive):
        """
        Sequence distance filter

        dist: sequence distance cut-off
        mode: 0/1. Switch filtering mode:
            0 filters residues WITHIN distance dist (abs(i-j)<dist), 
            1 filters residues OUTSIDE sequence distance dist (abs(i-j)>dist),
        """
        return self.services.provideServiceByName("sign_filter")(data, positive=ispositive)
# --------------
    def get_filter(self, data, reference, *args, **kwargs):
        """
        Generate a filtering matrix of type=np.bool, to filter
        *data*, a matrix referring to the *reference* object in PyMOL.
        """
        if not self.is_activated:
            return np.ones(data.shape, dtype=np.bool)

        self.positive = True if self.signf.getvalue() == "positive" else False

        return self.filter(data, self.positive)

# --------------
    def GUI(self, page):
        """
        Add the specific GUI elements for this plugin to the provided Tk widget.
        """
        G = Pmw.Group(page,tag_text='Sign')
        G.pack()
        F = Frame(G.interior())
        F.pack()

        # input file
        # distance cut-off
        # sign mode
        self.signf = Pmw.OptionMenu(F,
                                       labelpos = 'w',
                                       label_text = '',
                                       items = ["positive","negative"],
                                       initialitem = "positive"
                                       )        #self.balloon.bind(self.co_modef, """Sequence distance cutoff mode""")
        self.signf.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        return G


# --------------
    def check(self):
        """
        Check wether the plugin is ready for a *get_filter* call (e.g. all GUI inputs
        are valid, all information is available, ecc.).

        Inputs are checked on change, and have default values: return True
        """
        #chains=list(set([ i.chain for i in self.services.provideServiceByName("get_calphas")(reference) ]))
        return (True,"")

        

