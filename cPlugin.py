#!/usr/bin/env python
#  File: cPlugin.py 
#
#  Copyright (C) 2011 Marco Pasi <mf.pasi@gmail.com> 
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
# "cPlugin.py v0.1 (C) 2011 Marco Pasi"

from yapsy.IPlugin import IPlugin
import numpy as np
from Tkinter import *

# -------------------
class ServicePlugin(IPlugin):
    """
    A plugin that can provide services, and can use services
    provided by other plugins.
    """
    def __init__(self):                 # init must have no arguments
        IPlugin.__init__(self)
        self.services = None

    def setServices(self, services):
        """
        Provides the plugin with a reference to the Services class,
        where services are registered by all plugins.
        To be called by the PluginManager after loading all plugins.
        """
        self.services = services
    
# -------------------
class FilterPlugin(ServicePlugin):
    """
    A typical cPlot filter plugin; it must:
      - handle services;
      - be capable of providing cPlot with a GUI to allow
        user-defined parameters, to be added to the global GUI;
      - implement a *check()* function, that checks that
        all prerequisites for the calculation of the filter are
        satisfied;
      - implement a *get_filter()* function, that creates the
        filter that cPlot can use to filter the plot.
    """
    
    def __init__(self):
        ServicePlugin.__init__(self)
    
    def GUI(self, page):
        """
        Add the specific GUI elements for this plugin to the provided Tk widget.
        """
        pass

    def check(self):
        """
        Check wether the plugin is ready for a *get_filter* call (e.g. all GUI inputs
        are valid, all information is available, ecc.).
        """
        return False
        
    def get_filter(self, data, reference, *args, **kwargs):
        """
        Generate a filtering matrix of type=np.bool, to filter
        *data*, a matrix referring to the *reference* object in PyMOL.
        """
        if not self.is_activated:
            return np.ones(data.shape, type=np.bool)

        return np.zeros(data.shape, type=np.bool)

    def check(self, *args, **kwargs):
        if not self.is_activated: 
            return (True, '')
        else:
            return ('notready','')
    

cPlugin = FilterPlugin

