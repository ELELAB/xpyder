#!/usr/bin/env python

#  File: xPyder.py 
#
#  Copyright (C) 2008-2012 Marco Pasi, Tiberti Matteo, Arrigoni Alberto, Papaleo Elena <elena.papaleo@unimib.it>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#

import Pmw
from Tkinter import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.backends.backend_tkagg
import matplotlib.figure
    
# ------------------------------------------
class DataSummary():
    def __init__(self, data, parent=None, fg="orange", bg="#333333"):
        self.data = data
        self.parent = parent
        self.fg = fg
        self.bg = bg
        
        self.vars = [
            [StringVar() for i in (1,2,3)], # DpC
            [StringVar() for i in (1,2,3)], # DnC
            [StringVar()] ]                 # DaC

        self.zero()

        self.makeSummary()
        self._binder(parent)
        
    def _binder(self, a):
        a.bind("<ButtonRelease-1>", self.onclick, add=True)
        for child in a.winfo_children():
            self._binder(child)

    def update(self, data, posC=None, negC=None):
        if data is None:
            return True
        
        self.data = data
        self.setvars(posC, negC, data)
        return True

    def onclick(self, event):
        pass

    def setvars(self, posC, negC, data):
        self.vars[0][0].set("%5d"%np.sum(posC!=0))
        self.vars[0][1].set("%5.3f"%np.min(posC))
        self.vars[0][2].set("%5.3f"%np.max(posC))
        self.vars[1][0].set("%5d"%np.sum(negC!=0))
        self.vars[1][1].set("%5.3f"%np.min(negC))
        self.vars[1][2].set("%5.3f"%np.max(negC))
        self.vars[2][0].set("%5d"%data.shape[0])
        
    def zero(self):
        dummy = np.zeros((1,1))
        self.setvars(dummy,dummy,dummy)
    
    def makeSummary(self):
        h1 = Label(self.parent, text="Num")
        h2 = Label(self.parent, text="Min")
        h3 = Label(self.parent, text="Max")
        h1.grid(row=0, column=1)
        h2.grid(row=0, column=2)
        h3.grid(row=0, column=3)
        l1 = Label(self.parent, text="Values +")
        l1.grid(row=1, column=0)
        l2 = Label(self.parent, text="Values -")
        l2.grid(row=2, column=0)
        l3 = Label(self.parent, text="Residues")
        l3.grid(row=3, column=0)
        for r in range(len(self.vars)):
            for v in range(len(self.vars[r])):
                tmp = Label(self.parent, textvariable=self.vars[r][v], fg=self.fg, bg=self.bg)
                tmp.grid(row=r+1, column=v+1)

# ------------------------------------------
class DataSummaryH(DataSummary):
    def __init__(self, data, parent=None, root=None,
                 title="Histogram", fg="orange", bg="#eeeeee", xlabel="Matrix values", ylabel="Frequency",
                 normed=True):
        self.data = data
        self.parent = parent
        self.root = root
        self.title = title
        self.fg = fg
        self.bg = bg
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.normed = normed

        self.figure, self.plot = self.makeHistogram(window=False)
        self._binder(parent)
        
        
    def update(self, data, posC=None, negC=None):
        if data is None:
            return True
        #plt.close(self.figure)
        self.data = data
        self.plotHistogram(self.figure, tight=True)
        self.plot.show()
        return True

    def onclick(self, event):
        if self.data is None:
            return True
        ax = self.plotHistogram(self.makeHistogram(window=True), tight=True)
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        return True
    
    def plotHistogram(self, figure, tight=True):
        figure.clear()
        ax = figure.add_subplot(111, facecolor=self.bg)
        if self.data is not None:
            n, bins, patches = ax.hist(np.ravel(self.data), 50, normed=self.normed, facecolor=self.fg, alpha=0.75)
        else:
            return ax
        if tight:
            # print "hist: maxy:",max(n)," xrange:",bins[0],",",bins[-1]," int:",np.sum(n)
            ax.set_xlim(np.round(bins[0],2), np.round(bins[-1],2))
            ax.set_ylim(0.0, np.round(max(n),2))
            ax.spines['bottom'].set_position(('data',0))
            ax.spines['left'].set_smart_bounds(True)
            ax.spines['bottom'].set_smart_bounds(True)
            ax.tick_params(labelsize='xx-small')
            ax.yaxis.set_ticks_position('left')
            # ax.yaxis.set_tick_params(labelsize=0)  # hide y tick labels
            ax.xaxis.set_ticks_position('bottom')
            ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=6))
            figure.subplots_adjust(bottom=0.175)
            # plt.tight_layout() # for future versions
            ax.fmt_xdata = lambda x: "%5.2f"%x
            ax.fmt_ydata = lambda x: "%5.2f"%x
        return ax

    def makeHistogram(self, window=False):
        ######
        if window:
            dialog = Pmw.Dialog(
                self.root,
                buttons = ("Close plot",),
                title = self.title)
            plot = dialog.interior()
            dpi = 100
            figw= 600
            figh= 450
        else:
            plot = self.parent
            dpi = 100
            figw= 200.0
            figh= 130.0

        f = matplotlib.figure.Figure(figsize=(figw/dpi,figh/dpi), dpi=dpi, frameon=True, facecolor='#dbdbdb')
        f.set_clip_on(False)
        dataPlot = matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(f, master=plot)
        dataPlot.show()
        dataPlot.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1, padx=0, pady=0)
        dataPlot.get_tk_widget().config(selectborderwidth=0)

        if window:
            ## a toolbar
            # toolbar = xPTB( dataPlot, plot, width=figw)
            toolbar = matplotlib.backends.backend_tkagg.NavigationToolbar2TkAgg( dataPlot, plot)
            toolbar.update()
            toolbar.pack()
            ######
            return f
        return f, dataPlot
        
# ------------------------------------------
class DataSummaryHGraph(DataSummaryH):
    def __init__(self, data, parent=None, root=None, title="Histogram", fg="red", bg="#eeeeee"):
        DataSummaryH.__init__(self, data, parent=parent, root=root, title=title, fg=fg, bg=bg)
        
    def plotHistogram(self, figure, tight=True):
        figure.clear()
        ax = figure.add_subplot(111, facecolor=self.bg)
        if self.data is not None:
            n, bins, patches = ax.hist(np.ravel(self.data), bins=[ x+0.5 for x in range(0,max(self.data)+1) ], normed=False, facecolor=self.fg, alpha=0.75)
        else:
            n, bins, patches = ax.hist(0.2*np.random.randn(10000), 50, normed=True, facecolor=self.fg, alpha=0.75)
        if tight:
            # print "hist: maxy:",max(n)," xrange:",bins[0],",",bins[-1]," int:",np.sum(n)
            ax.set_xlim(np.round(bins[0],2), np.round(bins[-1],2))
            ax.set_ylim(0.0, np.round(max(n),2))
            ax.spines['bottom'].set_position(('data',0))
            ax.spines['left'].set_smart_bounds(True)
            ax.spines['bottom'].set_smart_bounds(True)
            ax.tick_params(labelsize='xx-small')
            ax.yaxis.set_ticks_position('left')
            # ax.yaxis.set_tick_params(labelsize=0)  # hide y tick labels
            ax.xaxis.set_ticks_position('bottom')
            ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True,steps=range(1,max(self.data)+1)))
            figure.subplots_adjust(bottom=0.175)
            # plt.tight_layout() # for future versions
            ax.fmt_xdata = lambda x: "%5.2f"%x
            ax.fmt_ydata = lambda x: "%5.2f"%x
        return ax

    def onclick(self, event):
        if self.data is None:
            return True
        ax = self.plotHistogram(self.makeHistogram(window=True), tight=True)
        ax.set_xlabel("Node degree")
        ax.set_ylabel("Number of nodes")
        return True

