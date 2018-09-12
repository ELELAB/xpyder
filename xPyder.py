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

VERSION=1.0

### Installation-dependent variables and initializations
import sys

INSTALLDIR=r"INSTALLDIR_PLACEHOLDER"

sys.path.append(INSTALLDIR)

#########################
import os,re,sys
import Tkinter
try:
    import ttk
except:
    print "xPyder: Failed to load system ttk. Defaulting to embedded fallback."
    import ttk_fallback as ttk
from Tkinter import *
#from ttk import *
import tkMessageBox
import Pmw
from pymol import cmd,selector
from pymol.cmd import _feedback,fb_module,fb_mask,is_list,_cmd
from pymol.cgo import *
from pymol import stored
import numpy as np
import tkColorChooser
from pymol.vfont import plain
import tkFileDialog
import networkx as nx
import time as t
import matplotlib as mpl
mpl.use('TkAgg')                        # XXX
import matplotlib.pyplot as plt
import matplotlib.backends.backend_tkagg
import matplotlib.figure
from random import random
from operator import itemgetter, attrgetter
from collections import deque
import tkFont
from cPluginManager import ServicePluginManager
from cPlugin import FilterPlugin
from plotting import DataSummary, DataSummaryH, DataSummaryHGraph

# Python backward-compatibility...
try:
    True
except:
    True = 1
try:
    False
except:
    False = 0


def __init__(self):
    # cPlotTk(self)
    self.menuBar.addmenuitem('Plugin', 'command',
                             'Launch xPyder',
                             label='xPyder',
                            command = lambda s=self: runCPlotTk(s))

# proxy function to intercept exceptions
def runCPlotTk(s):
    try:
        cPlotTk(s)
    except:
        import traceback
        traceback.print_exc()

# set the defaults
defaults = {
    "DCCM": "corr.dat",
    "p_mode":0,
    "matrix_offset":0,
    # plugins
    "seq-cutoff": 0,
    "seq-cutoff-mode": 1,
    "corr-cutoff": 0,
    "corr-cutoff2": 1.0,
    # display
    "minsize": .1,
    "maxsize": .5,
    "diagonal": 0,
    "dminsize": .2,
    "dmaxsize": 1.0,
    "colorp" : 'red',
    "colorn" : 'blue',
    "widthmode":0,
    # chains
    "colorchains":"green",
    "chainscolmode":0,
    "excludeneighscontrol":0,
    "minval":0.0,
    "maxval":1.0,
    "highest_avg_only_control":False,
    "exclude_leaves_control":False
    }

# ------------------------------------------
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

# ------------------------------------------
# ------------------------------------------
class cPlot:
    #############################
    #
    # 1)  choose dccm file [ & cutoffs, output files]
    # 1b) extract correlations -> self.data
    # 2)  choose selection, colors, sizes
    # 2b) plot correlations
    #
    #############################
# ---------------
    def __init__(self,
                 DCCM    = defaults['DCCM'],
                 object = 'corrs',
                 matrix_offset = defaults['matrix_offset'],
                 # display
                 p_mode = defaults['p_mode'],
                 minsize = defaults['minsize'],
                 maxsize = defaults['maxsize'],
                 diagonal= defaults['diagonal'],
                 dminsize = defaults['dminsize'],
                 dmaxsize = defaults['dmaxsize'],
                 colorp = defaults['colorp'],
                 colorn = defaults['colorn'],
                 widthmode = defaults['widthmode'],
                 # chains
                 objectchains = 'corr_chains',
                 colorchains = defaults['colorchains'],
                 selchains = "",
                 chainscolmode = defaults['chainscolmode'],
                 excludeneighscontrol=defaults['excludeneighscontrol'],
                 minval = defaults['minval'],
                 maxval = defaults['maxval'],
                 highest_avg_only_control = defaults['highest_avg_only_control'],
                 exclude_leaves_control = defaults['exclude_leaves_control']
                 ):

        """
        sco_mode     0 [lower than c.o.] 1 [higher...]
        """

       # self.memory_tracker = tracker.SummaryTracker()

        self.DCCM    = DCCM
        self.DCCM2   = "corr2.dat"         # second file, for delta calculation
        self.deltaA  = "align.dat"         # alignment file, for delta calc.
        self.object  = object
        self.moffset = int(matrix_offset)


        self.minval = minval
        self.maxval = maxval

        self.p_mode = p_mode        
        
        self.objectchains = objectchains
        self.reference= "all"           # the reference object

        self.chains_d=3
        self.chains_w=3
        
        self.minsize = float(minsize)
        self.maxsize = float(maxsize)
        self.diagonal= False
        self.dminsize = float(dminsize)
        self.dmaxsize = float(dmaxsize)
        
        self.chainscolmode=chainscolmode # 0: random; 1: user choice

        self.widthmode = widthmode # 1: Absolute, 0: Relative
        
        self.colorp = colorp
        self.colorn = colorn
        self.colorchains = colorchains
        
        self.colorRGBp = cmd.get_color_tuple(cmd._interpret_color(cmd, str(self.colorp)))
        self.colorRGBn = cmd.get_color_tuple(cmd._interpret_color(cmd, str(self.colorn)))
        self.colorRGBchains = cmd.get_color_tuple(cmd._interpret_color(cmd, str(self.colorchains)))
        
        self.exclude_neighs_control=excludeneighscontrol
        self.highest_avg_only_control=highest_avg_only_control
        self.exclude_leaves_control=exclude_leaves_control


        


        # ---- PLUGIN STUFF
        import logging
        #logging.basicConfig(level=logging.DEBUG)
        # ----
        #plugin_dirs.append(os.path.split(__file__)[0])
        #plugin_dirs.append(os.environ["PYMOL_PATH"])
        plugin_dirs = [INSTALLDIR]
        #plugin_dirs.append(os.path.expanduser("~/devel/cplot/package"))
        self.plugin_manager = ServicePluginManager(
            categories_filter = {"Filters": FilterPlugin},
            directories_list = plugin_dirs,
            plugin_info_ext = "xpyder-plugin")
        # register locally provided services
        self.plugin_manager.registerService(self.get_calphas, "get_calphas", VERSION)
        self.plugin_manager.registerService(self.parse_dat, "parse_dat", VERSION)
        self.plugin_manager.registerService(self.get_objects_names, "get_objects_names", VERSION)
        self.plugin_manager.registerService(self.resi2i, "resi2i", VERSION) # XXX change service name!
        self.plugin_manager.registerService(self.get_reference, "get_reference", VERSION) # XXX change service name!
        self.plugin_manager.registerService(self.filter_sign, "sign_filter", VERSION) # XXX change service name!
        self.plugin_manager.registerService(self.get_idxs, "get_idxs", VERSION) # XXX change service name!
        self.plugin_manager.registerService(self.parse_align, "parse_align", VERSION) # XXX change service name!
        

        # collect the plugins
        self.plugin_manager.collectPlugins()
        # XXX activate all plugins
        
        for plugin_info in self.plugin_manager.getPluginsOfCategory("Filters"):
            plugin_info.plugin_object.is_activated = False
        #    plugin_info.plugin_object.GUI(page)
        
        #---- data
        self.data =  None         # matrix
        self.data2=  None         # second matrix
        self.delta = None         # the delta matrix
        self.chains= None         # the chains filter?
        self.pdbdict=None         # ca.chain,ca.resi -> index in DCCM
        self.alignment=None       # the parsed alignment, for delta
        ##
        
        self.selchains = ""
        self.cor_res_selection="correlated_residue"

        self.logfile="correlations.log"
        self.graph=nx.Graph()

        self.inv_delta = False

        #---- cache
        self._cache_idx = {}
        self._cache_calphas = {}
        
# ------------------------------------------
    def msg(self, m):
        sys.stderr.write("cPlot: %s\n"%m)

    def ask(self, msg):
        self.msg(msg+": Proceeding...")
        return True
        
    def warn(self, msg):
        self.msg(msg)
        
    def error(self, msg):
        self.msg(msg)

#######################################################################################
#  UTILITIES
#######################################################################################

# -------------- provided
    def get_idxs(self, sel):
        #TTT
        t=time.time()

        cmd.select("tmpgetcalphas", sel+" and name CA and symbol C")
        caids = cmd.identify("tmpgetcalphas ")
        cmd.delete("tmpgetcalphas")

        #print "idx: %f"%(time.time()-t)
        return caids

# -------------- provided
    def get_calphas(self, sel):

        idxs = self.get_idxs(sel)
        if sel in self._cache_idx.keys() and idxs == self._cache_idx[sel]:
            return self._cache_calphas[sel]
        
        #TTT
        t=time.time()
        
        cas = []
        resi= []
        model = cmd.get_model(sel)
        res_pos = model.get_residues()
        i=0
        for res in res_pos:
            i+=1
            for a in range(res[0], res[1]):
                atom = model.atom[a]
                if atom.name == "CA" and atom.symbol == "C":
                   if not atom in cas:
                       #resi.append(int(atom.resi))
                       cas.append(atom)
                   else:
                       pass
            
        
        self._cache_idx[sel] = idxs
        self._cache_calphas[sel] = cas
        return cas

# --------------- provided
    def get_objects_names(self):
        return cmd.get_names("all")
    
# --------------- provided
    def parse_dat(self, f):
        DAT_COMMENTS_RE=re.compile('^([^#]+).*?$')
        ret=[]
        lines=file(f).readlines()
        for line in lines:
            line=line.strip()
            if line.startswith("#") or len(line) == 0: continue
            match=DAT_COMMENTS_RE.match(line)
            if not match: continue
            tmp=map(float, match.group(1).split())
            ret.append(tmp)
        return ret

# -------------- provided
    def resi2i(self, ca):
        """
        Given an atom *ca*, returns its index in the DCCM (self.data)
        """
        return self.pdbdict["%s%s"%(ca.chain,ca.resi)]

# -------------- provided
    def get_reference(self):
        return self.referencef.getvalue()    

# --------------- provided
    def parse_align(self, f):
        """
        Parse an alignment file to obtain an array of
        aligned sequences in alignment order.
        
        RETURN: alignment, labels
        """
        ALIGN_RE=re.compile('^([^\s\:\.\*\#]+)\s+([^\s]+).*$')
        align=[]
        lines=[]
        al_fd=file(f, 'r')
        line="/*cycle starter*/"
        labels=[]

        # PARSE FILE LINES TO GET ALIGNMENT LABELS
        # AND ALIGNMENT LINES
        while len(line) > 0:
            line=al_fd.readline()
            match = ALIGN_RE.match(line)
            if not match: continue
            if not match.group(1) in labels:
                # add label to labels
                labels.append(match.group(1))
            # add line to lines
            lines.append(match)

        # INITIALISE ARRAY WITH SPACES FOR EACH LABEL
        for i in range(len(labels)):
            align.append("")

        # FILL ARRAY WITH ALIGNED SEQUENCES (w/gaps)
        for line in lines:
            # concatenate alignment lines
            align[labels.index(line.group(1))] += line.group(2)
            #print "%s: %s -> %s"%(line, line.group(1),line.group(2))

        return align, labels

# --------------
    def align_to_indices(self, align, labels):
        """ Obtain an array [[idx_seq1 idx_seq2],...] of matching indices of the 2 sequences """
        # align = [ SEQ1 SEQ2 ]
        # labels= [ LAB1 LAB2 ]
        # SEQn  = [A-Z-]+ (- == gap)
        ret = []
        gaps= [0, 0]                    # num. of gaps encountered
        class isgap(Exception): pass    # a dummy exception
        for i in range(len(align[0])):  # cycle alignment sequence
            # keep track of gaps
            try:
                for a in range(2):
                    if align[a][i] == '-':
                        gaps[a] += 1
                        raise isgap
            except isgap:
                continue
            # it's not a gap: remember
            tmp = []
            for a in range(2):
                tmp.append(i-gaps[a])
            ret.append(tmp)
            
        return np.array(ret)
        
# --------------
    def sparse(self, data=None):
        """
        Return a sparse representation of data.
        """
        if data is None:
            data = self.data
        cas = self.get_calphas(self.reference)
        out = []
        for i in range(len(cas)):
            for j in range(i, len(cas)):
                if not data[i,j] == 0.0:
                    out.append((i,j,data[i,j]))
        return out

    def calc_graph(self, m):
        data = m.copy()
        filters = []
        cas = self.get_calphas(self.reference)
        if not data.shape[0] == data.shape[1]:
            self.error("Matrix must be square, it is %s.  Load aborted."%(",".join(["%d"%x for x in data.shape])))
            self.reset_net_graph()
		        # 2) check correspondence with reference object
        if not data.shape[0] == len(cas):
            self.error("Matrix size (%d) must match reference object %s (%d)."%(data.shape[0], self.reference, len(cas)))
            self.reset_net_graph()

        self.update_filters_status()
        for filter_plugin in self.plugin_manager.getPluginsOfCategory("Filters"):
            if filter_plugin.is_activated:
                print "Using filter: %s" % filter_plugin.name
                filters.append(filter_plugin.plugin_object.get_filter(data, self.reference))
        filters.append(self.filter_sign(data, positive=True))

        filtered_pm = self.apply_filters(data, filters)

        filtered_pm_checkm = filtered_pm.copy()
        for i in range(filtered_pm_checkm.shape[0]):
            filtered_pm_checkm[i,i] = 0.0
        if not filtered_pm_checkm.nonzero()[0].size:
            self.error("No matrix elements were left after filtering")
            self.reset_net_graph()
            
        self.graphmcache = filtered_pm.copy()

        self.graph = nx.Graph()

        for i in range(filtered_pm.shape[0]):
            for j in range(filtered_pm.shape[0]):
                if abs(filtered_pm[i][j]) != 0.0 and i!=j:
                    label1 = cas[i].chain + cas[i].resi
                    label2 = cas[j].chain + cas[j].resi
                    self.graph.add_edge(label1,label2,weight=filtered_pm[i][j])
                        
        components = [ self.graph.subgraph(c).copy() for c in nx.algorithms.connected_components(self.graph) ]

        connectivity_deg = {}
        deg = dict(self.graph.degree()).values()
        for i in range(m.shape[0]-1):
            if deg.count(i) != 0:
                connectivity_deg[i] = deg.count(i)
                        
        nodesn = len(self.graph.nodes())
        edgesn = len(self.graph.edges())
        
        return (nodesn,edgesn,components,deg,connectivity_deg)
                
    def select_hubs(self):
        if not self.graph or not self.graphdata:
            self.error("Please calculate the relationships graph first!")
            return
        deg = dict(self.graph.degree())
        try: 
            val = int(self.hubs_sel.getvalue())
        except:
            self.error("Please specify a correct value for minimum degree")
            return
        residues = [ i for i in deg.keys() if deg[i] >= val ]
        self.select_residues("hubs_deg_"+str(val),residues)
        
    def select_cluster(self):
        if not self.graph or not self.graphdata:
            self.error("Please calculate the relationships graph first!")
            return
        try:
            cln = int(self.n_clustersf.getcurselection()[0])# XXX GET CLUSTER
        except:
            self.error("Please specify a correct cluster number")
            return
        self.select_residues("cluster_residues_%d" % cln, self.graphdata[2][cln-1].nodes())
        
    def plot_cluster(self):
        if not self.graph or not self.graphdata:
            self.error("Please calculate the relationships graph first!")
            return
        try:
            cln = int(self.n_clustersf.getcurselection()[0])# XXX GET CLUSTER
        except:
            self.error("Please specify a correct cluster number")
            return
        self.plot_graph(self.graphdata[2][cln-1], "cluster_%d" % cln, self.graphmcache)

        
        
    def plot_graph(self, graph, objname, m, dofilter=False):
        if not self.graph or not self.graphdata:
            self.error("Please calculate the relationships graph first!")
            return
        try:
            cln = int(self.n_clustersf.getcurselection()[0])# XXX GET CLUSTER
        except:
            self.error("Please specify a correct cluster number")
            return
        self.reference = self.referencef.getvalue()
        cas = self.get_calphas(self.reference)
        data = m.copy()
        if not data.shape[0] == len(cas):
            self.error("Matrix size (%d) must match reference object %s (%d)."%(data.shape[0], self.reference, len(cas)))
            return

        self.minsize = float(self.minsizef.get())
        self.maxsize = float(self.maxsizef.get())
        self.dminsize = float(self.dminsizef.get())
        self.dmaxsize = float(self.dmaxsizef.get())
        self.minval  = float(self.minvalf.get())
        self.maxval  = float(self.maxvalf.get())
        self.object = self.objecti.getvalue()
        self.dolog = self.log_file_st.get()
        self.doimg = self.img_file_st.get()
        if self.dolog:
            self.log_file_name = self.log_file_var.get()
        if self.doimg:
            self.img_file_name = self.img_file_var.get()

        object = objname

        rgbp = self.colorRGBp
        # what to plot?
        if self.data is None:
            self.error("Please load a matrix first.")
            return False
        if self.graph is None or self.graphdata[0] is None:
            self.error("Please calculate the graph first.")
            return False
        maxval = self.maxval
        minval = self.minval

        if self.widthmode == 1:         # absolute
            above = zip(*np.where(abs(data)>maxval))
            below = zip(*np.where(abs(data)<minval))
            if len(below) > 0:
                self.msg("WARNING: some matrix values were below the specified minimum. Will be reduced to minimum value.")
            if len(above) > 0:
                self.msg("WARNING: some matrix values were above the specified maximum. Will be reduced to maximum value.")
            if self.minval != 0.0:
                for i in below:
                    if data[i[0],i[1]] > 0:
                        data[i[0],i[1]] = self.minval
                    else:
                        data[i[0],i[1]] = - self.minval
                for i in above:
                    if data[i[0],i[1]] > 0:
                        data[i[0],i[1]] = self.maxval
                    else:
                        data[i[0],i[1]] = - self.maxval
            else:                    
                for i in above:
                    if data[i[0],i[1]] > 0.0:
                        data[i[0],i[1]] = self.maxval
                    else:
                        data[i[0],i[1]] = - self.maxval

        # ------ compute filter list
        if dofilter:
            filters = []
            for filter_plugin in self.plugin_manager.getPluginsOfCategory("Filters"):
                if filter_plugin.is_activated:
                    # print "Using filter: %s" % filter_plugin.name
                    filters.append(filter_plugin.plugin_object.get_filter(data, self.reference))
            filters.append(self.filter_sign(data, positive=True))
    
            filtered_pm = self.apply_filters(data, filters)
        else:
            filtered_pm = data
        
        ncorrs = len(graph.edges())
        if ncorrs > 1000:
            if not self.ask("You're trying to plot %d relationships: do you want to continue?"%ncorrs):
                return
        
        if self.widthmode == 1: # Absolute
            mX=minval
            MX=maxval
        elif self.widthmode == 0: # Relative to cutoffs
            filtered_pm_range = filtered_pm.copy()
            for i in range(filtered_pm_range.shape[0]):
                filtered_pm_range[i,i] = 0.0
            clean_data = abs(filtered_pm_range[abs(filtered_pm_range)>0.0])
            MX = max(clean_data)
            mX = min(clean_data)
            if mX >= MX: 
                self.msg("ERROR: switching to absolute mode.\n")
                mX = minval
                MX = maxval
        sparse_m = []
        resis = {}
        for i in graph.nodes():
            resis[i] = [k for v,k in self.pdbdict.iteritems() if v == i][0]
        for i in graph.edges(data=True):
            sparse_m.append((resis[i[0]],resis[i[1]],i[2]['weight']))
        
        cgOb = []
        cgOb.extend(self.makeCGO4(cas, rgbp, sparse_m, mX, MX, self.diagonal))

        cmd.load_cgo(cgOb,object)
                
        if self.doimg:
            bidiG = nx.Graph()
            for i in range(filtered_pm.shape[0]):
                for j in range(filtered_pm.shape[0]):
                    if abs(filtered_pm[i][j]) > 0.0 and i!=j:
                        label1 = cas[i].chain + cas[i].resi 
                        label2 = cas[j].chain + cas[j].resi 
                        bidiG.add_edge(label1,label2,weight=filtered_pm[i][j]*5)
            self.save_graph_pic(bidiG, self.img_file_name, mX, MX, dpi=600)
            
        
        if self.dolog:
            self.save_log(self.sparse(filtered_pm), self.log_file_name)

    def get_clusters(self, cln):
        return nx.algorithms.components.connected_components(self.graph)
        
    def dfs_traversal(self,G,source,end,maxl):
        out = []
        stack = [ [source] ]
        while stack:
            curr = stack.pop()
            for child in G[curr[-1]]:
                if child not in curr:
                    if child == end:
                        out.append(curr+[child])
                        if len(out) == 1000:
                            if not self.ask("1000 paths have been identified already. Calculations might be long. Do you want to continue?"):
                                return out
                    elif len(curr) < maxl:
                        stack.append(curr+[child])
        return out

    def get_all_paths(self):
        if self.data is None:
            self.error("Please load a matrix first.")
            return
        if self.graphdata[0] is None:
            self.error("Please calculate the graph first.")
            return
        res1 = self.sp_res1f.getcurselection()[0]
        res2 = self.sp_res2f.getcurselection()[0]
        try:
            maxl = int(self.sp_maxhops.getvalue())
        except:
            self.error("You must specify an integer value for the maximum path length")
            return
        if maxl < 0:
            self.error("The maximum number of paths must be positive.")
            return
        if res1 == res2:
            self.error("Resiues 1 and 2 coincide.")
            return
        try:
            minl = nx.algorithms.shortest_paths.shortest_path_length(self.graph,res1,res2)            
        except nx.NetworkXNoPath:
            self.msg("No paths exist between the selected pair of residues")
            return
        if minl > maxl:            
            self.msg("No paths of length <= %d exist between the selected pair of residues; minimum path length is %d" % (maxl,minl))
            return
        try:
            self.graphpaths = ((res1,res2),self.dfs_traversal(self.graph, res1, res2, maxl))
        except:
            self.error("There was an error calculating the paths.")
            self.graphpaths = None
            raise

        map(self.table.delete, self.table.get_children('')) #`delete all lines

        pathlens = [len(i)-1 for i in self.graphpaths[1]]
        pathcumw = []
        pathavgw = []
        for i in self.graphpaths[1]:
            thiscumw = 0
            for j in range(len(i)-1):
                thiscumw += self.graph[ i[j] ][ i[j+1] ]['weight']
            pathcumw.append(thiscumw)
            pathavgw.append(thiscumw/(len(i)-1))
        tabledata = zip(map(lambda x: str(int(x+1)), range(len(self.graphpaths[1]))), 
                        map(lambda x: "%d"%x,pathlens), 
                        map(lambda x: "%3.2f"%x,pathcumw), 
                        map(lambda x: "%3.2f"%x,pathavgw) )
        for i in tabledata:
            self.table.insert('', 'end', values=i)
            # adjust columns lenghts if necessary
        for item in tabledata:
            for indx, val in enumerate(item):
                ilen = tkFont.Font().measure(val)
                if self.table.column(self.table_columns[indx], width=None) < ilen:
                    self.table.column(self.table_columns[indx], width=ilen)        

        
    def select_sp(self):
        if not self.graph or not self.graphdata:
            self.error("Please calculate the relationships graph first!")
            return
        for i in self.table.selection():
            spn = self.table.item(i)['values'][0]
            self.select_residues("sp_"+self.graphpaths[0][0]+"_"+self.graphpaths[0][1]+"_"+str(spn)+"_residues", self.graphpaths[1][spn-1])

    def plot_sp(self):
        if not self.graph or not self.graphdata:
            self.error("Please calculate the relationships graph first!")
            return
        spns = [ self.table.item(i)['values'][0] for i in self.table.selection() ]
        for spn in spns:
            sparse = []
            for i in range(len(self.graphpaths[1][spn-1])-1):
                n1 = self.graphpaths[1][spn-1][i]
                n2 = self.graphpaths[1][spn-1][i+1]
                sparse.append([n1, n2, self.graph.get_edge_data(n1, n2, 'weight')])
            tmpgraph = nx.Graph(sparse)
            #print sparse
            self.plot_graph(tmpgraph, "sp_"+self.graphpaths[0][0]+"_"+self.graphpaths[0][1]+"_"+str(spn)+"_residues", self.graphmcache)
        
    def get_mst(self):
        return nx.algorithm.minimum_spanning_edges(self.graph,data = True)

    def find_chained_corrs(self,root,d,w,G):
        Gout = nx.Graph()
        q = deque([[root]])
        out = []
        visited=[root]
        if self.exclude_neighs_control == 0:
            while len(q) > 0:
                tmp_path = q.popleft()
                last_node = tmp_path[-1]
                next_nodes = list(G.neighbors(last_node))
                next_nodes.sort(key = lambda x: G.get_edge_data(last_node, x, 'weight')['weight'], reverse=True)
                next_nodes_filtered = [i for i in next_nodes if i not in list(set(tmp_path+visited))]
                if len(next_nodes_filtered) == 0:
                    out.append(tmp_path)
                for link_node in next_nodes_filtered[:w]:
                    if len(tmp_path) <= d: # -1 because for k nodes there are k-1 edges in a walk, +1 because we're gonna add 1 if this passes. Total 0
                        q.append(tmp_path + [link_node])
                        visited.append(link_node)
                    else:
                        out.append(tmp_path)
        else:
            while len(q) > 0:
                tmp_path = q.popleft()
                last_node = tmp_path[-1]
                next_nodes = list(G.neighbors(last_node))
                next_nodes.sort(key = lambda x: G.get_edge_data(last_node, x, 'weight')['weight'], reverse=True)
                if len(tmp_path) >= 2:
                    next_nodes_filtered = [i for i in next_nodes if i not in list(set(tmp_path+visited+G.neighbors(tmp_path[-2])))]
                else:
                    next_nodes_filtered = [i for i in next_nodes if i not in list(set(tmp_path+visited))]
                if len(next_nodes_filtered) == 0:
                    out.append(tmp_path)
                for link_node in next_nodes_filtered[:w]:
                    if len(tmp_path) <= d: # -1 because for k nodes there are k-1 edges in a walk, +1 because we're gonna add 1 if this passes. Total 0
                        q.append(tmp_path + [link_node])
                        visited.append(link_node)
                    else:
                        out.append(tmp_path)
        for i in out:
            for j in range(len(i)-1):
                #print "%s %s", i[j],i[j+1]
                Gout.add_edge(i[j],i[j+1],weight=G[i[j]][i[j+1]]['weight'])
        return Gout

        


#######################################################################################
#  FILTERS
#######################################################################################

# --------------
# correlation sign filter
    def filter_sign(self, data, positive=True):
        """
        data: the DCCM (np.array)
        positive: True:  filter only positive correlations (and 0)
              False: filter only negative correlations
                  
        Filter on the basis of the sign of the correlation.
        
        """
        if positive:
            return data >= 0
        else:
            return data < 0

# --------------
    def apply_filters(self, data, filter_list):
        """
        Apply all the filters in *filter_list* to *data*, and return the result.
        """
        out = np.zeros(data.shape)
        if not filter_list: 
            return data
        out[:] = data                   # copy data
        total_filter = np.logical_and.reduce(np.array(filter_list), 0)
        # XXX:p check whether sizes match?
        out[np.logical_not(total_filter)] = 0.0
        #### Debug
        # self.msg(str(np.array(filter_list).shape))
        # self.msg(str(total_filter))
        # self.msg(str(data))
        # self.msg(str(out))
        ####
        return out

# ------------------------------------------
#  Main data-handling and plotting functions
    def doDCCM(self):
        """
        Parse the DCCM file
        """
        # check if file exists
        if self.reference not in cmd.get_names("all")+["all"]:
            self.error("The specified reference object is not valid.")
            return False
        DCCM = self.DCCM
        if DCCM == "":
            self.error('Please specify a valid matrix file name.')
            return False
        elif not os.path.isfile(DCCM):
            self.error("Matrix file not found.")
            return False
        # parse matrix
        self.msg("Parsing matrix file %s ..."%DCCM)
        data = None
        try: 
            data = self.parse_dat(DCCM)
        except:
            self.error("Couldn't read matrix from matrix file.")
            return False

        #---- the DCCM
        data = np.array(data)
        # consistency checks:
        # 1) check square
        if not data.shape[0] == data.shape[1]:
            self.error("Matrix must be square, it is %s.  Load aborted."%(",".join(["%d"%x for x in data.shape])))
            return False
        # 2) check correspondence with reference object
        cas = self.get_calphas(self.reference)
        if not data.shape[0] == len(cas):
            self.error("Matrix size (%d) must match reference object %s (%d)."%(data.shape[0], self.reference, len(cas)))
            return False
            
        #if self.data[self.data<self.minval] or self.data[self.data>self.maxval]:
        #    self.error("Matrix values must be within defined boundaries.")
        #       return False
        
        #-- apply matrix offset
        if self.moffset > 0:
            if np.any(data[:,-self.moffset:data.shape[0]]) != 0 or np.any(data[-self.moffset:data.shape[0],:]) != 0:
                self.warn("Applying matrix offset is removing non-zero matrix values!")
            offset_data = np.zeros(data.shape)
            offset_data[self.moffset:, self.moffset:] = data[:-self.moffset, :-self.moffset]

            data = offset_data

        #-- it's ok
        self.data = data
        
        #---- other 
        self.chains= None               # the chains filter?
        self.pdbdict= dict(zip(*[["%s%s"%(ca.chain,ca.resi) for ca in cas], range(self.data.shape[0])])) # atom -> index of DCCM

      
  # self.msg(self.pdbdict)
        # XXX:p postpone graph creation to doplot()
        # self.graph = nx.generators.classic.complete_graph(self.data.shape[0])
        ####
        # edges = []
        # for i in range(len(self.data)):
        #     self.graph.add_node(i+1)
        #     for k in range(len(self.data[i])):
        #         edges.append((i,j,self.data[i][k]))
        # self.graph.add_weighted_edges_from(edges)
        ####

        npos, nneg = self.DCCM1_params(self.data)
        self.msg("Found %d positive values and %d negative values for %d residues"%(npos,nneg,self.data.shape[0]))
        return True

# --------------
    def doDCCM2(self):
        """
        Parse a second matrix, check its consistency with the first one
        """
        # XXX:p consider moving parseing and consistency checks to a function
        # check if file exists
        DCCM = self.DCCM2
        if DCCM == "":
            self.error('Please specify a valid matrix file name.')
            return False
        elif not os.path.isfile(DCCM):
            self.error("Matrix file not found.")
            return False
        # parse matrix
        self.msg("Parsing matrix file %s ..."%DCCM)
        data = None
        try: 
            data = self.parse_dat(DCCM)
        except:
            self.error("Couldn't read matrix from matrix file.")
            return False

        #---- the DCCM
        data2 = np.array(data)
        # consistency checks:
        # 1) check square
        if not data2.shape[0] == data2.shape[1]:
            self.error("Matrix must be square, it is (%s) instead.  Load aborted."%(",".join(["%d"%x for x in data2.shape])))
            return False

        # apply matrix offset
        if self.moffset > 0:
            if np.any(data2[:,-self.moffset:data2.shape[0]]) != 0 or np.any(data2[-self.moffset:data2.shape[0],:]) != 0:
                self.warn("Applying matrix offset is removing non-zero matrix values!")
            offset_data = np.zeros(data2.shape)
            offset_data[self.moffset:, self.moffset:] = data2[:-self.moffset, :-self.moffset]

            data = offset_data

        self.data2 = data2
        npos, nneg = self.DCCM2_params(self.data2, foutp=None, foutn=None)
        self.msg("Second matrix: Found %d positive values and %d negative values for %d residues"%(npos,nneg,self.data2.shape[0]))
        return True

# ---------------        
    def dodelta(self):
        """Compute a deltaDCCM, and store it in self.delta"""
        
        if self.data is None:
            self.error("Please load a matrix first.")
            return False
        
        if self.data2 is None:
            self.error("Please load a second matrix first.")
            return False

        if self.alignment is not None:  # if alignment is specified
            # check that self.alignment is applicable to data & data2
            # XXX first sequence MUST be reference
            lens = np.max(self.alignment, 0)
            if lens[0] > self.data.shape[0]: # is this redundant, since we've checked reference in doAlignment?
                self.error("The first sequence in the alignment is larger than the first loaded matrix: %d elements instead of %d"%(
                    lens[0], self.data.shape[0]))
                return False
            if lens[1] > self.data2.shape[0]:
                self.error("The second sequence in the alignment is larger than the second loaded matrix: %d elements instead of %d"%(
                    lens[1], self.data2.shape[0]))
                return False

            idx1 = self.alignment[:,0]
            idx2 = self.alignment[:,1]
            self.delta = np.zeros(self.data.shape) # delta matrix MUST be the size of reference!
            self.delta[np.ix_(idx1,idx1)] = self.data[idx1][:,idx1] - self.data2[idx2][:,idx2]
        else:
            if not np.all(self.data.shape == self.data2.shape):
                self.error("The size of the second matrix doesn't match the first: (%s) instead of (%s).  An alignment is required!"%(
                    ",".join(["%d"%x for x in self.data2.shape]),
                    ",".join(["%d"%x for x in self.data.shape])))
                return False
            self.delta = self.data - self.data2
            
        if self.inv_delta:
            self.delta = -self.delta

        np.savetxt("xPyder_delta.dat", self.delta, fmt="%.2f")
                    
        npos, nneg = self.DCCMd_params(self.delta, foutp=None, foutn=None)
        self.msg("Delta: Calculated %d positive values and %d negative values for %d residues"%(npos,nneg,self.delta.shape[0]))
        return True


# ---------------        
    def doAlignment(self):
        """ Parse the alignment and obtain indices for delta matrices matching """
        
        # Check that alignment has been loaded
        if self.deltaA is None:
            self.error("Please load an alignment first.")
            return False

        align, labels = self.parse_align(self.deltaA)

        # Check that alignment contains only 2 sequences
        if len(align) > 2:
            self.warn("Alignment (%s) contains more than 2 sequences (%d)! Using the first two (labels %s and %s)..."%(
                self.deltaA, len(align), labels[0], labels[1]))
            align = align[:2]
            labels= labels[:2]
            # return False

        # XXX first sequence MUST be reference
        # Check that len1 refers to reference
        len1 = len(align[0])-align[0].count("-")
        len2 = len(align[1])-align[1].count("-")
        reflen = len(self.get_calphas(self.reference))
        if reflen != len1:
            self.error("First sequence alignment (label %s) does not match the length of the reference (%s): %d instead of %d."%(
                labels[0], len1, reflen))
            return False

        # XXX Check that seq1 corresponds to reference (aminoacidic sequence)
        # ...

        self.alignment = self.align_to_indices(align, labels) # np.array( [ [ idx_seq1, idx_seq2 ], ... ] )
        common = self.alignment.shape[0]
        self.msg("Alignment: Found %d aligned residues, skipped %d from first, %d from second"%(
            common,len1-common,len2-common))
        return True
        
# ---------------        
    def DCCM_params(self,
                    data,
                    foutp, foutn):
        """
        Write a summary of the loaded data in the provided structures, and
        optionally write out + and - matrices.
        """
        if data is None:
            return 0,0
            
        posC =  self.apply_filters(data, (self.filter_sign(data, positive=True),) )
        negC =  self.apply_filters(data, (self.filter_sign(data, positive=False),) )
        data = self.data[:]

        if foutp:
            posCmfh=open(foutp,'w')
            for i in posC:
                for j in i:
                    print >> posCmfh, j,
                print >> posCmfh, "\n"
            posCmfh.close()
                
        if foutn:
            negCmfh=open(foutn,'w')
            for i in negC:
                for j in i:
                    print >> negCmfh, j,
                print >> negCmfh, "\n"
            negCmfh.close()

        return np.sum(posC!=0), np.sum(negC!=0)

# ---------------        
    def DCCM1_params(self, data, foutp=None, foutn=None):
        "Dummy method to call original DCCM_params"
        return self.DCCM_params(data, foutp, foutn)
    
# ---------------        
    def DCCM2_params(self, data, foutp=None, foutn=None):
        "Dummy method to call original DCCM_params"
        return self.DCCM_params(data, foutp, foutn)

# ---------------        
    def DCCMd_params(self, data, foutp=None, foutn=None):
        "Dummy method to call original DCCM_params"
        return self.DCCM_params(data, foutp, foutn)

# ---------------    


    def select_sparse(self, selname, tuples, reference=0, divrange=50):

        if reference == 0:
            reference = self.reference

        selestring = reference + " and ( none "    

        if len(tuples) < 1:
            selestring += ")"
            cmd.select(selname, selestring)            
            return
	print "sostrunz"
        cmd.select(selname,"none")
	print "cicciaao"
        cas = self.get_calphas(self.reference)
        tuplecas=[]
        for i in tuples:
            tuplecas.append(cas[i[0]])
            tuplecas.append(cas[i[1]])
            tuplechains = [i.chain for i in tuplecas]
            tupleresi = [i.resi for i in tuplecas]
            allchains = list(set(tuplechains))
            chainresi={}
        print "tc", tuplecas
        for i in allchains:
            chainresi[i] = []
        for i in range(len(tuplecas)):
            chainresi[tuplechains[i]].append(tupleresi[i])
	print "cr", chainresi
        for i in allchains:
            chainresi[i] = list(set(chainresi[i]))
            divn = len(chainresi[i])/divrange
            divrest = len(chainresi[i]) % divrange
            divisions = []
            for j in range(divn):
                divisions.append(chainresi[i][j*divrange:(j+1)*divrange] )
            if divrest > 0:
                divisions.append(chainresi[i][-divrest:])
            print divisions
            for j in divisions:
                thistmpsel = selestring
                thistmpsel += " or ( chain %s and resi %s ) " % (i, ",".join(j))
                thistmpsel += " ) "
		print "selecting %s" % thistmpsel
                cmd.select("thistmpsel",thistmpsel)
		print "fatto1"
		print "selecting %s" % selname+" or thistmpsel "
                cmd.select(selname, selname+" or thistmpsel ")
		print "fatto2"
            cmd.delete("thistmpsel")
            
            
    def select_residues(self, selname, reslist, reference=0, divrange=50):

        if reference == 0:
            reference = self.reference

        selestring = reference + " and ( None "    

        if len(reslist) < 1:
            selestring += ")"
            cmd.select(selname, selestring)            
            return
        cmd.select(selname,None)
        #cas = self.get_calphas(self.reference)

        tuplecas = reslist
        for resi2ii in reslist:
            #tuplecas.append(cas[i])
            tuplechains = [i[0] for i in tuplecas]
            tupleresi = [i[1:] for i in tuplecas]
            allchains = list(set(tuplechains))
            chainresi={}
        for i in allchains:
            chainresi[i] = []
        for i in range(len(tuplecas)):
            chainresi[tuplechains[i]].append(tupleresi[i])
        for i in allchains:
            chainresi[i] = list(set(chainresi[i]))
            divn = len(chainresi[i])/divrange
            divrest = len(chainresi[i]) % divrange
            divisions = []
            for j in range(divn):
                divisions.append(chainresi[i][j*divrange:(j+1)*divrange] )
            if divrest > 0:
                divisions.append(chainresi[i][-divrest:])
            for j in divisions:
                thistmpsel = selestring
                thistmpsel += " or ( chain %s and resi %s ) " % (i, ",".join(j))
                thistmpsel += " ) "
                cmd.select("thistmpsel",thistmpsel)
                cmd.select(selname,selname+" or thistmpsel ")
            cmd.delete("thistmpsel")
            

    def doplot(self):
        """
        do the actual plotting
        """
        # get object name
        object = self.object

        # get display parameters
        rgbp = self.colorRGBp
        rgbn = self.colorRGBn
        self.msg("Getting alpha-carbons")
        cas = self.get_calphas(self.reference)

        self.msg("Found %d alpha-carbons in selection"%len(cas))

        # what to plot?
        if self.p_mode == 0:
            if self.data is None:
                self.error("Please load a matrix first.")
                return False
            data = self.data.copy()
            maxval = self.maxval
            minval = self.minval
        elif self.p_mode == 2:
            if self.delta is None:
                self.error("Please compute delta matrix first.")
                return False
            data = self.delta.copy()
            maxval = 2*self.maxval
            minval = 0.0
        if not data.shape[0] == data.shape[1]:
            self.error("Matrix must be square, it is %s.  Load aborted."%(",".join(["%d"%x for x in data.shape])))
            return False
	        # 2) check correspondence with reference object
        if not data.shape[0] == len(cas):
            self.error("Matrix size (%d) must match reference object %s (%d)."%(data.shape[0], self.reference, len(cas)))
            return False

        if self.widthmode == 1:         # absolute
            above = zip(*np.where(abs(data)>maxval))
            below = zip(*np.where(abs(data)<minval))
            if len(below) > 0:
                self.msg("WARNING: some matrix values were below the specified minimum. Will be reduced to minimum value.")
            if len(above) > 0:
                self.msg("WARNING: some matrix values were above the specified maximum. Will be reduced to maximum value.")
            if self.minval != 0.0:
                for i in below:
                    if data[i[0],i[1]] > 0:
                        data[i[0],i[1]] = self.minval
                    else:
                        data[i[0],i[1]] = - self.minval
                for i in above:
                    if data[i[0],i[1]] > 0:
                        data[i[0],i[1]] = self.maxval
                    else:
                        data[i[0],i[1]] = - self.maxval
            else:                    
                for i in above:
                    if data[i[0],i[1]] > 0.0:
                        data[i[0],i[1]] = self.maxval
                    else:
                        data[i[0],i[1]] = - self.maxval

        # ------ compute filter list
        filters = []
        filter_names = []
        for filter_plugin in self.plugin_manager.getPluginsOfCategory("Filters"):
            if filter_plugin.is_activated:
                # print "Using filter: %s" % filter_plugin.name
                filters.append(filter_plugin.plugin_object.get_filter(data, self.reference))
                filter_names.append(filter_plugin.name)

        ### DEBUG
        # for fi,f in enumerate(filters):
        #     f_s = self.sparse(self.apply_filters(data, [f]))
        #     self.msg("SINGLE FILTER:   %2d (%15s) %5d"%(fi, filter_names[fi], len(f_s)))
        #     f_s = self.sparse(self.apply_filters(data, filters[:fi+1]))
        #     self.msg("CUMLTV FILTER: 0-%2d (%15s) %5d"%(fi, filter_names[fi], len(f_s)))
        #########
        
        #TTT
        t=0
        # print "init: %f"%(time.time()-t)
        t=time.time()
        filtered_pm = self.apply_filters(data, filters)
        filtered_pm_checkm = filtered_pm.copy()
        for i in range(filtered_pm_checkm.shape[0]):
            filtered_pm_checkm[i,i] = 0.0
        if not filtered_pm_checkm.nonzero()[0].size:
            self.error("No matrix elements were left after filtering")
            return False

        
        # print "applied: %f"%(time.time()-t)
        
        # ------ warn if operation is very large
        ncorrs = len(filtered_pm[np.nonzero(filtered_pm)])-filtered_pm.shape[0]
        if ncorrs > 1000:
            if not self.ask("You're trying to plot %d relationships: do you want to continue?"%ncorrs):
                return
        
        if self.widthmode == 1: # Absolute
            mX=minval
            MX=maxval
        elif self.widthmode == 0: # Relative to cutoffs
            #xs = zip(*self.sparse(dataabs(filtered_pm))[2]
            filtered_pm_range = filtered_pm.copy()
            for i in range(filtered_pm_range.shape[0]):
                filtered_pm_range[i,i] = 0.0
            clean_data = abs(filtered_pm_range[abs(filtered_pm_range)>0.0])
            MX = max(clean_data)
            mX = min(clean_data)
            if mX >= MX: 
                self.msg("ERROR: switching to absolute mode.\n")
                mX = minval
                MX = maxval
        #TTT
        t=time.time()
        self.msg("filterpos")
        #-- PLUS matrix: apply filters
        filtered = self.apply_filters( filtered_pm,
                                       [self.filter_sign(filtered_pm, positive=True)] )
        # get sparse representation
        self.msg("got positive")
        filtered_sparse = self.sparse(filtered)
        
        # print "uno: %d"%(len(filtered_sparse)),
        self.select_sparse("tmpselplus",filtered_sparse)
        self.msg("got sparse")

    
        cgOb = []
        # define plotted range of correlations
        cgOb.extend(self.makeCGO4(cas, rgbp, filtered_sparse, mX, MX, self.diagonal))
        self.msg("got cgos")
        
        #-- MINUS matrix: apply filters
        filtered = abs(self.apply_filters( filtered_pm,
                                       [self.filter_sign(filtered_pm, positive=False)] ) )
        self.msg("gotfiltered")

        # get sparse representation
        filtered_sparse = self.sparse(abs(filtered))
        self.msg("gotsparse")
        self.select_sparse("tmpselminus",filtered_sparse)

        #TTT
        t=time.time()
        cgOb.extend(self.makeCGO4(cas, rgbn, filtered_sparse, mX, MX, self.diagonal))
        self.msg("gotcgo")

        cmd.select("%s_residues" % (self.object), "tmpselplus or tmpselminus")
        cmd.delete("tmpselplus")
        cmd.delete("tmpselminus")

        #TTT
        # print "cged: %f"%(time.time()-t)
        t=time.time()
        # print cgOb,
        # print "Finished"
        # print cgOb
        cmd.load_cgo(cgOb,object)
                
        bidiG = nx.Graph()
        edgewidth = []
        condegs = []
        sparsecas=[]

        if self.doimg:
            for i in range(filtered_pm.shape[0]):
                for j in range(filtered_pm.shape[0]):
                    if abs(filtered_pm[i][j]) > 0.0 and i!=j:
                        label1 = cas[i].chain + cas[i].resi 
                        label2 = cas[j].chain + cas[j].resi 
                        bidiG.add_edge(label1,label2,weight=filtered_pm[i][j]*5)
            self.save_graph_pic(bidiG, self.img_file_name, mX, MX, dpi=600)
            
        
        if self.dolog:
            self.save_log(self.sparse(filtered_pm), self.log_file_name)
        
    def save_log(self,tuples,fname,extension="log",suffix=""):
        try:
            fh = open(fname+suffix+"."+extension,'w')
        except:
            self.error("ERROR: Could not write file %s\n" % fname)
            return False
        for i in tuples:
            id1 = [k for k, v in self.pdbdict.iteritems() if v == i[0]][0]
            id2 = [k for k, v in self.pdbdict.iteritems() if v == i[1]][0]
            fh.write("%s\t%s\t%s\t%s\t%1.5f\n" % (id1[0],id1[1:],id2[0],id2[1:],i[2]))
        fh.close()
        
    def save_graph_pic(self,bidiG, fname, mX, MX, dpi, extension="eps", suffix=""):
        # print "mX: %s" %mX
        # print "MX: %s" %MX
        try:
            fh = open(fname+suffix+"."+extension,'w')
        except:
            self.error("ERROR: Could not write file %s\n" % fname) # XXX move to doDCCM?
            return False
        fh.close()
        edgewidth=[]
        edgecol=[]
        condegs=[]
        bidiGdegs = bidiG.degree()
        for (n1, n2, w) in bidiG.edges(data=True):
            edgewidth.append(abs(w['weight']))
            edgecol.append(w['weight'])
        for n in bidiG.nodes():
            condegs.append(20*bidiGdegs[n])
        gl=nx.spring_layout(bidiG)
        nx.draw_networkx_nodes(bidiG, gl, node_size=condegs, node_color='0.5', width=0.0, linewidths=0.0)
        nx.draw_networkx_edges(bidiG, gl, width=edgewidth, alpha=0.3, edge_vmin=-MX, edge_vmax=MX, edge_color=edgecol, edge_map=plt.cm.magma)
        nx.draw_networkx_labels(bidiG, gl, font_size=8)
        plt.savefig(fname+suffix+"."+extension,dpi=dpi)
        plt.clf()
                
        

# ---------------
    def dochains(self):
        """
        performs checks and call findchainedcorrs on all residues of the target selection.
        """
        selchains = self.selchains
        rgbs=[]
        cgOb=[]
        if selchains == "":
            selchains = "all"
            self.error('Please specify a valid selection.')
            return False
        if not self.selchains in cmd.get_names("all")+["all"]:
            self.error("Chained selection (%s) is not a valid selection."%sel1)
            return False
        else: 
            selchainsids=set(cmd.identify(selchains,1))
            for j in cmd.identify(selchains):
                if not j in cmd.identify(self.reference):
                    self.error("Residues must belong to the reference object.")
                    return False
        if self.data is None:
            self.error("Please load a matrix first.")
            return False
            
        data = self.data.copy()

        if not data.shape[0] == data.shape[1]:
            self.error("Matrix must be square, it is %s.  Load aborted."%(",".join(["%d"%x for x in data.shape])))
            return False
		        # 2) check correspondence with reference object        
        refcas = self.get_calphas(self.reference)
        if not data.shape[0] == len(refcas):
            self.error("Matrix size (%d) must match reference object %s (%d)."%(data.shape[0], self.reference, len(refcas)))
            return False
        cas = self.get_calphas(selchains)

        # get display parameters
        self.msg("Found %d alpha-carbons in selection"%len(cas))


        data = self.data.copy()

        if self.widthmode == 1:         # absolute
            above = zip(*np.where(abs(data)>maxval))
            below = zip(*np.where(abs(data)<minval))
            if len(below) > 0:
                self.msg("WARNING: some matrix values were below the specified minimum. Will be reduced to minimum value.")
            if len(above) > 0:
                self.msg("WARNING: some matrix values were above the specified maximum. Will be reduced to maximum value.")
            if self.minval != 0.0:
                for i in below:
                    if data[i[0],i[1]] > 0:
                        data[i[0],i[1]] = self.minval
                    else:
                        data[i[0],i[1]] = - self.minval
                for i in above:
                    if data[i[0],i[1]] > 0:
                        data[i[0],i[1]] = self.maxval  
                    else:
                        data[i[0],i[1]] = - self.maxval
            else:                    
                for i in above:
                    if data[i[0],i[1]] > 0.0:
                        data[i[0],i[1]] = self.maxval
                    else:
                        data[i[0],i[1]] = - self.maxval                    

        cgOb = []
        #-- PLUS matrix: apply filters
        filters = []
        filter_names = []
        for filter_plugin in self.plugin_manager.getPluginsOfCategory("Filters"):
            if filter_plugin.is_activated:
                print "Using filter: %s" % filter_plugin.name
                filters.append(filter_plugin.plugin_object.get_filter(data, self.reference))
                filter_names.append(filter_plugin.name)

        filtered = self.apply_filters( data,
                                  filters + [self.filter_sign(self.data, positive=True)] )
    
        minval = self.minval
        maxval = self.maxval
    
        if self.widthmode == 1: # Absolute
            mX=minval
            MX=maxval
        elif self.widthmode == 0: # Relative to cutoffs
            filtered_pm_range = filtered.copy()
            for i in range(filtered_pm_range.shape[0]):
                filtered_pm_range[i,i] = 0.0
            clean_data = abs(filtered_pm_range[abs(filtered_pm_range)>0.0])
	    if len(clean_data) == 0:
                self.error("The final matrix contained no valid elements. Are you filtering too much?")
                return 
    
            MX = max(clean_data)
            mX = min(clean_data)
            if mX >= MX: 
                self.msg("ERROR: switching to absolute mode.\n")
                mX = minval
                MX = maxval
                                
        graph = nx.Graph()
        for i in range(data.shape[0]):
            for j in range(data.shape[0]):
                if filtered[i][j] > 0.0 and i != j:
                    graph.add_edge(i,j,weight=filtered[i][j])
  
        #-- MINUS matrix: apply filters
        if self.chainscolmode == 0:
            rgbs = [ [random(),random(),random()]  for i in cas ]
        else:
            rgbs = [ self.colorRGBchains for i in cas ]

        graphnodes = graph.nodes()
        for ca in range(len(cas)):
            thistuple = []

            if self.resi2i(cas[ca]) not in graphnodes:
                self.error("ERROR: residue %s of chain %s has no values in the filtered matrix. Are you filtering too much?\n" % (cas[ca].resi,cas[ca].chain))
                continue
            chainedG = self.find_chained_corrs(self.resi2i(cas[ca]),self.chains_d,self.chains_w,graph)
            if self.highest_avg_only_control: # FIND ONLY HIGHEST WEIGHT PATH
                maxw = [0,[]]                
                for i in [x for x in chainedG.nodes() if len(chainedG.edge[x]) < 2 and x != self.resi2i(cas[ca]) ]:                    
                    wsum = 0
                    path= nx.algorithms.shortest_path(chainedG,self.resi2i(cas[ca]),i)
                    for j in range(len(path)-1):
                        wsum += chainedG.edge[path[j]][path[j+1]]['weight']
                    #print wsum/len(path)
                    if wsum/len(path) > maxw[0]:
                        maxw[0] = wsum/len(path)
                        maxw[1] = path
                path = []
                for j in range(len(maxw[1])-1):
                    path.append((maxw[1][j],maxw[1][j+1],{ 'weight':chainedG.edge[maxw[1][j]][maxw[1][j+1]]['weight'] }))
                chainedG = nx.Graph(path)
                
            if self.exclude_leaves_control: # DO NOT PLOT LEAF NODES
                for i in [x for x in chainedG.nodes() if len(chainedG.edge[x]) < 2 and x != self.resi2i(cas[ca]) ]:
                    chainedG.remove_node(i)

            thistuple = [ (i,j,data['weight']) for i,j,data in chainedG.edges(data=True) ]
            cmd.load_cgo(self.makeCGO4(refcas, rgbs[ca], thistuple, mX, MX), "%s_%s" % (self.objectchains, cas[ca].chain+cas[ca].resi))
            self.select_sparse("%s_%s_residues" % (self.objectchains, cas[ca].chain+cas[ca].resi), thistuple)

            if self.doimg:
                bidiG = nx.Graph()
                for i in thistuple:
                    label1 = refcas[i[0]].chain + refcas[i[0]].resi 
                    label2 = refcas[i[1]].chain + refcas[i[1]].resi 
                    bidiG.add_edge(label1, label2, weight=i[2]*5)
                self.save_graph_pic(bidiG, self.img_file_name, mX, MX, dpi=600, suffix="_"+cas[ca].chain+cas[ca].resi)
            
            if self.dolog:
                self.save_log(thistuple, self.log_file_name, suffix="_"+cas[ca].chain+cas[ca].resi)
                

# --------------- 
    def makechainsCGO(self,tuples,cas,rgbs,data,selchains):
        tmpdata=[]
        if self.widthmode == 1: # Absolute
            mX=0.0
            MX=1.0
        elif self.widthmode == 0: # Relative to cutoffs
            for i in data:
                tmpdata.extend([abs(j) for j in i])
            mX=min(tmpdata)
            MX=max(tmpdata)

        if self.widthmode == 1: # Absolute
            mX=self.minval
            MX=self.maxval
        elif self.widthmode == 0: # Relative to cutoffs
            #xs = zip(*self.sparse(dataabs(filtered_pm))[2]
            filtered_pm_range = filtered_pm.copy()
            for i in range(filtered_pm_range.shape[0]):
                filtered_pm_range[i,i] = 0.0
            clean_data = abs(filtered_pm_range[abs(filtered_pm_range)>0.0])
            MX = max(clean_data)
            mX = min(clean_data)
            if mX >= MX: 
                self.msg("ERROR: switching to absolute mode.\n")
                mX = self.minval
                MX = self.maxval

        cgOb = []
        for i in range(len(tuples)):
            tmpcgOb=[]
            for j in tuples[i]:
                selca1 = cmd.select("chainedtmpsel1", "chain "+str(j[0][0])+" and resi "+j[0][1:]+" and "+self.reference)
                selca2 = cmd.select("chainedtmpsel2", "chain "+str(j[1][0])+" and resi "+j[1][1:]+" and "+self.reference)
                ca1 = self.get_calphas("chainedtmpsel1")[0]
                ca2 = self.get_calphas("chainedtmpsel2")[0]
                tmpcgOb.append(CYLINDER)
                tmpcgOb.extend(ca1.coord)
                tmpcgOb.extend(ca2.coord)
                tmpcgOb.append( (j[2]-mX)/(MX-mX)*Dx+ms )
                tmpcgOb.extend(rgbs[i])
                tmpcgOb.extend(rgbs[i])
            cgOb.append(tmpcgOb)
        cmd.delete("chainedtmpsel1")
        cmd.delete("chainedtmpsel2")
        return cgOb


# ---------------
    def makeCGO4(self, cas, rgb, cors, mX, MX, diagonal=False):
        """
        Make CGObject
        """
        # check that DX is not 0
        if not len(cors): return []
        cgOb = []
        DX = MX-mX
        if DX == 0:
            DX = 1
        ms = self.minsize
        Dx = self.maxsize - ms
        # diagonal
        dms = self.dminsize
        dDx = self.dmaxsize - dms
        # build CGO object
        # pre-define color for spheres
        cgOb.append(COLOR)
        cgOb.extend(rgb)
        for [a,b,x] in cors:
            caa = cab = None
            caa = cas[a]
            cab = cas[b]
            if a==b and diagonal:
                cgOb.append(SPHERE)
                cgOb.extend(caa.coord)
                cgOb.append( (x-mX)/DX*dDx+dms )
            cgOb.append(CYLINDER)
            cgOb.extend(caa.coord)
            cgOb.extend(cab.coord)
            cgOb.append( (x-mX)/DX*Dx+ms )
            cgOb.extend(rgb)
            cgOb.extend(rgb)
        return cgOb

# ------------------------------------------
    def do(self):
        self.doDCCM()
        self.doplot()

#    def docc(self):
#        self.doDCCM()
#        self.dochains()

# ------------------------------------------
class cVariable(Variable):
    def __init__(self, variable, parent, label, callback, width=5, tooltip=None, balloon=None):
        Variable.__init__(self)
        self.variable = variable
        self.frame = Frame(parent)
        self.entry = Entry(self.frame,textvariable=self,bg='black',fg='green',width=width)
        self.scroll= Scrollbar(self.frame,orient="horizontal",command=callback)
        self.label = Label(self.frame, text=label)
        self.scroll.pack(side=LEFT)
        self.entry.pack(side=LEFT)
        self.label.pack(side=LEFT)
        self.frame.pack(fill='x',padx=4,pady=1) # vertical
        if tooltip and balloon:
            balloon.bind(self.frame, tooltip)

    def config(self, *args, **kwargs):
        all = [self.entry]
        for a in all:
            try: 
                a.config(*args, **kwargs)
            except:
                sys.stderr.write("couldn't config in %s\n"%(a.__class__))
                pass

    def __getattr__(self, attr):
        if attr in self.__dict__:
            # this object has it
            return getattr(self, attr)
        # proxy to the wrapped object
        return getattr(self.variable, attr)        
            
class unpackedcVariable(cVariable):
    def __init__(self, variable, parent, label, callback, width=5, tooltip=None, balloon=None):
        Variable.__init__(self)
        self.variable = variable
        self.frame = Frame(parent)
        self.entry = Entry(self.frame,textvariable=self,bg='black',fg='green',width=width)
        self.scroll= Scrollbar(self.frame,orient="horizontal",command=callback)
        self.label = Label(self.frame, text=label)
        self.scroll.pack(side=LEFT)
        self.entry.pack(side=LEFT)
        self.label.pack(side=LEFT)
        #self.frame.pack(fill='x',padx=4,pady=1) # vertical
        if tooltip and balloon:
            balloon.bind(self.frame, tooltip)    
# ------------------------------------------
class cPlotTk(cPlot):
    #############################
    #
    # 1)  choose dccm file [ & cutoffs, output files]
    # 1b) extract correlations -> self.data
    # 2)  choose selection, colors, sizes
    # 2b) plot correlations
    #
    #############################

    ############# Functions for GUI outside the notebook
    
    def sel_plot(self,var):
        if var == "Direct":
            self.p_mode=0
        elif var == "Chained":
            self.p_mode=1
        elif var == "Delta":
            self.p_mode=2
           
    def widthmode_selection(self,var):
        if var == "Absolute":
            self.widthmode=1
            self.value_optframe.config(relief=FLAT)
            self.minvalf.config(state="normal")
            self.maxvalf.config(state="normal")
        elif var == "Relative":
            self.widthmode=0
            self.value_optframe.config(relief=SUNKEN)
            self.minvalf.config(state=DISABLED)
            self.maxvalf.config(state=DISABLED)
            
    def chainscolmode_selection(self,var):
        if var == "Random":
            self.chainscolmode=0
            self.Lab_col_ch.config(text="         \n                    ")
            self.colorFchains.configure(state=DISABLED)
        else:
            self.chainscolmode=1
            self.Lab_col_ch.config(text="Select a \n custom color")
            self.colorFchains.configure(state="normal")
        
    def show_matrix_plot(self):
        #print self.data
        if self.data is None:
            self.error("Please load a matrix first.")
            return False

        m = self.data.copy()
        size = m.shape[0]
        if self.widthmode == 0: # Relative
            for i in range(size):
                m[i,i] = 0.0
            clean_data = abs(m[abs(m)>0.0])
            maxval = max(clean_data)
            minval = min(clean_data)
        elif self.widthmode == 1: # Absolute
            minval = float(self.minvalf.get())
            maxval = float(self.maxvalf.get())            
            above = zip(*np.where(abs(m)>maxval))
            below = zip(*np.where(abs(m)<minval))
            if len(below) > 0:
                self.msg("WARNING: some matrix values were below the specified minimum. Will be reduced to minimum value.")
            if len(above) > 0:
                self.msg("WARNING: some matrix values were above the specified maximum. Will be reduced to maximum value.")
            if self.minval != 0.0:
                for i in below:
                    if m[i[0],i[1]] > 0:
                        m[i[0],i[1]] = minval
                    else:
                        m[i[0],i[1]] = - minval
                for i in above:
                    if m[i[0],i[1]] > 0:
                        m[i[0],i[1]] = maxval
                    else:
                        m[i[0],i[1]] = - maxval
            else:                    
                for i in above:
                    if m[i[0],i[1]] > 0.0:
                        m[i[0],i[1]] = maxval
                    else:
                        m[i[0],i[1]] = - maxval
        plt.pcolor(m, cmap=plt.get_cmap("magma"), vmax=maxval, vmin=-maxval)
        plt.colorbar()
        plt.axis([0,size,0,size])
        plt.show()


        
        
        
    def update_filters_status(self, plugin_source=None, modified_state=None):
        """ Update the status and GUI of active filters, or of the specified one """
        #TTT
        t = time.time()
        plugins = []
        if plugin_source:
            plugins = [plugin_source]
        else:
            plugins = [p.name for p in self.plugin_manager.getPluginsOfCategory("Filters")]
        # run through plugins
        for k in plugins:
            if self.states[k].get() == 1: # active
                self.checks[k] = self.plugin_manager.checkPluginStatus(k)
                if self.checks[k][0]:   # check passed!
                    self.labels[k].config(background="green")
                    self.plugin_manager.activatePluginByName(k,category="Filters")
                    #self.states_labels[k].config(text=self.checks[k][1])
                    if modified_state == k:
                        self.msg(k+": Activated ("+self.checks[k][1]+")", ["filter_ok"])
                else:                   # check failed :(
                    self.labels[k].config(background="yellow")
                    self.plugin_manager.deactivatePluginByName(k,category="Filters")
                    #self.states_labels[k].config(text=self.checks[k][1])
                    self.msg(k+": ERROR: "+self.checks[k][1], ["filter_alert"])
            else:                       # inactive
                self.labels[k].config(background="red")
                self.plugin_manager.deactivatePluginByName(k,category="Filters")
                #self.states_labels[k].config(text="")
                if modified_state == k:
                    self.msg(k+": Inactivated")
        
        # print "update_filters: %f"%(time.time()-t)
        # for plugin_info in self.plugin_manager.getPluginsOfCategory("Filters"):
        #    print "%s: %s" % (plugin_info.name, plugin_info.is_activated)
        return True
            
    def __init__(self,app):
        
        parent = app.root
        self.parent = parent
        
        self.radio = True
        self.checks= {}                 # check status of plugins
        #parent.plugin_manager.registerService(self.get_reference, "get_reference", VERSION)


        cPlot.__init__(self)
        self.colorRGBp = cmd.get_color_tuple(cmd._interpret_color(cmd, str(self.colorp)))
        self.colorRGBn = cmd.get_color_tuple(cmd._interpret_color(cmd, str(self.colorn)))
        self.status=StringVar(value="Ok.")

        # tooltip handler
        self.balloon = Pmw.Balloon(parent)
        # main window
        self.closeString = 'Close xPyder'
        self.dialog = Pmw.Dialog(parent,
                                 buttons = (self.closeString,),
                                 title = 'xPyder',
                                 command = self.buttonPressed)
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))
        self.dialog.bind('<Return>',self.buttonPressed)
        
        ###################### FUNCTIONS "OUTSIDE" THE NOTEBOOK
        
        right_box=Frame(self.dialog.interior())
        right_box.pack(side=RIGHT,padx=10)
                    
        try:
            photo1 = PhotoImage(file=INSTALLDIR+"/xPyder_logo.gif")
            panel1 = Label(right_box, image=photo1)
            panel1.pack(side='top')
            panel1.image=photo1
        except:
            print "WARNING: xPyder Logo image could not be loaded."
     
     
        self.p_group = Pmw.Group(right_box,tag_text='xPyder plugin')
        self.p_group.pack(fill = 'both')
        
#        self.c_group = Pmw.Group(right_box,tag_text='Output settings')
#        self.c_group.pack(fill = 'both')
        
        self.f_group = Pmw.Group(right_box,tag_text='Plotting mode')
        self.f_group.pack(fill = 'both')
        
#        Pmw.Color.changecolor(self.p_group,
#        background = 'red3', foreground = 'white')
        
        ###############################PLUGIN RELATED CHECKBOXES
        self.frames={}
        self.chkbuts={}
        self.states=[]
        self.labels={}
        self.cbfuncs={}
        self.vars={}            
        self.plug_frame=Frame(self.p_group.interior())
        self.plug_frame.pack(side=TOP)
        # Create the ScrolledFrame.
        self.scroll = Pmw.ScrolledFrame(self.plug_frame,
                labelpos = 'n', label_text = 'Filters control',
                usehullsize = 1,
                hull_width = 260,
                hull_height = 300
        )
        self.scroll.pack(side="top",padx=5,pady=5)

        self.states = {}
        self.states_labels = {}
        self.filters=[ plugin_info.name for plugin_info in self.plugin_manager.getPluginsOfCategory("Filters") ]
        k=1

        def __acommand(plugin_name):
            return lambda: self.update_filters_status(modified_state=plugin_name)
        
        for i in self.filters:
            var = IntVar()
      #      self.frames[i]=Frame(self.p_group.interior())
            self.frames[i]=Frame(self.scroll.interior())
            self.frames[i].pack(fill="both",padx=5,pady=10)
      #      self.frames[i].grid(row=k,column=0)
            self.labels[i] = Label(self.frames[i],text = "    ",background="red",justify=LEFT)
      #      self.labels[i].pack(side=LEFT)
            self.labels[i].grid(row=0,column=0)
            self.states_labels[i] = Label(self.frames[i],text = "",justify=LEFT)
      #      self.states_labels[i].pack(side=RIGHT)
            self.states_labels[i].grid(row=0,column=2)
            chk = Checkbutton(self.frames[i], text=str(i), variable=var, 
                              command=__acommand(i), justify=LEFT)
            chk.grid(row=0, column=1)
      #      chk.pack(side=RIGHT,padx=10,pady=10)
            self.states[i] = var
            k=k+1
        
        frame_3=Frame(self.p_group.interior())
        frame_3.pack(side=BOTTOM)
        s=Button(frame_3,text="Refresh",command=self.update_filters_status)
#        s.pack(padx=10,pady=10,side="right")
#        root.mainloop()                                       

        # Plot selection
        self.plot_mode= Pmw.RadioSelect(self.f_group.interior(),
            buttontype = 'radiobutton',
            orient = 'vertical',
            labelpos = 'n',
            command = self.sel_plot,
            label_text = 'Mode selection',
            hull_borderwidth = 2,
            hull_relief = 'ridge'
            )

        self.plot_mode.pack(side="top",expand = 1, padx = 10, pady = 10)
        self.plot_mode.add('Direct')
        self.plot_mode.add('Chained')
        self.plot_mode.add('Delta')
        self.plot_mode.setvalue('Direct')
 
 
      
        # buttons
        self.pltB = Pmw.ButtonBox(self.f_group.interior(), padx=0)
        self.pltB.pack(side=BOTTOM,padx=10,pady=10)
        #self.pltB.add('Change Color',command=self.tk_color_dialog)
        self.pltB.add('Plot',command=self.doplot)
        
        # self.statusL = Label(self.dialog.interior(), textvariable=self.status, bg='black', fg='#cceeee', anchor='sw', height=2, justify='left')
        # self.statusL.pack(fill=BOTH,pady=10,side=BOTTOM)
        self.statusL = Pmw.ScrolledText(self.dialog.interior(), text_bg='black', text_fg='#cceeee', text_height=3, text_state=DISABLED)
        self.statusL.tag_config("msg", foreground="#cceeee")
        self.statusL.tag_config("alert", foreground="red")
        self.statusL.tag_config("filter_ok", foreground="green")
        self.statusL.tag_config("filter_alert", foreground="yellow")
        self.statusL.pack(fill=BOTH,pady=5,anchor='sw',side=BOTTOM)
        self.msg("Ready...")
        self.balloon.bind(self.statusL, 'Status bar')
                
        ######################

        # the title
        w = Label(self.dialog.interior(),
                  text = 'xPyder \n Marco Pasi, Matteo Tiberti, Alberto Arrigoni, Elena Papaleo \n <xpyder.pymol@gmail.com>',
                  background = 'black',
                  foreground = 'green', 
                  #pady = 20,
                  )
        w.pack(fill='x', padx = 4, pady = 4, side=TOP)
        self.balloon.bind(w, 'Thank you for using xPyder!')

        # the basic notebook
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both',padx=3,pady=3,expand=1)
        page = self.notebook.add('xPyder')

############################################################# log and graph file saving
        self.back_x=Pmw.Group(page, tag_text="Graph and Log")
        self.back_x.pack(side=BOTTOM,fill='both')

        self.img_frame=Frame(self.back_x.interior())
        self.img_frame.grid(row=0,column=1)
        self.img_file_var=StringVar()    
        self.img_file_var.set(os.getcwd()+"/correlations_graph")
        self.img_file=Entry(self.img_frame,textvariable=self.img_file_var, state=DISABLED, width=40)
        self.img_file.pack(side=BOTTOM)
        self.img_file_butt = Button(self.img_frame, text="Save img file as...", command=self.save_img_file,state=DISABLED)    
        self.img_file_butt.pack(side=BOTTOM)
        self.balloon.bind(self.img_file_butt, 'Save graph image')
        self.img_file_st=IntVar()
        self.img_file_opt = Checkbutton(self.img_frame, text="Save graph image", variable=self.img_file_st, command=self.img_file_mod)
        self.img_file_opt.pack(side=BOTTOM)
    
        self.log_frame=Frame(self.back_x.interior())
        self.log_frame.grid(row=0,column=0,padx=5)
        self.log_file_var=StringVar()    
        self.log_file_var.set(os.getcwd()+"/correlations")
        self.log_file=Entry(self.log_frame,textvariable=self.log_file_var, state=DISABLED, width=40)
        self.log_file.pack(side=BOTTOM)
        self.log_file_butt = Button(self.log_frame, text="Save log file as...", command=self.save_log_file,state=DISABLED)    
        self.log_file_butt.pack(side=BOTTOM)
        self.balloon.bind(self.log_file_butt, 'Save log file')
        self.log_file_st=IntVar()
        self.log_file_opt = Checkbutton(self.log_frame, text="Save log file", variable=self.log_file_st, command=self.log_file_mod)
        self.log_file_opt.pack(side=BOTTOM)    

    #############################################################
    ########################### Matrix value

        #### DCCM: first TAB
        ## handle extraction of data from Correlation Matrix
        corrG = Pmw.Group(page,tag_text='Matrix')
        corrG.pack(fill = 'both')
        # input file
        left = Frame(corrG.interior())
        left.pack(side=LEFT, fill=X, expand=1)
        self.order=Frame(left)
        self.order.pack(side=LEFT,expand=1)
        self.DCCMf = Pmw.EntryField(self.order,
                                    labelpos='w',
                                #    label_pyclass = FileDialogButtonClassFactory.get(self.setDCCM,filter=("*.dat")),
                                    validate = {'validator': self.quickFileValidation,},
                                    value = self.DCCM)
                                 #   label_text = 'Choose matrix file:')
        def open_DCCM_matrix():
            self.DCCM = tkFileDialog.askopenfilename()
            self.DCCMf.setentry(self.DCCM)
        self.DCCMf.grid(column=1,row=0,columnspan=2)

        trial=Button(self.order,command=open_DCCM_matrix,text="Choose matrix file")
    #    trial.pack(side=TOP)
        trial.grid(row=0,column=0)

        corB = Pmw.ButtonBox(self.order, padx=0)
     #   corB.pack(side=RIGHT)
        corB.grid(column=1,row=1)
        corB.add('Load matrix',command=self.doDCCM)
        #showMatrix = Pmw.ButtonBox(self.order, padx=0)
     #   corB.pack(side=RIGHT)
        corB.add('Show matrix',command=self.show_matrix_plot)
        #self.offset_opt = Checkbutton(self.dsize_frame, text="Matrix offset", variable=self.matrix_offset, command=matrix_offset_cb)
        #self.diagonal_opt.pack(side=TOP)
        #ds_optframe = Frame(self.dsize_frame)
        #ds_optframe.pack(side=BOTTOM)
        #self.dminsizef = self.makeDoubleVar(ds_optframe, "Minimum Size", self.changeDmSize,
        #                                    tooltip="""Minimum size of plotted spheres,


        #self.moffset_frame.grid(row=1,column=0)
        #self.moffset_frame.grid(column=0,row=1,sticky=S,padx=13,pady=13)


        self.moffset_frame = Frame(self.order, relief=SUNKEN, bd=2, padx=5, pady=5)
        self.moffset_frame.grid(row=1, column=0, sticky= S,padx=13, pady=13)

        self.moffset_b = IntVar()

        def moffset_cb():
            self.moffset = self.moffset_b.get()
            if self.moffset:
                self.moffset_frame.config(relief=FLAT)
                self.moffsetf.config(state="normal")
            else:
                self.moffset_frame.config(relief=SUNKEN)
                self.moffsetf.config(state=DISABLED)
                self.moffsetf.set(0)

        self.moffset_opt = Checkbutton(self.moffset_frame, text="Matrix offset", variable=self.moffset_b, command=moffset_cb)
        #self.offset_opt.grid(row=1, column=0)
        self.moffset_opt.pack()
        #moffset_optframe = Frame(self.dsize_frame)
        #moffset_optframe.pack(side=BOTTOM)
        self.moffsetf = self.makeIntVar(self.moffset_frame, "Offset amount", self.changeMatrixOffset,
                                            tooltip="""Number of lines and rows for the matrix to be offset""",
                                            packed=False)
        self.moffsetf.frame.pack()
    #    self.DCCMf.pack(fill=X, expand=1, padx = 5, pady = 5)
        self.moffsetf.set(self.moffset)
        moffset_cb()

   



        
        ## RIGHT-side data
        right = Frame(corrG.interior())
        self.balloon.bind(right, 'Information about the values\nloaded from the matrix file.')
        #self.sum1 = DataSummary(self.data, right, "orange", "#333333")
        self.sum1 = DataSummaryH(self.data, parent=right, root=parent, title="Matrix", fg="orange", bg="#eeeeee")
        # DCCM parameters
        self.DCCM1_params(self.data)
        right.pack(side=BOTTOM, padx=0, pady=0)

        #### STATUS
        groupP=Pmw.Group(page,tag_text='Pymol objects')
        groupP.pack(side=TOP, fill=X)
        self.referencef = Pmw.EntryField(groupP.interior(),
                                      labelpos='w',
                                      label_text='Input object name',
                                      value=self.reference,
                                      validate = {'validator':self.referenceValidation,},
                                      )
        self.referencef.pack(side="left",pady=10,padx=10) # vertical 
        
        self.objecti = Pmw.EntryField(groupP.interior(),
                                      labelpos='w',
                                      label_text='Output object name',
                                      value=self.object,
                                      )
        self.objecti.pack(side="top",pady=10,padx=10) # vertical 
    #    self.balloon.bind(self.objecti, """Name of the PyMOL object containing Plot.""")

        ##### FILTERS PAGE (second page)
        
        self.page_filt = self.notebook.add('Filters')
        # XXX activate all plugins
        self.labels_list = []
        
        def bind_children(father, function, kwargs):
            """ temp function to bind to update function relevant events of ALL CHILDREN """
            EVENTS = "<ButtonRelease-1> <KeyRelease>".split()
            for event in EVENTS:
                father.bind(event, lambda ev: function(**kwargs), '+')
            for child in father.winfo_children():
                bind_children(child, function, kwargs)
        
        # Create the ScrolledFrame.
        self.filt_scroll = Pmw.ScrolledFrame(self.page_filt,
                labelpos = 'n', label_text = 'Filters',
                usehullsize = 0.1,
                hull_width = 800,
                hull_height = 650,
        )
        self.filt_scroll.pack(side=TOP,padx=10,pady=10,fill='both',expand=1)
        for plugin_info in self.plugin_manager.getPluginsOfCategory("Filters"):
            outstring = "\n"            
            # print "Initializing plugin   ",plugin_info.name
            pGroup = plugin_info.plugin_object.GUI(self.filt_scroll.interior())
            #pGroup.pack(side=TOP,padx=10,pady=5,fill=BOTH,expand=1)
            bind_children(pGroup, self.update_filters_status, {"plugin_source": plugin_info.name})
            for i in plugin_info.depends:
                try:
                    outstring = outstring+"Will use %s capability from plugin %s\n" % (i, plugin_info.depends_plugins[i])
                except:
                    pass
            if outstring == "\n":
                outstring = ""
            self.labels_list.append(Label(pGroup.interior(),text = outstring))
            self.labels_list[-1].pack()
        
        ##### PLOT SETTING (third page)
        
#        self.page = self.notebook.add('Plot settings')
        self.colors = Pmw.Group(page,tag_text='Plot settings')
        self.colors.pack(fill = 'both')
        # self.cols_frame=Frame(self.colors.interior())
        # self.cols_frame.pack(side=TOP,padx=5,pady=5)
        col_label = Frame(self.colors.interior())
        # col_label.pack(side=TOP)
        col_label.grid(column=1,row=0,columnspan=2)
        self.colorFp = Button(col_label, relief=SUNKEN, bd=2, height=1, width=1, text="+", font="Helvetica 16", command=self.tk_color_dialogp)
        self.colorFp.grid(row=0,column=0)
        self.balloon.bind(self.colorFp, "Color for positive values")
        #self.colorFp.create_text(12.5,12.5,text='+',font="Helvetica 16")
        self.colorFn   = Button(col_label, relief=SUNKEN, bd=2, height=1, width=1, text="-", font="Helvetica 16", command=self.tk_color_dialogn)
        self.colorFn.grid(row=0,column=3)
        self.balloon.bind(self.colorFn, "Color for negative values")
        #self.colorFn.create_text(12.5,12.5,text='-',font="Helvetica 16")
        self.update_color_frames()
        lb_col_l = Label(col_label,text="Color for positive values")
        lb_col_l.grid(row=0,column=1)
        lb_col_r = Label(col_label,text="Color for negative values")
        lb_col_r.grid(row=0,column=4)
        
        
        ######### CGO related
        widthmode = Frame(self.colors.interior(), borderwidth = 2, relief = 'ridge')
        self.widthmodef= Pmw.RadioSelect(widthmode,
        buttontype = 'radiobutton',
        orient = 'vertical',
        labelpos = 'n',
        command = self.widthmode_selection,
        label_text = 'Width normalization mode',
        Button_state = 'active'
        )

        self.value_optframe = Frame(widthmode, relief=SUNKEN, bd=2)
        self.minvalf = self.makeDoubleVar(self.value_optframe, "Min value", self.changemVal,
                                           tooltip="""Minimum value,
i.e. size corresponding to the lowest
matrix values.""")
        self.minvalf.set(0.0)
        self.maxvalf = self.makeDoubleVar(self.value_optframe, "Max value", self.changeMVal,
                                          tooltip="""Maximum value,
i.e. size corresponding to the highest
matrix values.""")
        self.maxvalf.set(1.0)

        self.widthmodef.add('Relative')
        self.widthmodef.add('Absolute')
        # initialise opt frame
        self.widthmodef.setvalue('Relative')
        self.widthmode_selection('Relative')
        self.widthmodef.pack(side = 'top', expand = 0, padx = 10, pady = 5)
        self.value_optframe.pack(side='bottom',padx=3,pady=3)

        widthmode.grid(column=0,row=0,rowspan=2,padx=10,pady=10)

        # size
        size_frame=Frame(self.colors.interior())
        size_frame.grid(column=1,row=1,sticky=S)
        bottom_frame=Frame(size_frame)
        bottom_frame.pack(side=BOTTOM,padx=15,pady=15)
        self.minsizef = self.makeDoubleVar(bottom_frame, "Minimum Size", self.changemSize,
                                           tooltip="""Minimum size of plotted sticks,
i.e. size corresponding to the lowest
matrix values.""")
        self.minsizef.set(self.minsize)
        self.maxsizef = self.makeDoubleVar(bottom_frame, "Maximum Size", self.changeMSize,
                                           tooltip="""Maximum size of plotted sticks,
i.e. size corresponding to the highest
matrix values.""")
        self.maxsizef.set(self.maxsize)
        self.balloon.bind(self.referencef, """Name of the PyMOL object containing Plot.""")

        # size 2
        self.dsize_frame=Frame(self.colors.interior(), relief=SUNKEN, bd=2)
        self.dsize_frame.grid(column=2,row=1,sticky=S,padx=13,pady=13)

        self.diagonal_b = IntVar()
        def diagonal_cb():
            self.diagonal = self.diagonal_b.get()
            if self.diagonal:
                self.dsize_frame.config(relief=FLAT)
                self.dmaxsizef.config(state="normal")
                self.dminsizef.config(state="normal")
            else:
                self.dsize_frame.config(relief=SUNKEN)
                self.dmaxsizef.config(state=DISABLED)
                self.dminsizef.config(state=DISABLED)

        self.diagonal_opt = Checkbutton(self.dsize_frame, text="Plot diagonal", variable=self.diagonal_b, command=diagonal_cb)
        self.diagonal_opt.pack(side=TOP)
        ds_optframe = Frame(self.dsize_frame)
        ds_optframe.pack(side=BOTTOM)
        self.dminsizef = self.makeDoubleVar(ds_optframe, "Minimum Size", self.changeDmSize,
                                            tooltip="""Minimum size of plotted spheres,
i.e. size corresponding to the lowest
matrix values.""")
        self.dminsizef.set(self.dminsize)
        self.dmaxsizef = self.makeDoubleVar(ds_optframe, "Maximum Size", self.changeDMSize,
                                            tooltip="""Maximum size of plooted spheres,
i.e. size corresponding to the highest
matrix values.""")
        self.dmaxsizef.set(self.dmaxsize)
        # initialise entries
        diagonal_cb()
        self.empty_l = Label(size_frame)
        self.empty_l.pack(side=BOTTOM)
        #all_sizes.grid(column=1,row=1)
        #XXX checkbox: add state to object / replace object
        
        ##### CHAINED CORRELATIONS PAGE (FOURTH PAGE) #####

        self.page = self.notebook.add('Chained')
        # chainsF = Frame(page)
        self.chainsG = Pmw.Group(self.page,tag_text='Chained parameters')
        self.colors_ch = Pmw.Group(self.page,tag_text='Color settings')
        #self.chainsK = Pmw.Group(self.page,tag_text='Cylinder settings')
        self.chainsG.pack(fill = 'both',padx=10,pady=10)
        #self.chainsK.pack(fill = 'both',padx=10,pady=10)
        self.colors_ch.pack(fill = 'both',padx=10,pady=10)
        self.chainsF = Frame(self.chainsG.interior())
        self.chainsF.pack(fill = 'both',padx=10,pady=10, side=TOP)
        self.chainsH = Frame(self.chainsG.interior())
        self.chainsH.pack(fill = 'both',padx=10,pady=10)
        self.selchainsf = Pmw.EntryField(self.chainsF,
                                            labelpos='w',
                                            label_text='Selection: ',
                                            value=self.selchains,
                                            validate = {'validator':self.selectionValidation,},
                                            )

        self.selchainsf.pack(fill='both',padx=10,pady=30,expand=0) # vertical
        self.balloon.bind(self.selchainsf, """PYMOL selection to be used \n for chained plot""")
        
        self.wf = self.makeDoubleVar(self.chainsF, "Width", self.changew, tooltip="Width")
        self.wf.set(3)
        self.df = self.makeDoubleVar(self.chainsF, "Depth", self.changed, tooltip="Depth")
        self.df.set(3)

        #self.objectchainsf = Pmw.EntryField(self.chainsF,
        #                                labelpos='w',
        #                                label_text='Object name',
        #                                value=self.objectchains,
        #                                validate = {'validator':self.objectValidation,},
        #                                )
        #self.objectchainsf.pack(fill='x',padx=4,pady=1,expand=0) # vertical
        #self.balloon.bind(self.referencef, """Name of the PyMOL object containing Plot.""")

        def exclude_neighs_cb():
            self.exclude_neighs_control = self.exclude_neighs.get()
        def highest_avg_only_cb():
            self.highest_avg_only_control = self.highest_avg_only.get()
        def exclude_leaves_cb():
            self.exclude_leaves_control = self.exclude_leaves.get()

        self.exclude_neighs = IntVar()
        self.exclude_neighsc = Checkbutton(self.chainsF, text="Exclude neighbors of previous root node", variable = self.exclude_neighs, 
                                command=exclude_neighs_cb)
        self.exclude_neighsc.pack(pady=3,side=BOTTOM,anchor=W)
        self.highest_avg_only = IntVar()        
        self.highest_avg_onlyc = Checkbutton(self.chainsF, text="Plot path with highest average weight only", variable = self.highest_avg_only, 
                                command=highest_avg_only_cb)     
        self.highest_avg_onlyc.pack(pady=3,side=BOTTOM,anchor=W)
        self.exclude_leaves = IntVar()        
        self.exclude_leavesc = Checkbutton(self.chainsF, text="Do not plot terminal nodes", variable = self.exclude_leaves, 
                                command=exclude_leaves_cb)     
        self.exclude_leavesc.pack(pady=3,side=BOTTOM,anchor=W)



#        self.colors_ch=Frame(self.chainsG.interior())
#        self.colors_ch.pack(side=BOTTOM,padx=20,pady=20)
        
        self.custom_col=0
        self.chainscolmodef= Pmw.RadioSelect(self.colors_ch.interior(),
                                        buttontype = 'radiobutton',
                                        orient = 'vertical',
                                        labelpos = 'n',
                                        command = self.chainscolmode_selection,
                                        label_text = 'Color mode',
                                        hull_borderwidth = 2,
                                        hull_relief = 'ridge',
                                        )
        self.chainscolmodef.pack()
    #    self.chainscolmodef.grid(row=1,column=3)
        self.chainscolmodef.add('Random')
        self.chainscolmodef.add('Choose color')
        self.chainscolmodef.setvalue('Random')
        
#       self.plot_mode.getvalue('Correlations')
        
        self.Lab_col_ch=Label(self.colors_ch.interior(),text="         \n                    ")
        self.Lab_col_ch.pack(side=BOTTOM)    
        
        self.Lab_col_c=Label(self.colors_ch.interior(),text="       ")
    #   self.Lab_col_c.grid(row=1,column=4) 
        
        self.colorFchains = Button(self.colors_ch.interior(), relief=SUNKEN, bd=2, height=1, width=1, text="", font="Helvetica 16", command=self.tk_color_dialogchains, state=DISABLED)
    #    self.colorFchains = Button(self.colors_ch.interior(), relief=SUNKEN, bd=2, height=1, width=1, text="", font="Helvetica 16", command=self.tk_color_dialogchains, state=DISABLED)
        self.colorFchains.pack(padx=10,pady=10)
    #   self.colorFchains.grid(row=1,column=2)
        self.balloon.bind(self.colorFchains, "Color for chained")
        self.update_color_frames_chains()
        

        self.pltBchains = Pmw.ButtonBox(self.chainsH, padx=0)
        self.pltBchains.pack(side=BOTTOM)

############# NETWORK ANALYSES
        self.graphdata = [None,None,None,None,None]
        self.graph = None
        page = self.notebook.add('Network Analysis')
        NET_FIRST = Pmw.Group(page,tag_text='Generate graph')
        NET_FIRST.pack(fill='x',expand=0)
        self.left_net=Frame(NET_FIRST.interior())
        self.left_net.pack(side=LEFT)
       
        def create_net_graph():
            if self.data is None:
                self.error("Please load a matrix first!")
                self.graph = None
                self.graphdata = [None,None,None,None,None]
                return
            try:
                self.graphdata = self.calc_graph(self.data)                
                self.lab_one.config(text="Number of nodes: %d" %self.graphdata[0])
                self.lab_two.config(text="Number of edges: %d" %self.graphdata[1])
                self.lab_three.config(text="Number of components: %d" %len(self.graphdata[2]))  
                self.n_clustersf.setlist(range(1,len(self.graphdata[2])+1))
                self.n_clustersf.selectitem(0,setentry=1)
                nodes = sorted(self.graph.nodes, key=lambda x: (x[0], int(x[1:])))
                self.sp_res1f.setlist(nodes)
                self.sp_res2f.setlist(nodes)
                self.sp_res1f.selectitem(0,setentry=1)
                self.sp_res2f.selectitem(0,setentry=1)
                self.sumgraph.update(self.graphdata[3])
            except:
                self.error("An error occurred while creating the graph.")
                self.reset_net_graph()
                raise
                
                
                
        butt_graph=Button(self.left_net,command=create_net_graph,text="Generate graph")
        butt_graph.pack(padx=10,pady=10)

        self.lab_one=Label(self.left_net,text="")
        self.lab_one.pack()
        self.lab_two=Label(self.left_net,text="")
        self.lab_two.pack()
        self.lab_three=Label(self.left_net,text="")
        self.lab_three.pack()
        
        self.right_net=Frame(NET_FIRST.interior())
        self.right_net.pack(side=RIGHT)
        self.sumgraph = DataSummaryHGraph(None, parent=self.right_net, root=parent, title="Degree distribution", fg='lightblue', bg='#eeeeee')
        self.sumgraph.update(self.graphdata[3])

        NET_SECOND = Pmw.Group(page,tag_text='Hubs')
        NET_SECOND.pack(fill='x',expand=0)
        self.hubs_n=1
        self.hubs_sel = Pmw.EntryField(NET_SECOND.interior(),
                                     labelpos='w',
                                     label_text = 'Minimum degree: ',
                                     value = self.hubs_n)
        self.hubs_sel.pack(side=LEFT,padx=10,pady=10)
        butt_hubs=Button(NET_SECOND.interior(),command=self.select_hubs,text="Select hubs")
        butt_hubs.pack(padx=10,pady=10,side=LEFT) 
        #butt_hubs_rep=Button(NET_SECOND.interior(),command=self.select_hubs,text="Hubs report")
        #butt_hubs_rep.pack(padx=10,pady=10,side=LEFT) 

        NET_THIRD = Pmw.Group(page,tag_text='Components')
        NET_THIRD.pack(fill='x',expand=0)
        self.left_clust=Frame(NET_THIRD.interior())
        self.left_clust.pack(side=LEFT,pady=10)

        def n_clustersf_callback(x=None):
            try:
                cursel = int(self.n_clustersf.getcurselection()[0])
                self.lab_clusters_data.config(text="Number of nodes: %d" %len(self.graphdata[2][cursel-1].nodes()) )
                self.lab_clusters_data2.config(text="Number of edges: %d" %len(self.graphdata[2][cursel-1].edges()) )

            except:
                pass


        self.n_clustersf = Pmw.ComboBox(self.left_clust,
              label_text = 'Component #',
              labelpos = 'nw',
              selectioncommand = n_clustersf_callback,
              scrolledlist_items = [],
     
              )
        self.n_clustersf.pack(side=LEFT)
        self.n_clustersf.configure(entry_state='readonly')

        
        self.lab_clusters_data = Label(self.left_clust,text="")
        self.lab_clusters_data.pack(side=TOP)
        self.lab_clusters_data2 = Label(self.left_clust,text="")
        self.lab_clusters_data2.pack(side=BOTTOM)


        def changeText():
            pass

        self.buttons_frame=Frame(NET_THIRD.interior())
        self.buttons_frame.pack(side=RIGHT)
        butt_sel_clu=Button(self.buttons_frame,command=self.select_cluster,text="Select cluster")
        butt_sel_clu.pack(padx=10,pady=10,side=LEFT) 
        
        butt_plot_clusters=Button(self.buttons_frame,command=self.plot_cluster,text="Plot Cluster")
        butt_plot_clusters.pack(padx=10,pady=10,side=LEFT)
        
        NET_FOURTH = Pmw.Group(page,tag_text='Paths')
        NET_FOURTH.pack(fill='x',expand=0)
        self.left_short=Frame(NET_FOURTH.interior())
        self.left_short.pack(side=LEFT)
        self.big_right=Frame(self.left_short)
        self.big_right.pack(side=RIGHT)
        self.center_frame=Frame(self.big_right)
        self.center_frame.pack(side=LEFT)
        def changeText1(x=None):
            pass
        self.sp_res1f = Pmw.ComboBox(self.left_short,
                label_text = 'Residue 1',
                labelpos = 'nw',
                scrolledlist_items = [],
                entry_width = 3
        )
        self.sp_res1f.pack(side=TOP)
        self.sp_res1f.configure(entry_state='readonly')

        self.sp_res2f = Pmw.ComboBox(self.left_short,
                label_text = 'Residue 2',
                labelpos = 'nw',
                scrolledlist_items = [],
                entry_width = 3
        )
        self.sp_res2f.pack(side=TOP)
        self.sp_res2f.configure(entry_state='readonly')
        
        self.sp_maxhops = Pmw.EntryField(self.left_short,
                label_text = 'Maximum length',
                labelpos = 'nw',
                entry_width = 3
        )
        self.sp_maxhops.pack()
        
        self.right_short=Frame(NET_FOURTH.interior())
        self.right_short.pack(side=RIGHT)
        
##################################
        def selectionCommand():
            self.sele_clust = self.box_clust.getcurselection()
            if len(sels) == 0:
                print 'No selection'
            else:
                print 'Selection:', sels[0]
        self.sf = Pmw.ScrolledFrame(self.right_short,
                labelpos = 'n', label_text = '',
                usehullsize = 1,
                hull_width = 200,
                hull_height = 100,
        )
     #   self.sf.pack(side=LEFT)
        self.box_clust = Pmw.ScrolledListBox(self.sf.interior(),
                items=[],
                labelpos='nw',
                label_text='Selected cluster',
            #    listbox_height = 6,
                selectioncommand=selectionCommand,
            #    usehullsize = 1,
                hull_width = 200,
                hull_height = 100,
        )
#        self.box_clust.pack(side=LEFT)
        
####################################

        butt_sp_calcspf = Button(self.left_short,command=self.get_all_paths,text="Calculate paths")
        butt_sp_calcspf.pack(padx=10,pady=10,side=BOTTOM)
      #  def treeview_sort_column(tv, col, reverse):
     #       l = [(tv.set(k, col), k) for k in tv.get_children('')]
    #        l.sort(reverse=reverse)
    # rearrange items in sorted positions
   #         for index, (val, k) in enumerate(l):
  #              tv.move(k, '', index)
    # reverse sort next time
#            tv.heading(col, command=lambda: \
 #               treeview_sort_column(tv, col,not reverse))
#########################################################
        def sortby(tree, col, descending):
            """Sort tree contents when a column is clicked on."""
    # grab values to sort
            data = [(tree.set(child, col), child) for child in tree.get_children('')]

    # reorder data
            #print data
            data.sort(reverse=descending,key=lambda x: float(itemgetter(0)(x)))
            for indx, item in enumerate(data):
                tree.move(item[1], '', indx)

    # switch the heading so that it will sort in the opposite direction
            tree.heading(col,command=lambda col=col: sortby(tree, col,int(not descending)))        
#########################################################
        
        self.table_columns=("#","length","sum of weights","average weight")
        self.table=ttk.Treeview(self.center_frame,column=self.table_columns,show="headings",height=7)
   
        vsb = ttk.Scrollbar(self.center_frame, orient="vertical", command=self.table.yview)
        hsb = ttk.Scrollbar(self.center_frame, orient="horizontal", command=self.table.xview)
        #vsb.pack(side=RIGHT, expand=YES, fill=Y)
        self.table.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.table.grid(column=0, row=0, sticky='nsew')
        vsb.grid(column=1, row=0, sticky='ns',in_=self.center_frame)
#        hsb.grid(column=0, row=1, sticky='ew',in_=self.center_frame)

        self.center_frame.grid_columnconfigure(0, weight=1)
        self.center_frame.grid_rowconfigure(0, weight=1)       
        column_widths=[20,90,120,120]
        for i in range(len(self.table_columns)):
            col = self.table_columns[i]
            self.table.heading(col, text=col.title(),
                command=lambda c=col: sortby(self.table, c, 0))
            # XXX tkFont.Font().measure expected args are incorrect according
            #     to the Tk docs
            #self.table.column(col,width=tkFont.Font().measure(col.title()))
            self.table.column(col,width=column_widths[i])

        tree_data = []
        #self.table.pack()
        butt_sp_selectf=Button(self.right_short,command=self.select_sp,text="Select path")
        butt_sp_selectf.pack(padx=10,pady=10,side=TOP) 
        butt_sp_plotf=Button(self.right_short,command=self.plot_sp,text="Plot path")
        butt_sp_plotf.pack(padx=10,pady=10,side=BOTTOM)

#### DELTA CORRELATIONS (fifth page)
        page = self.notebook.add('Delta')
        DCCM2G = Pmw.Group(page,tag_text='Load another matrix')
        DCCM2G.pack(fill = 'both')
        lDCCM2F = Frame(DCCM2G.interior())
        lDCCM2F.pack(side=LEFT, fill=X, expand=1)
        self.matrix_1=Frame(lDCCM2F)
        self.matrix_1.pack(side=LEFT)
        self.DCCM2f = Pmw.EntryField(self.matrix_1,
                                     labelpos='w',
                                #     label_text = 'Choose matrix file:',
                                #     label_pyclass = FileDialogButtonClassFactory.get(self.setDCCM2,filter=("*.dat")),
                                     validate = {'validator': self.quickFileValidation,},
                                     value = self.DCCM2)
        # self.balloon.bind(self.DCCM, 'Choose a DCCM file.')
    ####################################################
        def open_delta_matrix1():
            self.DCCM2 = tkFileDialog.askopenfilename()
            self.DCCM2f.setentry(self.DCCM2)
        matrix_1=Button(self.matrix_1,command=open_delta_matrix1,text="Choose second matrix file")
    #    trial.pack(side=TOP)
        matrix_1.grid(row=0,column=0)
    #    self.DCCMf.pack(fill=X, expand=1, padx = 5, pady = 5)
        self.DCCM2f.grid(column=1,row=0,columnspan=2)
    ####################################################
    
    #    self.DCCM2f.pack(fill=X, expand=1, padx = 5, pady = 5)
        cor2B = Pmw.ButtonBox(self.matrix_1, padx=0)
        cor2B.grid(column=1,row=1)
        cor2B.add('Load second matrix',command=self.doDCCM2)

        ### summary
        # rightG = Pmw.Group(DCCM2G.interior(), tag_text='Second DCCM:')
        # self.balloon.bind(rightG, 'Information about the correlations\nloaded from the second DCCM file.')
        # right = Frame(rightG.interior())
        # self.makeSummary(right, [self.DpC2, self.DnC2, self.DaC2], self.data2,
        #                  fg='lightblue', bg='$333333')
        # right.pack(fill=X,padx=0, pady=0)
        # rightG.pack(side=BOTTOM, padx=10, pady=5)
        ### histogram w/ popup
        right = Frame(DCCM2G.interior())
        self.balloon.bind(right, 'Information about the values\nloaded from the second matrix file.')
        self.sum2 = DataSummaryH(self.data2, parent=right, root=parent, title="Second matrix histogram", 
                                fg='lightblue', bg='#eeeeee')
        # DCCM2 parameters
        self.DCCM2_params(self.data2)
        right.pack(side=BOTTOM,padx=0, pady=0)
        ###
        
        #---
        deltaG = Pmw.Group(page,tag_text='Delta')
        deltaG.pack(fill = 'both')
        ldeltaF = Frame(deltaG.interior())
        ldeltaF.pack(side=LEFT, fill = X, expand=1)

        self.matrix_2=Frame(ldeltaF)
        self.matrix_2.pack(side=TOP)

        self.deltaAf = Pmw.EntryField(self.matrix_2,
                                     labelpos='w',
                                #     label_text = 'Choose alignment file:',
                                #     label_pyclass = FileDialogButtonClassFactory.get(self.setDCCM3,filter=("*.*")),
                                     validate = {'validator': self.quickFileValidation,},
                                     value = self.deltaA)
        self.balloon.bind(self.deltaAf, 'The first sequence in the aligment must be relative to the reference.')
    #    self.deltaAf.pack(fill=X, expand = 1, padx = 10, pady = 5)
        
       
    ####################################################
        def open_alignment():
            self.deltaA = tkFileDialog.askopenfilename()
            self.deltaAf.setentry(self.deltaA)
        alignment=Button(self.matrix_2,command=open_alignment,text="Choose alignment file")
    #    trial.pack(side=TOP)
        alignment.grid(row=0,column=0,pady=10)
    #    self.DCCMf.pack(fill=X, expand=1, padx = 5, pady = 5)
        self.deltaAf.grid(column=1,row=0,columnspan=2)
    ####################################################        
        
        #XXX move out
        def set_inv_delta(v):
            self.inv_delta = v == 'Second matrix - first matrix'

        ldeltaFF = Frame(ldeltaF)
        ldeltaFF.pack(side=RIGHT, expand=1, anchor=SE)
        
        # aliB = Pmw.ButtonBox(ldeltaFF, padx=0, pady=0)
        # aliB.pack()
        # aliB.add('Load alignment',command=self.doAlignment)
        
        
        self.deltadirf= Pmw.RadioSelect(ldeltaF,
                                        buttontype = 'radiobutton',
                                        orient = 'vertical',
                                        labelpos = 'n',
                                        command = set_inv_delta,
                                        label_text = 'Delta direction',
                                        hull_borderwidth = 2,
                                        hull_relief = 'ridge',
                                        )
        self.deltadirf.pack(side = RIGHT, padx = 5, pady = 5, anchor=SE)
        self.deltadirf.add('First matrix - second matrix')
        self.deltadirf.add('Second matrix - first matrix')
        self.deltadirf.setvalue('First matrix - second matrix')
        
        self.delta_bAlign = IntVar()
        self.delta_bA_opt = Checkbutton(ldeltaFF, text="Use alignment", variable=self.delta_bAlign)
        self.delta_bA_opt.pack(side=LEFT,pady=10,padx=10)

        ### summary
        # rightG = Pmw.Group(deltaG.interior(), tag_text='Delta DCCM:')
        # self.balloon.bind(rightG, 'Information about the delta correlations.')
        # right = Frame(rightG.interior())
        # self.makeSummary(right, [self.dDpC, self.dDnC, self.dDaC], self.delta,
        #                  fg='white', bg='$333333')
        # right.pack(fill=X,padx=0, pady=0)
        # rightG.pack(side=BOTTOM, padx=10, pady=5)
        ### histogram w/ popup
        right = Frame(deltaG.interior())
        self.balloon.bind(right, 'Information about the delta matrix.')
        self.sumd = DataSummaryH(self.delta, parent=right, root=parent, title="Delta matrix histogram", 
                                fg='white')
        # DCCMdelta parameters
        self.DCCMd_params(self.delta)
        right.pack(side=BOTTOM, padx=0, pady=0)
        ###
        
        #---
        dummyF = Frame(page)
        dummyF.pack(side = RIGHT, anchor=NE, pady=10)
        delB = Pmw.ButtonBox(dummyF, padx = 0, pady=0)
        delB.pack(side=LEFT)
        delB.add('Compute delta matrix',command=self.dodelta)


        dpltB = Pmw.ButtonBox(dummyF, padx = 10, pady = 5,
                              hull_borderwidth = 2,
                              hull_relief = 'groove')
       # dpltB.pack(side=LEFT)
        dpltB.add('Plot',command=self.delta_do_all)
        


#### ABOUT TAB
        page = self.notebook.add('About')
        group = Pmw.Group(page, tag_text='About xPyder')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        text = """
  xPyder creates a 3D Plot using a matrix and a selection.

  1) First  xPyder loads a  matrix in
     ASCII format from a file.

  2) After  the matrix  is loaded, xPyder filters and  plots it 3D  on a
     specified selection in the PyMOL workspace.

  Created by Marco Pasi, Matteo Tiberti, Alberto Arrigoni, Elena Papaleo \n <xpyder.pymol@gmail.com>'
"""
        lfre=Frame(group.interior())
        bar=Scrollbar(lfre,)
        ll=Text(lfre,yscrollcommand=bar.set,background="#ffffaa",font="Times 15",height=10)
        bar.config(command=ll.yview)    
        ll.insert(END,text)
        ll.pack(side=LEFT,expand="yes",fill="both")
        bar.pack(side=LEFT,expand="yes",fill="y")
        lfre.pack(expand="yes",fill="both")
            
        self.notebook.setnaturalsize()
        self.showAppModal()


                
    def reset_net_graph(self):
        self.graph = None
        self.graphdata = [None,None,None,None,None]
        self.graphpaths = None
        self.graphmcache = None
        self.lab_one.config(text="")
        self.lab_two.config(text="")
        self.lab_three.config(text="")  
        self.n_clustersf.setlist([])
        self.sp_res1f.setlist([])
        self.sp_res2f.setlist([])
        self.sp_res1f.clear()
        self.sp_res2f.clear()
        self.n_clustersf.clear()
        self.lab_clusters_data.config(text="")
        self.lab_clusters_data2.config(text="")
        map(self.table.delete, self.table.get_children(''))
        #self.sumgraph.update(None)
        self.sumgraph.figure.clear()
        self.sumgraph.plot.show()

    def reset_delta_graph(self):
        #self.sumgraph.update(None)
        self.sumd.figure.clear()
        self.sumd.plot.show()

    def save_log_file(self):
         self.log_file_name=tkFileDialog.asksaveasfilename(initialdir=self.log_file_var.get())
         self.log_file_var.set(self.log_file_name)
                  
    def save_img_file(self):        
         self.img_file_name=tkFileDialog.asksaveasfilename(initialdir=self.log_file_var.get())
         self.img_file_var.set(self.img_file_name)

    def log_file_mod(self):
        if int(self.log_file_st.get())==1:
            self.log_file_butt.config(state=ACTIVE)
            self.log_file.config(state="normal")
        if int(self.log_file_st.get())==0:
            self.log_file_butt.config(state=DISABLED)
            self.log_file.config(state="disabled")

    def img_file_mod(self):
        if int(self.img_file_st.get())==1:
            self.img_file_butt.config(state=ACTIVE)
            self.img_file.config(state="normal")
        if int(self.img_file_st.get())==0:
            self.img_file_butt.config(state=DISABLED)
            self.img_file.config(state="disabled")


# ------------------------------------------
# GUI building functions
    def makeVar(self, var, parent, label, callback, width=5, tooltip=None, packed=True):
        if packed:
            return cVariable(var, parent, label, callback, width=width, tooltip=tooltip, balloon=self.balloon)
        else:
            return unpackedcVariable(var, parent, label, callback, width=width, tooltip=tooltip, balloon=self.balloon)
        
    def makeIntVar(self, parent, label, callback, tooltip=None, packed=True):
        var = IntVar()
        return self.makeVar(var, parent, label, callback, tooltip=tooltip, packed=packed)
        
    def makeDoubleVar(self, parent, label, callback, tooltip=None):
        var = DoubleVar()
        return self.makeVar(var, parent, label, callback, tooltip=tooltip)

# ------------------------------------------
    status_buffer = ""
    def msg(self, m, tags=[]):
        #XXX 
        cPlot.msg(self, m)
        try:
            if len(self.status_buffer):
                self.statusL.insert('end', self.status_buffer)
                self.status_buffer=""
            #self.status.set( "\n".join(self.status.get().split("\n")[-5:])  + "\n" + m)
            self.statusL.component('text').config(state='normal')
            self.statusL.insert('end', '\n'+m, *tags)
            self.statusL.see('end')
            self.statusL.component('text').config(state=DISABLED)
            #self.statusL.config(fg='#cceeee')
            #self.statusL.update()
        except:
            # label is not loaded yet
            self.status_buffer+='\n'+m
            # import traceback
            # traceback.print_exc()

    def alertmsg(self, m):
        self.msg(m, tags=["alert"])

    def ask(self, msg):
        return tkMessageBox.askokcancel('Warning',msg)
        
    def warn(self, msg):
        tkMessageBox.showwarning('Warning',msg)
        self.msg("Ok.")
        
    def error(self, msg):
        tkMessageBox.showerror('Error',msg)
        self.msg("Ok.")

    def finish(self):
        self.msg("Ok.")
        
    def DCCM_params(self, data,
                    sumobj,
                    foutp, foutn):
        """
        Write a summary of the loaded data in the provided structures, and
        optionally write out + and - matrices.
        """
        if data is None:
            sumobj.update(None, None, None)
            return 0,0
            
        posC =  self.apply_filters(data, (self.filter_sign(data, positive=True),) )
        negC =  self.apply_filters(data, (self.filter_sign(data, positive=False),) )
        data = data[:]
        sumobj.update(data, posC, negC)

        if foutp:
            posCmfh=open(foutp,'w')
            for i in posC:
                for j in i:
                    print >> posCmfh, j,
                print >> posCmfh, "\n"
            posCmfh.close()
                
        if foutn:
            negCmfh=open(foutn,'w')
            for i in negC:
                for j in i:
                    print >> negCmfh, j,
                print >> negCmfh, "\n"
            negCmfh.close()

        return np.sum(posC!=0), np.sum(negC!=0)

# ---------------        
    def DCCM1_params(self, data, foutp=None, foutn=None):
        "Dummy method to call original DCCM_params"
        return self.DCCM_params(
            data, self.sum1,
            foutp, foutn)
        
# ---------------        
    def DCCM2_params(self,data, foutp=None, foutn=None):
        "Dummy method to call original DCCM_params"
        return self.DCCM_params(
            data, self.sum2,
            foutp, foutn)
                                
# ---------------        
    def DCCMd_params(self,data, foutp=None, foutn=None):
        "Dummy method to call original DCCM_params"
        return self.DCCM_params(
            data, self.sumd,
            foutp, foutn)

# ------------------------------------------
# Callback functions
    def buttonPressed(self,result):
        if result == self.closeString or result is None:
            self.dialog.withdraw()

    def update_color_frames(self):
        rgbSp = "#%s" % "".join(["%02x"%(x*255.) for x in self.colorRGBp])
        rgbSn = "#%s" % "".join(["%02x"%(x*255.) for x in self.colorRGBn])
        self.colorFp.configure(background=rgbSp)
        self.colorFn.configure(background=rgbSn)

    def update_color_frames_chains(self):
        rgbcSchains = "#%s" % "".join(["%02x"%(x*255.) for x in self.colorRGBchains])
        self.colorFchains.configure(background=rgbcSchains)
        
    def tk_color_dialogp(self):
        self.colorRGBp = self.tk_color_dialog(self.colorRGBp)
        self.msg("Positive color changed to [%5.3f,%5.3f,%5.3f]"%(self.colorRGBp[0],
                                                                  self.colorRGBp[1],
                                                                  self.colorRGBp[2]))
        self.update_color_frames()
            
    def tk_color_dialogchains(self):
        self.colorRGBchains = self.tk_color_dialog(self.colorRGBchains)
        self.msg("Chained color changed to [%5.3f,%5.3f,%5.3f]"%(self.colorRGBchains[0],
                                                                              self.colorRGBchains[1],
                                                                              self.colorRGBchains[2]))
        self.update_color_frames_chains()

    def tk_color_dialogn(self):
        self.colorRGBn = self.tk_color_dialog(self.colorRGBn)
        self.msg("Negative color changed to [%5.3f,%5.3f,%5.3f]"%(self.colorRGBn[0],
                                                                  self.colorRGBn[1],
                                                                  self.colorRGBn[2]))
        self.update_color_frames()
        
    def tk_color_dialog(self, col):
        rgb = tuple([x*255 for x in col])
        color = tkColorChooser.Chooser(
            initialcolor=rgb,title='Choose color').show()
        if color[0] is not None:
            return [color[0][0]/255.,
                    color[0][1]/255.,
                    color[0][2]/255.]
        else:
            return col

    def dummy(self, a, *args):
        pass

    def changeMatrixOffset(self, a, *args):
        if len(args)!=0:
            a=args[0]
        val = max(0, int(round(float(self.moffsetf.get()))) + int(a))
        self.moffsetf.set(val)
        self.msg("Matrix offset value changed to %d" %val)
    def changemSize(self, a, *args):
        if len(args)!=0:
            a=args[0]
        val = max(0.0, float(self.minsizef.get()) + float(a)*0.05)
        self.minsizef.set(val)
        self.msg("Minimum Size changed to %f"%val)
    def changeMSize(self, a, *args):
        if len(args)!=0:
            a=args[0]
        val = max(0.0, float(self.maxsizef.get()) + float(a)*0.05)
        self.maxsizef.set(val)
        self.msg("Maximum Size changed to %f"%val)
    def changeDmSize(self, a, *args):
        if len(args)!=0:
            a=args[0]
        val = max(0.0, float(self.dminsizef.get()) + float(a)*0.05)
        self.dminsizef.set(val)
        self.msg("Minimum Size (spheres) changed to %f"%val)
    def changeDMSize(self, a, *args):
        if len(args)!=0:
            a=args[0]
        val = max(0.0, float(self.dmaxsizef.get()) + float(a)*0.05)
        self.dmaxsizef.set(val)
        self.msg("Maximum Size (spheres) changed to %f"%val)
    def changemSizechains(self, a, *args):
        if len(args)!=0:
            a=args[0]
        val = max(0.0, float(self.minsizechainsf.get()) + float(a)*0.05)
        self.minsizechainsf.set(val)
        self.msg("Minimum Size changed to %f"%val)
    def changeMSizechains(self, a, *args):
        if len(args)!=0:
            a=args[0]
        val = max(0.0, float(self.maxsizechainsf.get()) + float(a)*0.05)
        self.maxsizechainsf.set(val)
        self.msg("Maximum Size changed to %f"%val)
    def changemVal(self, a, *args):
        if len(args)!=0:
            a=args[0]
        val = max(0.0, float(self.minvalf.get()) + float(a)*0.05)
        self.minvalf.set(val)
        self.msg("Minimum Value changed to %f"%val)
    def changeMVal(self, a, *args):
        if len(args)!=0:
            a=args[0]
        val = max(0.0, float(self.maxvalf.get()) + float(a)*0.05)
        self.maxvalf.set(val)
        self.msg("Maximum Value changed to %f"%val)
    def setDCCM(self, f):
        self.DCCMf.setentry(f)
    def setDCCM2(self, f):
        self.DCCM2f.setentry(f) 
    def setDCCM3(self, f):
        self.deltaAf.setentry(f)
    def changew(self, a, *args):
        if len(args)!=0:
            a=args[0]
        val = max(1, int(self.wf.get()) + int(a))
        self.wf.set(val)
        self.msg("asd")
    def changed(self, a, *args):
        if len(args)!=0:
            a=args[0]
        val = max(1, int(self.df.get()) + int(a))
        self.df.set(val)
        self.msg("premium")
    
# ------------------------------------------
# Validator functions
    def quickFileValidation(self, s):
        # print "checking %s: %s"%(s,os.path.isfile(s))
        if s == '':
            return Pmw.PARTIAL
        elif os.path.isfile(s):
            self.msg("Matrix set to %s"%s)
            return Pmw.OK
        elif os.path.exists(s):
            return Pmw.PARTIAL
        else:
            return Pmw.PARTIAL

    def selectionValidation(self, s):
        if s in cmd.get_names("all")+["all"]:
            self.msg("Selection changed to %s"%s)
            return Pmw.OK
        return Pmw.PARTIAL

    def objectValidation(self, s):
        #XXX if exists, ask if append (state) or substitute
        if s in cmd.get_names("all"):
            self.msg("Selection changed to %s"%s)
            return Pmw.OK
        return Pmw.PARTIAL
#         return Pmw.OK

    def referenceValidation(self, s):
        try:
            self.reference = self.referencef.getvalue()
        except:
            # print "fuori"
            pass
        if s in cmd.get_names("all")+["all"]:
            self.msg("Selection changed to %s"%s)
            return Pmw.OK
        return Pmw.PARTIAL

        
# ------------------------------------------
# Main data-handling and plotting functions
    def doDCCM(self):
        # ------- Get parameters from GUI
        oldDCCM = self.DCCM
        self.DCCM = self.DCCMf.getvalue()
        oldreference = self.reference
        self.reference = self.referencef.getvalue()
        self.moffset = self.moffsetf.get()
        # -------
        if not cPlot.doDCCM(self):
            # restore
            self.reference = oldreference
            self.DCCM = oldDCCM
            return False
        self.reset_net_graph()
        self.reset_delta_graph()
        self.finish()
        return True

# ---------------
    def doplot(self):
        self.minsize = float(self.minsizef.get())
        self.maxsize = float(self.maxsizef.get())
        self.dminsize = float(self.dminsizef.get())
        self.dmaxsize = float(self.dmaxsizef.get())
        self.minval  = float(self.minvalf.get())
        self.maxval  = float(self.maxvalf.get())
        self.object = self.objecti.getvalue()
        self.dolog = self.log_file_st.get()
        self.doimg = self.img_file_st.get()
        if self.dolog:
            self.log_file_name = self.log_file_var.get()
        if self.doimg:
            self.img_file_name = self.img_file_var.get()
        self.reference = self.referencef.getvalue()
        self.update_filters_status()
        if self.p_mode == 0:
            cPlot.doplot(self)
        elif self.p_mode == 1:
            self.dochains()
        elif self.p_mode == 2:
            cPlot.doplot(self)          # check of data is done inside doplot
        self.finish()

# --------------- 
    def dochains(self):
        self.minsizechains = float(self.minsizef.get())
        self.maxsizechains = float(self.maxsizef.get())
        self.objectchains = self.objecti.getvalue()
        self.selchains = self.selchainsf.getvalue()
        self.chains_w =	int(self.wf.get())
        self.chains_d =	int(self.df.get())
        #self.update_filters_status()
        cPlot.dochains(self)

# --------------- 
    def dodelta(self):
        self.alignment = None
        if self.delta_bAlign.get() == 1:
            if not self.doAlignment():
                return False
        if not cPlot.dodelta(self):
            return False
        self.finish()
        return True
        
# --------------- 
    def doDCCM2(self):
        # ------- Get parameters from GUI
        oldDCCM2 = self.DCCM2
        self.DCCM2 = self.DCCM2f.getvalue()
        self.moffset = self.moffsetf.get()
        self.reset_delta_graph()
        # -------
        if not cPlot.doDCCM2(self):
            self.DCCM2 = oldDCCM2
            return False
        self.finish()
        return True

# --------------- 
    def doAlignment(self):
        # ------- Get parameters from GUI
        oldA = self.deltaA
        self.deltaA = self.deltaAf.getvalue()
        # -------
        if self.deltaA != "":
            if not cPlot.doAlignment(self):
                self.deltaA = oldA
                return False
        self.reset_delta_graph()
        return True

# --------------- 
    def delta_do_all(self):
        if not self.doDCCM2():
            return False
        if not self.doAlignment():
            return False
        if not self.dodelta():
            return False
        self.plot_mode.invoke("Delta")
        self.p_mode = 2
        if not self.doplot():
            return False
        self.finish()
        return True

# ------------------------------------------
    def showAppModal(self):
        self.dialog.geometry('1000x700')
        self.dialog.show()
             
#
# The classes PmwFileDialog and PmwExistingFileDialog and the _errorpop function
# are taken from the Pmw contrib directory.  The attribution given in that file
# is:
################################################################################
# Filename dialogs using Pmw
#
# (C) Rob W.W. Hooft, Nonius BV, 1998
#
# Modifications:
#
# J. Willem M. Nissink, Cambridge Crystallographic Data Centre, 8/2002
#    Added optional information pane at top of dialog; if option
#    'info' is specified, the text given will be shown (in blue).
#    Modified example to show both file and directory-type dialog
#
# No Guarantees. Distribute Freely. 
# Please send bug-fixes/patches/features to <r.hooft@euromail.com>
#
################################################################################
import os,fnmatch,time
import Tkinter,Pmw
#Pmw.setversion("0.8.5")

def _errorpop(master,text):
    d=Pmw.MessageDialog(master,
                        title="Error", 
                        message_text=text,
                        buttons=("OK",))
    d.component('message').pack(ipadx=15,ipady=15)
    d.activate()
    d.destroy()
    
class PmwFileDialog(Pmw.Dialog):
    """File Dialog using Pmw"""
    def __init__(self, parent = None, **kw):
    # Define the megawidget options.
        optiondefs = (
            ('filter',    '*',              self.newfilter),
            ('directory', os.getcwd(),      self.newdir),
            ('filename',  '',               self.newfilename),
            ('historylen',10,               None),
            ('command',   None,             None),
            ('info',      None,             None),
            )
        self.defineoptions(kw, optiondefs)
        # Initialise base class (after defining options).
        Pmw.Dialog.__init__(self, parent)

        self.withdraw()

         # Create the components.
        interior = self.interior()
    
        if self['info'] is not None:
            rowoffset=1
            dn = self.infotxt()
            dn.grid(row=0,column=0,columnspan=2,padx=3,pady=3)
        else:
            rowoffset=0

        dn = self.mkdn()
        dn.grid(row=0+rowoffset,column=0,columnspan=2,padx=3,pady=3)
        del dn

    # Create the directory list component.
        dnb = self.mkdnb()
        dnb.grid(row=1+rowoffset,column=0,sticky='news',padx=3,pady=3)
        del dnb

    # Create the filename list component.
        fnb = self.mkfnb()
        fnb.grid(row=1+rowoffset,column=1,sticky='news',padx=3,pady=3)
        del fnb
  
    # Create the filter entry
        ft = self.mkft()
        ft.grid(row=2+rowoffset,column=0,columnspan=2,padx=3,pady=3)
        del ft

    # Create the filename entry
        fn = self.mkfn()
        fn.grid(row=3+rowoffset,column=0,columnspan=2,padx=3,pady=3)
        fn.bind('<Return>',self.okbutton)
        del fn

    # Buttonbox already exists
        bb=self.component('buttonbox')
        bb.add('OK',command=self.okbutton)
        bb.add('Cancel',command=self.cancelbutton)
        del bb

        Pmw.alignlabels([self.component('filename'),
                         self.component('filter'),
                         self.component('dirname')])

    def infotxt(self):
        """ Make information block component at the top """
        return self.createcomponent(
                'infobox',
                (), None,
                Tkinter.Label, (self.interior(),),
                width=51,
                relief='groove',
                foreground='darkblue',
                justify='left',
                text=self['info']
            )

    def mkdn(self):
        """Make directory name component"""
        return self.createcomponent(
        'dirname',
        (), None,
        Pmw.ComboBox, (self.interior(),),
        entryfield_value=self['directory'],
        entryfield_entry_width=40,
            entryfield_validate=self.dirvalidate,
        selectioncommand=self.setdir,
        labelpos='w',
        label_text='Directory:')

    def mkdnb(self):
        """Make directory name box"""
        return self.createcomponent(
        'dirnamebox',
        (), None,
        Pmw.ScrolledListBox, (self.interior(),),
        label_text='directories',
        labelpos='n',
        hscrollmode='none',
        dblclickcommand=self.selectdir)

    def mkft(self):
        """Make filter"""
        return self.createcomponent(
        'filter',
        (), None,
        Pmw.ComboBox, (self.interior(),),
        entryfield_value=self['filter'],
        entryfield_entry_width=40,
        selectioncommand=self.setfilter,
        labelpos='w',
        label_text='Filter:')

    def mkfnb(self):
        """Make filename list box"""
        return self.createcomponent(
        'filenamebox',
        (), None,
        Pmw.ScrolledListBox, (self.interior(),),
        label_text='files',
        labelpos='n',
        hscrollmode='none',
        selectioncommand=self.singleselectfile,
        dblclickcommand=self.selectfile)

    def mkfn(self):
        """Make file name entry"""
        return self.createcomponent(
        'filename',
        (), None,
        Pmw.ComboBox, (self.interior(),),
        entryfield_value=self['filename'],
        entryfield_entry_width=40,
            entryfield_validate=self.filevalidate,
        selectioncommand=self.setfilename,
        labelpos='w',
        label_text='Filename:')
    
    def dirvalidate(self,string):
        if os.path.isdir(string):
            return Pmw.OK
        else:
            return Pmw.PARTIAL
        
    def filevalidate(self,string):
        if string=='':
            return Pmw.PARTIAL
        elif os.path.isfile(string):
            return Pmw.OK
        elif os.path.exists(string):
            return Pmw.PARTIAL
        else:
            return Pmw.OK
        
    def okbutton(self):
        """OK action: user thinks he has input valid data and wants to
           proceed. This is also called by <Return> in the filename entry"""
        fn=self.component('filename').get()
        self.setfilename(fn)
        if self.validate(fn):
            self.canceled=0
            self.deactivate()

    def cancelbutton(self):
        """Cancel the operation"""
        self.canceled=1
        self.deactivate()

    def tidy(self,w,v):
        """Insert text v into the entry and at the top of the list of 
        the combobox w, remove duplicates"""
        if not v:
            return
        entry=w.component('entry')
        entry.delete(0,'end')
        entry.insert(0,v)
        list=w.component('scrolledlist')
        list.insert(0,v)
        index=1
        while index<list.index('end'):
            k=list.get(index)
            if k==v or index>self['historylen']:
                list.delete(index)
            else:
                index=index+1
        w.checkentry()

    def setfilename(self,value):
        if not value:
            return
        value=os.path.join(self['directory'],value)
        dir,fil=os.path.split(value)
        self.configure(directory=dir,filename=value)
        
        c=self['command']
        if callable(c):
            c()

    def newfilename(self):
        """Make sure a newly set filename makes it into the combobox list"""
        self.tidy(self.component('filename'),self['filename'])
    
    def setfilter(self,value):
        self.configure(filter=value)

    def newfilter(self):
        """Make sure a newly set filter makes it into the combobox list"""
        self.tidy(self.component('filter'),self['filter'])
        self.fillit()

    def setdir(self,value):
        self.configure(directory=value)

    def newdir(self):
        """Make sure a newly set dirname makes it into the combobox list"""
        self.tidy(self.component('dirname'),self['directory'])
        self.fillit()

    def singleselectfile(self):
        """Single click in file listbox. Move file to "filename" combobox"""
        cs=self.component('filenamebox').curselection()
        if cs!=():
            value=self.component('filenamebox').get(cs)
            self.setfilename(value)

    def selectfile(self):
        """Take the selected file from the filename, normalize it, and OK"""
        self.singleselectfile()
        value=self.component('filename').get()
        self.setfilename(value)
        if value:
            self.okbutton()

    def selectdir(self):
        """Take selected directory from the dirnamebox into the dirname"""
        cs=self.component('dirnamebox').curselection()
        if cs!=():
            value=self.component('dirnamebox').get(cs)
            dir=self['directory']
            if not dir:
                dir=os.getcwd()
            if value:
                if value=='..':
                    dir=os.path.split(dir)[0]
                else:
                    dir=os.path.join(dir,value)
            self.configure(directory=dir)
            self.fillit()

    def askfilename(self,directory=None,filter=None):
        """The actual client function. Activates the dialog, and
       returns only after a valid filename has been entered 
           (return value is that filename) or when canceled (return 
           value is None)"""
        if directory!=None:
            self.configure(directory=directory)
        if filter!=None:
            self.configure(filter=filter)
        self.fillit()
        self.canceled=1 # Needed for when user kills dialog window
        self.activate()
        if self.canceled:
            return None
        else:
            return self.component('filename').get()

    lastdir=""
    lastfilter=None
    lasttime=0
    def fillit(self):
        """Get the directory list and show it in the two listboxes"""
        # Do not run unnecesarily
        if self.lastdir==self['directory'] and self.lastfilter==self['filter'] and self.lasttime>os.stat(self.lastdir)[8]:
            return
        self.lastdir=self['directory']
        self.lastfilter=self['filter']
        self.lasttime=time.time()
        dir=self['directory']
        if not dir:
            dir=os.getcwd()
        dirs=['..']
        files=[]
        try:
            fl=os.listdir(dir)
            fl.sort()
        except os.error,arg:
            if arg[0] in (2,20):
                return
            raise
        for f in fl:
            if os.path.isdir(os.path.join(dir,f)):
                dirs.append(f)
            else:
                filter=self['filter']
                if not filter:
                    filter='*'
                if fnmatch.fnmatch(f,filter):
                    files.append(f)
        self.component('filenamebox').setlist(files)
        self.component('dirnamebox').setlist(dirs)
    
    def validate(self,filename):
        """Validation function. Should return 1 if the filename is valid, 
           0 if invalid. May pop up dialogs to tell user why. Especially 
           suited to subclasses: i.e. only return 1 if the file does/doesn't 
           exist"""
        return 1

class PmwExistingFileDialog(PmwFileDialog):
    def filevalidate(self,string):
        if os.path.isfile(string):
            return Pmw.OK
        else:
            return Pmw.PARTIAL
        
    def validate(self,filename):
        if os.path.isfile(filename):
            return 1
        elif os.path.exists(filename):
            _errorpop(self.interior(),"This is not a plain file")
            return 0
        else:
            _errorpop(self.interior(),"Please select an existing file")
            return 0

# Create demo in root window for testing.
if __name__ == '__main__':
    class App:
        def my_show(self,*args,**kwargs):
            pass
        def Cmd(self, do=None):
            if do:
                cmd.do("@%s"%do)
        def destroy(self):
            app.root.destroy()
            cmd.quit()

    app = App()
    app.root = Tkinter.Tk()
    Pmw.initialise(app.root)
    app.root.title('Some Title')

    # import pymol, thread
    # pymol.glutThread = thread.get_ident()
    widget = cPlotTk(app)
    exitButton = Tkinter.Button(app.root, text = 'Exit', command = app.destroy)
    exitButton.pack()
    if len(sys.argv) > 1:
        app.Cmd(sys.argv[1])
    else:
        app.Cmd()
    app.root.mainloop()


def cPlotfun(_self, **kwargs):
    """
    syntax: cPlot [options]
    """
    kwargs.update({
        'sel1': 'a3',
        'sel2': '',
        })
    c = cPlot(**kwargs)
    c.do()
    #c.docc()

cmd.extend('cPlot', cPlotfun)
