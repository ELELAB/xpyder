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
import tkMessageBox
import tkFileDialog
import Pmw
from pymol import cmd
import networkx as nx
import copy 


class SecondaryStructureFilterPlugin(FilterPlugin):
    """
    Depends:  get_calphas
    Provides: distance_matrix
    
    The "get_dm" method will be provided using the capability name "distance_matrix".
    """
    def __init__(self):                 # init must have no arguments
        FilterPlugin.__init__(self)
        # data
        self.ss       = None            # secondary structure definition (see parse_ss)
        self.doalpha  = True
        self.dobeta   = True
        self.betaco   = 0.7
        # parameters
        self.ss_fname = ""              # the ss file, either DSSP or PDSSP
        #self.ss_state 
        # parameters used to compute current ss, for caching
        self._ss_fname = ""
        self._reference= None
        self.ss_mode = "file"

    def parse_ss(self, filen):
        """
        filen: filename

        parse_ss outputs a dictionary composed of:
        - a bool (key "complete")
        - two lists of lists (keys "ah" and "be")
        {
          "be":[[],...,[]],
          "ah":[[],...,[]],
          "complete":bool
        }

        "complete" is True if complete beta-sheet information was available (dssp), False otherwise.

        For each secondary structure definition of a protein/complex of N residues, an
        alpha-helices (alpha) and b-sheets (beta) list is generated. Each list
        contains N elements, which are ordered as the lines in the parsed file (i.e. element 0
        corresponds to the first "useful" line of the file, element 1 to the second and so on).
        Each element of the list has '' as a value if the given residue wasn't in the secondary
        structure of that list (e.g. if residue 12 is '' in alpha, it means that it wasn't found
        in an alpha-helix), or a value otherwise. All the elements having the same value (except '')
        in one of the alpha or the beta lists belong to the same structural element (e.g. all the
        elements with "3" in the alpha list belong to the same alpha-helix.)

        Both the "ah" and the "be" lists contain k alpha and beta lists, respectively, where k is the
        number of secondary structure definitions found in the parsed file. For dssp this will always
        be one, while it may be one or more for pdss.
        
        """
        mx=[]
        dd=[]
        b_sheet=[]
        b_sheet_out=[]
        dssp = open(filen,'r')

        p = dssp.readlines()

        dels=[]   

        if p[0].startswith('='):
            mx.append([])
            b_sheet.append([])
            for lines in p:
                if not ((lines.strip().endswith('.') or lines.strip().startswith('#') or '*!' in lines)):
                    l = lines[16]
                    if l == 'E':
                        mx[-1].append('')
                        b_sheet[-1].append(lines.rstrip()[33])
                    elif l == 'H':
                        mx[-1].append('H')
                        b_sheet[-1].append('')
                    else:
                        mx[-1].append('')
                        b_sheet[-1].append('')
            complete=True    
        else:
            for line in p:
                if line.strip().startswith('#'):
                    dels.append(line)
            for deletion in dels:
                p.remove(deletion)

            k=len(p[0].strip().split())-1
            for j in range(0,k/2):
                mx.append([])
                b_sheet.append([])
            for lines in p:
                if not lines.startswith('#'):
                    words = lines.strip().split()
                for j in range(0,k/2):
                    l=words[1+j*2]
                    if l == 'E':
                        mx[j].append('')
                        b_sheet[j].append('B')
                    elif l == 'H':
                        mx[j].append('H')
                        b_sheet[j].append('')
                    else:
                        mx[j].append('')
                        b_sheet[j].append('')
            complete=False

        for i in range(len(mx)):
            b_sheet_out.append([])
 
        #if complete:
        b_sheet_out=b_sheet
        return {'ah':mx,'be':b_sheet_out,'complete':complete}

    def parse_ss_from_fname(self, fname):
        if not fname:
            return None
        try: 
            open(fname,'r')
        except:
            tkMessageBox.showerror('Error',"Could not open the secondary structure file.")
            self.ss = ''            
            return None            
        return self.parse_ss(fname)
        


    def compute_ss(self, reference, fromscratch=False):
        """
        Generate secondary structure definition for a given PyMol object
        starting from the built-in object definition (default behaviour) or 
        using util.ss (fromscratch=True). The default behaviour could be
        particularly useful after using the DSSP plugin to redefine
        secondary structures.
        
        reference:   name of the object from which ss has to be extracted (usually the reference)
        fromscratch: switch for chosing ss generation mode 
        """
        ss={"ah":[[]],"be":[[]],"complete":False}        
        if fromscratch:            
            cmd.util.ss(reference)
        cas = self.services.provideServiceByName("get_calphas")(reference)
        for ca in cas:
            ss["be"][0].append('B' if ca.ss=='S' else '')
            ss["ah"][0].append('H' if ca.ss=='H' else '')
        return ss

    def get_ss(self, reference, ss_fname=None, fromscratch=False):
        """
        Secondary structure filter
        #doalpha, dobeta,betaco,ss
        dm:      atomic distance matrix
        metass:  secondary structure as parsed by parse_ss
        idx:     number of reference ss in metass
        dobeta:  beta filtering switch (True = On)
        doalpha: alpha filtering switch (True = On)
        betaco:  distance cut-off for beta-sheets identification
        """
        #ss={}
        #ss["complete"]=metass["complete"]
        #ss["ah"]=metass["ah"][idx]
        #ss["be"]=metass["be"][idx]
        if fromscratch and ss_fname:
            raise "Error"
        if self.ss_mode == 'file':
            return self.parse_ss(ss_fname)
        return self.compute_ss(reference, fromscratch=fromscratch)

    def get_filter(self, data, reference, *args, **kwargs):
        # get from GUI
        self.betaco = float(self.cof.get())
        if self.ss_mode == 'file':
            ss_choice = self.cho_str.curselection()
            if len(ss_choice) > 0:            
                ss_sele = int(self.cho_str.curselection()[0])
            else:
                ss_sele=0
            self.ss = self.get_ss(reference, ss_fname=self.ss_fname)
            ss = {'ah':self.ss['ah'][ss_sele],'be':self.ss['be'][ss_sele],'complete':self.ss['complete']}
        elif self.ss_mode == 'pdb': 
            self.ss = self.get_ss(reference)
            ss = {'ah':self.ss['ah'][0],'be':self.ss['be'][0],'complete':self.ss['complete']}
        elif self.ss_mode == 'pymol':
            self.ss = self.get_ss(reference, fromscratch=True)
            ss = {'ah':self.ss['ah'][0],'be':self.ss['be'][0],'complete':self.ss['complete']}
         
        dm = self.services.provideServiceByName("distance_matrix")(reference)
        return self.filter(reference, ss, dm, self.betaco, self.doalpha, self.dobeta)
    
            
            
            
        #self.ss = self.get_ss(reference)
        #self.dm = # matrix loading through service?
        #self.doalpha = ..
        #self.dobeta = ..
        #self.betaco = .. # get from GUI
        #self.idx = .. 
        #return filter(self, self.ss, self.betaco, self.doalpha, self.dobeta)      
        #return True            
                   
    def filter(self, reference, ss, dm, betaco, doalpha, dobeta):
    # doalpha 
    # dobeta
    # chainz
        cas = self.services.provideServiceByName("get_calphas")(reference)
        lencas=len(cas)
        chains = [ i.chain for i in cas ]
	if dobeta:
            if ss["complete"]:
                betas = np.array(ss['be'])
                sheets = list(set(ss['be']))
                filter_matrix_be = np.ones((lencas,lencas), dtype=np.bool)
                try:
                    sheets.remove('')
                except:
                    pass
                try:
                    sheets.remove(' ')
                except:
                    pass
                for i in sheets:
                    beta_positions = np.where(betas == i)[0]
                    for j in beta_positions:
                        for k in beta_positions:
                            filter_matrix_be[j][k] = False
                            filter_matrix_be[k][j] = False

            else:
                k=0
	        betaposs=[]
                for j in range(len(ss['be'])):
                    if ss['be'][j] == 'B':
                        betaposs.append(k)
                    elif j!=0:
                        if ss['be'][j-1] == '':
                            k=k+1
                        betaposs.append('')
                    else:
                        betaposs.append('')
                betas = np.array(betaposs)
                sheets = list(set(betaposs))
                try:
                    sheets.remove('')
                except:
                    pass
                try:
                    sheets.remove(' ')
                except:
                    pass

                beta_residues=np.array([],dtype=np.int)
                filter_matrix_pdssp = np.ones((lencas,lencas), dtype=np.bool)
                for i in sheets:                    
                    beta_positions = np.where(betas == str(i))[0]
                    
                    beta_residues = np.concatenate((beta_residues,beta_positions))
                    for j in beta_positions:
                        for k in beta_positions:
                            filter_matrix_pdssp[int(j)][int(k)] = False
                            filter_matrix_pdssp[int(k)][int(j)] = False

                dm_graph_matrix = dm < betaco
                #dm_graph_matrix = np.ma.masked_where(dm<betaco,dm).mask
                beta_graph_matrix = np.zeros((lencas,lencas), dtype=np.bool)
		for i in beta_residues:
                    for j in beta_residues:
                        beta_graph_matrix[i][j] = True
			beta_graph_matrix[j][i] = True
                G = nx.Graph(beta_graph_matrix & dm_graph_matrix)
                filter_matrix_bsheets = np.ones((lencas,lencas), dtype=np.bool)
                bsheets = nx.connected_components(G)
                for i in bsheets:
                    if len(i) > 1:
                        for j in i:
                            for k in i:
                                filter_matrix_bsheets[j][k] = False
                                filter_matrix_bsheets[k][j] = False

                    #nx.algorithms.traversal.breadth_first_search.bfs_edges(G, i)                    
                    #if np.any(beta_residues == i):
                    #    filter_matrix_dist[i]=True
                    #    filter_matrix_dist[:,i]=True
                filter_matrix_be = filter_matrix_pdssp & filter_matrix_bsheets
        else:
            filter_matrix_be = np.ones((lencas,lencas), dtype=np.bool)
            
        if doalpha:
            dd = []
            k = 0
            for j in range(len(ss['ah'])):
                if ss['ah'][j] == 'H':
                    dd.append(k)
                elif j!=0:
                    if dd[j-1] != 'H':
                        k=k+1
                    dd.append('')
                else:
                    dd.append('')
            filter_matrix_ah = np.ones((lencas,lencas), dtype=np.bool)
            alphas = np.array(dd)
            helices = list(set(alphas))
            try:
                helices.remove('')
            except:
                pass
            for i in range(len(helices)):
                alpha_positions = np.where(alphas == helices[i])[0]
                for j in alpha_positions:
                    for k in alpha_positions:
                        filter_matrix_ah[j][k] = False
        else:
            filter_matrix_ah=np.ones((lencas,lencas), dtype=np.bool)

        chain_total=0
        chain_ids=[]
        chain_filter=np.ones((lencas,lencas), dtype=np.bool)
        chains=[ i.chain for i in cas ]
        for chain in chains: 
            if chain in chain_ids:
                continue 
            chain_ids.append(chain)
        j=1
        for i in chain_ids:
            j=chains.count(i)
            chain_filter[chain_total:chain_total+j,chain_total:chain_total+j]=False
            chain_total+=j
        
        return filter_matrix_ah & filter_matrix_be | chain_filter

# --------------

    def ss_selection(self,var):
        if var == "alpha":
            self.doalpha=True
            self.dobeta=False
        elif var == "beta":
            self.doalpha=False
            self.dobeta=True
        elif var == "both":
            self.doalpha=True
            self.dobeta=True

    def mode_selection(self,var):
        if var == 'read from file':
            self.ss_mode='file'
        elif var == 'read from loaded structure':
            self.ss_mode='pdb'
        elif var == 'calculate with pymol':
            self.ss_mode='pymol'            

            
    def GUI(self, page):
        """
        Add the specific GUI elements for this plugin to the provided Tk widget.
        """
        def select_file():
            self.cho_list.set("")
            self.ss_fname=tkFileDialog.askopenfilename(title="DSSP/PDSSP selection")
            #label["text"]= str(self.ss_fname)+" selected"
	    try:
                tmp_ss=self.parse_ss_from_fname(self.ss_fname)
	    except:
                tkMessageBox.showerror('Error','An error occured while parsing your secondary structure file.')
		self.ss = None
                return
            if tmp_ss:
                self.ss = tmp_ss
                ss_len=len(tmp_ss['ah'])
                for t in range(ss_len,0,-1):
                    self.cho_str.insert(0, t)
                self.cho_str.selection_set(0)

            
            
        group = Pmw.Group(page,tag_text='Secondary Structure')
        group.pack(padx=40,pady=10)
        frame = Frame(group.interior())
        frame.pack(side=TOP)

        self.ss_modef= Pmw.RadioSelect(frame,
        buttontype = 'radiobutton',
        orient = 'vertical',
        labelpos = 'n',
        command = self.mode_selection,
       # label_wtext = 'Use a ss',
       # hull_borderwidth = 2,
       # hull_relief = 'ridge'
        )
        self.ss_modef.pack(side = 'left', expand = 1, padx = 10, pady = 10)
        self.ss_modef.add('read from file')
        self.ss_modef.add('read from loaded structure')
        self.ss_modef.add('calculate with pymol')
        self.ss_modef.setvalue('read from file')



        self.ss_state = IntVar()
        self.cho_list = StringVar()
#       self.cho_list.set("1 2 3 4")
  #      filelist=["dssp","pdssp"]
 
        self.loadssbutton = Button(frame, text="Select file", command=select_file)
        self.loadssbutton.pack(side=LEFT)

        self.cho_str=Listbox(frame,selectmode=BROWSE,listvariable=self.cho_list,height=2)
        self.scroll_ss=Scrollbar(frame, orient=VERTICAL,command=self.cho_str.yview)
        self.scroll_ss.pack(side=LEFT, fill=Y)
        self.cho_str.pack(side=LEFT, fill=BOTH, expand=1)
        #load_ss = Button(frame, text="Load the selected\n structure", command=load_structure, state=DISABLED)
        #load_ss.pack(side=LEFT)
        

        self.ss_type= Pmw.RadioSelect(frame,
        buttontype = 'radiobutton',
        orient = 'vertical',
        labelpos = 'n',
        command = self.ss_selection,
       # label_wtext = 'Use a ss',
       # hull_borderwidth = 2,
       # hull_relief = 'ridge'
        )
        self.ss_type.pack(side = 'left', expand = 1, padx = 10, pady = 10)
        self.ss_type.add('alpha')
        self.ss_type.add('beta')
        self.ss_type.add('both')
        self.ss_type.setvalue('both')


            # self.ss_type.config(Button_state = "active"
                # self.ss_type.config(Button_state = ACTIVE)
        
        self.ss_state.set(0)
        self.balloon=Pmw.Balloon(frame)
        self.balloon.bind(self.loadssbutton, 'Select a DSSP or PDSSP file \n which contains many structures as you wish\n and let\'s hope for the best')
        label = Label(frame, text=" ")
        label.pack(side=BOTTOM)
        self.prova=Frame(group.interior())
        self.prova.pack(side=BOTTOM)
        
        self.cof = self.makeDoubleVar(self.prova, "Distance cut-off (nm)", self.change_co,
                                      tooltip="""\
Only correlations between residues with
sequence distance lower/higher than this cutoff
will be plotted.""")
        self.cof.set(self.betaco)
        return group



    def change_co(self, a, *args):
        if len(args)!=0: 
            a=args[0]
        val = int(self.cof.get()) + int(a) 
        self.cof.set(val)
    
        
        # input file

#----
    def change_co(self, a, *args):
        if len(args)!=0: 
            a=args[0]
        val = float(self.cof.get()) + 0.05*float(a)
        self.cof.set(val)
        
    def validate_file(self, s):
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
        
        if self.ss_mode == 'file':
            if not self.ss_fname:
                return (False,"Load a ss file first")
            if not self.ss:
                return (False,"There was a problem parseing your ss file")
            ss_choice = self.cho_str.curselection()
            if len(ss_choice) > 0:
                ss_sele = int(ss_choice[0])
            else:
                ss_sele = 0
            reference = self.services.provideServiceByName("get_reference")()
            cas = self.services.provideServiceByName("get_idxs")(reference)
            if len(cas) != len(self.ss['ah'][ss_sele]):
                return (False, "Lengths of loaded ss and of reference object to not match")
            return (True,"ss loaded from file")
        if self.ss_mode == 'pdb':
            return (True,"ss loaded from pdb")
        if self.ss_mode == 'pymol':
            return (True,"ss calculated by PyMol")
        return (True, "")
