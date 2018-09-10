#!/usr/bin/env python
#  File: cPluginManager.py 
#
#  Copyright (C) 2011 Marco Pasi <mf.pasi@gmail.com> 
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
# "cPluginManager.py v0.1 (C) 2011 Marco Pasi"
import logging
import ConfigParser
import Pmw
from distutils.version import StrictVersion
from yapsy.VersionedPluginManager import VersionedPluginInfo
from yapsy.IPlugin import IPlugin
from yapsy.PluginManagerDecorator import PluginManagerDecorator
import networkx as nx
import time

# -------------------
class ServicePluginInfo(VersionedPluginInfo):
    """
    Gather some info about a plugin such as its name, author,
    description...
    """
    
    def __init__(self, plugin_name, plugin_path):
        """
        Set the name and path of the plugin as well as the default
        values for other usefull variables.
        """
        VersionedPluginInfo.__init__(self, plugin_name, plugin_path)
        # dict of dependencies: { capability_name: ( version or None ) }
        self.depends    = {}
        # dict of depending plugins: { capability_name: provider_plugin_name } 
        # It is filled at check time
        self.depends_plugins = {}        
        # dict of provides: { capability_name: method_name }
        self.provides   = {}
        # dict of 
        self.depending_plugins = {} 

    def addDepends(self, capability, vstring=""):
        if vstring == "":
            version = None
        else:
            version = StrictVersion(vstring)
        self.depends[capability] = version

    def addProvides(self, capability, method_name):
        self.provides[capability] = method_name

# -------------------
class ServicePluginManager(PluginManagerDecorator):
    """
    A plugin manager that handles the possibility for plugins to
    provide services, and to use services provided by other plugins.

    The plugins must define 2 new sections in their configuration file:

      [Provides]
      service_name = property_name

    where the provided services are listed along with the corresponding
    property of the providing plugin; and

      [Depends]
      service_name = [ required_version ]

    where the dependencies of the plugin are listed, with optionally
    the required version of the service.
    
    """
    
    def __init__(self, 
                 decorated_manager=None,
                 categories_filter={"Default":IPlugin}, 
                 directories_list=None, 
                 plugin_info_ext="yapsy-plugin"):
        """
        """
        # Create the base decorator class
        PluginManagerDecorator.__init__(self,decorated_manager,
                                        categories_filter,
                                        directories_list,
                                        plugin_info_ext)
        self.setPluginInfoClass(ServicePluginInfo)
        # prepare the Services class and the dependency list.
        self._provides = {} # { capability : version } (kept for checking purposes)
        self._unresolved = []
        self.services = VersionedServices()

    def addServicePluginInfo(self, candidate_infofile, plugin_info):
        """
        Add information on dependencies/provides from
        the candidate plugins, and return the modified plugin_info.
        """
        # we need to parse the configuration file
        config_parser = ConfigParser.SafeConfigParser()
        try:
            config_parser.read(candidate_infofile)
        except:
            logging.debug("Could not parse the plugin file %s" % candidate_infofile)					
            return {}
        
        # collect additional information
        if config_parser.has_section("Provides"):
            plugin_info.provides.update(config_parser.items("Provides"))
        if config_parser.has_section("Depends"):
            for capability, vstring in config_parser.items("Depends"):
                plugin_info.addDepends(capability, vstring)

        return plugin_info

    def locatePlugins(self):
        """
        Once the candidates have been defined, set up and solve the
        dependency chain, and place plugins that cannot be loaded in
        the *_unresolved* group.
        """
	def allprede(G,i,nodes):
            """
            Recursively find all predecessors of a given node of a directed nx graph
            up to root(s). Also works with cycles.
            """
            tmp=G.predecessors(i)
            for j in tmp:
                if j not in nodes:
                    nodes.append(j)
                    allprede(G,j,nodes)
        # first let the decorated manager do the locating
        self._component.locatePlugins()
        # create dependencies tree graph 
        depG=nx.MultiDiGraph()
        
        # Add services provided directly to the Services class
        # self._provides.update(self.services._list)
        # logging.debug("PluginManager gathered %d services: "%(len(self._provides))+(",".join(self._provides.keys())))
        
        pcs=self.getPluginCandidates()
        removenda=[]
        
        for finfo, path, info in pcs:
            self.addServicePluginInfo(finfo, info) # Add additional informations to the plugin infos
            depG.add_node((finfo,path,info)) # Add the plugin to the dep tree

        for finfo1, path1, info1 in pcs:
            for dep1, ver1 in info1.depends.iteritems():
                if dep1 in self.services._list:
                    logging.debug("Capability %s, which is needed for plugin %s, is provided by core cPlot functionalities."%(dep1,info1.name))
                    continue
                providers=[[],[]]
                for finfo2, path2, info2 in pcs:
                    if path1 == path2:
                        continue
                    if dep1 in info2.provides.keys():
                        providers[0].append((finfo2,path2,info2))
                        providers[1].append(info2.version)
                if len(providers[0]) == 0: 
                    removenda.append((finfo1, path1, info1))  
                    logging.debug("A dependency for plugin %s cannot be satisfied: capability %s could not be found"%(info1.name,dep1)) 
                    continue
                maxver=max(providers[1])
                if ver1:
                    if maxver < ver1:
                        removenda.append((finfo1, path1, info1))
                        logging.debug("Dependency %s for plugin %s cannot be satisfied: version %s required, while %s was available"%(dep1,info1.name,ver1,maxver))
                        continue    
                # if a dep is not satisfied, keep track the incriminated plugin; it will be removed later
                depG.add_edge((finfo1,path1,info1),providers[0][providers[1].index(maxver)]) # only connect as dep the capability with the highest version. Disregard the others.     
                info1.depends_plugins[dep1] = providers[0][providers[1].index(maxver)][2].name                
                
        unmetdeps=len(removenda)        
                                             
        for finfo, path, info in list(set(removenda)):
            deldeps=[]
            allprede(depG,(finfo,path,info),deldeps)
            removenda.extend(deldeps) 
        if len(removenda) == 0:
            logging.debug("All the plugins were accepted. Hell yeah!")
        else:
            if len(removenda) != unmetdeps:
                logging.debug("%d plugins had to be removed due to additional unment capabilities requirements"%(len(removenda)-unmetdeps))
            for i in removenda:
                self.removePluginCandidate(i)
                self._unresolved.append(i)
            logging.debug("\n-----\nThe following plugins have been removed due to unmet dependencies: "+", ".join(x[2].name for x in removenda))

        for finfo, path, info in self._candidates:
            # add information about provides/depdends
            # update self._provides            
            for capability in info.provides.keys():
                if self._provides.has_key(capability):
                    if info.version <= self._provides[capability]:
                        logging.debug("Capability %s of %s is already provided at a newer version (new: %s <= current: %s)"%(capability, info.name, info.version, self._provides[capability]))
                        continue
                    else:
                        logging.debug("Capability %s of %s overrides present copy (new: %s > current: %s)"%(capability, info.name, info.version, self._provides[capability]))                        
                self._provides[capability] = info.version
   
        logging.debug("PluginManager gathered %d services: "%(len(self._provides))+(",".join(self._provides.keys())))
        self._provides.update(self.services._list)
        return len(self._candidates)

    
    def loadPlugins(self, callback=None):
        """
        After loading the plugins, register the provided capabilities
        in the Services() class.
        """
        # first let the decorated manager do the loading
        self._component.loadPlugins(callback)
        
        # then, register all capabilities to self.services 
        # and provide the plugins with a reference to self.services
        for categ in self.getCategories():
            allPlugins = self.getPluginsOfCategory(categ)
            for plugin in allPlugins:
                for capability, property_name in plugin.provides.iteritems():
                    # some obvious considerations
                    assert capability in self._provides.keys()
                    # assert self._provides[capability] == plugin.version
                    # one less obvious one
                    assert hasattr(plugin.plugin_object, property_name)
                    
                    # register the service
                    self.registerService(getattr(plugin.plugin_object, property_name),
                                         capability,
                                         version = plugin.version)
                plugin.plugin_object.setServices(self.services)
        return True

    def registerService(self, property, capability, version = None):
        """
        Provide a utility function to register a capability without
        handling the Services class directly.
        """
        return self.services.registerService(property, capability, version)

    def checkPluginsStatus(self):
        """
        Return a dictionary of the plugin statuses
        """
        outdict = {}
        allPlugins = self.getPluginsOfCategory("Filters")
        for plugin in allPlugins:
            outdict[plugin.name] = plugin.plugin_object.check()
        return outdict

    def checkPluginStatus(self, plugin_name):
        """
        Return the status of a plugin
        """
        plugin = self.getPluginByName(plugin_name, category="Filters")
        status = (False, "")
        if plugin:
            status = plugin.plugin_object.check()
        return status
        
# -------------------
class Service(object):
    """
    A single service.
    """

    def __init__(self, property, capability):
        self.property = property
        self.name = capability

    def __repr__(self):
        return "Service %s->%s"%(self.name, self.property)
        
# -------------------
class Services(object):
    """
    The container for all registered capabilities.
    """

    service_class = Service
    
    def __init__(self):
        self._list  = {}  # { capability: service }, for latest version

    def registerService(self, property, capability):
        """
        Register a service, that can be accessed by its name *capability*,
        and consists in the *property* (an object or method).  
        """
        service = self.service_class(property, capability)
        
        if self._list.has_key(capability):
            return None
        self._list[capability] = service
        return service

    def provideServiceByName(self, name):
        """
        Provide a service by capability name.
        """
        if not self._list.has_key(name):
            raise AttributeError
        return self._list[name].property
        
    def provideServiceByProperty(self, property_name):
        """
        Should we provide this quick & dirty way to access services?
        """
        # check in latest versions
        for capability, service in self._list.iteritems():
            if service.property.__name__ == property_name:
                return self.provideServiceByName(capability)
        logging.debug("No service providing property %s was found"%property_name)
        raise AttributeError

    def __getattr__(self, name):        # fallback method
        try:
            return self.provideServiceByProperty(name)
        except AttributeError:
            logging.debug("Couldn't provide service %s."%name)
            return None

# -------------------
class VersionedService(Service):
    """
    A service with StrictVersion definition
    """

    def __init__(self, property, capability, version=None):
        Service.__init__(self, property, capability)
        if not version:
            version = StrictVersion("0.0")
        self.version = version
        
    def __repr__(self):
        return Service.__repr__(self)+"@v%s"%self.version

# -------------------
class VersionedServices(Services):
    """
    Services that handle versions.  Older versions of a same capability are kept in the attic.
    """

    service_class = VersionedService
    
    def __init__(self):
        Services.__init__(self)
        self._attic = {}  # { capability: [service, service, ...] }, for old versions

    def registerService(self, property, capability, version=None):
        """
        Register a new service, or another version of a known service.  The newest version
        of the service is placed in *self._list*, while all older versions are kept in
        *self._attic*.
        """
        service = self.service_class(property, capability, version)
        
        if self._list.has_key(capability):
            # already got this capability
            if self._list[capability].version == version:
                # same version: keep the one we have
                logging.debug("The capability %s is already registered at version %s.  Not registering."%(capability, version))
                return None
            elif self._list[capability].version > version:
                # older: add this version to attic
                logging.debug("The capability %s is already registered at version %s.  Adding older version %s."%(capability, self._list[capability].version, version))
                self._attic[capability].append(service)
                return service
            else:
                # newer: put current version in attic
                logging.debug("The capability %s is already registered at version %s.  Adding newer version %s [default]."%(capability, self._list[capability].version, version))
                self._attic[capability].append(self._list[capability])
        else:
            # new capability: prepare attic
            self._attic[capability] = []

        self._list[capability] = service
        return service

    def provideServiceByName(self, name, version=None):
        """
        Provide a service by capability name.  The newest version is used
        if none specified.  Otherwise the right version is searched in the
        attic.
        """
        if not self._list.has_key(name):
            raise AttributeError

        if not version or self._list[name].version == version:
            # provide latest version
            return self._list[name].property
        else:
            # look for the right version
            available_versions = []
            for service in self._attic[name]:
                available_versions.append(service.version)
                if service.version == version:
                    return service.property
            logging.debug("Couldn't provide service %s at version %s (only have %s)"%
                          (name, version, ",".join(available_versions)))
            raise AttributeError
        
    def provideServiceByProperty(self, property_name):
        """
        Dirty, but kept for backward compatibility with the Services class.
        """
        try:
            # check in latest versions
            return Services.provideServiceByProperty(self, property_name)
        except AttributeError:
            # check in attic
            for capability, services in self._attic.iteritems():
                for service in services:
                    if service.property.__name__ == property_name:
                        return self.provideServiceByName(capability, version=service.version)

        logging.debug("No service providing property %s was found in the attic"%property_name)
        raise AttributeError
