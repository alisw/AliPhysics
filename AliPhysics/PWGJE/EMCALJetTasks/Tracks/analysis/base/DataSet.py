#**************************************************************************
#* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
#*                                                                        *
#* Author: The ALICE Off-line Project.                                    *
#* Contributors are mentioned in the code where appropriate.              *
#*                                                                        *
#* Permission to use, copy, modify and distribute this software and its   *
#* documentation strictly for non-commercial purposes is hereby granted   *
#* without fee, provided that the above copyright notice appears in all   *
#* copies and that both the copyright notice and this permission notice   *
#* appear in the supporting documentation. The authors make no claims     *
#* about the suitability of this software for any purpose. It is          *
#* provided "as is" without express or implied warranty.                  *
#**************************************************************************
from ROOT import TIter, TList
from copy import copy,deepcopy
from PWGJE.EMCALJetTasks.Tracks.analysis.base.MergeException import MergeException
from PWGJE.EMCALJetTasks.Tracks.analysis.base.JetContainer import JetContainer
from PWGJE.EMCALJetTasks.Tracks.analysis.base.SpectrumContainer import TrackContainer, ClusterContainer

class DataSet(object):
    """
    Data set for a given trigger class. A data set contains a set of cluster containers and an set of track containers
    """
    
    class ContentException(Exception):
        
        def __init__(self, searchkey, container):
            self.__searchkey = searchkey
            self.__container = container
            
        def __str__(self):
            return "%s already present in container %s" %(self.__searchkey, self.__container)
    
    def __init__(self):
        """
        Constructor
        """
        self.__trackContainers = {}
        self.__clusterContainers = {}
        self.__jetContainer = JetContainer()
        
    def __copy__(self):
        """
        shallow copy constructor
        """
        print "Simple copy called from %s" %(self.__class__)
        newobject = DataSet()
        for name,tc in self.__trackContainers.iteritems():
            newobject.AddTrackContainer(name, copy(tc))
        for name,cc in self.__clusterContainers.iteritems():
            newobject.AddClusterContainer(name, copy(cc))
        return newobject
         
    def __deepcopy__(self, memo):
        """
        deep copy constructor
        """
        print "deep copy called from %s" %(self.__class__)
        newobject = DataSet()
        for name,tc in self.__trackContainers.iteritems():
            newobject.AddTrackContainer(name, deepcopy(tc, memo))
        for name,cc in self.__clusterContainers.iteritems():
            newobject.AddClusterContainer(name, deepcopy(cc, memo))
        return newobject
        
    def AddTrackContainer(self, name, data):
        """
        Add a new track container to the dataset
        """
        if name in self.__trackContainers.keys():
            raise DataSet.ContentException(name, "TrackContainer")
        self.__trackContainers[name] = data
        
    def AddJetSpectrum(self, spectrum, jetpt, isMCkine):
        self.__jetContainer.SetJetPtHist(jetpt, spectrum, isMCkine)
        
    def AddEventHistForJets(self, hist):
        self.__jetContainer.SetEventHist(hist)
        
    def AddClusterContainer(self, name, data):
        """
        Add a new cluster container to the dataset
        """
        if name in self.__clusterContainers.keys():
            raise DataSet.ContentException(name, "ClusterContainer")
        self.__clusterContainers[name] = data
        
    def FindTrackContainer(self, name):
        """
        Find a track container within the dataset
        """
        if not name in self.__trackContainers.keys():
            return None
        return self.__trackContainers[name]

    def FindClusterContainer(self, name):
        """
        Find a cluster container within the dataset
        """
        if not name in self.__clusterContainers:
            return None
        return self.__clusterContainers[name]
    
    def GetJetContainer(self):
        return self.__jetContainer
    
    def GetListOfTrackContainers(self):
        """
        Get a list of track container names
        """
        return self.__trackContainers.keys()
    
    def GetListOfClusterContainers(self):
        """
        Get a list of cluster container names
        """
        return self.__clusterContainers.keys()
    
    def Add(self, other):
        """
        Add other data set to this one
        """
        if not isinstance(other, DataSet):
            raise MergeException("Incompatible types: this(Dataset), other(%s)" %(str(other.__class__)))
        nfailure = 0
        for cont in self.GetListOfTrackContainers():
            othercont = other.FindTrackContainer(cont)
            if othercont:
                self.__trackContainers[cont].Add(othercont)
            else:
                nfailure += 1
        for cont in self.GetListOfClusterContainers():
            othercont = other.FindClusterContainer(cont)
            if othercont:
                self.__clusterContainers[cont].Add(othercont)
            else:
                nfailure += 1
        if nfailure > 0:
            raise MergeException("Several containers have not been found inside the other datase")
         
    def Scale(self, scalefactor):
        """
        Scale all track or cluster containers with the underlying scale factor
        """
        for cont in self.__trackContainers.values():
            cont.Scale(scalefactor)
        for cont in self.__clusterContainers.values():
            cont.Scale(scalefactor)
     
    def GetRootPrimitive(self, listname):
        """
        Make root primitives (for root IO)
        """
        result = TList()
        result.SetName(listname)
        tracklist = TList()
        tracklist.SetName("trackContainers")
        for name,tc in self.__trackContainers.iteritems():
            tracklist.Add(tc.GetRootPrimitive("trackcontainer_%s" %(name)))
        clusterlist = TList()
        clusterlist.SetName("clusterContainers")
        for name,cc in self.__clusterContainers.iteritems():
            clusterlist.Add(cc.GetRootPrimitive("clustercontainer_%s" %(name)))
        result.Add(tracklist)
        result.Add(clusterlist)
        return result
     
    @staticmethod
    def BuildFromRootPrimitive(rootprimitive):
        """
        Build dataset from a root primitive
        """
        result = DataSet()
        trackContainers = rootprimitive.FindObject("trackContainers")
        clusterContainers = rootprimitive.FindObject("clusterContainers")
        trackIter = TIter(trackContainers)
        clusterIter = TIter(clusterContainers)
        currententry = trackIter.Next()
        while currententry:
            entryname = currententry.GetName()
            if "trackcontainer" in entryname:
                contname = entryname.replace("trackcontainer_","")
                result.AddTrackContainer(contname, TrackContainer.BuildFromRootPrimitive(currententry))
            currententry = trackIter.Next()
        currententry = clusterIter.Next()
        while currententry:
            entryname = currententry.GetName()
            if "clustercontainer" in entryname:
                contname = entryname.replace("clustercontainer_","")
                result.AddClusterContainer(contname, ClusterContainer.BuildFromRootPrimitive(currententry))
            currententry = clusterIter.Next()
        return result
     
    def Print(self):
        """
        Print content of the data set
        """
        print "Dataset content:"
        print "============================================"
        print "   Track Containers:"
        for cont in self.__trackContainers.keys():
            print "      %s" %(cont)
        print "   Cluster Containers:"
        for cont in self.__clusterContainers.keys():
            print "      %s" %(cont)
        print "--------------------------------------------"
        print "Status of the different containers:"
        for contname, container in self.__trackContainers.iteritems():
            print "   %s:" %(contname)
            container.Print()
        for contname, container in self.__clusterContainers.iteritems():
            print "   %s:" %(cont)
            container.Print()

 