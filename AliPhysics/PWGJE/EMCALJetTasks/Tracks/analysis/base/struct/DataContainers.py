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
from PWGJE.EMCALJetTasks.Tracks.analysis.base.Helper import NormaliseBinWidth
from PWGJE.EMCALJetTasks.Tracks.analysis.base.struct.ClusterTHnSparse import ClusterTHnSparseNew, ClusterTHnSparseOld
from PWGJE.EMCALJetTasks.Tracks.analysis.base.struct.TrackTHnSparse import TrackTHnSparseNew, TrackTHnSparseOld
from PWGJE.EMCALJetTasks.Tracks.analysis.base.MergeException import MergeException
from copy import copy,deepcopy
from ROOT import TList
 
class DataContainer(object):
    """
    Data container (event and spectum container)
    """
    
    class DataException(Exception):
        """
        Exception handling incomplete data
        """
        def __init__(self, raiseobject):
            """
            constructor
            """
            self.__object = raiseobject
            
        def __str__(self):
            """
            Produce error string
            """
            return "Container %s missing" 
    
    
    def __init__(self, eventHist = None, dataHist = None, histotype=None):
        """
        Construct spectrum container
        """
        self._events = eventHist
        self._spectrum = dataHist
        self.__doNormBW = True
        self._datahistname = "DataHist"
        
    def __copy__(self):
        """
        Shallow copy constructor
        """
        print "Simple copy called from %s" %(self.__class__)
        newobject = DataContainer(self.__events, self._spectrum, None)
        newobject._Copy(self)
        return newobject
             
    def __deepcopy__(self, memo):
        """
        Deep copy constructor
        """
        print "deep copy called from %s" %(self.__class__)
        newobject = DataContainer(None, None, None)
        newobject._Deepcopy(self, memo)
        return newobject
    
    def GetValues(self):
        """
        For copy constructors
        """
        return {"dn": self.__doNormBW, "hn": self._datahistname}
    
    def __SetValues(self, values):
        """
        For copy constructors
        """
        self.__doNormBW = values["dn"]
        self._datahistname = values["hn"]
                        
    def SetEventHist(self, eventHist):
        """
        Set the event hist
        """
        self._events = eventHist
        
    def SetSpectrumContainer(self, cont):
        """
        Initialise spectrum with container
        """
        self._spectrum = cont
        
    def GetEventHist(self):
        """
        Access event counter histogram
        """
        return self._events
    
    def GetSpectrumContainer(self):
        """
        Access underlying spectrum container
        """
        return self.__spectrum
    
    # Property definitions
    EventHist = property(GetEventHist, fset=SetEventHist)
    SpectrumHist = property(GetSpectrumContainer, fset=SetSpectrumContainer)

    
    def _Copy(self, other):
        """
        Make a shallow copy of the other object into this
        """
        evhist = other.EventHist
        speccont = other.SpectrumHist
        if evhist:
            self.EventHist = copy(evhist)
        if speccont:
            self.SpectrumHist = copy(speccont)
        self.__SetValues(other.GetValues())
     
    def _Deepcopy(self, other, memo):
        """
        Make a deep copy of the other object into this
        """
        evhist = other.EventHist
        speccont = other.SpectrumHist
        if evhist:
            self.EventHist = deepcopy(evhist, memo)
        if speccont:
            self.SpectrumHist = deepcopy(speccont, memo)
        self.__SetValues(other.GetValues())
         
    def Add(self, other):
        """
        Add other data container to this data container
        """
        if not isinstance(other, DataContainer):
            raise MergeException("Incompatible types: this(DataContainer), other(%s)" %(str(other.__class__)))
        if self._events and other.EventHist:
            self._events.Add(other.EventHist)
        else:
            if other.EventHist:
                self._events = deepcopy(other.EventHist)
        if self._spectrum and other.SpectrumHist:
            self._spectrum.Add(other.SpectrumHist)
        else:
            if other.SpectrumHist:
                self._spectrum = other.SpectrumHist
         
    def Scale(self, scalefactor):
        """
        Scale the underlying spectrum container with the scale factor
        """
        self.__spectrum.Scale(scalefactor)
         
    def GetRootPrimitive(self, name):
        """
        Convert object to a root primitive (so that it can be wirtten to a rootfile)
        """
        result = TList()
        result.SetName(name)
        events = self._events.GetROOTHisto()
        events.SetName("events")
        result.Add(events)
        spectrum = self.__spectrum.GetHistogram()
        spectrum.SetName("spectrum")
        result.Add(spectrum)
        return result
     
    @staticmethod
    def GetHistoPair(rootprimitive):
        """
        Convert to TList to python dictionary
        """
        return {"events":rootprimitive.FindObject("events"), "spectrum":rootprimitive.FindObject("spectrum")}
        
    def GetEventCount(self):
        """
        Get the number of selected events
        """
        return self._events.GetEventCount()
            
    def MakeProjection(self, dim, histname = None, xtitle = None, ytitle = None, doNorm = True):
        """
        Make event-normalised projection to 1D
        """
        if not self._spectrum:
            raise DataContainer.DataException(self._datahistname)
        if not self._events:
            raise DataContainer.DataException("EventHist")
        
        if not histname:
            histname = "%sProjected" %(self._spectrum.GetHistogramName())
        if not xtitle:
            xtitle = ""
        if not ytitle:
            ytitle = ""
        projected = self._spectrum.Projection1D(histname, self._spectrum.GetAxisDefinition().GetAxisName(dim))
        projected.GetXaxis().SetTitle(xtitle)
        projected.GetYaxis().SetTitle(ytitle)
        projected.Sumw2()
        if doNorm:
            NormaliseBinWidth(projected)
            # Normalise by number of events
            projected.Scale(1./self.GetEventCount())
        return projected         
          
    def Print(self):
        """
        Print status of the data container
        """
        statusEvents = "no"
        if self._events:
            statusEvents = "yes"
        statusSpectrum = "no"
        if self._spectrum:
            statusSpectrum = "yes"
        print "Datacontainer status: events: %s , Spectrum: %s" %(statusEvents, statusSpectrum)
        if self._spectrum:
            self._spectrum.Print() 
            
class TrackContainer(DataContainer):
    """
    Data representation of track specific histograms
    """
    
    def __init__(self, eventHist = None, trackHist = None, dataformat = "new"):
        """
        Constructor, initialising base class and additional data members
        """
        trackwrapper = None
        if trackHist:
            if dataformat == "new":
                trackwrapper = TrackTHnSparseNew(trackHist)
            else:
                trackwrapper = TrackTHnSparseOld(trackHist)
        DataContainer.__init__(self, eventHist, trackwrapper, "tracks")
        self._datahistname = "TrackHist"
    
    def SetVertexRange(self, minv, maxv):
        """
        Apply vertex selection both to the event counter and the track hist
        """
        if self._spectrum:
            self._spectrum.SetVertexCut(minv, maxv)
        if self._events:
            self._events.SetVertexRange(minv, maxv)
        
    def SetEtaRange(self, etamin, etamax):
        """
        Select tracks in a given eta range
        """
        if self._spectrum:
            self._spectrum.SetEtaCut(etamin, etamax)
            
    def SetPhiRange(self, phimin, phimax):
        """
        Select tracks in a given eta range
        """       
        if self._spectrum:
            self._spectrum.SetPhiCut(phimin, phimax)
        
    def SetPileupRejection(self, on):
        """
        Apply pileup rejection (yes or no)
        """
        self._events.SetUsePileupRejected(on)
        if self._spectrum:
            self._spectrum.SetPileupRejection(on)
            
    def SelectTrackCuts(self, cutID):
        """
        Select a set of track cuts
        """
        if self._spectrum:
            self._spectrum.SelectTrackCuts(cutID)
    
    def RequestSeenInMinBias(self):
        if self._spectrum:
            self._spectrum.SetRequestSeenInMB()
        
    def __copy__(self):
        """
        Shallow copy constructor
        """
        print "Simple copy called from %s" %(self.__class__)
        newobject = TrackContainer(self._events, self._spectrum, None)
        newobject._Copy(self)
        return newobject
             
    def __deepcopy__(self, memo):
        """
        Deep copy constructor
        """
        print "deep copy called from %s" %(self.__class__)
        newobject = TrackContainer(None, None, None)
        newobject._Deepcopy(self, memo)
        return newobject
   
    @staticmethod
    def BuildFromRootPrimitive(rootprimitive):
        """
        Build a cluster container from a root primitive
        """
        infopair = DataContainer.GetHistoPair(rootprimitive)
        return TrackContainer(infopair["events"], infopair["spectrum"])
         
class ClusterContainer(DataContainer):
    """
    Data representation of cluster specific histograms
    """
    
    def __init__(self, eventHist = None, clusterHist = None, dataformat = "new"):
        """
        Constructor, initialising base class and additional data members
        """
        clusterwrapper = None
        if dataformat == "new":
            clusterwrapper = ClusterTHnSparseNew(clusterHist)
        else:
            clusterwrapper = ClusterTHnSparseOld(clusterHist)
        DataContainer.__init__(self, eventHist, clusterwrapper, "clusters")
        self._datahistname = "ClusterHist"

    def SetVertexRange(self, minv, maxv):
        """
        Apply vertex selection both to the event counter and the track hist
        """
        if self._spectrum:
            self._spectrum.SetVertexCut(minv, maxv)
        if self._events:
            self._events.SetVertexRange(minv, maxv)
        
    def SetEtaRange(self, etamin, etamax):
        """
        Select tracks in a given eta range
        """
        if self._spectrum:
            self._spectrum.SetEtaCut(etamin, etamax)
            
    def SePhiRange(self, phimin, phimax):
        """
        Select tracks in a given eta range
        """       
        if self._spectrum:
            self._spectrum.SetPhiCut(phimin, phimax)
        
    def SetPileupRejection(self, on):
        """
        Apply pileup rejection (yes or no)
        """
        self._events.SetUsePileupRejected(on)
        self._spectrum.SetPileupRejection(on)
            
    def RequestSeenInMinBias(self):
        if self._spectrum:
            self._spectrum.SetRequestSeenInMB()

    def __copy__(self):
        """
        Shallow copy constructor
        """
        print "Simple copy called from %s" %(self.__class__)
        newobject = ClusterContainer(self._events, self._spectrum, None)
        newobject._Copy(self)
        return newobject
             
    def __deepcopy__(self, memo):
        """
        Deep copy constructor
        """
        print "deep copy called from %s" %(self.__class__)
        newobject = ClusterContainer(None, None, None)
        newobject._Deepcopy(self, memo)
        return newobject
       
    @staticmethod
    def BuildFromRootPrimitive(rootprimitive):
        """
        Build a cluster container from a root primitive
        """
        infopair = DataContainer.GetHistoPair(rootprimitive)
        return ClusterContainer(infopair["events"], infopair["spectrum"])

