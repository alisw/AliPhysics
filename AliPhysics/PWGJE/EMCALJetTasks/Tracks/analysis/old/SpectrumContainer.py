#! /usr/bin/env python

from base.Helper import NormaliseBinWidth
from base.struct.ClusterTHnSparse import AxisFormatClustersOld, AxisFormatClustersNew
from base.struct.TrackTHnSparse import AxisFormatTracksOld, AxisFormatTracksNew
from base.struct.EventHistogram import EventHistogramOld, EventHistogramNew
from base.MergeException import MergeException
from copy import copy,deepcopy
from ROOT import TList
 
class DataContainer(object):
    """
    Data container (event and spectum container)
    """
    
    class SpectrumCut(object):
        """
        Helper structure storing a  cut definition for a given dimension
        """
        def __init__(self, dimension, minv, maxv):
            """
            Constructor
            """
            self.__dimension = dimension
            self.__min = minv
            self.__max = maxv
            
        def SetDimension(self, dimension):
            """
            Set the dimension where we want to apply the cut
            """
            self.__dimension = dimension
            
        def GetDimension(self):
            """
            Return the value of the cut
            """
            return self.__dimension
        
        def SetLimits(self, minv, maxv):
            """
            Set the cut range
            """
            self.__min = minv
            self.__max = maxv
            
        def SetMinimum(self, minv):
            """
            Set the minimum of the range
            """
            self.__min = minv
            
        def SetMaximum(self, maxv):
            """
            Set the maximum of the range
            """
            self.__max = maxv
            
        def GetMinimum(self):
            """
            Return minimum value
            """
            return self.__min
        
        def GetMaximum(self):
            """
            Return maximum value
            """
            return self.__max
        
        def __copy__(self):
            print "Simple copy called from %s" %(self.__class__)
            newobject = DataContainer.SpectrumCut(self.Dimension, self.Minimumm, self.Maximum)
            return newobject
         
        def __deepcopy__(self, memo):
            print "deep copy called from %s" %(self.__class__)
            return self.__copy__()
         
        # Properties
        Minimum = property(GetMinimum, SetMinimum)
        Maximum = property(GetMaximum, SetMaximum)
        Dimension = property(GetDimension, SetDimension)
    
    
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
    
    
    def __init__(self, eventHist = None, dataHist = None, histotype=None, dataformat = None):
        """
        Construct spectrum container
        """
        self._events = None
        if eventHist:
            self._events = EventHistFactory.CreateEventHist(deepcopy(eventHist), dataformat)
        self.__spectrum = None
        if dataHist:
            self.__spectrum = SpectrumContainer(dataHist)
        self.__cutList = []
        self.__doNormBW = True
        self._AxisDefinition = None
        if dataformat:
            self._AxisDefinition = AxisFactory.GetAxisFormat(histotype, True if dataformat == "old" else False)
        self._datahistname = "DataHist"
        
    def __copy__(self):
        """
        Shallow copy constructor
        """
        print "Simple copy called from %s" %(self.__class__)
        newobject = DataContainer(self.__events, self.__spectrum, None, None)
        newobject._Copy(self)
        return newobject
             
    def __deepcopy__(self, memo):
        """
        Deep copy constructor
        """
        print "deep copy called from %s" %(self.__class__)
        newobject = DataContainer(None, None, None, None)
        newobject._Deepcopy(self, memo)
        return newobject
    
    def GetAxisDefinition(self):
        return self._AxisDefinition
        
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
        
    def SetTrackHist(self, trackhist):
        """
        Set the track hist
        """
        self.__spectrum = SpectrumContainer(trackhist)
        
    def SetSpectrumContainer(self, cont):
        """
        Initialise spectrum with container
        """
        self.__spectrum = cont
        
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

    
    def GetCutList(self):
        return self.__cutList
    
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
        self._AxisDefinition = copy(other.GetAxisDefinition())
        othercuts = other.GetCutList()
        for cut in othercuts:
            self.__cutList.append(copy(cut))
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
        self._AxisDefinition = deepcopy(other.GetAxisDefinition(), memo)
        othercuts = other.GetCutList()
        for cut in othercuts:
            self.__cutList.append(deepcopy(cut, memo))
        self.__SetValues(other.GetValues())
         
    def Add(self, other):
        """
        Add other data container to this data container
        """
        if not isinstance(other, DataContainer):
            raise MergeException("Incompatible types: this(DataContainer), other(%s)" %(str(other.__class__)))
        if self.__events and other.EventHist:
            self.__events.Add(other.EventHist)
        else:
            if other.EventHist:
                self.__events = deepcopy(other.EventHist)
        if self.__spectrum and other.SpectrumHist:
            self.__spectrum.Add(other.SpectrumHist)
        else:
            if other.SpectrumHist:
                self.__spectrum = other.SpectrumHist
         
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
        spectrum = self.__spectrum.GetRootPrimitive()
        spectrum.SetName("spectrum")
        result.Add(spectrum)
        return result
     
    @staticmethod
    def GetHistoPair(rootprimitive):
        """
        Convert to TList to python dictionary
        """
        return {"events":rootprimitive.FindObject("events"), "spectrum":rootprimitive.FindObject("spectrum")}
        
    def _AddCut(self, dimension, minv, maxv):
        """
        Add cut for a given dimension
        """
        cutFound = False
        for cut in self.__cutList:
            if cut.GetDimension() == dimension:
                #cut needs to be changed
                cut.SetLimits(minv, maxv)
                cutFound = True
                break
        if not cutFound:
            self.__cutList.append(DataContainer.SpectrumCut(dimension, minv, maxv))

    def GetEventCount(self):
        """
        Get the number of selected events
        """
        return self._events.GetEventCount()
            
    def MakeProjection(self, dim, histname = None, xtitle = None, ytitle = None, doNorm = True):
        """
        Make event-normalised projection to 1D
        """
        if not self.__spectrum:
            raise DataContainer.DataException(self._datahistname)
        if not self.__events:
            raise DataContainer.DataException("EventHist")
        # Apply cuts
        for cut in self.__cutList:
            #print "Processing next cut"
            self.__spectrum.ApplyCut(cut.GetDimension(), cut.GetMinimum(), cut.GetMaximum())
        if not histname:
            histname = "%s/" %(self.__spectrum.GetData().GetName())
        if not xtitle:
            xtitle = ""
        if not ytitle:
            ytitle = ""
        projected = self.__spectrum.ProjectToDimension(dim, histname, xtitle, ytitle)
        projected.Sumw2()
        if doNorm:
            NormaliseBinWidth(projected)
            # Normalise by number of events
            projected.Scale(1./self.GetEventCount())
        return projected         
          
    def Reset(self):
        """ 
        Reset underlying spectrum
        """
        self.__spectrum.Reset()
        for clhist in self.__clusters.values():
            clhist.Reset()
            
    def Print(self):
        """
        Print status of the data container
        """
        statusEvents = "no"
        if self.__events:
            statusEvents = "yes"
        statusSpectrum = "no"
        if self.__spectrum:
            statusSpectrum = "yes"
        print "Datacontainer status: events: %s , Spectrum: %s" %(statusEvents, statusSpectrum)
        if self.__spectrum:
            self.__spectrum.Print() 
            
class TrackContainer(DataContainer):
    """
    Data representation of track specific histograms
    """
    
    def __init__(self, eventHist = None, trackHist = None, dataformat = "new"):
        """
        Constructor, initialising base class and additional data members
        """
        DataContainer.__init__(self, eventHist, trackHist, "tracks", dataformat)
        self._datahistname = "TrackHist"
    
    def SetVertexRange(self, minv, maxv):
        """
        Apply vertex selection both to the event counter and the track hist
        """
        vzaxis = self._AxisDefinition.FindAxis("vertexz")
        if vzaxis > -1:
            self._events.SetVertexRange(minv, maxv)
            self._AddCut(vzaxis, minv, maxv)
        
    def SetEtaRange(self, etamin, etamax):
        """
        Select tracks in a given eta range
        """
        etaaxis = self._AxisDefinition.FindAxis("eta")
        if etaaxis > -1:
            self._AddCut(etaaxis, etamin, etamax)
            
    def SePhiRange(self, phimin, phimax):
        """
        Select tracks in a given eta range
        """
        phiaxis = self._AxisDefinition.FindAxis("phi")
        if phiaxis > -1:
            self._AddCut(phiaxis, phimin, phimax)
        
    def SetPileupRejection(self, on):
        """
        Apply pileup rejection (yes or no)
        """
        self._events.SetUsePileupRejected(on)
        pileupaxis = self._AxisDefinition.FindAxis("pileup")
        if pileupaxis > -1:
            self._AddCut(pileupaxis, 1., 1.)
            
    def SelectTrackCuts(self, cutID):
        """
        Select a set of track cuts
        """
        traxis = self._AxisDefinition.FindAxis("trackcuts")
        if traxis > -1:
            self._AddCut(traxis, cutID, cutID)
    
    def RequestSeenInMinBias(self):
        mbaxis = self._AxisDefinition.FindAxis("MBtrigger")
        if mbaxis > -1:
            self._AddCut(mbaxis, 1., 1.)
        
    def __copy__(self):
        """
        Shallow copy constructor
        """
        print "Simple copy called from %s" %(self.__class__)
        newobject = TrackContainer(self.__events, self.__spectrum, None)
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
        DataContainer.__init__(self, eventHist, clusterHist, "clusters", dataformat)
        self._datahistname = "ClusterHist"

    def SetVertexRange(self, minv, maxv):
        """
        Apply vertex selection both to the event counter and the track hist
        """
        vzaxis = self._AxisDefinition.FindAxis("vertexz")
        if vzaxis > -1:
            self._events.SetVertexRange(minv, maxv)
            self._AddCut(vzaxis, minv, maxv)
        
    def SetEtaRange(self, etamin, etamax):
        """
        Select tracks in a given eta range
        """
        etaaxis = self._AxisDefinition.FindAxis("eta")
        if etaaxis > -1:
            self._AddCut(etaaxis, etamin, etamax)
            
    def SePhiRange(self, phimin, phimax):
        """
        Select tracks in a given eta range
        """
        phiaxis = self._AxisDefinition.FindAxis("phi")
        if phiaxis > -1:
            self._AddCut(phiaxis, phimin, phimax)
        
    def SetPileupRejection(self, on):
        """
        Apply pileup rejection (yes or no)
        """
        self._events.SetUsePileupRejected(on)
        pileupaxis = self._AxisDefinition.FindAxis("pileup")
        if pileupaxis > -1:
            self._AddCut(pileupaxis, 1., 1.)
           
    def RequestSeenInMinBias(self):
        mbaxis = self._AxisDefinition.FindAxis("MBtrigger")
        if mbaxis > -1:
            self._AddCut(mbaxis, 1., 1.)

    def __copy__(self):
        """
        Shallow copy constructor
        """
        print "Simple copy called from %s" %(self.__class__)
        newobject = ClusterContainer(self.__events, self.__spectrum, None)
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
            
class SpectrumContainer(object):
    """
    Container class for combined spectrum
    """
    
    class DimensionException(Exception):
        """
        Exception class handling user requests to dimensions which do not exist
        """
        
        def __init__(self, dim, maxv):
            """
            Constructor, initialising max. and requested dimension
            """
            self.__dim = dim
            self.__max = maxv
        
        def __str__(self):
            """
            Make string representation of the range exception
            """
            return "Dimension outside range: max %d, requested %d" %(self.__max, self.__dim)
 
    class RangeException(Exception):
        """
        Exception class handling user request which are out of range
        """
        
        def __init__(self, dim, value, minv, maxv):
            """
            Constructor, initialising basic information about the range exception
            """
            self.__dimension = dim
            self.__value = value
            self.__minimum = minv
            self.__maximum = maxv
            
        def __str__(self):
            """
            Make string representation of the range exception
            """
            return "Range exceeded for dimension %d: %f not in [%f,%f]" %(self.__dimension, self.__value, self.__minimum, self.__maximum)
           
    def __init__(self, hsparse):
        """
        Constructor, defining underlying THnSparse
        """
        self.__hsparse = hsparse
#        self.__hsparse.Sumw2()
        
    def __copy__(self):
        """
        Shallow copy constructor
        """
        print "Simple copy called from %s" %(self.__class__)
        newobject = SpectrumContainer(self.__hsparse)
        return newobject
 
    def __deepcopy__(self, memo):
        """
        Deep copy constructor
        """
        print "deep copy called from %s" %(self.__class__)
        newobject = SpectrumContainer(deepcopy(self.__hsparse, memo))
        return newobject
        
    def GetData(self):
        """
        Access to underlying histogram
        """
        return self.__hsparse
    
    def SetData(self, data):
        """
        Setter for data
        """
        self.__hsparse = data
        
    Data = property(GetData, SetData)
    
    def Add(self, other):
        """
        Implement add method for the spectrum container, adding up the internal content
        """
        if not isinstance(other, SpectrumContainer):
            raise MergeException("Incompatible types: this(SpectrumContainer), other(%s)" %(str(other.__class__)))
        self.__hsparse.Add(other.GetData())
         
    def Scale(self, scalefactor):
        """
        Scale the underlying THnSparse with the scale factor
        """
        self.__hsparse.Scale(scalefactor)
        
    def GetRootPrimitive(self, name = None):
        """
        Return the root primitive
        """
        return self.__hsparse
    
    @staticmethod
    def BuildFromRootPrimitive(rootprimitive):
        """
        Create spectrum container from root primitive
        """
        return SpectrumContainer(rootprimitive)
        
    def ApplyCut(self, dim, minv, maxv):
        """
        Apply restrictions in one axis
        """
        kVerySmall = 1e-7
        cutaxis = None
        if isinstance(dim, int):
            cutaxis = self.__FindAxisByNumber(dim)
        elif isinstance(dim,str):
            cutaxis = self.__FindAxisByName(dim)
        if not self.__IsInRange(minv, cutaxis):
            raise self.RangeException(dim, minv, cutaxis.GetXmin(), cutaxis.GetXmax())
        if not self.__IsInRange(maxv, cutaxis):
            raise self.RangeException(maxv, cutaxis.GetXmin(), cutaxis.GetXmax())
        binmin = cutaxis.FindBin(minv + kVerySmall)
        binmax = cutaxis.FindBin(maxv - kVerySmall)
        #print "Setting range in axis %d from %.f [%d] to %.f [%d]" %(dim, min, binmin, max, binmax)
        self.__hsparse.GetAxis(dim).SetRange(binmin, binmax)
        
    def Reset(self):
        """
        Remove all cuts
        """
        for iaxis in range(0,self.__hsparse.GetNdimensions()):
            self.__hsparse.GetAxis(iaxis).SetRange()
            
    def Print(self):
        """
        Print content of the spectrum container
        """
        status = "no"
        if self.__hsparse:
            status = "yes"
        print "Spectrum, content set: %s" %(status)
            
    def ProjectToDimension(self, dimension, histname, xtitle = "", ytitle = ""):
        """
        Make 1D projection of the multi-dimensional histogram
        """
        if dimension >= self.__hsparse.GetNdimensions():
            raise self.DimensionException(self.__hsparse.GetNdimensions(), dimension)
        result = self.__hsparse.Projection(dimension)
        result.SetName(histname)
        if len(xtitle):
            result.GetXaxis().SetTitle(xtitle)
        if len(ytitle):
            result.GetYaxis().SetTitle(ytitle)
        return result
    
    def __IsInRange(self, value, axis):
        """
        Check whether value is in the range
        """
        if value < axis.GetXmin() or value > axis.GetXmax():
            return False
        return True
    
    def __FindAxisByNumber(self, dim):
        """
        Find axis for a given dimension number
        """
        if dim >= self.__hsparse.GetNdimensions():
            raise self.DimensionException(dim, self.__hsparse.GetNdimensions())
        return self.__hsparse.GetAxis(dim)
    
    def __FindAxisByName(self, name):
        """
        Find axis for a given dimension name
        """
        result = None
        for r in range (0, self.__hsparse.GetNdimensions()):
            axis = self.__hsparse.GetAxis(r)
            if name == axis.GetName():
                result = axis
                break
        return result
    
    def GetDimension(self, axisname):
        """
        Find dimension for a given name
        """
        result = -1
        for r in range (0, self.__hsparse.GetNdimensions()):
            if axisname == self.__hsparse.GetAxis(r).GetName():
                result = r
                break
        return result
    
class AxisFactory(object):
    
    def __init__(self):
        pass
    
    @staticmethod
    def GetAxisFormat(histotype, useold):
        result = None
        if histotype == "tracks":
            if useold:
                result = AxisFormatTracksOld()
            else:
                result = AxisFormatTracksNew()
        elif histotype == "clusters":
            if useold:
                result = AxisFormatClustersOld()
            else:
                result = AxisFormatClustersNew()
        return result
    


class EventHistFactory(object):
    
    def __init__(self):
        pass
    
    @staticmethod
    def CreateEventHist(histo, dataformat):
        if dataformat == "new":
            return EventHistogramNew(histo)
        elif dataformat == "old":
            return EventHistogramOld(histo)

