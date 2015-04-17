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
from copy import copy, deepcopy

class EventHistogram(object):
    
    def __init__(self, histo):
        self._histo = histo
        self._vertexrange = {}
    
    def GetROOTHisto(self):
        return self._histo
    
    def GetVertexRange(self):
        return self._vertexrange
    
    def SetVertexRange(self, vtxmin, vtxmax):
        self._vertexrange["min"] = vtxmin
        self._vertexrange["max"] = vtxmax

    def GetEventCount(self):
        print "Method virtual - to be implemented by inheriting classes"
        
    def _Deepcopy(self, other, memo):
        underlyinghist = other.GetROOTHisto()
        self._histo = deepcopy(underlyinghist, memo)
        self._vertexrange = deepcopy(other.GetVertexRange(), memo)
    
    def _Copy(self, other):
        underlyinghist = other.GetROOTHisto()
        self._histo = copy(underlyinghist)
        self._vertexrange = copy(other.GetVertexRange())
        
    def Add(self, otherhisto):
        self._histo.Add(otherhisto)
    
    def Scale(self, scalefactor):
        self._histo.Scale(scalefactor)
    
        
class EventHistogramOld(EventHistogram):
    
    def __init__(self, histo):
        EventHistogram.__init__(self, histo)
        self.__usePileupRejected = True
        
    def SetUsePileupRejected(self, doUse = True):
        self.__usePileupRejected = doUse
        
    def IsUsingPileupRejected(self):
        return self.__usePileupRejected
        
    def GetEventCount(self):
        if len(self._vertexrange):
            binMin = self._histo.GetYaxis().FindBin(self._vertexrange["min"])
            binMax = self._histo.GetYaxis().FindBin(self._vertexrange["max"])
            eventcounter = self._histo.ProjectionX("eventCounter", binMin, binMax)
        else:
            eventcounter = self._histo.ProjectionX("eventcounter")
        pileupbin = 1
        if self.__usePileupRejected:
            pileupbin = 2
        return eventcounter.GetBinContent(pileupbin)
    
    def __copy__(self, other):
        newobj = EventHistogramOld(None)
        newobj._Copy(other)
        newobj.SetUsePileupRejected(other.IsUsingPileupRejected())
        return newobj
        
    def __deepcopy__(self, other, memo):
        newobj = EventHistogramOld(None)
        newobj._Deepcopy(other, memo)
        newobj.SetUsePileupRejected(other.IsUsingPileupRejected())
        return newobj
    
    
class EventHistogramNew(EventHistogram):

    def __init__(self, histo):
        EventHistogram.__init__(self, histo)
        
    def SetUsePileupRejected(self, doUse = True):
        pass
        
    def GetEventCount(self):
        if not len(self._vertexrange):
            return self._histo.Integral()
        else:
            binmin = self._histo.GetXaxis().FindBin(self._vertexrange["min"])
            binmax = self._histo.GetXaxis().FindBin(self._vertexrange["max"])
            return self._histo.Integral(binmin, binmax)
    
    def __copy__(self, other):
        newobj = EventHistogramNew(None)
        newobj._Copy(other)
        return newobj
    
    def __deepcopy__(self, other, memo):
        newobj = EventHistogramNew(None)
        newobj._Deepcopy(other, memo)
        return newobj