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
"""
Handling of jet-dependenet histograms in the analysis framework

@author: Markus Fasel
"""

from PWGJE.EMCALJetTasks.Tracks.analysis.base.struct.JetTHnSparse import JetTHnSparse
from PWGJE.EMCALJetTasks.Tracks.analysis.base.struct.EventHistogram import EventHistogramNew
from PWGJE.EMCALJetTasks.Tracks.analysis.base.Helper import NormaliseBinWidth
from copy import deepcopy

class JetPtBin():
    
    def __init__(self, jetpt):
        self.__jetpt = jetpt
        self.__spectrum = None
        self.__spectrumTrue = None
        
    def SetRecSpectrum(self, spec):
        self.__spectrum = JetTHnSparse(spec)
        
    def SetMCKineSpectrum(self, spec):
        self.__spectrumTrue = JetTHnSparse(spec)
         
    def GetJetPt(self):
        return self.__jetpt
    
    def GetRecSpectrum(self):
        return self.__spectrum
    
    def GetMCKineSpectrum(self):
        return self.__spectrumTrue
    
    def SetVertexRange(self, vtxmin, vtxmax):
        if self.__spectrum:
            self.__spectrum.SetVertexCut(vtxmin, vtxmax)
        if self.__spectrumTrue:
            self.__spectrumTrue.SetVertexCut(vtxmin, vtxmax)
    
    def SetMinJetPt(self, minpt):
        if self.__spectrum:
            self.__spectrum.ApplyCut("jetpt", minpt, 1000.)
        if self.__spectrumTrue:
            self.__spectrumTrue.ApplyCut("jetpt", minpt, 1000.)
    
    def SetEtaRange(self, etamin, etamax):
        if self.__spectrum:
            self.__spectrum.SetEtaCut(etamin, etamax)
        if self.__spectrumTrue:
            self.__spectrumTrue.SetEtaCut(etamin, etamax)
            
    def SetPhiRange(self, phimin, phimax):
        if self.__spectrum:
            self.__spectrum.SetPhiCut(phimin, phimax)
        if self.__spectrumTrue:
            self.__spectrumTrue.SetPhiCut(phimin, phimax)
    
    def SetSeenInMinBias(self):
        if self.__spectrum:
            self.__spectrum.ApplyCut(5, 1, 1)
        if self.__spectrumTrue:
            self.__spectrumTrue.ApplyCut(5, 1, 1)
            
    def GetProjectedRecSpectrum(self, dim, name, doNormBW = False):
        projected = self.__spectrum.Projection1D(name,self.__spectrum.GetAxisDefinition().GetAxisName(dim))
        if doNormBW:
            NormaliseBinWidth(projected)
        return projected
        
    def GetProjectedMCKineSpectrum(self, dim, name, doNormBW = False):
        projected = self.__spectrumTrue.Projection1D(name, self.__spectrumTrue.GetAxisDefinition().GetAxisName(dim))
        if doNormBW:
            NormaliseBinWidth(projected)
        return projected
    
class JetContainer(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.__jetpts =  []
        self.__vertexrange = {"min":None, "max":None}
        self.__eventhist = None
        
    def SetEventHist(self, hist):
        self.__eventhist = EventHistogramNew(deepcopy(hist))
        
    def AddJetPt(self, jetpt):
        existing = self.FindJetPt(jetpt)
        if not existing:
            self.__jetpts.append(JetPtBin(jetpt))
        
    def SetJetPtHist(self, jetpt, spec, isTrueKine):
        existing = self.FindJetPt(jetpt)
        if not existing:
            existing = JetPtBin(jetpt)
            self.__jetpts.append(existing)
        if isTrueKine:
            existing.SetMCKineSpectrum(spec)
        else:
            existing.SetRecSpectrum(spec)
            
    def FindJetPt(self, jetpt):
        found =  None 
        for entry in self.__jetpts:
            if entry.GetJetPt() == jetpt:
                found = entry
                break
        return found

    def SetVertexRange(self, vtxmin, vtxmax):
        self.__vertexrange["min"] = vtxmin
        self.__vertexrange["max"] = vtxmax
        for entry in self.__jetpts:
            entry.SetVertexRange(vtxmin, vtxmax)
        
    def SetEtaRange(self, etamin, etamax):
        for entry in self.__jetpts:
            entry.SetEtaRange(etamin, etamax)
            
    def SetPhiRange(self, phimin, phimax):
        for entry in self.__jetpts:
            entry.SetPhiRange(phimin, phimax)
            
    def SetRequestSeenInMinBias(self):
        for entry in self.__jetpts:
            entry.SetRequestSeenInMinBias()
            
    def GetListOfJetPts(self):
        listpts = []
        for entry in self.__jetpts:
            listpts.append(entry.GetJetPt())
        return listpts
            
    def MakeProjectionRecKine(self, jetpt, dimension, name, doNorm = False):
        existing = self.FindJetPt(jetpt)
        if not existing:
            return None
        projected = existing.GetProjectedRecSpectrum(dimension, name, doNorm)
        if projected and doNorm:
            projected.Scale(1./self.__eventhist.GetEventCount())
        return projected

    def MakeProjectionMCKine(self, jetpt, dimension, name, doNorm = False):
        existing = self.FindJetPt(jetpt)
        if not existing:
            return None
        projected = existing.GetProjectedMCKineSpectrum(dimension, name, doNorm)
        if projected and doNorm:
            projected.Scale(1./self.__eventhist.GetEventCount())
        return projected
