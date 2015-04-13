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
Representation of the track-based THnSparse

@author: Markus Fasel
"""
from PWGJE.EMCALJetTasks.Tracks.analysis.base.struct.THnSparseWrapper import THnSparseWrapper,AxisFormat
from copy import copy, deepcopy

class AxisFormatTracksOld(AxisFormat):
    '''
    Axis format for old track container
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        AxisFormat.__init__(self, "tracksold")
        self._axes["pt"] = 0
        self._axes["eta"] = 1
        self._axes["phi"] = 2
        self._axes["vertexz"] = 3
        self._axes["pileup"] = 4
        self._axes["trackcuts"] = 5
        self._axes["MBtrigger"] = 6
        
    def __deepcopy__(self, other, memo):
        '''
        deep copy constructor
        '''
        newobj = AxisFormatTracksOld()
        newobj._Deepcopy(other, memo)
        return newobj
    
    def __copy__(self, other):
        '''
        shallow copy constructor
        '''
        newobj = AxisFormatTracksOld()
        newobj._Copy()
        return newobj
        
class AxisFormatTracksNew(AxisFormat):
    '''
    Axis format for new track container
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        AxisFormat.__init__(self, "tracksnew")        
        self._axes["pt"] = 0
        self._axes["eta"] = 1
        self._axes["phi"] = 2
        self._axes["vertexz"] = 3
        self._axes["MBtrigger"] = 4
        
    def __deepcopy__(self, other, memo):
        '''
        deep copy constructor
        '''
        newobj = AxisFormatTracksNew()
        newobj._Deepcopy(other, memo)
        return newobj
    
    def __copy__(self, other):
        '''
        shallow copy constructor
        '''
        newobj = AxisFormatTracksNew()
        newobj._Copy()
        return newobj

class TrackTHnSparse(THnSparseWrapper):
    '''
    Base class for THnSparse in track format
    '''

    def __init__(self, roothist):
        '''
        Constructor
        '''
        THnSparseWrapper.__init__(self, roothist)
        
    def SetEtaCut(self, etamin, etamax):
        '''
        Apply cut in eta
        '''
        self.ApplyCut("eta",etamin,etamax)
    
    def SetPhiCut(self, phimin, phimax):
        '''
        Apply cut in phi
        '''
        self.ApplyCut("phi", phimin, phimax)
    
    def SetVertexCut(self, vzmin, vzmax):
        '''
        Apply cut on vertex-z
        '''
        self.ApplyCut("vertexz", vzmin, vzmax)

    def SetRequestSeenInMB(self, vzmin, vzmax):
        '''
        Request that the track was also in a min. bias event
        '''
        self.ApplyCut("mbtrigger", 1., 1.)
        
    def SelectTrackCuts(self, trackcuts):
        '''
        Select track cuts (old format only)
        '''
        if self._axisdefinition.FindAxis("trackcuts"):
            self.ApplyCut("trackcuts", trackcuts, trackcuts)
            
    def SetPileupRejection(self, on):
        if on and self._axisdefinition.FindAxis("pileup"):
            self.ApplyCut("pileup", 1., 1.)
            
    def Print(self):
        pass
        
class TrackTHnSparseOld(TrackTHnSparse):
    '''
    Class for old format track container
    '''
    
    def __init__(self, roothist):
        '''
        Constructor
        '''
        TrackTHnSparse.__init__(self, roothist)
        self._axisdefinition = AxisFormatTracksOld()

    def __deepcopy__(self, memo):
        '''
        Deep copy constructor
        '''
        result = TrackTHnSparseOld(deepcopy(self._rootthnsparse))
        result.CopyCuts(self._cutlist, True)
        return result
        
    def __copy__(self):
        '''
        Shallow copy constructor
        '''
        result = TrackTHnSparseOld(copy(self._rootthnsparse))
        result.CopyCuts(self._cutlist, False)
        return result

class TrackTHnSparseNew(TrackTHnSparse):
    '''
    Class for new format track container
    '''
    
    def __init__(self, roothist):
        '''
        Constructor
        '''
        TrackTHnSparse.__init__(self, roothist)
        self._axisdefinition = AxisFormatTracksNew()

    def __deepcopy__(self, memo):
        '''
        Deep copy constructor
        '''
        result = TrackTHnSparseNew(deepcopy(self._rootthnsparse))
        result.CopyCuts(self._cutlist, True)
        return result
        
    def __copy__(self):
        '''
        Shallow copy constructor
        '''
        result = TrackTHnSparseNew(copy(self._rootthnsparse))
        result.CopyCuts(self._cutlist, False)
        return result

