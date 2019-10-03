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
from PWGJE.EMCALJetTasks.Tracks.analysis.base.struct.THnSparseWrapper import THnSparseWrapper, AxisFormat
from copy import copy, deepcopy

class AxisDefinitinonTrueParticles(AxisFormat):
    
    def __init__(self):
        AxisFormat.__init__(self, "trueparticles")
        self._axes["pt"] = 0
        self._axes["eta"] = 1
        self._axes["phi"] = 2
        self._axes["vertexz"] = 3
    
    def __deepcopy__(self, other, memo):
        '''
        deep copy constructor
        '''
        newobj = AxisDefinitinonTrueParticles()
        newobj._Deepcopy(other, memo)
        return newobj
    
    def __copy__(self, other):
        '''
        shallow copy constructor
        '''
        newobj = AxisDefinitinonTrueParticles()
        newobj._Copy()
        return newobj
    
class AxisDefinitionTrueParticlesOld(AxisDefinitinonTrueParticles):
    
    def __init__(self):
        AxisDefinitinonTrueParticles.__init__(self)
        self._axes["isPileup"] = 4
        
class AxisDefinitionTrueParticlesNew(AxisDefinitinonTrueParticles):
    
    def __init__(self):
        AxisDefinitinonTrueParticles.__init__(self)
            
class ParticleTHnSparse(THnSparseWrapper):
    '''
    THnSparse representation for MC truth
    '''

    def __init__(self, roothist):
        '''
        Constructor
        '''
        THnSparseWrapper.__init__(self, roothist)
        
    def __deepcopy__(self, memo):
        '''
        Deep copy constructor
        '''
        result = ParticleTHnSparse(deepcopy(self._rootthnsparse))
        result.CopyCuts(self._cutlist, True)
        return result
        
    def __copy__(self):
        '''
        Shallow copy constructor
        '''
        result = ParticleTHnSparse(copy(self._rootthnsparse))
        result.CopyCuts(self._cutlist, False)
        return result 

    def SetEtaCut(self, etamin, etamax):
        '''
        Apply cut in eta
        '''
        self.ApplyCut(self,etamin,etamax)
    
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
        
    def SetPileupCut(self):
        '''
        Apply Cut on pileup events(old format only)
        '''
        self.ApplyCut("isPileup", 1, 1)
        
class ParticleTHnSparseNew(ParticleTHnSparse):
    
    def __init__(self, roothist):
        ParticleTHnSparse.__init__(self, roothist)
        self._axisdefinition = AxisDefinitionTrueParticlesNew()
    
class ParticleTHnSparseOld(ParticleTHnSparse):
    
    def __init__(self, roothist):
        ParticleTHnSparse.__init__(self, roothist)
        self._axisdefinition = AxisDefinitionTrueParticlesOld()
        
def CreatePartcileTHnSparse(roothist, isNew):
    if isNew:
        return ParticleTHnSparseNew(roothist)
    else:
        return ParticleTHnSparseOld(roothist)
        