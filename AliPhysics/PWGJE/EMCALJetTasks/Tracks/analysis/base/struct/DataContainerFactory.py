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
from PWGJE.EMCALJetTasks.Tracks.analysis.base.struct.DataContainers import TrackContainer, ClusterContainer
from PWGJE.EMCALJetTasks.Tracks.analysis.base.struct.EventHistogram import EventHistogramOld, EventHistogramNew
from PWGJE.EMCALJetTasks.Tracks.analysis.base.struct.ParticleTHnSparse import CreatePartcileTHnSparse
from copy import deepcopy

class DataContainerFactory(object):
    '''
    classdocs
    '''


    def __init__(self, dataformat):
        '''
        Constructor
        '''
        self.__dataformat = dataformat
        
    def SetDataFormat(self, df):
        self.__dataformat = df
       
    def CreateTrackContainer(self, eventhist, trackhist):
        return TrackContainer(self.MakeEventHist(eventhist), trackhist, self.__dataformat)
    
    def CreateClusterContainer(self, eventhist, clusterhist):
        return ClusterContainer(self.MakeEventHist(eventhist), clusterhist, self.__dataformat)
    
    def CreateParticleContainer(self, particlehist):
        return CreatePartcileTHnSparse(particlehist, True if self.__dataformat == "new" else False)
    
    def MakeEventHist(self, eventhist):
        if self.__dataformat == "new":
            return EventHistogramNew(deepcopy(eventhist))
        else:
            return EventHistogramOld(deepcopy(eventhist))