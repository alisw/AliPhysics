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
from PWGJE.EMCALJetTasks.Tracks.analysis.base.SpectrumFitter import MinBiasFitter, TriggeredSpectrumFitter
import numpy as np

class TriggeredSpectrumScaler(object):
    '''
    classdocs
    '''


    def __init__(self, minbiasspectrum, triggeredSpectrum):
        """
        Constructor
        """
        self.__triggeredspectrum = triggeredSpectrum
        self.__mbfitter = MinBiasFitter("minbiasfitter", minbiasspectrum)
        self.__trfitter = TriggeredSpectrumFitter("triggeredfittter", self.__triggeredspectrum)
        self.__isScaled = False
        
    def ScaleDownTriggeredSpectrum(self):
        """
        Scale down the triggered spectrum by a mean scale factor
        evaluated from the two fits above 60 GeV
        """
        if not self.__isScaled:
            pointlist = [60,65,70,75,80,85]
            scalingfactors = []
            for point in pointlist:
                scalingfactors.append(self.__mbfitter.GetParameterisedValueAt(point)/self.__trfitter.GetParameterisedValueAt(point))
            scalingfactor = np.array(scalingfactors, np.float64).mean()
            print "Using scaling factor %.2f" %(scalingfactor)
            self.__triggeredspectrum.Scale(scalingfactor)
            self.__isScaled = True
        
    def GetScaledTriggeredSpectrum(self):
        if not self.__isScaled:
            self.ScaleDownTriggeredSpectrum()
        return self.__triggeredspectrum