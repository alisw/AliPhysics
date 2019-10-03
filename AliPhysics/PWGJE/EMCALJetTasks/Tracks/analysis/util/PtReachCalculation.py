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
from scipy.optimize import fsolve

class PtReachCalculator(object):
    """
    classdocs
    """

    def __init__(self, name, data, isMinBias, limit):
        '''
        Constructor
        '''
        self.__fitter = None
        if isMinBias:
            self.__fitter = MinBiasFitter(name, data)
        else:
            self.__fitter = TriggeredSpectrumFitter(name, data)
        self.__limit = limit
        
    def GetPtReach(self, numberOfEvents):
        """
        Get the Pt reach for a given number of events
        """
        model = lambda p : numberOfEvents * self.__fitter.GetParameterisedValueAt(p) - self.__limit
        initialGuess = 10.
        result = fsolve(model, initialGuess)
        return result
    
    def GetPtReachForIntegral(self, numberOfEvents):
        """
        Get the Pt reach for a given number of events using integrated yield above 
        """
        model = lambda p : numberOfEvents * self.__fitter.GetNormalisedIntegralAbove(p) - self.__limit
        initialGuess = 10.
        result = fsolve(model, initialGuess)
        return result
