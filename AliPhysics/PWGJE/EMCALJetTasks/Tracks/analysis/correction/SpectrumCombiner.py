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
from copy import deepcopy
from PWGJE.EMCALJetTasks.Tracks.analysis.correction.TriggeredSpectrumScaler import TriggeredSpectrumScaler

class SpectrumCombiner(object):
    """
    Class combining the min. bias spectrum and the scaled-down triggered spectrum
    """

    def __init__(self, minbiasspectrum, triggeredspectrum):
        """
        Constructor
        """
        self.__minbiasspectrum = minbiasspectrum
        scaler = TriggeredSpectrumScaler(minbiasspectrum, triggeredspectrum)
        self.__triggeredspectrum = scaler.GetScaledTriggeredSpectrum()
        
    def MakeCombinedSpectrum(self, swappt):
        """
        Create a combined spectrum from the min bias spectrum and the triggered spectrum, using
        the points from the min bias spectrum up to a given pt, and the points from the triggered
        spectrum from that pt on
        """
        result = deepcopy(self.__minbiasspectrum)
        result.Sumw2()
        for mybin in range(1, result.GetXaxis().GetNbins()+1):
            inputspectrum = None
            if result.GetXaxis().GetBinUpEdge(mybin) <= swappt:
                inputspectrum = self.__minbiasspectrum
            else:
                inputspectrum = self.__triggeredspectrum
            result.SetBinContent(mybin, inputspectrum.GetBinContent(mybin))
            result.SetBinError(mybin, inputspectrum.GetBinError(mybin))
        return result