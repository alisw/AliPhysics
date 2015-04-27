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

class TriggerEfficiency:
    """
    Class calculating the trigger efficiency from a given min. bias container and a given triggered container
    """

    def __init__(self, triggername, minbiascontainer, triggeredcontainer):
        """
        Constructor
        """
        self.__triggername = triggername
        self.__minbiascontainer = minbiascontainer
        self.__triggeredcontainer = triggeredcontainer
        self.__triggerefficiency = None
        self.__CalculateTriggerEfficiency()
        
    def __MakeNormalisedSpectrum(self, container, name):
        container.SetVertexRange(-10., 10.)
        container.SetPileupRejection(True)
        if container.__class__ == "TrackContainer":
            container.SelectTrackCuts(1)
            container.RequestSeenInMinBias()
        return container.MakeProjection(0, "ptSpectrum%s" %(name), "p_{#rm{t}} (GeV/c)", "1/N_{event} 1/(#Delta p_{#rm t}) dN/dp_{#rm{t}} ((GeV/c)^{-2}", doNorm = False)
        
    def __CalculateTriggerEfficiency(self):
        minbiasspectrum = self.__MakeNormalisedSpectrum(self.__minbiascontainer, "minbias")
        self.__triggerefficiency = self.__MakeNormalisedSpectrum(self.__triggeredcontainer, self.__triggername)
        self.__triggerefficiency.Divide(self.__triggerefficiency, minbiasspectrum, 1., 1., "b")
        self.__triggerefficiency.SetName("triggerEff%s" %(self.__triggername))
        
    def GetEfficiencyCurve(self):
        return self.__triggerefficiency
        