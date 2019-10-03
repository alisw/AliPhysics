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
from PWGJE.EMCALJetTasks.Tracks.analysis.base.SpectrumFitter import MinBiasFitter
from PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import GraphicsObject
from PWGJE.EMCALJetTasks.Tracks.analysis.base.DataCollection import DataCollection, Datapoint

class TriggerTurnonCurve:
    
    def __init__(self, name, emcaldata, minbiasdata, fitmin = 15):
        self.__name = name
        self.__mbfitter = MinBiasFitter("mbfitter", minbiasdata)
        self.__values = self.__Create(emcaldata)
        
    def __Create(self, emcaldata):
        result = DataCollection("turnonCurve%s" %(self.__name));
        for mybin in range(1, emcaldata.GetXaxis().GetNbins()+1):
            minval = emcaldata.GetXaxis().GetBinLowEdge(mybin)
            if minval < 15:
                continue
            maxval = emcaldata.GetXaxis().GetBinUpEdge(mybin)
            binnedMb = self.__mbfitter.CalculateBinMean(minval, maxval)
            statError = emcaldata.GetBinError(mybin)/binnedMb
            datapoint = Datapoint(emcaldata.GetXaxis().GetBinCenter(mybin), emcaldata.GetBinContent(mybin)/binnedMb, emcaldata.GetXaxis().GetBinWidth(mybin)/2.)
            datapoint.AddErrorSource("stat", statError, statError)
            result.AddDataPoint(datapoint)
        return result;
    
    def GetPoints(self):
        return self.__values.MakeErrorGraphForSource("stat")
    
    def GetName(self):
        return self.__name
    
    def MakeGraphicsObject(self, style):
        return GraphicsObject(self.GetPoints(), style)
         
    def WriteData(self, name):
        self.GetPoints().Write("turnonCurve%s" %(name))
