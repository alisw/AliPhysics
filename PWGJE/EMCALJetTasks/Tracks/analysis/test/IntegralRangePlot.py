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
from PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import SinglePanelPlot, GraphicsObject, Style, Frame
from PWGJE.EMCALJetTasks.Tracks.analysis.base.SpectrumFitter import SpectrumFitter
from PWGJE.EMCALJetTasks.Tracks.analysis.base.DataCollection import DataCollection, Datapoint
from ROOT import kRed

class IntegralRangeChecker:
    
    def __init__(self, minbiasspectrum):
        self.__spectrumfitter = SpectrumFitter("minbiasfitter", minbiasspectrum)
        self.__spectrumfitter.DoFit(15., 50.)
        self.__upperrange = DataCollection("IntMax")
        self.__intmin = 0.
        
    def MakeIntMaxDependence(self, intmin = 30):
        for intmax in range(50, 1050, 50):
            value = self.__spectrumfitter.CalculateIntegral(float(intmin), float(intmax))
            print "Max.: %f, value: %e" %(intmax, value)
            self.__upperrange.AddDataPoint(Datapoint(float(intmax), value, 0.))
        self.__intmin = intmin
            
    def GetGraphics(self):
        return GraphicsObject(self.__upperrange.MakeLimitCurve(None, direction="central"), Style(kRed, 24))
    
    def GetLabelText(self):
        return "Min. p_{t}: %1.f GeV/c" %(self.__intmin)

class IntegralRangePlot(SinglePanelPlot):
    
    def __init__(self, minbiasspectrum):
        SinglePanelPlot.__init__(self)
        self.__integralChecker = IntegralRangeChecker(minbiasspectrum)
        
    def Create(self, intmin = 30.):
        self._OpenCanvas("IntegralRangeCheck", "Check of the upper limit of the integral")
        pad = self._GetFramedPad()
        frame = Frame("irframe", 0., 1100., 1e-6, 1e-4)
        frame.SetXtitle("Max p_{t} for integral")
        frame.SetYtitle("Integrated yield above p_{t}")
        pad.DrawFrame(frame)
        self.__integralChecker.MakeIntMaxDependence(intmin)
        pad.DrawGraphicsObject(self.__integralChecker.GetGraphics(), False, "", "p")
        pad.DrawLabel(0.15, 0.8, 0.45, 0.84, self.__integralChecker.GetLabelText())
        
        