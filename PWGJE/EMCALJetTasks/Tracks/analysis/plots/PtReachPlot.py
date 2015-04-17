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
from PWGJE.EMCALJetTasks.Tracks.analysis.util.PtReachCalculation import PtReachCalculator
from PWGJE.EMCALJetTasks.Tracks.analysis.base.DataCollection import DataCollection, Datapoint
from PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import SinglePanelPlot, Frame, GraphicsObject
from ROOT import TFile

class PtReachData:    
    
    def __init__(self, name, data, isMinBias, doIntegral = True):
        self.__calculator = PtReachCalculator(name, data, isMinBias, 50)
        self.__data = DataCollection("data%s" %(name))
        self.__isIntegral = doIntegral
        self.__rootobject = None
        self.__CreateData()
        
    def __CreateData(self):
        events = [500000, 1000000, 20000000, 50000000, 100000000, 200000000, 300000000, 400000000, 500000000, 750000000, 1000000000]
        PtReachCalculation = lambda p : self.__calculator.GetPtReach(p)
        if self.__isIntegral:
            PtReachCalculation = lambda p : self.__calculator.GetPtReachForIntegral(p)
        for nevents in events:
            self.__data.AddDataPoint(Datapoint(nevents, PtReachCalculation(nevents), 0.))
            
    def MakeGraphics(self):
        if not self.__rootobject:
            self.__MakeRootObject()
        return self.__rootobject
    
    def IsIntegral(self):
        return self.__isIntegral
    
    def __MakeRootObject(self):
        self.__rootobject = self.__data.MakeLimitCurve(None, direction = "central")
    
class PtReachPlot(SinglePanelPlot):
    
    def __init__(self):
        SinglePanelPlot.__init__(self)
        self.__ptreachdata = {}
        
    def AddData(self, name, data, style):
        self.__ptreachdata[name] ={"data": data, "style": style}
        
    def Create(self):
        self._OpenCanvas("PtReachPlot", "Pt-reach depending on the number of events")
        
        pad = self._GetFramedPad()
        frame = Frame("framePtReach", 0., 1100000000, 0., 150)
        frame.SetXtitle("Number of events")
        frame.SetYtitle("Max. p_{t} (GeV/c)")
        pad.DrawFrame(frame)
        isIntegral = False
        for key,data in self.__ptreachdata.iteritems():
            graphics = GraphicsObject(data["data"].MakeGraphics(), data["style"])
            pad.DrawGraphicsObject(graphics, True, key , "p")
            isIntegral = data["data"].IsIntegral()
        pad.CreateLegend(0.5, 0.15, 0.89, 0.35)
        labelText = "At least 50 tracks at p_{t}"
        if isIntegral:
            labelText = "At least 50 tracks above p_{t}"
        pad.DrawLabel(0.5, 0.37, 0.7, 0.42, labelText)
        
    def WriteData(self, filename):
        """
        Write root data to rootfiles
        """
        result = TFile(filename, "RECREATE")
        for key,data in self.__ptreachdata.iteritems():
            graph = data["data"].MakeGraphics()
            graph.Write(key)
        result.Close()
        
            