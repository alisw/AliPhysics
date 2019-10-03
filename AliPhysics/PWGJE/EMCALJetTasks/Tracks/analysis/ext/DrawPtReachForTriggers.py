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
from PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import SinglePanelPlot, Frame, GraphicsObject, Style
from ROOT import gROOT, TF1, TFile, kBlack, kBlue, kGreen, kRed, kOrange

class PtReachData:

    def __init__(self, data, trigger):
        self.__data = data
        self.__trigger = trigger
        self.__param = self.__CreateParameterisation()
        
    def __CreateParameterisation(self):
        param = TF1("powerlaw%s" %(self.__trigger), "[0] * TMath::Power(x, [1])", 0, 1e12)
        param.SetParLimits(0, 1e-5, 500)
        param.SetParLimits(1, 1e-1, 1)
        param.SetParameter(0, 2)
        param.SetParameter(1, 0.2)
        self.__data.Fit(param, "N", "", 1e8, 1e11)
        return param
    
    def Draw(self, pad, style):
        pad.DrawGraphicsObject(GraphicsObject(self.__data, style), True, self.__trigger, "epsame")
        pad.DrawGraphicsObject(GraphicsObject(self.__param, style), False, "param%s" %(self.__trigger), "lsame")

    
class PtReachPlot(SinglePanelPlot):
    
    class StyledData:
        
        def __init__(self, data, style):
            self.__data = data
            self.__style = style
            
        def GetData(self):
            return self.__data
        
        def GetStyle(self):
            return self.__style
        
        def Draw(self, pad):
            return self.__data.Draw(pad, self.__style)
    
    def __init__(self):
        SinglePanelPlot.__init__(self)
        self.__data = []
        
    def AddTriggerData(self, data, style):
        self.__data.append(self.StyledData(data, style))
        
    def Create(self):
        self._OpenCanvas("PtReachPlot", "PtReach for different triggers")
        pad = self._GetFramedPad()
        pad.GetPad().SetLogx()
        frame = Frame("ptreach", 1e5, 1e12, 0., 500.)
        frame.SetXtitle("Number of events")
        frame.SetYtitle("Max. p_{t} (GeV/c)")
        pad.DrawFrame(frame)
        for mydat in self.__data:
            mydat.Draw(pad)
        pad.CreateLegend(0.15, 0.7, 0.45, 0.89)
        
class PtReachDrawer:
    
    def __init__(self, limit):
        self.__plot = PtReachPlot()
        self.__limit = limit
        
    def Create(self):
        styles = {"MinBias":Style(kBlack, 20), "EMCJHigh":Style(kRed, 24), "EMCJLow":Style(kBlue, 25), "EMCGHigh":Style(kOrange, 26), "EMCGLow":Style(kGreen, 27)}
        for trigger in styles.keys():
            mytrigger = trigger
            if trigger == "MinBias":
                mytrigger = "MinBiasCombined"
            self.__plot.AddTriggerData(PtReachData(self.__ReadFile("PtReach%s.root" %(mytrigger)), trigger), styles[trigger])
        self.__plot.Create()
        
    def __ReadFile(self, filename):
        inputfile = TFile.Open(filename)
        gROOT.cd()
        graphlist = inputfile.Get("Results%d" %(self.__limit))
        ptreachgraph = graphlist.FindObject("ptreach")
        inputfile.Close()
        return ptreachgraph
        
    def GetPlot(self):
        return self.__plot
    
def MakePlot():
    drawer = PtReachDrawer(50)
    drawer.Create()
    return drawer.GetPlot()