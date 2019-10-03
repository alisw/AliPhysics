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
from PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import TwoPanelPlot, GraphicsObject, Style, Frame
from PWGJE.EMCALJetTasks.Tracks.analysis.base.SpectrumFitter import SpectrumFitter
from PWGJE.EMCALJetTasks.Tracks.analysis.base.DataCollection import DataCollection, Datapoint
from copy import deepcopy
from ROOT import kRed, kBlue, kBlack, TF1, TList, TFile, TObject

class YieldCalculator:
    
    def __init__(self, rawspectrum, Model = None):
        self.__rawspectrum = rawspectrum
        self.__fitter = SpectrumFitter("fitter", rawspectrum, Model)
            
    def GetRawSpectrum(self):
        return self.__rawspectrum
    
    def GetFitCurve(self):
        return self.__fitter.GetFitFunction()
    
    def Initialise(self, fitrange):
        self.__fitter.DoFit(fitrange[0], fitrange[1])
            
    def GetRawYieldsForNumberOfEvents(self, nevents):
        histo = self.__fitter.MakeBinnedParameterisation(100, 0., 1000.)
        histo.Scale(nevents)
        return histo
    
    def FindMaxPt(self, nevents, limit):
        dndpthist = self.GetRawYieldsForNumberOfEvents(nevents)
        binabove = -1
        for mybin in range(2, dndpthist.GetXaxis().GetNbins() + 1):
            if dndpthist.GetBinContent(mybin) >= limit:
                binabove = mybin
            else:
                # we are below the limit, so break
                break
        if binabove == -1:
            return 0
        return dndpthist.GetBinCenter(mybin)
    
class EvDepPtReachData:
    
    def __init__(self, rawspectrum, yieldCalculator, limit):
        self.__limit = limit
        self.__points = DataCollection("ptreachdata%s" %(limit))
        self.__parameterisation = TF1("powerlaw%d" %(self.__limit), "[0] * TMath::Power(x, [1])", 0, 1e12)
        self.__datamodel = yieldCalculator
        self.__CreateCurve()
        
    def GetRawSpectrum(self):
        return self.__datamodel.GetRawSpectrum()
    
    def GetFit(self):
        return self.__datamodel.GetFitCurve()
    
    def GetReachForNumberOfEvents(self, nevents):
        return self.__parameterisation.Eval(nevents)
    
    def GetParameterisation(self, style):
        return GraphicsObject(self.__parameterisation, style)
        
    def __CreateCurve(self):
        eventstocheck = [1e6, 2e6, 5e6, 1e7, 5e7, 1e8, 5e8, 1e9, 5e9, 1e10, 5e10, 1e11, 5e11]
        for point in eventstocheck:
            self.__points.AddDataPoint(Datapoint(point, self.__datamodel.FindMaxPt(point, self.__limit), 0))
        self.__CreatePowerLawParameterisation()

    def GetGraphics(self, style):
        return GraphicsObject(self.__points.MakeLimitCurve(None, direction="central"), style)
    
    def __CreatePowerLawParameterisation(self):
        points = self.__points.MakeLimitCurve(None, direction="central")
        self.__parameterisation.SetParLimits(0, 1e-5, 500)
        self.__parameterisation.SetParLimits(1, 1e-1, 1)
        self.__parameterisation.SetParameter(0, 2)
        self.__parameterisation.SetParameter(1, 0.2)
        points.Fit(self.__parameterisation, "N", "", 1e8, 1e11)
  
    def MakeOutputList(self):
        result = TList()
        result.SetName("Results%d" %(self.__limit))
        rawspectrum = deepcopy(self.__datamodel.GetRawSpectrum())
        rawspectrum.SetName("rawspectrum")
        fitcurve = deepcopy(self.__datamodel.GetFitCurve())
        fitcurve.SetName("parameterisation")
        ptreach = self.__points.MakeLimitCurve(None, direction="central")
        ptreach.SetName("ptreach")
        result.Add(rawspectrum)
        result.Add(fitcurve)
        result.Add(ptreach)
        return result
    
class TrackYieldEventPlot(TwoPanelPlot):
    '''
    classdocs
    '''

    def __init__(self, rawspectrum, fitrange, model = None):
        '''
        Constructor
        '''
        TwoPanelPlot.__init__(self)
        self.__yieldCalculator = YieldCalculator(rawspectrum, model)
        self.__yieldCalculator.Initialise(fitrange)
        self.__models = {10: EvDepPtReachData(rawspectrum, self.__yieldCalculator, 10), \
                          50:EvDepPtReachData(rawspectrum,self.__yieldCalculator, 50), \
                          100:EvDepPtReachData(rawspectrum, self.__yieldCalculator, 100)} 
        self.__labeltext = None
        
    def SetLabel(self, labeltext):
        self.__labeltext = labeltext
        
    def GetData(self, limit):
        return self.__models[limit]
       
    def Create(self):
        self._OpenCanvas("trackyield", "Pt reach plot", 1000, 500)
        pad = self._OpenPad(1)
        pad.GetPad().SetLogx()
        pad.GetPad().SetLogy()
        sframe = Frame("sframe", 0., 100, 1e-8, 100)
        sframe.SetXtitle("p_{t} (GeV/c)")
        sframe.SetYtitle("1/N_{ev} dN/dp_{t} ((GeV/c)^{-1})")
        pad.DrawFrame(sframe)
        pad.DrawGraphicsObject(GraphicsObject(self.__models[10].GetRawSpectrum(),Style(kRed, 24)), True, "Data", "epsame")
        pad.DrawGraphicsObject(GraphicsObject(self.__models[10].GetFit(), Style(kBlue, 25)), True, "Fit", "lsame")
        pad.CreateLegend(0.5, 0.75, 0.89, 0.89)
        if self.__labeltext:
            pad.DrawLabel(0.15, 0.15, 0.75, 0.23, self.__labeltext)
        
        pad = self._OpenPad(2)
        pad.GetPad().SetLogx()
        pframe = Frame("pframe", 1e5, 1e12, 0., 300.)
        pframe.SetXtitle("Number of events")
        pframe.SetYtitle("Max p_{t} (GeV/c)")
        pad.DrawFrame(pframe)
        styles = {10:Style(kBlack, 24), 50:Style(kRed, 25), 100:Style(kBlue,26)}
        for limit, model in self.__models.iteritems():
            pad.DrawGraphicsObject(model.GetGraphics(styles[limit]), True, "At least %d tracks" %(limit), "epsame")
            pad.DrawGraphicsObject(model.GetParameterisation(styles[limit]), False, "Fit%d" %(limit), "lsame")
        pad.CreateLegend(0.15, 0.65, 0.55, 0.85)
        pad.DrawLabel(0.15, 0.55, 0.45, 0.6, "Tracks / 10 GeV/c")
        self._canvas.cd()
        
    def WriteData(self, filename):
        result = TFile(filename, "Recreate")
        for model in self.__models.itervalues():
            olist = model.MakeOutputList()
            olist.Write(olist.GetName(), TObject.kSingleKey)
        result.Close()
        