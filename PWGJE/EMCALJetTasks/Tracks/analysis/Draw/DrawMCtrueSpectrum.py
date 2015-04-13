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
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.MonteCarloFileHandler import MonteCarloFileHandler
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import Style
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.Helper import NormaliseBinWidth
from os import getcwd
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.plots.MCTruthPlot import MCTrueSpectrumPlot, MCRecSpectrumPlot, MCWeightPlot
from ROOT import kRed, kBlue, kBlack, kGreen, kMagenta, kViolet, kOrange, kTeal, kYellow, kGray, kCyan

class PtHatBinDrawer(object):

    def __init__(self, plotter = None):
        self.__handler = MonteCarloFileHandler(True)
        self.__plotter = plotter
        self.__nbins = 1
        self._datacol = None

    def SetNumberOfPtHatBins(self, nbins):
        self.__nbins = nbins
        
    def SetBaseDirectory(self, path):
        for i in range(1, self.__nbins + 1):
            self.__handler.AddFile("%s/%02d/AnalysisResults.root" %(path, i), i)

    def CreatePlot(self):
        self._datacol = self.__handler.GetCollection()
        styles = [Style(kGreen-2, 24), Style(kBlue, 25), Style(kRed, 26), Style(kGreen, 27), Style(kMagenta, 28), Style(kOrange, 29), \
                  Style(kTeal, 30), Style(kViolet, 31), Style(kGray, 32), Style(kYellow + 2, 33), Style(kCyan+3, 34), Style(kRed-9, 35)]
        for i in range(1, self.__nbins + 1):
            spectrum = self._SpectrumForPtHatBin(i)
            self._datacol.GetWeigthHandler().ReweightSpectrum(i, spectrum)
            self.__plotter.AddMCSpectrum(i, spectrum, styles[i])
        self.__plotter.Create()
        return self.__plotter
    
    def _SpectrumForPtHatBin(self, pthatbin):
        return None

class MCTrueDrawer(PtHatBinDrawer):
    
    def __init__(self):
        PtHatBinDrawer.__init__(self, MCTrueSpectrumPlot())
    
    def _SpectrumForPtHatBin(self, pthatbin):
        cont = self._datacol.GetData(pthatbin).GetMCTruth()
        spectrum = cont.ProjectToDimension(0,"MCtruth%d" %(pthatbin))
        spectrum.Sumw2()
        NormaliseBinWidth(spectrum)
        trackCont = self._datacol.GetData(pthatbin).GetData("MinBias").FindTrackContainer("tracksAll")
        spectrum.Scale(1./trackCont.GetEventCount())
        return spectrum
   
class MCRecDrawer(PtHatBinDrawer):
    
    def __init__(self, trigger):
        PtHatBinDrawer.__init__(self, MCRecSpectrumPlot(trigger))
        self.__trigger = trigger
        self.__trackcontname = "tracksAll"
        
    def SelectContainer(self, contname):
        self.__trackcontname = contname
        
    def __MakeNormalisedSpectrum(self, container, pthatbin):
        container.SetVertexRange(-10., 10.)
        container.SetPileupRejection(True)
        if container.__class__ == "TrackContainer":
            container.SelectTrackCuts(1)
        return container.MakeProjection(0, "ptSpectrum%d%s" %(pthatbin, self.__trigger), "p_{#rm{t}} (GeV/c)", "1/N_{event} 1/(#Delta p_{#rm t}) dN/dp_{#rm{t}} ((GeV/c)^{-2}", doNorm = True)


    def _SpectrumForPtHatBin(self, pthatbin):
        trackCont = self._datacol.GetData(pthatbin).GetData(self.__trigger).FindTrackContainer(self.__trackcontname)
        return self.__MakeNormalisedSpectrum(trackCont, pthatbin)
    
class WeightPlotter():
    
    def __init__(self):
        self.__handler = MonteCarloFileHandler(True)
        self.__plotter = None
        self.__nbins = 1
    
    def SetNumberOfPtHatBins(self, nbins):
        self.__nbins = nbins
        
    def SetBaseDirectory(self, path):
        for i in range(1, self.__nbins + 1):
            self.__handler.AddFile("%s/%02d/AnalysisResults.root" %(path, i), i)
            
    def CreatePlot(self):
        self.__plotter = MCWeightPlot(self.__handler.GetCollection().GetWeigthHandler())
        self.__plotter.Create()
        return self.__plotter

def DrawWeights(basedir = None):
    drawer = WeightPlotter()
    drawer.SetNumberOfPtHatBins(9)
    if not basedir:
        basedir = getcwd()
    print "Using results from directory %s" %(basedir)
    drawer.SetBaseDirectory(basedir)
    return drawer.CreatePlot()    

def DrawMCtrue(basedir = None):
    drawer = MCTrueDrawer()
    #drawer.SetNumberOfPtHatBins(9)
    drawer.SetNumberOfPtHatBins(2)
    if not basedir:
        basedir = getcwd()
    print "Using results from directory %s" %(basedir)
    drawer.SetBaseDirectory(basedir)
    return drawer.CreatePlot()

def DrawMCrec(trigger, trackcont = "tracksAll", basedir = None):
    drawer = MCRecDrawer(trigger)
    drawer.SetNumberOfPtHatBins(9)
    drawer.SelectContainer(trackcont)
    if not basedir:
        basedir = getcwd()
    print "Using results from directory %s" %(basedir)
    drawer.SetBaseDirectory(basedir)
    return drawer.CreatePlot()