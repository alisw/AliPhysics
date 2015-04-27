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
from ROOT import kRed,kBlue,kGreen,kOrange,kBlack
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.FileHandler import ResultDataBuilder
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import Style
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.plots.TriggeredSpectrumComparisonPlot import PtTriggeredSpectrumComparisonPlot, EnergyTriggeredSpectrumComparisonPlot

class RawSpectrumDrawer(object):
    styles = {"MinBias": Style(kBlack, 20), "EMCJHigh":Style(kRed,24), "EMCJLow":Style(kOrange,26), "EMCGHigh":Style(kBlue,25), "EMCGLow":Style(kGreen, 27)}
    
    def __init__(self, filename, filetype = "lego"):
        reader = ResultDataBuilder(filetype, filename)
        self.__data = reader.GetResults()
        self.__plots = {}
        for tracks in ["tracksAll", "tracksWithClusters"]:
            self.__plots[tracks] = self.MakeTrackComparisonPlot(tracks)
        for clusters in ["Calib", "Uncalib"]:
            self.__plots[clusters] = self.MakeClusterComparisonPlot(clusters)
        
    def MakeTrackComparisonPlot(self, contname):
        plot = PtTriggeredSpectrumComparisonPlot(contname)
        for trg in self.styles.keys():
            plot.AddSpectrum(trg, self.MakeNormalisedSpectrum(self.__data.GetData(trg).FindTrackContainer(contname), "%s%s" %(trg, contname), True), self.styles[trg])
        plot.Create()
        return plot
    
    def MakeClusterComparisonPlot(self, contname):
        plot = EnergyTriggeredSpectrumComparisonPlot(contname)
        nminbias = self.__data.GetData("MinBias").FindClusterContainer(contname).GetEventCount()
        for trg in self.styles.keys():
            plot.AddSpectrum(trg, self.MakeNormalisedSpectrum(self.__data.GetData(trg).FindClusterContainer(contname), "%s%s" %(trg, contname), False), self.styles[trg])
        plot.Create()
        return plot
        
    def MakeNormalisedSpectrum(self, spectrum, name, istrack):
        spectrum.SetVertexRange(-10., 10.)
        spectrum.SetPileupRejection(True)
        if istrack:
            spectrum.SelectTrackCuts(1)
        projected = spectrum.MakeProjection(0, "rawspectrum%s" %(name), "p_{#rm{t}} (GeV/c)", "1/N_{event} 1/(#Delta p_{#rm t}) dN/dp_{#rm{t}} ((GeV/c)^{-2}")
        return projected

    def MakeNormalisedSpectrumV1(self, spectrum, name, nminbias, istrack):
        spectrum.SetVertexRange(-10., 10.)
        spectrum.SetPileupRejection(True)
        if istrack:
            spectrum.SelectTrackCuts(1)
        spectrum.RequestSeenInMinBias()
        projected =  spectrum.MakeProjection(0, "rawspectrum%s" %(name), "p_{#rm{t}} (GeV/c)", "1/N_{event} 1/(#Delta p_{#rm t}) dN/dp_{#rm{t}} ((GeV/c)^{-2}", doNorm = False)
        projected.Scale(nminbias)
        return projected
    
    def FindPlots(self, plotname):
        return self.__plots[plotname]
    
    def GetListOfPlots(self):
        return self.__plots


def CreatePlot(filename, filetype = "lego"):
    plotter = RawSpectrumDrawer(filename, filetype)
    return plotter

