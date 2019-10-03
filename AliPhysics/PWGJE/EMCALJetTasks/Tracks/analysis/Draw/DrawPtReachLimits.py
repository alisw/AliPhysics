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
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.struct.DataContainers import SpectrumContainer
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.FileHandler import LegoTrainFileReader
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.plots.TrackYieldEventPlot import TrackYieldEventPlot
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.SpectrumCombiner import SpectrumCombiner
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.SpectrumFitter import PowerLawModel, ModifiedHagedornModel

def MakeNormalisedSpectrum(container):
    inputcontainer = container.FindTrackContainer("tracksAll")
    try:
        inputcontainer.SetVertexRange(-10., 10.)
        inputcontainer.SetPileupRejection(True)
        inputcontainer.SelectTrackCuts(1)
    except SpectrumContainer.RangeException as e:
        print str(e)
    return inputcontainer.MakeProjection(0, "ptSpectrum%s", "p_{#rm{t}} (GeV/c)", "1/N_{event} 1/(#Delta p_{#rm t}) dN/dp_{#rm{t}} ((GeV/c)^{-2}")

def Draw(filename, trigger = "MinBias"):
    reader = LegoTrainFileReader(filename)
    content = reader.ReadFile()
    content.SetName("RawSpectra")
    isMinBias = trigger == "MinBias"
    fitrange = [15., 50]
    if not isMinBias:
        fitrange = [50., 100.]
    plot = TrackYieldEventPlot(MakeNormalisedSpectrum(content.GetData(trigger)), fitrange, PowerLawModel())
    plot.SetLabel("Fit model: Power Law")
    plot.Create()
    return plot

def DrawCombined(filename, trigger = "EMCJHigh"):
    reader = LegoTrainFileReader(filename)
    content = reader.ReadFile()
    spectrumCombiner = SpectrumCombiner(MakeNormalisedSpectrum(content.GetData("MinBias")), \
                                        MakeNormalisedSpectrum(content.GetData(trigger)))
    
    plot = TrackYieldEventPlot(spectrumCombiner.MakeCombinedSpectrum(50.), [2., 90.], ModifiedHagedornModel())
    plot.SetLabel("Fit model: Modified Hagedorn Function")
    plot.Create()
    return plot
