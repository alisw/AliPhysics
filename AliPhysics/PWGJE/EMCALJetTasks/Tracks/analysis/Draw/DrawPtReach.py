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
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import Style
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.plots.PtReachPlot import PtReachData, PtReachPlot
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.FileHandler import LegoTrainFileReader
from ROOT import kBlack, kRed, kBlue, kGreen, kOrange

def MakeNormalisedSpectrum(container):
    inputcontainer = container.FindTrackContainer("tracksAll")
    try:
        inputcontainer.SetVertexRange(-10., 10.)
        inputcontainer.SetPileupRejection(True)
        inputcontainer.SelectTrackCuts(1)
    except SpectrumContainer.RangeException as e:
        print str(e)
    return inputcontainer.MakeProjection(0, "ptSpectrum%s", "p_{#rm{t}} (GeV/c)", "1/N_{event} 1/(#Delta p_{#rm t}) dN/dp_{#rm{t}} ((GeV/c)^{-2}")


def DrawPtReach(filename, doIntegral = False):
    reader = LegoTrainFileReader(filename)
    content = reader.ReadFile()
    content.SetName("RawSpectra")
    
    plot = PtReachPlot()
    triggers = ["MinBias", "EMCJHigh", "EMCJLow", "EMCGHigh", "EMCGLow"]
    styles = {"MinBias":Style(kBlack, 20), "EMCJHigh":Style(kRed,24), "EMCJLow":Style(kOrange, 25), "EMCGHigh":Style(kBlue, 26), "EMCGLow":Style(kGreen, 27)}
    isMinBias = False
    for trigger in triggers:
        if trigger == "MinBias":
            isMinBias = True
        else:
            isMinBias = False
        plot.AddData(trigger, PtReachData(trigger, MakeNormalisedSpectrum(content.GetData(trigger)), isMinBias, doIntegral), styles[trigger])
    plot.Create()
    return plot
