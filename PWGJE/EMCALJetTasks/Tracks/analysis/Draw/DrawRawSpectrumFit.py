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
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.FileHandler import ResultDataBuilder
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.struct.DataContainers import SpectrumContainer

from plots.RawDataFittingPlot import RawDataFittingPlot

def MakeNormalisedRawSpectrum(data):
    try:
        data.SetVertexRange(-10., 10.)
        data.SetPileupRejection(True)
        data.SelectTrackCuts(1)
    except SpectrumContainer.RangeException as e:
        print str(e)
    return data.MakeProjection(0, "ptSpectrum%s", "p_{#rm{t}} (GeV/c)", "1/N_{event} 1/(#Delta p_{#rm t}) dN/dp_{#rm{t}} ((GeV/c)^{-2}")

def DrawRawSpectrumIntegral(filename, trigger = "MinBias"):
    reader = ResultDataBuilder("lego", filename)
    content = reader.ReadFile()
    isMinBias = (trigger == "MinBias")
    
    plot = RawDataFittingPlot(MakeNormalisedRawSpectrum(content.GetData(trigger).FindTrackContainer("tracksAll")), isMinBias)
    plot.Create()
    return plot
