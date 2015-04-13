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
from base.FileHandler import LegoTrainFileReader

from ROOT import TCanvas

def CountTracks(data, ptmin, ptmax):
    binmin = data.GetXaxis().FindBin(ptmin)
    binmax = data.GetXaxis().FindBin(ptmax)
    trackCount = 0
    for mybin in range(binmin, binmax + 1):
        trackCount = trackCount + data.GetBinContent(mybin)
    return trackCount

def PrintTrackCountPerBin(data, ptmin = 2.):
    kVerySmall = 1e-5
    binmin = data.GetXaxis().FindBin(ptmin + kVerySmall)
    binmax = data.GetXaxis().FindBin(100.)
    for mybin in range(binmin, binmax + 1):
        print "x = %f +-  %f GeV/c: %d counts" %(data.GetXaxis().GetBinCenter(mybin), data.GetXaxis().GetBinWidth(mybin)/2., data.GetBinContent(mybin)) 


def MakeRawSpectrum(inputdata, name):
    """
    Normalise spectrum by the number of events and by the bin width
    """
    inputdata.SetVertexRange(-10., 10.)
    inputdata.SetPileupRejection(True)
    inputdata.SelectTrackCuts(1)
    return inputdata.MakeProjection(0, "ptSpectrum%s" %(name), "p_{t} (GeV/c)", "1/N_{event} 1/(#Delta p_{t}) dN/dp_{t} ((GeV/c)^{-2})", doNorm = False)

def Run():
    fhandler = LegoTrainFileReader("AnalysisResults.root")
    filecontent = fhandler.ReadFile()
    canvas = TCanvas("spectra", "Spectra", 1000, 800)
    canvas.Divide(3,2)
    counter =1
    results = []
    mainclasses = ["EMCGHigh", "EMCGLow", "EMCJHigh", "EMCJLow", "MinBias"]
    for trigger in sorted(filecontent.GetListOfTriggers()):
        rawspectrum = MakeRawSpectrum(filecontent.GetData(trigger).FindTrackContainer("tracksAll"), trigger)
        print "%s: 50 - G0 GeV/c: %d, 80 - 100 GeV/c: %d" %(trigger, CountTracks(rawspectrum, 51, 59), CountTracks(rawspectrum, 82, 100))
        PrintTrackCountPerBin(rawspectrum, 2.)
        if trigger in mainclasses:
            canvas.cd(counter)
            rawspectrum.Draw("ep")
            results.append(rawspectrum)
            counter = counter + 1
    results.append(canvas)
    return results

if __name__ == "__main__":
    Run()