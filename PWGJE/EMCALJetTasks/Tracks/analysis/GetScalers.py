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
from ROOT import TFile, TH1D
from ROOT import gDirectory, gROOT

triggerlookup = {"MinBias":0, "EMCJHigh":1, "EMCJLow":2, "EMCGHigh":3, "EMCGLow":4}

def GetTriggerScalers(filename):
    inputfile = TFile.Open(filename)
    inputfile.cd("PtEMCalTriggerTask")
    gDirectory.ls()
    tasklist = gDirectory.Get("results")
    histlist = tasklist.FindObject("histosPtEMCalTriggerHistograms")
    gROOT.cd()
    triggerhist = histlist.FindObject("hEventTriggers")
    inputfile.Close()
    
    # Get number of Min. Bias counts
    mbcounts = GetCounts(triggerhist, "MinBias")
    print "MinBias counts: %d" %(mbcounts)

    triggerhist.GetAxis(0).SetRange(2,2)
    triggercounts = {}
    for trigger in triggerlookup.keys():
        if trigger == "MinBias":
            continue
        triggercounts[trigger] = GetCounts(triggerhist, trigger)
        print "Number of events for trigger %s: %d" %(trigger, triggercounts[trigger])
    
    hScalers = TH1D("triggerScalers", "trigger scalers", len(triggercounts), -0.5, len(triggercounts) - 0.5)
    counter = 1
    for trigger in triggercounts.keys():
        scaler = float(mbcounts)/float(triggercounts[trigger])
        print "Scaler for trigger %s: %f" %(trigger, scaler)
        hScalers.GetXaxis().SetBinLabel(counter, trigger)
        hScalers.SetBinContent(counter, scaler)
        counter += 1
        
    outputfile = TFile("TriggerScalers.root", "RECREATE")
    outputfile.cd()
    hScalers.Write()
    outputfile.Close()
    
def GetCounts(triggerhist, triggerclass):
    projection = triggerhist.Projection(triggerlookup[triggerclass])
    return projection.GetBinContent(2)