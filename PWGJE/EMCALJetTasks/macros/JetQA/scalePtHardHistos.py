#! /usr/bin/env python

# Macro to scale histograms of all Pt-hard bins, using xsec and ntrials from AliAnalysisTaskPWGJEQA.
# This script expects files X/AnalysisResultsPtHardX.root, and will output scaled histograms
# to the same file, in a new output list with suffix "Scaled". The script will automatically loop over
# all output lists, subject to some simple criteria that covers basic use cases (can be adapted as needed).
#
# To get scale factors from a reference file, use the option "-f".
#
# Author: James Mulligan (james.mulligan@yale.edu)
#

import ROOT
import argparse

###################################################################################
# Main function
def scalePtHardHistos(referenceFile):

  PtHardBins = 20
  ptHardLo = [ 5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235 ]
  ptHardHi = [ 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, -1 ]
  
  # Get a list of all the output list names, in order to scale all of them
  f = ROOT.TFile("1/AnalysisResultsPtHard1.root", "READ")
  qaListKeys = f.GetListOfKeys()
  qaListNames = []
  eventList = ""
  for key in qaListKeys:
    name = key.GetName()
    if "PWGJEQA" in name or "Jet" in name or "Emcal" in name:
      qaListNames.append(name)
    if "PWGJEQA" in name or "Jet" in name:
      eventList = name
      print "Using " + name + " for event list."
  f.Close()

  # Create histogram of NEvents per pT-hard bin
  hNEvents = ROOT.TH1F("hNEvents", "hNEventsAccepted", PtHardBins+1, 0, PtHardBins+1)
  nEventsTotal = 0.
  for bin in range(0,PtHardBins):
    hNEvents.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))
  for bin in range(0,PtHardBins):
    inputFile = "%d/AnalysisResultsPtHard%d.root" % (bin+1, bin+1)
    f = ROOT.TFile(inputFile, "UPDATE")
    qaList = f.Get(eventList)
    hNEventsPtHard = qaList.FindObject("fHistEventCount")
    if hNEventsPtHard:
      nEvents = hNEventsPtHard.GetBinContent(1)
      if bin is 1:
        print "Getting nEvents from fHistEventCount."
    else:
      eventCutList = qaList.FindObject("EventCutOutput")
      hNEventsPtHard = eventCutList.FindObject("fCutStats")
      nEvents = hNEventsPtHard.GetBinContent(16)
      if bin is 1:
        print "Getting nEvents from EventCutOutput."

    nEventsTotal += nEvents
    hNEvents.Fill(bin, nEvents)
    hNEvents.Write()
    hNEvents.Fill(bin,-nEvents) # undo the entry, since we will later hadd the histograms
    f.Close()
  nEventsAvg = nEventsTotal/PtHardBins
  
  # If a reference file is provided, get the scale factors from there, and scale histos
  if referenceFile:
    print "Scaling from reference file: " + referenceFile
    for bin in range(0,PtHardBins):
      
      # Open ref file and get scale factor for given bin
      refFile = ROOT.TFile(referenceFile, "READ")
      scaleFactorHist = refFile.Get("hScaleFactor")
      scaleFactor = scaleFactorHist.GetBinContent(bin+1)
      refFile.Close()
      
      # Open input file and get relevant lists
      inputFile = "%d/AnalysisResultsPtHard%d.root" % (bin+1, bin+1)
      print("Scaling Pt-hard bin %d" % (bin+1))
      f = ROOT.TFile(inputFile, "UPDATE")
      
      hNEvents = f.Get("hNEvents")
      nEvents = hNEvents.GetBinContent(bin+1)
      eventScaleFactor = nEventsAvg/nEvents # also scale to account that there are different number of events in each Pt-hard bin
      
      for qaListName in qaListNames:
        qaList = f.Get(qaListName)
      
        # Now, scale all the histograms
        print "Scaling list: " + qaList.GetName()
        for obj in qaList:
          ScaleAllHistograms(obj, scaleFactor * eventScaleFactor, f)
        
        # Write the histograms to file
        qaList.Write("%sScaled" % qaListName, ROOT.TObject.kSingleKey)
    
      f.Close()

  # If no reference file is provided, compute the scale factors and write them, and scale the histos
  else:
    hNEvents = ROOT.TH1F("hNEvents", "hNEventsAccepted", PtHardBins+1, 0, PtHardBins+1)
    hXSecPerEvent = ROOT.TH1F("hXSecPerEvent", "hXSecPerEvent", PtHardBins+1, 0, PtHardBins+1)
    hNTrialsPerEvent = ROOT.TH1F("hNTrialsPerEvent", "hNTrialsPerEvent", PtHardBins+1, 0, PtHardBins+1)
    hScaleFactor = ROOT.TH1F("hScaleFactor", "hScaleFactor", PtHardBins+1, 0, PtHardBins+1)
    
    nEventsTotal = 0.
    for bin in range(0,PtHardBins):
      # Label histograms
      hNEvents.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))
      hXSecPerEvent.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))
      hNTrialsPerEvent.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))
      hScaleFactor.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))
      
    for bin in range(0,PtHardBins):
      # Fill nEvents histogram, for all Pt-hard bins
      inputFile = "%d/AnalysisResultsPtHard%d.root" % (bin+1, bin+1)
      f = ROOT.TFile(inputFile, "UPDATE")
      qaList = f.Get(eventList)
      hNEventsPtHard = qaList.FindObject("fHistEventCount")
      nEvents = hNEventsPtHard.GetBinContent(1)
      nEventsTotal += nEvents
      hNEvents.Fill(bin, nEvents)
      hNEvents.Write()
      hNEvents.Fill(bin,-nEvents) # undo the entry, since we will later hadd the histograms
      f.Close()
    nEventsAvg = nEventsTotal/PtHardBins

    for bin in range(0,PtHardBins):
      
      # Open input file and get relevant lists
      inputFile = "%d/AnalysisResultsPtHard%d.root" % (bin+1, bin+1)
      print("Scaling Pt-hard bin %d" % (bin+1))
      f = ROOT.TFile(inputFile, "UPDATE")

      qaList = f.Get(qaListName[0])
      
      hXsecPtHard = qaList.FindObject("hXsec")
      hTrialsPtHard = qaList.FindObject("hNtrials")
      hNEvents = f.Get("hNEvents")

      # Compute: scale factor = xsec per event / trials per event
      nEvents = hNEvents.GetBinContent(bin+1)
      xsec = hXsecPtHard.GetBinContent(1) / hXsecPtHard.GetEntries()
      trials = 1.*hTrialsPtHard.GetBinContent(1) / nEvents
      scaleFactor = xsec/trials
      eventScaleFactor = nEventsAvg/nEvents # also scale to account that there are different number of events in each Pt-hard bin
      
      hXSecPerEvent.Fill(bin, xsec)
      hNTrialsPerEvent.Fill(bin, trials)
      hScaleFactor.Fill(bin, scaleFactor*eventScaleFactor)

      # Now, scale all the histograms
      for obj in qaList:
        ScaleAllHistograms(obj, scaleFactor,f)

      # Write the histograms to file
      hXSecPerEvent.Write()
      hXSecPerEvent.Reset()
      hNTrialsPerEvent.Write()
      hNTrialsPerEvent.Reset()
      hScaleFactor.Write()
      hScaleFactor.Reset()
      qaList.Write("%sScaled" % qaListName[0], ROOT.TObject.kSingleKey)
      f.Close()

###################################################################################
# Function to iterate recursively through an object to scale all TH1/TH2/THnSparse
def ScaleAllHistograms(obj, scaleFactor,f):
  if obj.InheritsFrom(ROOT.TProfile.Class()):
    print("TProfile %s not scaled..." % obj.GetName())
  elif obj.InheritsFrom(ROOT.TH2.Class()):
    obj.Sumw2()
    obj.Scale(scaleFactor)
    print("TH2 %s was scaled..." % obj.GetName())
  elif obj.InheritsFrom(ROOT.TH1.Class()):
    obj.Sumw2()
    obj.Scale(scaleFactor)
    print("TH1 %s was scaled..." % obj.GetName())
  elif obj.InheritsFrom(ROOT.THnSparse.Class()):
    obj.Sumw2()
    obj.Scale(scaleFactor)
    print("THnSparse %s was scaled..." % obj.GetName())
  else:
    print("Not a histogram!")
    print obj.GetName()
    for subobj in obj:
      ScaleAllHistograms(subobj, scaleFactor,f)

#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  print("Executing scalePtHardHistos.py...")
  
  # Define arguments
  parser = argparse.ArgumentParser(description="Scale pT-hard bins")
  parser.add_argument("-f", "--referenceFile", action="store",
                      type=str, metavar="referenceFile",
                      default="",
                      help="Reference file to get pT-hard scale factors from")
  
  # Parse the arguments
  args = parser.parse_args()
  
  scalePtHardHistos(args.referenceFile)
