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
    
    # Get all relevant lists for which scaling should be done
    if referenceFile:
      if "PWGJEQA" in name or "Jet" in name or "Emcal" in name: # For the case of using reference scale factors, we want essentially all user tasks
        qaListNames.append(name)
    else:
      if "PWGJEQA" in name: # For the case of computing the scale factors, we want only the PWGJEQA task
        qaListNames.append(name)

    # Get a list that has the event histograms
    if "PWGJEQA" in name or "Jet" in name:
        eventList = name
    # Note: In the case of embedding, we assume that we only need the number of accepted events, and that internal event selection
    # is activated, in which case the event count can be read from any task that has the internal event selection applied

  print "Using " + eventList + " for event list."
  f.Close()

  # Create histogram of NEvents accepted and NEvents acc+rej, as a function of pT-hard bin
  hNEventsAcc = ROOT.TH1F("hNEventsAcc", "hNEventsAccepted", PtHardBins+1, 0, PtHardBins+1)
  hNEventsTot = ROOT.TH1F("hNEventsTot", "hNEventsTotal", PtHardBins+1, 0, PtHardBins+1)
  nEventsAccSum = 0
  for bin in range(0,PtHardBins):
    hNEventsAcc.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))
    hNEventsTot.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))
  for bin in range(0,PtHardBins):
    GetNEvents(eventList, bin, hNEventsTot, bAcceptedEventsOnly=False)
    nEventsAcc = GetNEvents(eventList, bin, hNEventsAcc, bAcceptedEventsOnly=True)
    nEventsAccSum += nEventsAcc
  nEventsAccAvg = nEventsAccSum/PtHardBins
  
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
      inputFile = "{0}/AnalysisResultsPtHard{0}.root".format(bin+1)
      print("Scaling Pt-hard bin %d" % (bin+1))
      f = ROOT.TFile(inputFile, "UPDATE")
      
      # Scale further to account for a different number of events in each pT-hard bin
      hNEventsAccepted = f.Get("hNEventsAcc")  # Get NEvents histos from file, since we undo the entry after writing it
      nEventsAcc = hNEventsAccepted.GetBinContent(bin+1)
      eventScaleFactor = nEventsAccAvg/nEventsAcc
      
      print "nEventsAcc: {0}".format(nEventsAcc)
      print "scaleFactor: {0}".format(scaleFactor)
      print "eventScaleFactor: {0}".format(eventScaleFactor)
      
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
    hXSecPerEvent = ROOT.TH1F("hXSecPerEvent", "hXSecPerEvent", PtHardBins+1, 0, PtHardBins+1)
    hNTrialsPerEvent = ROOT.TH1F("hNTrialsPerEvent", "hNTrialsPerEvent", PtHardBins+1, 0, PtHardBins+1)
    hScaleFactor = ROOT.TH1F("hScaleFactor", "hScaleFactor", PtHardBins+1, 0, PtHardBins+1)
    
    for bin in range(0,PtHardBins):
      # Label histograms
      hXSecPerEvent.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))
      hNTrialsPerEvent.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))
      hScaleFactor.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))

    for bin in range(0,PtHardBins):
      
      # Open input file and get relevant lists
      inputFile = "{0}/AnalysisResultsPtHard{0}.root".format(bin+1)
      print("Scaling Pt-hard bin %d" % (bin+1))
      f = ROOT.TFile(inputFile, "UPDATE")

      qaList = f.Get(qaListNames[0])
      print "Computing scaling factors with list: " + qaList.GetName()
      
      hXsecPtHard = qaList.FindObject("hXsec")
      hTrialsPtHard = qaList.FindObject("hNtrials")
      hNEventsTotal = f.Get("hNEventsTot") # Get NEvents histos from file, since we undo the entry after writing it
      hNEventsAccepted = f.Get("hNEventsAcc")

      # Compute: scale factor = xsec per event / trials per event
      nEventsTot = hNEventsTotal.GetBinContent(bin+1)
      nEventsAcc = hNEventsAccepted.GetBinContent(bin+1)
      xsec = hXsecPtHard.GetBinContent(1) / hXsecPtHard.GetEntries()
      trials = 1.*hTrialsPtHard.GetBinContent(1) / nEventsTot
      scaleFactor = xsec/trials
      eventScaleFactor = nEventsAccAvg/nEventsAcc # also scale to account that there are different number of events in each Pt-hard bin

      print "nEventsAcc: {0}".format(nEventsAcc)
      print "nEventsTot: {0}".format(nEventsTot)
      print "scaleFactor: {0}".format(scaleFactor)
      print "eventScaleFactor: {0}".format(eventScaleFactor)
      
      hXSecPerEvent.Fill(bin+0.5, xsec)
      hNTrialsPerEvent.Fill(bin+0.5, trials)
      hScaleFactor.Fill(bin+0.5, scaleFactor)

      # Now, scale all the histograms
      print "Scaling list: " + qaList.GetName()
      for obj in qaList:
        ScaleAllHistograms(obj, scaleFactor * eventScaleFactor,f)

      # Write the histograms to file
      hXSecPerEvent.Write()
      hXSecPerEvent.Reset()
      hNTrialsPerEvent.Write()
      hNTrialsPerEvent.Reset()
      hScaleFactor.Write()
      hScaleFactor.Reset()
      qaList.Write("%sScaled" % qaListNames[0], ROOT.TObject.kSingleKey)
      f.Close()

###################################################################################
# Given event list name eventList, pT-hard bin number, and histogram hNEvents of appropriate form, this function fills
# the number of events (accepted events only if bAcceptedEventsOnly=True, otherwise all events)
def GetNEvents(eventList, bin, hNEvents, bAcceptedEventsOnly = True):
  
  if bin is 0:
    if bAcceptedEventsOnly:
      print "Getting accepted number of events..."
    else:
      print "Getting total (acc+rej) number of events..."
  
  inputFile = "{0}/AnalysisResultsPtHard{0}.root".format(bin+1)
  f = ROOT.TFile(inputFile, "UPDATE")
  
  qaList = f.Get(eventList)
  nEvents = 0

  # Look for the EventCutOutput from AliEventCuts, and if it doesn't exist, look for histo fHistEventCount
  eventCutList = qaList.FindObject("EventCutOutput")
  if eventCutList:
    hNEventsPtHard = eventCutList.FindObject("fCutStats")
    if bAcceptedEventsOnly:
      nEvents = hNEventsPtHard.GetBinContent(16)
    else:
      nEvents = hNEventsPtHard.GetBinContent(1)
    if bin is 0:
      print "from EventCutOutput."
  else:
    hNEventsPtHard = qaList.FindObject("fHistEventCount")
    if bAcceptedEventsOnly:
      nEvents = hNEventsPtHard.GetBinContent(1)
    else:
      nEvents = hNEventsPtHard.GetBinContent(1) + hNEventsPtHard.GetBinContent(2)
    if bin is 0:
      print "from fHistEventCount."

  hNEvents.Fill(bin+0.5, nEvents)
  hNEvents.Write()
  hNEvents.Fill(bin+0.5,-nEvents) # undo the entry, since we will later hadd the histograms
  f.Close()

  return nEvents

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
