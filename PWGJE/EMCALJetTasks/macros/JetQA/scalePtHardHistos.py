#! /usr/bin/env python

# Macro to scale histograms of all Pt-hard bins, using xsec and ntrials from AliAnalysisTaskPWGJEQA.
# This script expects files X/AnalysisResultsPtHardX.root, and will output scaled histograms
# to the same file, in a new output list with suffix "Scaled". The script will automatically loop over
# all output lists, subject to some simple criteria that covers basic use cases (can be adapted as needed).
#
# There is an option "bRemoveOutliers" (off by default) to remove outliers from certain histograms. The features are
# currently hard-coded below so you will need to modify the code as needed. This feature is adapted from code of Raymond Ehlers.
#
# To get scale factors from a reference file, use the option "-f".
#
# Author: James Mulligan (james.mulligan@yale.edu)
#

import ROOT
import argparse
import ctypes
import os
import gc

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

###################################################################################
# Main function
def scalePtHardHistos(referenceFile, StartPtHardBin, EndPtHardBin):

  PtHardBins = 20
  ptHardLo = [ 5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235 ]
  ptHardHi = [ 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, -1 ]
  
  # Option to remove outliers from specified histograms (see below)
  bRemoveOutliers = True

  outlierLimit=2          #(standard =2) average bin content of 5 consecutive bins. If the average bin content of this window drops below "outlierLimit",  the bin and all following bins are marked for a candiate for removal
  outlierNBinsThreshold=4 #if the average bin content of the window stays below the "outlierLimit" for "outlierNBinsThreshold" bins after the candidate was marked for removal it is finally declared an outlier
  
  # Option to print out detailed info about scaling and outlier removal
  verbose = False
  
  # Get a list of all the output list names, in order to scale all of them
  f = ROOT.TFile("{}/AnalysisResultsPtHard{}.root".format(StartPtHardBin,StartPtHardBin), "READ")
  qaListKeys = f.GetListOfKeys()
  qaListNames = []
  eventList = ""
  print("Scaling the following lists:")
  for key in qaListKeys:
    name = key.GetName()
    
    # Get all relevant lists for which scaling should be done
    if referenceFile:
      if "PWGJEQA" in name or "Jet" in name or "Emcal" in name: # For the case of using reference scale factors, we want essentially all user tasks
        qaListNames.append(name)
        print("    {}".format(name))
    else:
      if "PWGJEQA" in name: # For the case of computing the scale factors, we want only the PWGJEQA task
        qaListNames.append(name)
        print("    {}".format(name))
    # Get a list that has the event histograms
    if "PWGJEQA" in name or "JetPerformance" in name:
        eventList = name
    # Note: In the case of embedding, we assume that we only need the number of accepted events, and that internal event selection
    # is activated, in which case the event count can be read from any task that has the internal event selection applied

  print("")
  print("......... Get the Number of events .........")
  print("Using " + eventList + " for event list.")
  f.Close()

  # Create histogram of NEvents accepted and NEvents acc+rej, as a function of pT-hard bin
  hNEventsAcc = ROOT.TH1F("hNEventsAcc", "hNEventsAccepted", PtHardBins+1, 0, PtHardBins+1)
  hNEventsTot = ROOT.TH1F("hNEventsTot", "hNEventsTotal", PtHardBins+1, 0, PtHardBins+1)
  nEventsAccSum = 0
  nEventsTotSum=0
  for bin in range(StartPtHardBin-1, EndPtHardBin):
    hNEventsAcc.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))
    hNEventsTot.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))
  for bin in range(StartPtHardBin-1, EndPtHardBin):
    nNEventsTot= GetNEvents(eventList, bin, hNEventsTot, verbose, bAcceptedEventsOnly=False)
    nEventsAcc = GetNEvents(eventList, bin, hNEventsAcc, verbose, bAcceptedEventsOnly=True)
    nEventsAccSum += nEventsAcc
    nEventsTotSum += nNEventsTot
  nEventsAccAvg = nEventsAccSum/PtHardBins
  nNEventsTotAvg= nEventsTotSum/PtHardBins

  print("")
  print("......... Start the scaling .........")
  print("......... pT hard bins {} - {} ......".format(StartPtHardBin, EndPtHardBin))
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # If a reference file is provided, get the scale factors from there, and scale histos
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if referenceFile:
    print("Get Scale Factors from Reference file: {}".format(referenceFile))
    controlFile = open('ControlTable.txt', 'w') # This is a measure to check the factors used for the single bins
    controlFile.write('| *Bin* | *pT-hard range (!GeV)* | *Factor From File* | *Factor Applied (xEvtFactor)* |\n')
    controlFile.close()
    for bin in range(StartPtHardBin-1, EndPtHardBin):
    #for bin in range(0,PtHardBins):
      
      # Open ref file and get scale factor for given bin
      refFile = ROOT.TFile(referenceFile, "READ")
      scaleFactorHist = refFile.Get("hScaleFactor")
      scaleFactor = scaleFactorHist.GetBinContent(bin+1)
      refFile.Close()
      
      # Open input file and get relevant lists
      inputFile = "{0}/AnalysisResultsPtHard{0}.root".format(bin+1)
      print("................................................................")
      print("ooo Scaling Pt-hard bin %d" % (bin+1))
      f = ROOT.TFile(inputFile, "UPDATE")
 
      # Scale further to account for a different number of events in each pT-hard bin
      hNEventsAccepted = f.Get("hNEventsAcc")  # Get NEvents histos from file, since we undo the entry after writing it
      nEventsAcc = hNEventsAccepted.GetBinContent(bin+1)
      eventScaleFactorAcc = nEventsAccAvg/nEventsAcc
      
      print("ooo nEventsAcc: {0}".format(nEventsAcc))
      print("ooo scaleFactor: {0}".format(scaleFactor))
      print("ooo eventScaleFactor: {0}".format(eventScaleFactorAcc))
      print("ooo combined ScaleFactor: {0}".format(eventScaleFactorAcc*scaleFactor))
      controlFile = open('ControlTable.txt', 'a')
      controlFile.write("| *{0}* | {1} - {2} | {3} | {4} | \n".format((bin+1), ptHardLo[bin], ptHardHi[bin], scaleFactor, scaleFactor*eventScaleFactorAcc))
      controlFile.close()

      for qaListName in qaListNames:
        qaList = f.Get(qaListName)
        # Now, scale all the histograms
        print("")
        RadiusName=""
        for obj in qaList:
          name = obj.GetName()
          ScaleAllHistograms(obj, scaleFactor * eventScaleFactorAcc, f, verbose, bRemoveOutliers, outlierLimit, outlierNBinsThreshold, bin, EndPtHardBin, name, qaList.GetName())
        
        # Write the histograms to file
        qaList.Write("%sScaled" % qaListName, ROOT.TObject.kSingleKey)
        # Make sure we delete the lists and sublists, to prevent memory leak
        for obj in qaList:
          obj.SetOwner(True)
          obj.Delete()
        qaList.SetOwner(True)
        qaList.Delete()
      f.Close()
      f.Delete()
      gc.collect()
    print("ooo Save control table for current run")
    controlFile.close()  # you can omit in most cases as the destructor will call it

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # If no reference file is provided, compute the scale factors and write them, and scale the histos
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  else:
    print("ooo Determine Scale Factors from single pT Hard Bins")
    # Prepare histograms
    hXSecPerEvent = ROOT.TH1F("hXSecPerEvent", "hXSecPerEvent", PtHardBins+1, 0, PtHardBins+1)
    hNTrialsPerEvent = ROOT.TH1F("hNTrialsPerEvent", "hNTrialsPerEvent", PtHardBins+1, 0, PtHardBins+1)
    hScaleFactor = ROOT.TH1F("hScaleFactor", "hScaleFactor", PtHardBins+1, 0, PtHardBins+1)
    twikifile = open('twikiTableoutput.txt', 'w') # This is to be added to https://twiki.cern.ch/twiki/bin/viewauth/ALICE/JetMCProductionsCrossSections
    twikifile.write('This is to be added to https://twiki.cern.ch/twiki/bin/viewauth/ALICE/JetMCProductionsCrossSections \n')  # python will convert \n to os.linesep
    twikifile.write('| *Bin* | *pT-hard range (!GeV)* | *Factor xSec* | *Event Scale Factor (Acc)* | *Event Scale Factor (Tot)* |\n')
    twikifile.close()

    for bin in range(0,PtHardBins):
      # Label histograms
      hXSecPerEvent.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))
      hNTrialsPerEvent.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))
      hScaleFactor.GetXaxis().SetBinLabel(bin+1, "%d-%d" % (ptHardLo[bin],ptHardHi[bin]))

    # Extract cross sections from pT Hard bins
    for bin in range(0,PtHardBins):
      
      # Open input file and get relevant lists
      inputFile = "{0}/AnalysisResultsPtHard{0}.root".format(bin+1)
      print("ooo Scaling Pt-hard bin %d" % (bin+1))
      print("  o File: {}".format(inputFile))
      f = ROOT.TFile(inputFile, "UPDATE")

      qaList = f.Get(qaListNames[0])
      print("ooo Computing scaling factors with list: " + qaList.GetName())
      
      hXsecPtHard      = qaList.FindObject("hXsec")
      hTrialsPtHard    = qaList.FindObject("hNtrials")
      hNEventsTotal    = f.Get("hNEventsTot") # Get NEvents histos from file, since we undo the entry after writing it
      hNEventsAccepted = f.Get("hNEventsAcc")

      # Compute: scale factor = xsec per event / trials per event
      nEventsTot = hNEventsTotal.GetBinContent(bin+1)
      nEventsAcc = hNEventsAccepted.GetBinContent(bin+1)
      xsec = hXsecPtHard.GetBinContent(1) / hXsecPtHard.GetEntries() #entries in this case are number of files that were merged together
      #print("ooo Test entries: {0}".format(hXsecPtHard.GetEntries()))
      trials = 1.*hTrialsPtHard.GetBinContent(1) / nEventsTot
      scaleFactor = xsec/trials
      if nEventsAcc>0:
        eventScaleFactorAcc = nEventsAccAvg/nEventsAcc # also scale to account that there are different number of events in each Pt-hard bin
      else:
        eventScaleFactorAcc = 0
      if nEventsTot>0:
        eventScaleFactorTot = nNEventsTotAvg/nEventsTot
      else:
        eventScaleFactorTot = 0

      print("ooo nEventsAcc: {0}".format(nEventsAcc))
      print("ooo nEventsTot: {0}".format(nEventsTot))
      print("ooo nTrials: {0}".format(trials))
      print("ooo scaleFactor: {0}".format(scaleFactor))
      print("ooo eventScaleFactor Acc: {0}".format(eventScaleFactorAcc))
      print("ooo eventScaleFactor All: {0}".format(eventScaleFactorTot))
      print("ooo combined ScaleFactor: {0}".format(eventScaleFactorAcc*scaleFactor))
      print("")
      twikifile = open('twikiTableoutput.txt', 'a')
      twikifile.write("| *{0}* | {1} - {2} | {3} | {4} | {5} | \n".format((bin+1), ptHardLo[bin], ptHardHi[bin], scaleFactor, eventScaleFactorAcc, eventScaleFactorTot ))
      twikifile.close()

      hXSecPerEvent.Fill(bin+0.5, xsec)
      hNTrialsPerEvent.Fill(bin+0.5, trials)
      hScaleFactor.Fill(bin+0.5, scaleFactor)

      # Now, scale all the histograms
      print("ooo Scaling list: " + qaList.GetName())
      RadiusName=""
      for obj in qaList:
        name = obj.GetName()
        ScaleAllHistograms(obj, scaleFactor * eventScaleFactorAcc, f, verbose, bRemoveOutliers, outlierLimit, outlierNBinsThreshold, bin, EndPtHardBin, name, qaList.GetName())

      # Write the histograms to file
      hXSecPerEvent.Write()
      hXSecPerEvent.Reset()
      hNTrialsPerEvent.Write()
      hNTrialsPerEvent.Reset()
      hScaleFactor.Write()
      hScaleFactor.Reset()
      qaList.Write("%sScaled" % qaListNames[0], ROOT.TObject.kSingleKey)
      
      # Make sure we delete the lists and sublists, to prevent memory leak
      for obj in qaList:
        obj.SetOwner(True)
        obj.Delete()
        qaList.SetOwner(True)
        qaList.Delete()
      
      f.Close()
    twikifile.close()  # you can omit in most cases as the destructor will call it


###################################################################################
# Given event list name eventList, pT-hard bin number, and histogram hNEvents of appropriate form, this function fills
# the number of events (accepted events only if bAcceptedEventsOnly=True, otherwise all events)
def GetNEvents(eventList, bin, hNEvents, verbose, bAcceptedEventsOnly = True):
  
  inputFile = "{0}/AnalysisResultsPtHard{0}.root".format(bin+1)
  f = ROOT.TFile(inputFile, "UPDATE")
  
  qaList = f.Get(eventList)
  nEvents = 0
  if not qaList:
    print("ERROR no list found")

  # Look for the EventCutOutput from AliEventCuts, and if it doesn't exist, look for histo fHistEventCount
  eventCutList = qaList.FindObject("EventCutOutput")
  if eventCutList:
    hNEventsPtHard = eventCutList.FindObject("fCutStats")
    if bAcceptedEventsOnly:
      nEvents = hNEventsPtHard.GetBinContent(hNEventsPtHard.GetXaxis().FindBin("All cuts")) #Bin All cuts
    else:
      nEvents = hNEventsPtHard.GetBinContent(hNEventsPtHard.GetXaxis().FindBin("No cuts")) #Bin No cuts
    if bin is 0:
      if bAcceptedEventsOnly:
        print("Getting accepted number of events from EventCutOutput.")
      else:
        print("Getting total (acc+rej) number of events from EventCutOutput.")
  else:
    hNEventsPtHard = qaList.FindObject("fHistEventCount")
    if bAcceptedEventsOnly:
      nEvents = hNEventsPtHard.GetBinContent(1)
    else:
      nEvents = hNEventsPtHard.GetBinContent(1) + hNEventsPtHard.GetBinContent(2)
    if bin is 0:
      if bAcceptedEventsOnly:
        print("Getting accepted number of events from fHistEventCount.")
      else:
        print("Getting total (acc+rej) number of events from fHistEventCount.")

  if verbose:
    if bAcceptedEventsOnly:
      print("Acc. Events in bin {0} = {1}".format((bin+1), nEvents))
    else:
      print("All Events in bin {0} = {1}".format((bin+1), nEvents))

  hNEvents.Fill(bin+0.5, nEvents)
  hNEvents.Write()
  hNEvents.Fill(bin+0.5,-nEvents) # undo the entry, since we will later hadd the histograms
  f.Close()

  return nEvents

###################################################################################
# Function to iterate recursively through an object to scale all TH1/TH2/THnSparse
def ScaleAllHistograms(obj, scaleFactor, f, verbose, bRemoveOutliers=False, limit=2, nBinsThreshold=4, pTHardBin=0, EndPtHardBin=20, listName="", taskName=""):
  
  if obj.InheritsFrom(ROOT.TProfile.Class()):
    if verbose:
      print("TProfile %s not scaled..." % obj.GetName())
  # mostly for PbPb case
  elif obj.InheritsFrom(ROOT.TH3.Class()):
    obj.Sumw2()
    if bRemoveOutliers:
      name = obj.GetName()
      #only perform outlier removal on these couple histograms
      if "JESshiftEMCal" in name or "ResponseMatrixEMCal" in name or "hNEFVsPtEMCal" in name or "hMatchingDistance" in name:
        removeOutliers(pTHardBin, EndPtHardBin, obj, verbose, limit, nBinsThreshold, 3, listName, taskName)
    obj.Scale(scaleFactor)
    if verbose:
      print("TH3 %s was scaled..." % obj.GetName())

  # mostly for pp case
  elif obj.InheritsFrom(ROOT.TH2.Class()):
    obj.Sumw2()
    if bRemoveOutliers:
      name = obj.GetName()
      #only perform outlier removal on these couple histograms
      if "JESshiftEMCal" in name or "ResponseMatrixEMCal" in name or "hNEFVsPtEMCal" in name:
        removeOutliers(pTHardBin, EndPtHardBin, obj, verbose, limit, nBinsThreshold, 2, listName, taskName)
    obj.Scale(scaleFactor)
    if verbose:
      print("TH2 %s was scaled..." % obj.GetName())
  elif obj.InheritsFrom(ROOT.TH1.Class()):
    obj.Sumw2()
    obj.Scale(scaleFactor)
    if verbose:
      print("TH1 %s was scaled..." % obj.GetName())
  elif obj.InheritsFrom(ROOT.THnSparse.Class()):
    obj.Sumw2()
    obj.Scale(scaleFactor)
    if verbose:
      print("THnSparse %s was scaled..." % obj.GetName())
  else:
    if verbose:
      print("Not a histogram!")
      print(obj.GetName())
    for subobj in obj:
      ScaleAllHistograms(subobj, scaleFactor, f, verbose, bRemoveOutliers, limit, nBinsThreshold, pTHardBin, EndPtHardBin, listName, taskName)

###################################################################################
# Function to remove outliers from a TH3 (i.e. truncate the spectrum), based on projecting to the y-axis
# It truncates the 3D histogram based on when the 1D projection 4-bin moving average has been above
# "limit" for "nBinsThreshold" bins.
def removeOutliers(pTHardBin, EndPtHardBin, hist, verbose, limit=2, nBinsThreshold=4, dimension =3, listName="", taskName=""):

  #Project to the pT Truth axis
  #typically for PbPb case
  if dimension==3:
    if "hMatchingDistance" in hist.GetName():
      JetRadius = getRadiusFromlistName(taskName)
      print("list name {}, detectred Radius: {}".format(taskName,JetRadius))
      #restrict the matching distance first
      hist.GetZaxis().SetRangeUser(0,JetRadius)
      histToCheck = hist.ProjectionY("{}_projBefore_{}_{}".format(hist.GetName(), listName, taskName))
    else:
      histToCheck = hist.ProjectionY("{}_projBefore_{}_{}".format(hist.GetName(), listName, taskName))
  #plotResponsePDF(histToCheck, pTHardBin, EndPtHardBin, outputFilename, verbose, drawOptions = ""):
  #typically for pp case
  if dimension==2:
    if hist.GetName()=="hResponseMatrixEMCal":
      outputFilename = "./POST_Response_{}___{}.pdf".format(listName[:17], taskName[-20:-7])
      #plotResponsePDF(hist, pTHardBin, EndPtHardBin, outputFilename, verbose) #uses a lot of memory
      histToCheck = hist.ProjectionY("{}_projBefore_{}_{}".format(hist.GetName(), listName, taskName))
    else:
      histToCheck = hist.ProjectionX("{}_projBefore_{}_{}".format(hist.GetName(), listName, taskName))

  # Check with moving average
  foundAboveLimit = False
  cutLimitReached = False
  # The cut index is where we decided cut on that row
  cutIndex = -1
  nBinsBelowLimitAfterLimit = 0
  # nBinsThreshold= n bins that are below threshold before all bins are cut

  if verbose:
    (preMean, preMedian) = GetHistMeanAndMedian(histToCheck)
    
  for index in range(0, histToCheck.GetNcells()):
    if verbose:
      print("---------")
    avg = MovingAverage(histToCheck, index = index, numberOfCountsBelowIndex = 2, numberOfCountsAboveIndex = 2)
    if verbose:
      print("Index: {0}, Avg: {1}, BinContent: {5}, foundAboveLimit: {2}, cutIndex: {3}, cutLimitReached: {4}".format(index, avg, foundAboveLimit, cutIndex, cutLimitReached, histToCheck.GetBinContent(index)))
    if avg > limit:
      foundAboveLimit = True
        
    if not cutLimitReached:
      if foundAboveLimit and avg <= limit:
        if cutIndex == -1:
          cutIndex = index
        nBinsBelowLimitAfterLimit += 1
          
      if nBinsBelowLimitAfterLimit != 0 and avg > limit:
        # Reset
        cutIndex = -1
        nBinsBelowLimitAfterLimit = 0

      if nBinsBelowLimitAfterLimit > nBinsThreshold:
        cutLimitReached = True
        break #no need to continue the loop - we found our cut index

  # Do not perform removal here because then we miss values between the avg going below
  # the limit and crossing the nBinsThreshold
  if verbose:
    print("Hist checked: {0}, cut index: {1}".format(histToCheck.GetName(), cutIndex))
  
  # Use on both TH1 and TH2 since we don't start removing immediately, but instead only after the limit
  if cutLimitReached:
    if verbose:
      print("--> --> --> Removing outliers")
    # Check for values above which they should be removed by translating the global index
    x = ctypes.c_int(0)
    y = ctypes.c_int(0)
    z = ctypes.c_int(0)
    for index in range(0, hist.GetNcells()):
      # Get the bin x, y, z from the global bin
      hist.GetBinXYZ(index, x, y, z)
      #PbPb case
      if dimension==3:
        if y.value >= cutIndex:
          if hist.GetBinContent(index) > 1e-3:
            if verbose:
              print("Cutting for index {}. y bin {}. Cut index: {}".format(index, y, cutIndex))
            hist.SetBinContent(index, 0)
            hist.SetBinError(index, 0)
      #pp case
      if dimension==2:
        #for the response matrix the pT Truth is on the y-Axis
        if hist.GetName()=="hResponseMatrixEMCal":
          x.value=y.value
        if x.value >= cutIndex:
          if hist.GetBinContent(index) > 1e-3:
            if verbose:
              print("Cutting for index {}. x bin {}. Cut index: {}".format(index, x, cutIndex))
            hist.SetBinContent(index, 0)
            hist.SetBinError(index, 0)
  else:
    if verbose:
      print("Hist {} did not have any outliers to cut".format(hist.GetName()))

  # Check the mean and median
  # Use another temporary hist
  if dimension==3:
    histToCheckAfter = hist.ProjectionY("{}_projAfter_{}_{}".format(hist.GetName(), listName, taskName))
  if dimension==2:
    if hist.GetName()=="hResponseMatrixEMCal":
      histToCheckAfter = hist.ProjectionY("{}_projAfter_{}_{}".format(hist.GetName(), listName, taskName))
    else:
      histToCheckAfter = hist.ProjectionX("{}_projAfter_{}_{}".format(hist.GetName(), listName, taskName))

  if verbose:
    (postMean, postMedian) = GetHistMeanAndMedian(histToCheckAfter)
    print("Pre  outliers removal mean: {}, median: {}".format(preMean, preMedian))
    print("Post outliers removal mean: {}, median: {}".format(postMean, postMedian))
  outlierFilename = "./POST_{}_{}___{}.pdf".format(hist.GetName(),listName[:17], taskName[-20:-7])
  if "TrEff" not in taskName and "_0GeV" in taskName:
    if "hMatchingDistance" in hist.GetName() or "hResponseMatrixEMCal" in hist.GetName():
      plotOutlierPDF(histToCheck,histToCheckAfter, pTHardBin, EndPtHardBin, outlierFilename, verbose, "hist E", True)

########################################################################################################
def GetHistMeanAndMedian(hist):
  # Median
  # See: https://root-forum.cern.ch/t/median-of-histogram/7626/5
  x = ctypes.c_double(0)
  q = ctypes.c_double(0.5)
  # Apparently needed to be safe(?)
  hist.ComputeIntegral()
  hist.GetQuantiles(1, x, q)
    
  mean = hist.GetMean()
  return (mean, x.value)

########################################################################################################
def MovingAverage(hist, index, numberOfCountsBelowIndex = 0, numberOfCountsAboveIndex = 2):
  """
  # [-2, 2] includes -2, -1, 0, 1, 2
  """
  # Check inputs
  if numberOfCountsBelowIndex < 0 or numberOfCountsAboveIndex < 0:
    print("Moving average number of counts above or below must be >= 0. Please check the values!")
          
  count = 0.
  average = 0.
  for i in range(index - numberOfCountsBelowIndex, index + numberOfCountsAboveIndex + 1):
    # Avoid going over histogram limits
    if i < 0 or i >= hist.GetNcells():
      continue
    #print("Adding {}".format(hist.GetBinContent(i)))
    average += hist.GetBinContent(i)
    count += 1
    
  #if count != (numberOfCountsBelowIndex + numberOfCountsAboveIndex + 1):
  #    print("Count: {}, summed: {}".format(count, (numberOfCountsBelowIndex + numberOfCountsAboveIndex + 1)))
  #exit(0)

  return average / count

########################################################################################################
# Plot basic histogram    ##############################################################################
########################################################################################################
def plotOutlierPDF(h, hAfter, pTHardBin, EndPtHardBin, outputFilename, verbose, drawOptions = "", setLogy = False):
  
  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd()
  if setLogy:
    c.SetLogy()
  h.GetXaxis().SetRangeUser(0,250)
  h.Draw("hist")
  h.SetLineColor(616)
  hAfter.SetLineColor(820)
  hAfter.Draw("same hist")
  
  leg1 = ROOT.TLegend(0.17,0.7,0.83,0.85,"outlier removal of Bin {}".format(pTHardBin))
  leg1.SetFillColor(10)
  leg1.SetBorderSize(0)
  leg1.SetFillStyle(0)
  leg1.SetTextSize(0.04)
  leg1.AddEntry(h, "before", "l")
  leg1.AddEntry(hAfter, "after", "l")
  leg1.Draw("same")

  #if file does not yet exist create it
  if not os.path.exists(outputFilename):
    if verbose:
      print("Add first pT Hard bin to pdf with name: {0}".format(outputFilename))
    c.Print("{}(".format(outputFilename))#this opens a .pdf file with all figures
  #otherwise add pages to the file or close it
  else:
    if pTHardBin==EndPtHardBin-1:
      if verbose:
        print("Add last pT Hard bin to pdf with name: {0}".format(outputFilename))
      c.Print("{})".format(outputFilename)) #this closes the .pdf file with all figures
    else:
      if verbose:
        print("Add further pT Hard bin to pdf with name: {0}".format(outputFilename))
      c.Print("{}".format(outputFilename)) #this adds to the .pdf file with all figures

  c.Close()
########################################################################################################
# Plot response single pT hard    ######################################################################
########################################################################################################
def plotResponsePDF(hResponse, pTHardBin, EndPtHardBin, outputFilename, verbose):
  
  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd()
  c.SetLogz()
  hResponse.GetXaxis().SetRangeUser(0,250)
  hResponse.GetYaxis().SetRangeUser(0,250)
  hResponse.Draw("colz")

  leg1 = ROOT.TLegend(0.17,0.7,0.83,0.93,"outlier removal")
  leg1.SetFillColor(10)
  leg1.SetBorderSize(0)
  leg1.SetFillStyle(0)
  leg1.SetTextSize(0.04)
  leg1.AddEntry(hResponse, "pT bin {}".format(pTHardBin), "l")
  leg1.Draw("same")

  #if file does not yet exist create it
  if not os.path.exists(outputFilename):
    if verbose:
      print("Add first pT Hard bin to pdf with name: {0}".format(outputFilename))
    c.Print("{}(".format(outputFilename))#this opens a .pdf file with all figures
  #otherwise add pages to the file or close it
  else:
    if pTHardBin==EndPtHardBin-1:
      if verbose:
        print("Add last pT Hard bin to pdf with name: {0}".format(outputFilename))
      c.Print("{})".format(outputFilename)) #this closes the .pdf file with all figures
    else:
      if verbose:
        print("Add further pT Hard bin to pdf with name: {0}".format(outputFilename))
      c.Print("{}".format(outputFilename)) #this adds to the .pdf file with all figures
  c.Close()
########################################################################################################
# Get Jet radius from list analysis label                       #############################################
########################################################################################################
def getRadiusFromlistName(listName):

  radius = 0.0
  if "01" in listName:
    radius = 0.1
  if "02" in listName:
    radius = 0.2
  elif "03" in listName:
    radius = 0.3
  elif "04" in listName:
    radius = 0.4
  elif "05" in listName:
    radius = 0.5
  elif "06" in listName:
    radius = 0.6
  return radius

#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  print("Executing scalePtHardHistos.py...")
  print("")
  
  # Define arguments
  parser = argparse.ArgumentParser(description="Scale pT-hard bins")
  parser.add_argument("-f", "--referenceFile", action="store",
                      type=str, metavar="referenceFile",
                      default="",
                      help="Reference file to get pT-hard scale factors from")
  parser.add_argument("-s", "--StartPtHardBin", action="store",
                      type=int, metavar="StartPtHardBin",
                      default=1,
                      help="start pT bin for scaling")
  parser.add_argument("-e", "--EndPtHardBin", action="store",
                      type=int, metavar="EndPtHardBin",
                      default=20,
                      help="end pT bin for scaling")
  # Parse the arguments
  args = parser.parse_args()
  
  scalePtHardHistos(referenceFile = args.referenceFile,StartPtHardBin = args.StartPtHardBin,EndPtHardBin = args.EndPtHardBin)
