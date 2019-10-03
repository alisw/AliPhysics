#! /usr/bin/env python

# Macro to plot PWGJE QA histograms, using AliAnalysisTaskPWGJEQA.
#
# It automatically detects what to plot, based on the content of your analysis output file:
# whether to do track/calo/jet/event QA, as well as MC vs. data, PbPb vs. pp, Run1 vs. Run2, Phos vs. no Phos.
#
# To run:
#     python plotPWGJEQA.py -f "/my/dir/AnalysisResults.root" -o "/my/dir/outputQA/" -i ".png"
#
#     (or, run without options: defaults are "AnalysisResults.root" and "./outputQA/" and ".pdf")
#
# If not using standard AOD collections, you need to set the list names in the config below.
# You may need to set some of the displayed ranges on the plots.
#
# Note: It is possible you will have to change the scaling on a couple plots, to give them reasonable ranges.
#
# Note: AliAnalysisTaskPWGJEQA uses variable binning for centrality, track pT, track pT-res, and cluster E.
#       Relevant histograms are plotted below using "width" scaling option to divide by bin width, when applicable.
#
# Note: Changing the binning in the analysis task may break some functionality here.
#
# Author: James Mulligan (james.mulligan@yale.edu)
# Track plotting based in part on code from plotJETrackQA.C

# General
import os
import sys
import argparse
import itertools
import math
# ROOT
import ROOT

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

def plotPWGJEQA(inputFile, outputDir, referenceFile, fileFormat):

  # Open input file and get relevant lists
  f = ROOT.TFile(inputFile)
  
  # Set directory for QA output
  if not outputDir.endswith("/"):
      outputDir = outputDir + "/"
  if not os.path.exists(outputDir):
    os.makedirs(outputDir)
  
  # Detect whether this is a Pt-hard production (only returns true if the histos have been scaled, with scalePtHardHistos.py)
  isPtHard = False
  for key in f.GetListOfKeys():
    if "Scaled" in key.GetName():
      isPtHard = True
  print("Is Pt-hard: %s" % isPtHard)

  # Configure the plotting macro
  qaTaskBaseName = "AliAnalysisTaskPWGJEQA"

  # Input names
  tracksListName = "tracks"
  generatorTrackThnName = "tracks_PhysPrim"
  matchedTrackThnName = "tracks_Matched"

  # Handles older QA task
  if "EmcalTrackingQA" in qaTaskBaseName:
      tracksListName = "fTracks"
      generatorTrackThnName = "fParticlesPhysPrim"
      matchedTrackThnName = "fParticlesMatched"

  # Get the output list of AliAnalysisTaskPWGJEQA
  qaTaskName = determineQATaskName(qaTaskBaseName, f, isPtHard)
  print("Found qaTaskName \"{0}\"".format(qaTaskName))
  qaList = f.Get(qaTaskName)

  # If not a Pt-hard production (since it is done already), we need to set Sumw2 since we will scale and divide histograms
  if not isPtHard:
    print("Setting Sumw2 on histograms.")
    for obj in qaList:
      SetSumw2(obj)

  # Get the lists for tracks, cells, clusters, full jets, charged jets, and event QA
  trackTHnSparse = qaList.FindObject(tracksListName)
  cellQAList = qaList.FindObject("emcalCells")
  clusterQAList = qaList.FindObject("caloClusters")
  chargedJetList = qaList.FindObject("Jet_AKTChargedR020_tracks_pT0150_pt_scheme")
  fullJetList = qaList.FindObject("Jet_AKTFullR020_tracks_pT0150_caloClusters_E0300_pt_scheme")

  nEventsRef = 0
  # If reference file provided, get its analysis lists
  qaListRef = ""
  trackTHnSparseRef = ""
  clusterQAListRef = ""
  cellQAListRef = ""
  chargedJetListRef = ""
  fullJetListRef = ""
  if referenceFile:
    fRef = ROOT.TFile(referenceFile)
    qaListRef = fRef.Get(qaTaskName)
    if not isPtHard:
      print("Setting Sumw2 on reference histograms.")
      for obj in qaListRef:
        SetSumw2(obj)
    trackTHnSparseRef = qaListRef.FindObject(tracksListName)
    trackTHnSparseRef.SetName("trackRef")
    clusterQAListRef = qaListRef.FindObject("caloClusters")
    clusterQAListRef.SetName("caloClustersRef")
    cellQAListRef = qaListRef.FindObject("emcalCells")
    cellQAListRef.SetName("emcalCellsRef")
    chargedJetListRef = qaListRef.FindObject("Jet_AKTChargedR020_tracks_pT0150_pt_scheme")
    chargedJetListRef.SetName("chargedJetListRef")
    fullJetListRef = qaListRef.FindObject("Jet_AKTFullR020_tracks_pT0150_caloClusters_E0300_pt_scheme")
    fullJetListRef.SetName("fullJetListRef")

    histNEventRef = qaListRef.FindObject("fHistEventCount")
    nEventsRef = histNEventRef.GetBinContent(1)
    print("N events ref: %d" % nEventsRef)

  # Get number of events
  histNEvent = qaList.FindObject("fHistEventCount")
  nEvents = histNEvent.GetBinContent(1)
  print("N events: %d" % nEvents)

  # Set config: ispp, isMC, isRun2, includePhos
  if qaList.FindObject("fHistCentrality"):
    ispp = False
  else:
    ispp = True
  print("Is pp: %s" % ispp)

  if qaList.FindObject(generatorTrackThnName):
    isMC = True
  else:
    isMC = False
  print("Is MC: %s" % isMC)

  if clusterQAList:
      clusterTHnSparse = clusterQAList.FindObject("clusterObservables")
      if ispp:
        hClusterType = clusterTHnSparse.Projection(3)
      else:
        hClusterType = clusterTHnSparse.Projection(4)
      isRun2 = hClusterType.GetBinContent(2) > 0
      includePhos = hClusterType.GetBinContent(3) > 0
      print("Is Run 2: %s" % isRun2)
      print("Include Phos: %s" % includePhos)
  else:
      isRun2 = False
      includePhos = False

  # Plotting options
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)

  # Plot QA
  print("Plotting QA...")
  if trackTHnSparse:
    plotTrackQA(ispp, isMC, trackTHnSparse, generatorTrackThnName, matchedTrackThnName, qaList, nEvents, outputDir, qaListRef, trackTHnSparseRef, nEventsRef, fileFormat)
  if clusterQAList:
    plotCaloQA(ispp, isRun2, includePhos, clusterQAList, cellQAList, nEvents, outputDir, clusterQAListRef, cellQAListRef, nEventsRef, fileFormat)
  if chargedJetList:
    plotChargedJetQA(ispp, isPtHard, chargedJetList, outputDir, chargedJetListRef, nEvents, nEventsRef, fileFormat)
  if fullJetList:
    plotFullJetQA(ispp, isPtHard, isRun2, includePhos, fullJetList, outputDir, fullJetListRef, nEvents, nEventsRef, fileFormat)
  if qaList.FindObject("eventQA"):
    plotEventQA(ispp, isRun2, includePhos, qaList, outputDir, fileFormat)
  if isPtHard:
    plotPtHard(f, qaList,nEvents, qaListRef ,nEventsRef, outputDir, fileFormat)

def determineQATaskName(qaTaskBaseName, f, isPtHard):
    """ Determine the task name based on a wide variety of possible names.

    Since the task name varies depending on what input objects are included,
    we need to guess the name.

    Args:
        qaTaskBaseName (str): Base name of the QA task without any of the input object names
        f (TFile): Root file containing the QA task

    """
    # Get all task names stored in the input file
    possibleTaskNames = [key.GetName() for key in f.GetListOfKeys()]

    # Possible input object names
    tracksName = "tracks"
    mcTracksName = "mcparticles"
    cellsName = "emcalCells"
    clustersName = "caloClusters"

    # Compile into a list for easy processing
    possibleNames = [tracksName, mcTracksName, cellsName, clustersName]
    suffix = "histos"
    if isPtHard:
      suffix = "histosScaled"

    for length in range(0, len(possibleNames)+1):
        for elements in itertools.permutations(possibleNames, length):
            joined = "_".join(elements)
            testTaskName = qaTaskBaseName
            if joined:
                testTaskName += "_" + joined
            # Also Try ESD
            testTaskNameESD = testTaskName.replace("emcalCells", "EMCALCells").replace("caloClusters", "CaloClusters").replace("tracks", "Tracks").replace("mcparticles", "MCParticles")

            for taskName in [testTaskName, testTaskNameESD]:
              taskName = "{0}_{1}".format(taskName, suffix)
              if taskName in possibleTaskNames:
                  return taskName

    print("Could not determine QA task name! Please check your spelling!")
    exit(1)

########################################################################################################
# Plot track histograms   ##############################################################################
########################################################################################################
def plotTrackQA(ispp, isMC, trackTHnSparse, generatorTrackThnName, matchedTrackThnName, qaList, nEvents, outputDir, qaListRef, trackTHnSparseRef, nEventsRef, fileFormat):
  
  # Create subdirectory for Tracks
  outputDirTracks = outputDir + "Tracks/"
  if not os.path.exists(outputDirTracks):
    os.makedirs(outputDirTracks)

  # trackTHnSparse consists of (Centrality, Pt, Eta, Phi, Track type, sigma(pT)/pT)
  if isMC:
    generatorTHnSparse = qaList.FindObject(generatorTrackThnName) # (Centrality, Pt, Eta, Phi, findable)
    matchedTHnSparse = qaList.FindObject(matchedTrackThnName) # (Pt-gen, Eta-gen, Phi-gen, Pt-det, Eta-det, Phi-det, (pT-gen - pT-det)/pT-det, Track type)

  #---------------------------------------------------------------------------------------------------
  #                        phi distribution of hybrid tracks
  #---------------------------------------------------------------------------------------------------

  c1 = ROOT.TCanvas("c1","c1: Phi",600,450)
  c1.cd()

  # Project to (Phi, Track type)
  if ispp:
    hPhiTracktype = trackTHnSparse.Projection(2,3)
  else:
    hPhiTracktype = trackTHnSparse.Projection(3,4)

  hPhiGlobal = hPhiTracktype.ProjectionY("PhiGlobal", 1, 1)
  hPhiComplementary = hPhiTracktype.ProjectionY("PhiComplementary", 2, 2)

  hPhiGlobal.SetLineColor(2)
  hPhiGlobal.SetLineWidth(3)
  hPhiGlobal.SetLineStyle(1)
  hPhiComplementary.SetLineStyle(1)
  hPhiComplementary.SetLineColor(4)
  hPhiComplementary.SetLineWidth(3)

  hPhiSum = hPhiGlobal.Clone()
  hPhiSum.Add(hPhiComplementary)
  hPhiSum.SetTitle("hPhiSum")
  hPhiSum.SetName("hPhiSum")
  hPhiSum.SetLineColor(1)
  hPhiSum.SetMarkerColor(1)
  hPhiSum.SetLineStyle(1)

  hPhiGlobal.Scale(1./nEvents)
  hPhiComplementary.Scale(1./nEvents)
  hPhiSum.Scale(1./nEvents)

  hPhiGlobal.SetTitle("#phi Distribution of Hybrid Tracks")
  hPhiGlobal.GetYaxis().SetTitle("#frac{1}{N_{evts}} #frac{dN}{d#phi}")
  hPhiGlobal.GetYaxis().SetTitleSize(0.06)
  hPhiGlobal.GetXaxis().SetTitleSize(0.06)
  hPhiGlobal.GetXaxis().SetTitleOffset(0.5)
  hPhiGlobal.GetYaxis().SetRangeUser(0,15.)
  if ispp:
    hPhiGlobal.GetYaxis().SetRangeUser(0,0.2)
    if isMC:
      hPhiGlobal.GetYaxis().SetRangeUser(0,0.25)
  ROOT.gPad.SetLeftMargin(0.15)
  ROOT.gPad.SetRightMargin(0.05)
  ROOT.gPad.SetBottomMargin(0.13)
  ROOT.gPad.SetTopMargin(0.05)

  hPhiGlobal.Draw("hist")
  hPhiComplementary.Draw("hist same")
  hPhiSum.Draw("hist same")

  leg1 = ROOT.TLegend(0.17,0.7,0.83,0.93,"Hybrid tracks")
  leg1.SetFillColor(10)
  leg1.SetBorderSize(0)
  leg1.SetFillStyle(0)
  leg1.SetTextSize(0.04)
  leg1.AddEntry(hPhiGlobal, "w/ SPD & ITSrefit", "l")
  leg1.AddEntry(hPhiComplementary, "w/o SPD & w/ ITSrefit", "l")
  leg1.AddEntry(hPhiSum, "sum", "l")
  leg1.Draw("same")

  textNEvents = ROOT.TLatex()
  textNEvents.SetNDC()
  c1.cd()
  textNEvents.DrawLatex(0.52,0.68,"#it{N}_{events} = %d" % nEvents)

  outputFilename = os.path.join(outputDirTracks, "hTrackPhi" + fileFormat)
  c1.SaveAs(outputFilename)

  # Also plot the TH2 phi vs. pT -- make sure that phi is uniform at all pT
  # Project to (Pt, Phi)
  if ispp:
    hPhiPtSum = trackTHnSparse.Projection(2,0)
  else:
    hPhiPtSum = trackTHnSparse.Projection(3,1)
  hPhiPtSum.Scale(1.,"width")
  hPhiPtSum.GetZaxis().SetRangeUser(1e-7, 3e5)
  outputFilename = os.path.join(outputDirTracks, "hTrackPhiPt" + fileFormat)
  plotHist(hPhiPtSum, outputFilename, "colz", False, True)

  #---------------------------------------------------------------------------------------------------
  #                        pT distribution of hybrid tracks
  #---------------------------------------------------------------------------------------------------
  
  # Project to (Pt, Track type)
  if ispp:
    hPtTracktype = trackTHnSparse.Projection(0,3)
  else:
    hPtTracktype = trackTHnSparse.Projection(1,4)

  hPtGlobal = hPtTracktype.ProjectionY("PtGlobal", 1, 1)
  hPtComplementary = hPtTracktype.ProjectionY("PtComplementary", 2, 2)

  hPtSum = hPtGlobal.Clone()
  hPtSum.Add(hPtComplementary)

  # If reference distribution supplied, project to (Pt, Track type)
  hPtSumRef = ""
  if trackTHnSparseRef and qaListRef:
    if ispp:
      hPtTracktypeRef = trackTHnSparseRef.Projection(0,3)
    else:
      hPtTracktypeRef = trackTHnSparseRef.Projection(1,4)
      
    hPtGlobalRef = hPtTracktypeRef.ProjectionY("PtGlobalRef", 1, 1)
    hPtComplementaryRef = hPtTracktypeRef.ProjectionY("PtComplementaryRef", 2, 2)
    hPtSumRef = hPtGlobalRef.Clone()
    hPtSumRef.Add(hPtComplementaryRef)
      
  outputFilename = os.path.join(outputDirTracks, "hTrackPt" + fileFormat)
  xRangeMax = 100
  yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [GeV^{-1}]"
  legendTitle = "Track selection"
  legendRunLabel = "Hybrid tracks"
  legendRefLabel = "Hybrid tracks, all runs"
  ratioYAxisTitle = "Ratio: run / all runs"
    
  plotSpectra(hPtSum, hPtSumRef, hPtGlobal, hPtComplementary, nEvents, nEventsRef, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, "width", "w/ SPD & ITSrefit", "w/o SPD & w/ ITSrefit")
  
  # Plot also ratio of central track spectrum to peripheral track spectrum
  trackTHnSparse.GetAxis(0).SetRangeUser(0,10)
  hPt010 = trackTHnSparse.Projection(1)
  hPt010.SetName("hPt010")
  
  trackTHnSparse.GetAxis(0).SetRangeUser(50,90)
  hPt5090 = trackTHnSparse.Projection(1)
  hPt5090.SetName("hPt5090")
  
  outputFilename = os.path.join(outputDirTracks, "hTrackPtRatio" + fileFormat)
  xRangeMax = 75
  yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dE_{T}} [GeV^{-1}]"
  legendTitle = "Tracks"
  h1legendLabel = "50-90 %"
  h2legendLabel = "0-10 %"
  ratioYAxisTitle = "Central / Peripheral"
  yRatioMax = 12
  plotSpectraCent(hPt5090, hPt010, "", nEvents, ispp, outputFilename, xRangeMax, yAxisTitle, ratioYAxisTitle, legendTitle, h1legendLabel, h2legendLabel, "", "width", yRatioMax)
  
  trackTHnSparse.GetAxis(0).SetRangeUser(0,90)

  #---------------------------------------------------------------------------------------------------
  #                        pT resolution of hybrid tracks -- from track fitting
  #---------------------------------------------------------------------------------------------------

  c5 = ROOT.TCanvas("c5","c5: pT resolution",600,450)
  c5.cd()

  # Project to (Pt, Track type, pT resolution)
  if ispp:
    hPtTracktypePtSigma1Pt = trackTHnSparse.Projection(0,3,4)
  else:
    hPtTracktypePtSigma1Pt = trackTHnSparse.Projection(1,4,5)

  # Project to global tracks and take profile, to get the pT resolution as a function of pT (Profile of pT vs pT*sigma(1/pT), i.e. pT vs sigma(pT)/pT)
  # Note: No need to scale by bin width (despite pt-res having variable binning), since we take a profile (rather than e.g. plot a spectrum).
  hPtTracktypePtSigma1Pt.GetYaxis().SetRange(1,1)
  hPtPtSigma1PtGlobal = hPtTracktypePtSigma1Pt.Project3D("zx")
  hPtPtSigma1PtGlobal.SetName("hPtPtSigma1PtGlobal")
  profPtPtSigma1PtGlobal = hPtPtSigma1PtGlobal.ProfileX()
  profPtPtSigma1PtGlobal.SetName("profPtPtSigma1PtGlobal")
  profPtPtSigma1PtGlobal.SetLineColor(2)
  profPtPtSigma1PtGlobal.SetLineWidth(3)
  profPtPtSigma1PtGlobal.SetMarkerStyle(21)
  profPtPtSigma1PtGlobal.SetMarkerColor(2)
  profPtPtSigma1PtGlobal.SetMaximum(0.3)
  #profPtPtSigma1PtGlobal.GetYaxis().SetTitle("#it{p}_{T} #times #sigma(1/#it{p}_{T})")
  profPtPtSigma1PtGlobal.GetYaxis().SetTitle(" #sigma(#it{p}_{T}) / #it{p}_{T}")
  profPtPtSigma1PtGlobal.GetXaxis().SetTitleSize(0.06)
  profPtPtSigma1PtGlobal.GetYaxis().SetTitleSize(0.06)
  profPtPtSigma1PtGlobal.GetYaxis().SetRangeUser(0,0.15)
  ROOT.gPad.SetLeftMargin(0.15)
  ROOT.gPad.SetRightMargin(0.05)
  ROOT.gPad.SetBottomMargin(0.14)
  ROOT.gPad.SetTopMargin(0.05)
  #profPtPtSigma1PtGlobal.GetXaxis().SetRangeUser(0, 100)
  profPtPtSigma1PtGlobal.Draw()

  # Project to complementary tracks and take profile
  hPtTracktypePtSigma1Pt.GetYaxis().SetRange(2,2)
  hPtPtSigma1PtComplementary = hPtTracktypePtSigma1Pt.Project3D("zx")
  hPtPtSigma1PtComplementary.SetName("hPtPtSigma1PtComplementary")
  profPtPtSigma1PtComplementary = hPtPtSigma1PtComplementary.ProfileX()
  profPtPtSigma1PtComplementary.SetName("profPtPtSigma1PtComplementary")
  profPtPtSigma1PtComplementary.SetLineColor(4)
  profPtPtSigma1PtComplementary.SetLineWidth(3)
  profPtPtSigma1PtComplementary.SetMarkerStyle(24)
  profPtPtSigma1PtComplementary.SetMarkerColor(4)
  profPtPtSigma1PtComplementary.Draw("same")

  leg3 = ROOT.TLegend(0.21,0.6,0.88,0.88,"Hybrid tracks")
  leg3.SetFillColor(10)
  leg3.SetBorderSize(0)
  leg3.SetFillStyle(0)
  leg3.SetTextSize(0.04)
  leg3.AddEntry(profPtPtSigma1PtGlobal, "w/ SPD & ITSrefit", "lp")
  leg3.AddEntry(profPtPtSigma1PtComplementary, "w/o SPD & w/ ITSrefit", "lp")
  leg3.Draw("same")

  outputFilename = os.path.join(outputDirTracks, "profTrackPtResolution" + fileFormat)
  c5.SaveAs(outputFilename)

  #---------------------------------------------------------------------------------------------------
  #                        pT resolution of hybrid tracks -- from MC
  #---------------------------------------------------------------------------------------------------
  # (the error bars on this histogram, which denote the resolution, are not working at present...)
  if isMC:
    # Plot distribution (pT-gen - pT-det)/pT-det
    c25 = ROOT.TCanvas("c25","c25: pT Res Dist MC",600,450)
    c25.cd()
    c25.SetLogy()
    if ispp:
      hPtRes = matchedTHnSparse.Projection(6)
    else:
      hPtRes = matchedTHnSparse.Projection(7)
    hPtRes.GetYaxis().SetTitle("counts")
    hPtRes.Draw("hist E")
    outputFilename = os.path.join(outputDirTracks, "hTrackPtResolutionMC" + fileFormat)
    c25.SaveAs(outputFilename)

    # Plot mean of the distribution as a function of pT, with error bars as the standard deviation of the distribution
    c24 = ROOT.TCanvas("c24","c24: pT Resolution MC",600,450)
    c24.cd()
    # Project to (Pt, pT resolution, Track type)
    if ispp:
      hPtTracktypePtRes = matchedTHnSparse.Projection(3,7,6)
    else:
      hPtTracktypePtRes = matchedTHnSparse.Projection(4,8,7)
    # Project to global tracks and take profile, to get the pT resolution as a function of pT
    hPtTracktypePtRes.GetYaxis().SetRange(1,1)
    hPtPtResGlobal = hPtTracktypePtRes.Project3D("zx")
    hPtPtResGlobal.SetName("hPtPtResGlobal")
    profPtPtResGlobal = hPtPtResGlobal.ProfileX("prof",1,-1,"s") # set errors to standard deviation (rather than standard error on mean)
    profPtPtResGlobal.SetName("profPtPtResGlobal")
    profPtPtResGlobal.SetLineColor(2)
    profPtPtResGlobal.SetLineWidth(3)
    profPtPtResGlobal.SetMarkerStyle(21)
    profPtPtResGlobal.SetMarkerColor(2)
    profPtPtResGlobal.SetMaximum(0.3)
    profPtPtResGlobal.GetYaxis().SetTitle("(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{det}")
    profPtPtResGlobal.GetXaxis().SetTitleSize(0.06)
    profPtPtResGlobal.GetYaxis().SetTitleSize(0.06)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTopMargin(0.05)
    profPtPtResGlobal.GetYaxis().SetRangeUser(-0.5, 1)
    profPtPtResGlobal.GetXaxis().SetRangeUser(0,100)
    profPtPtResGlobal.Draw("E")
    # Project to complementary tracks and take profile
    hPtTracktypePtRes.GetYaxis().SetRange(2,2)
    hPtPtResComplementary = hPtTracktypePtRes.Project3D("zx")
    hPtPtResComplementary.SetName("hPtPtResComplementary")
    profPtPtResComplementary = hPtPtResComplementary.ProfileX()
    profPtPtResComplementary.SetName("profPtPtResComplementary")
    profPtPtResComplementary.SetLineColor(4)
    profPtPtResComplementary.SetLineWidth(3)
    profPtPtResComplementary.SetMarkerStyle(24)
    profPtPtResComplementary.SetMarkerColor(4)
    profPtPtResComplementary.Draw("same E")
    leg3 = ROOT.TLegend(0.21,0.6,0.88,0.88,"Hybrid tracks")
    leg3.SetFillColor(10)
    leg3.SetBorderSize(0)
    leg3.SetFillStyle(0)
    leg3.SetTextSize(0.04)
    leg3.AddEntry(profPtPtResGlobal, "w/ SPD & ITSrefit", "lp")
    leg3.AddEntry(profPtPtResComplementary, "w/o SPD & w/ ITSrefit", "lp")
    leg3.Draw("hist same")
    textPtRes = ROOT.TLatex()
    textPtRes.SetNDC()
    textPtRes.DrawLatex(0.45,0.9,"Data points: mean value")
    textPtRes.DrawLatex(0.45, 0.8,"Error bars: stdev (resolution)")
    outputFilename = os.path.join(outputDirTracks, "profTrackPtResolutionMC" + fileFormat)
    c24.SaveAs(outputFilename)

  #---------------------------------------------------------------------------------------------------
  #                        Tracking efficiency
  #---------------------------------------------------------------------------------------------------

  if isMC:
    # Plot ratio of pT-gen-matched to pT-gen
    c26 = ROOT.TCanvas("c26","c26: TrackingEfficiency",600,450)
    c26.cd()

    for dim in ["1D", "2D"]:
      if dim == "1D":
        # 1D case
        if ispp:
          hPtGenMatched = matchedTHnSparse.Projection(0)
          hPtGen1D = generatorTHnSparse.Projection(0, 3)
        else:
          hPtGenMatched = matchedTHnSparse.Projection(1)
          hPtGen1D = generatorTHnSparse.Projection(1, 4)
        hPtGenFindable = hPtGen1D.ProjectionY("trackEff", 2, 2)
      elif dim == "2D":
        # 2D case
        if ispp:
          hPtGenMatched = matchedTHnSparse.Projection(1, 0)
          hPtGen2D = generatorTHnSparse.Projection(0, 1, 3)
        else:
          hPtGenMatched = matchedTHnSparse.Projection(2, 1)
          hPtGen2D = generatorTHnSparse.Projection(1, 2, 4)
        hPtGen2D.GetZaxis().SetRange(2, 2)
        hPtGenFindable = hPtGen2D.Project3D("yx")

      hTrackingEfficiency = hPtGenMatched.Clone()
      hTrackingEfficiency.Divide(hPtGenMatched, hPtGenFindable, 1., 1., "B")
      hTrackingEfficiency.SetMarkerStyle(21)
      hTrackingEfficiency.SetMarkerColor(2)
      
      if hTrackingEfficiency.InheritsFrom(ROOT.TH2.Class()):
        hTrackingEfficiency.Draw("colz")
      else:
        hTrackingEfficiency.GetYaxis().SetTitle("Tracking Efficiency")
        hTrackingEfficiency.GetYaxis().SetRangeUser(0.6,1)
        hTrackingEfficiency.GetXaxis().SetRangeUser(0,50)
        hTrackingEfficiency.Draw("P")

      outputFilename = os.path.join(outputDirTracks, "hTrackingEfficiency{0}".format(dim) + fileFormat)
      c26.SaveAs(outputFilename)

  #---------------------------------------------------------------------------------------------------
  #                        eta distribution of hybrid tracks
  #---------------------------------------------------------------------------------------------------

  c6 = ROOT.TCanvas("c6","c6: Eta",600,450)
  c6.cd()

  # Project to (Eta, Track type)
  if ispp:
    hEtaTracktype = trackTHnSparse.Projection(1,3)
  else:
    hEtaTracktype = trackTHnSparse.Projection(2,4)

  hEtaGlobal = hEtaTracktype.ProjectionY("EtaGlobal", 1, 1)
  hEtaComplementary = hEtaTracktype.ProjectionY("EtaComplementary", 2, 2)

  hEtaGlobal.SetLineColor(2)
  hEtaGlobal.SetLineWidth(3)
  hEtaGlobal.SetLineStyle(1)
  hEtaComplementary.SetLineStyle(1)
  hEtaComplementary.SetLineColor(4)
  hEtaComplementary.SetLineWidth(3)

  hEtaSum = hEtaGlobal.Clone()
  hEtaSum.Add(hEtaComplementary)
  hEtaSum.SetTitle("hEtaSum")
  hEtaSum.SetName("hEtaSum")
  hEtaSum.SetLineColor(1)
  hEtaSum.SetMarkerColor(1)
  hEtaSum.SetLineStyle(1)

  hEtaGlobal.Scale(1./nEvents)
  hEtaComplementary.Scale(1./nEvents)
  hEtaSum.Scale(1./nEvents)

  hEtaGlobal.SetTitle("#eta Distribution of Hybrid Tracks")
  hEtaGlobal.GetYaxis().SetTitle("#frac{1}{N_{evts}} #frac{dN}{d#eta}")
  hEtaGlobal.GetYaxis().SetTitleSize(0.06)
  hEtaGlobal.GetXaxis().SetTitleSize(0.06)
  hEtaGlobal.GetXaxis().SetTitleOffset(0.7)
  hEtaGlobal.GetYaxis().SetRangeUser(0,20.)
  if ispp:
    hEtaGlobal.GetYaxis().SetRangeUser(0,0.2)
    if isMC:
      hEtaGlobal.GetYaxis().SetRangeUser(0,0.3)
  ROOT.gPad.SetLeftMargin(0.15)
  ROOT.gPad.SetRightMargin(0.05)
  ROOT.gPad.SetBottomMargin(0.13)
  ROOT.gPad.SetTopMargin(0.05)

  hEtaGlobal.Draw("hist")
  hEtaComplementary.Draw("hist same")
  hEtaSum.Draw("hist same")

  leg1 = ROOT.TLegend(0.17,0.7,0.83,0.93,"Hybrid tracks")
  leg1.SetFillColor(10)
  leg1.SetBorderSize(0)
  leg1.SetFillStyle(0)
  leg1.SetTextSize(0.04)
  leg1.AddEntry(hEtaGlobal, "w/ SPD & ITSrefit", "l")
  leg1.AddEntry(hEtaComplementary, "w/o SPD & w/ ITSrefit", "l")
  leg1.AddEntry(hEtaSum, "sum", "l")
  leg1.Draw("same")

  textNEvents = ROOT.TLatex()
  textNEvents.SetNDC()
  textNEvents.DrawLatex(0.65,0.87,"#it{N}_{events} = %d" % nEvents)

  outputFilename = os.path.join(outputDirTracks, "hTrackEta" + fileFormat)
  c6.SaveAs(outputFilename)

  # Also plot the TH2 eta vs. pT -- make sure that eta is uniform at all pT
  # Project to (Pt, Eta)
  if ispp:
    hEtaPtSum = trackTHnSparse.Projection(1,0)
  else:
    hEtaPtSum = trackTHnSparse.Projection(2,1)
  hEtaPtSum.Scale(1.,"width")
  hEtaPtSum.GetZaxis().SetRangeUser(1e-7, 3e5)
  outputFilename = os.path.join(outputDirTracks, "hTrackEtaPt" + fileFormat)
  plotHist(hEtaPtSum, outputFilename, "colz", False, True)

  #---------------------------------------------------------------------------------------------------
  #                        eta-phi distribution of hybrid tracks
  #---------------------------------------------------------------------------------------------------

  # Project to (Eta, Phi)
  if ispp:
    hEtaPhiSum = trackTHnSparse.Projection(1,2)
  else:
    hEtaPhiSum = trackTHnSparse.Projection(2,3)
  hEtaPhiSum.SetName("hEtaPhiSum")
  outputFilename = os.path.join(outputDirTracks, "hTrackEtaPhi" + fileFormat)
  plotHist(hEtaPhiSum, outputFilename, "colz")

  # And plot the eta-phi distribution for high-pT tracks
  ROOT.gStyle.SetOptTitle(1)
  if ispp:
    trackTHnSparse.GetAxis(0).SetRangeUser(10,150)
    hTrackEtaPhiHighPt = trackTHnSparse.Projection(1,2)
  else:
    trackTHnSparse.GetAxis(1).SetRangeUser(10,150)
    hTrackEtaPhiHighPt = trackTHnSparse.Projection(2,3)
  hTrackEtaPhiHighPt.SetTitle("Track Occupancy, p_{T} > 10 GeV")
  outputFilename = os.path.join(outputDirTracks, "hTrackEtaPhiHighPt" + fileFormat)
  plotHist(hTrackEtaPhiHighPt, outputFilename, "colz")
  if ispp:
    trackTHnSparse.GetAxis(0).SetRangeUser(0,150)
  else:
    trackTHnSparse.GetAxis(1).SetRangeUser(0,150)
  ROOT.gStyle.SetOptTitle(0)

########################################################################################################
# Plot cluster histograms ##############################################################################
########################################################################################################
def plotCaloQA(ispp, isRun2, includePhos, clusterQAList, cellQAList, nEvents, outputDir, clusterQAListRef, cellQAListRef, nEventsRef, fileFormat):
  
  # Create subdirectory for Cells, Clusters
  outputDirCells = outputDir + "Cells/"
  if not os.path.exists(outputDirCells):
    os.makedirs(outputDirCells)
  outputDirClusters = outputDir + "Clusters/"
  if not os.path.exists(outputDirClusters):
    os.makedirs(outputDirClusters)
  
  clusterTHnSparse = clusterQAList.FindObject("clusterObservables")
  # (Centrality, E_clus, eta, phi, clusterType)

  if clusterQAListRef:
    clusterTHnSparseRef = clusterQAListRef.FindObject("clusterObservables")

  # Plot Eta-Phi of ALL CLUSTERS -----------------------------------------------------

  # Project to (Eta, Phi)
  if ispp:
    hClusPhiEta = clusterTHnSparse.Projection(2,1)
  else:
    hClusPhiEta = clusterTHnSparse.Projection(3,2)
  hClusPhiEta.SetName("clusterEMCalObservables_proj_eta_phi")
  outputFilename = os.path.join(outputDirClusters, "hClusPhiEta" + fileFormat)
  hClusPhiEta.GetXaxis().SetRangeUser(-1.5,0.8)#ELIANE -0.8,0.8
  hClusPhiEta.GetYaxis().SetRangeUser(1.2,5.8)
  plotHist(hClusPhiEta, outputFilename, "colz")

  # Plot ratio to reference run, if supplied
  if clusterQAListRef:
    if ispp:
      hClusPhiEtaRef = clusterTHnSparseRef.Projection(2,1)
    else:
      hClusPhiEtaRef = clusterTHnSparseRef.Projection(3,2)

    hClusPhiEta.Scale(1./nEvents)
    hClusPhiEtaRef.Scale(1./nEventsRef)
    hClusPhiEtaRatio = hClusPhiEta.Clone()
    hClusPhiEtaRatio.Divide(hClusPhiEtaRef)
    ROOT.gStyle.SetOptTitle(1)
    hClusPhiEtaRatio.SetTitle("Cluster Occupancy (per event): Current Run / All Runs")
    outputFilename = os.path.join(outputDirClusters, "hClusPhiEtaRatio" + fileFormat)

    plotHist(hClusPhiEtaRatio, outputFilename, "colz", False, True)
    ROOT.gStyle.SetOptTitle(0)

  # Plot EMCAL CLUSTERS --------------------------------------------------------------

  # Project to (Energy, Eta, Phi, EMCal Cluster type)
  if ispp:
    clusterTHnSparse.GetAxis(3).SetRange(1,1)
    hClusEMCalEta = clusterTHnSparse.Projection(1)
    hClusEMCalPhi = clusterTHnSparse.Projection(2)
    hClusEMCalEnergy = clusterTHnSparse.Projection(0)
  else:
    clusterTHnSparse.GetAxis(4).SetRange(1,1)
    hClusEMCalEta = clusterTHnSparse.Projection(2)
    hClusEMCalPhi = clusterTHnSparse.Projection(3)
    hClusEMCalEnergy = clusterTHnSparse.Projection(1)
  hClusEMCalEta.SetName("ClusEtaEmcal")
  hClusEMCalPhi.SetName("ClusPhiEmcal")
  hClusEMCalEnergy.SetName("ClusEnergyEmcal")

  # Plot phi distribution
  outputFilename = os.path.join(outputDirClusters, "hClusEMCalPhi" + fileFormat)
  plotHist(hClusEMCalPhi, outputFilename, "hist E")

  # Plot eta distribution
  outputFilename = os.path.join(outputDirClusters, "hClusEMCalEta" + fileFormat)
  plotHist(hClusEMCalEta, outputFilename, "hist E")

  # Plot energy distribution
  hClusEMCalEnergy.SetName("hClusEMCalEnergy")

  hClusEMCalEnergyRef = ""
  if clusterQAListRef:
    if ispp:
      clusterTHnSparseRef.GetAxis(3).SetRange(1,1)
      hClusEMCalEnergyRef = clusterTHnSparseRef.Projection(0)
    else:
      clusterTHnSparseRef.GetAxis(4).SetRange(1,1)
      hClusEMCalEnergyRef = clusterTHnSparseRef.Projection(1)
    hClusEMCalEnergyRef.SetName("clusterEMCalObservablesRef_proj_energy")

  outputFilename = os.path.join(outputDirClusters, "hClusEMCalEnergy" + fileFormat)
  xRangeMax = 100
  if ispp:
    xRangeMax = 80
  yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dE_{T}} [GeV^{-1}]"
  legendTitle = "EMCal Clusters"
  legendRunLabel = "Current run"
  legendRefLabel = "All runs"
  ratioYAxisTitle = "Ratio: run / all runs"

  plotSpectra(hClusEMCalEnergy, hClusEMCalEnergyRef, 0, 0, nEvents, nEventsRef, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, "width")

  # Plot DCAL CLUSTERS (if isRun2) ----------------------------------------------------

  if isRun2:
    # Project to (Energy, Eta, Phi, DCal Cluster type)
    if ispp:
      clusterTHnSparse.GetAxis(3).SetRange(2,2)
      hClusDCalEta = clusterTHnSparse.Projection(1)
      hClusDCalPhi = clusterTHnSparse.Projection(2)
      hClusDCalEnergy = clusterTHnSparse.Projection(0)
    else:
      clusterTHnSparse.GetAxis(4).SetRange(2,2)
      hClusDCalEta = clusterTHnSparse.Projection(2)
      hClusDCalPhi = clusterTHnSparse.Projection(3)
      hClusDCalEnergy = clusterTHnSparse.Projection(1)
    hClusDCalEta.SetName("ClusEtaDcal")
    hClusDCalPhi.SetName("ClusPhiDcal")
    hClusDCalEnergy.SetName("ClusEnergyDcal")

    # Plot phi distribution
    outputFilename = os.path.join(outputDirClusters, "hClusDCalPhi" + fileFormat)
    plotHist(hClusDCalPhi, outputFilename, "hist E")

    # Plot eta distribution
    outputFilename = os.path.join(outputDirClusters, "hClusDCalEta" + fileFormat)
    plotHist(hClusDCalEta, outputFilename, "hist E")

    # Plot energy distribution
    hClusDCalEnergy.SetName("hClusDCalEnergy")
    
    hClusDCalEnergyRef = ""
    if clusterQAListRef:
      if ispp:
        clusterTHnSparseRef.GetAxis(3).SetRange(2,2)
        hClusDCalEnergyRef = clusterTHnSparseRef.Projection(0)
      else:
        clusterTHnSparseRef.GetAxis(4).SetRange(2,2)
        hClusDCalEnergyRef = clusterTHnSparseRef.Projection(1)
      hClusDCalEnergyRef.SetName("clusterDCalObservablesRef_proj_energy")

    outputFilename = os.path.join(outputDirClusters, "hClusDCalEnergy" + fileFormat)
    xRangeMax = 100
    if ispp:
      xRangeMax = 50
    yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dE_{T}} [GeV^{-1}]"
    legendTitle = "DCal Clusters"
    legendRunLabel = "Current run"
    legendRefLabel = "All runs"
    ratioYAxisTitle = "Ratio: run / all runs"
    
    plotSpectra(hClusDCalEnergy, hClusDCalEnergyRef, 0, 0, nEvents, nEventsRef, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, "width")

  # Plot PHOS CLUSTERS (if includePhos) -----------------------------------------------

  if includePhos:
    # Project to (Energy, Eta, Phi, PHOS Cluster type)
    if ispp:
      clusterTHnSparse.GetAxis(3).SetRange(3,3)
      hClusPHOSEta = clusterTHnSparse.Projection(1)
      hClusPHOSPhi = clusterTHnSparse.Projection(2)
      hClusPHOSEnergy = clusterTHnSparse.Projection(0)
    else:
      clusterTHnSparse.GetAxis(4).SetRange(3,3)
      hClusPHOSEta = clusterTHnSparse.Projection(2)
      hClusPHOSPhi = clusterTHnSparse.Projection(3)
      hClusPHOSEnergy = clusterTHnSparse.Projection(1)
    hClusPHOSEta.SetName("ClusEtaPHOS")
    hClusPHOSPhi.SetName("ClusPhiPHOS")
    hClusPHOSEnergy.SetName("ClusEnergyPHOS")

    # Plot phi distribution
    outputFilename = os.path.join(outputDirClusters, "hClusPHOSPhi" + fileFormat)
    plotHist(hClusPHOSPhi, outputFilename, "hist E")

    # Plot eta distribution
    outputFilename = os.path.join(outputDirClusters, "hClusPHOSEta" + fileFormat)
    plotHist(hClusPHOSEta, outputFilename, "hist E")

    # Plot energy distribution
    hClusPHOSEnergy.SetName("hClusPHOSEnergy")
    
    hClusPHOSEnergyRef = ""
    if clusterQAListRef:
      if ispp:
        clusterTHnSparseRef.GetAxis(3).SetRange(3,3)
        hClusPHOSEnergyRef = clusterTHnSparseRef.Projection(0)
      else:
        clusterTHnSparseRef.GetAxis(4).SetRange(3,3)
        hClusPHOSEnergyRef = clusterTHnSparseRef.Projection(1)
      hClusPHOSEnergyRef.SetName("clusterPHOSObservablesRef_proj_energy")

    outputFilename = os.path.join(outputDirClusters, "hClusPHOSEnergy" + fileFormat)
    xRangeMax = 100
    yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dE_{T}} [GeV^{-1}]"
    legendTitle = "PHOS Clusters"
    legendRunLabel = "Current run"
    legendRefLabel = "All runs"
    ratioYAxisTitle = "Ratio: run / all runs"
    
    plotSpectra(hClusPHOSEnergy, hClusPHOSEnergyRef, 0, 0, nEvents, nEventsRef, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, "width")

  # Plot the ratio of cluster spectra in EMCal/DCal/PHOS
  if isRun2 and includePhos:
    outputFilename = os.path.join(outputDirClusters, "hClusEnergyRatio" + fileFormat)
    xRangeMax = 250
    yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dE_{T}} [GeV^{-1}]"
    legendTitle = "Calo clusters"
    legendRunLabel = "EMCal clusters"
    legendRefLabel = "PHOS clusters"
    ratioYAxisTitle = "Ratio to PHOS"
    h2LegendLabel = "DCal clusters"
    # Note: the spectra already have been scaled by nEvents, bin width
    plotSpectra(hClusEMCalEnergy, hClusPHOSEnergy, hClusDCalEnergy, "", 1., 1., ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, "", h2LegendLabel)

  # Plot also the ratio of DCal to EMCal
  if isRun2:
    outputFilename = os.path.join(outputDirClusters, "hClusEnergyRatioEMC" + fileFormat)
    xRangeMax = 250
    if ispp:
      xRangeMax = 80
    yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dE_{T}} [GeV^{-1}]"
    legendTitle = "Calo clusters"
    legendRunLabel = "DCal clusters"
    legendRefLabel = "EMCal clusters"
    ratioYAxisTitle = "DCal / EMCal"
    # Note: the spectra already have been scaled by nEvents, bin width
    plotSpectra(hClusDCalEnergy, hClusEMCalEnergy, "", "", 1., 1., ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename)
  # Plot some PHOS QA plots
  if includePhos:
    # Plot also PHOS SM spectra
    SMlist = clusterQAList.FindObject("BySM")
    c2 = ROOT.TCanvas("c2","c2: hist",600,450)
    c2.cd()
    c2.SetLogy()

    leg = ROOT.TLegend(0.3,0.6,0.88,0.83,"PHOS SM")
    leg.SetFillColor(10)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)

    for sm in range(1,5):
      hSM = SMlist.FindObject("hPhosClusEnergy_SM" + str(sm))
      hSM.SetLineColor(sm)
      hSM.SetLineStyle(1)
      hSM.GetXaxis().SetRangeUser(0,100)
      if sm is 1:
        hSM.Draw("hist E")
      else:
        hSM.Draw("hist E same")
      leg.AddEntry(hSM, "SM " + str(sm), "l")
    leg.Draw("same")

    outputFilename = os.path.join(outputDirClusters, "hClusPHOSEnergyBySM" + fileFormat)
    c2.SaveAs(outputFilename)
    c2.Close()

  # Plot some PHOS QA plots
  if includePhos:
    hPhosNCellsVsEnergy = clusterQAList.FindObject("hPhosNCellsVsEnergy")
    outputFilename = os.path.join(outputDirClusters, "hClusPHOSNCellsVsEnergy" + fileFormat)
    plotHist(hPhosNCellsVsEnergy, outputFilename, "colz", True, True)

    hPhosM02VsEnergy = clusterQAList.FindObject("hPhosM02VsEnergy")
    outputFilename = os.path.join(outputDirClusters, "hClusPHOSM02VsEnergy" + fileFormat)
    plotHist(hPhosM02VsEnergy, outputFilename, "colz", True, True)

    hPhosCellIdVsEnergy = clusterQAList.FindObject("hPhosCellIdVsEnergy")
    outputFilename = os.path.join(outputDirClusters, "hClusPHOSCellIdVsEnergy" + fileFormat)
    plotHist(hPhosCellIdVsEnergy, outputFilename, "colz", True, True)

  # Plot EMCAL CELLS --------------------------------------------------------------

  hCellEnergy = cellQAList.FindObject("fHistCellEnergy")
  outputFilename = os.path.join(outputDirCells, "hCellEnergy" + fileFormat)
  plotHist(hCellEnergy, outputFilename, "hist E", True)

  profCellAbsIdEnergy = cellQAList.FindObject("fProfCellAbsIdEnergy")
  outputFilename = os.path.join(outputDirCells, "profCellAbsIdEnergy" + fileFormat)
  plotHist(profCellAbsIdEnergy, outputFilename)

  hCellTime = cellQAList.FindObject("fHistCellTime")
  outputFilename = os.path.join(outputDirCells, "hCellTime" + fileFormat)
  plotHist(hCellTime, outputFilename, "hist E")

  profCellAbsIdTime = cellQAList.FindObject("fProfCellAbsIdTime")
  outputFilename = os.path.join(outputDirCells, "profCellAbsIdTime" + fileFormat)
  profCellAbsIdTime.GetYaxis().SetRangeUser(-0.2e-6,0.2e-6)
  plotHist(profCellAbsIdTime, outputFilename)

  # Plot the CELL energy spectrum with and without timing cuts
  hCellEnergyTall = cellQAList.FindObject("fHistCellEvsTime")
  hCellEnergyTall = hCellEnergyTall.ProjectionY()
  hCellEnergyTall.SetName("cell_Allproj_energy")
  hCellEnergyTall.GetXaxis().SetTitle("E_{Cell} [GeV]")
  outputFilename = os.path.join(outputDirCells, "hCellEnergyTall" + fileFormat)
  plotHist(hCellEnergyTall, outputFilename, "hist E", True)

  hCellEnergyTsel = cellQAList.FindObject("fHistCellEvsTime")
  hCellEnergyTsel.GetXaxis().SetRangeUser(-50e-9,50e-9) #recomended time cut
  hCellEnergyTsel = hCellEnergyTsel.ProjectionY()
  hCellEnergyTsel.SetName("cell_Selproj_energy")
  hCellEnergyTsel.GetXaxis().SetTitle("E_{Cell} |t_{cell}|<50ns [GeV]")
  outputFilename = os.path.join(outputDirCells, "hCellEnergyTsel" + fileFormat)
  plotHist(hCellEnergyTsel, outputFilename, "hist E", True)

  #refernce histograms
  if cellQAListRef:
    hCellEnergyTallRef = cellQAListRef.FindObject("fHistCellEvsTime")
    hCellEnergyTallRef = hCellEnergyTallRef.ProjectionY()
    hCellEnergyTallRef.SetName("cellRef_Allproj_energy")

    hCellEnergyTselRef = cellQAListRef.FindObject("fHistCellEvsTime")
    hCellEnergyTselRef.GetXaxis().SetRangeUser(-50e-9,50e-9)
    hCellEnergyTselRef = hCellEnergyTselRef.ProjectionY()
    hCellEnergyTselRef.SetName("cellRef_Selproj_energy")

    xRangeMax = 100
    if ispp:
      xRangeMax = 80
    yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dE_{Cell}} [GeV^{-1}]"
    legendTitle = "EMCal Cells"
    legendRunLabel = "Current run"
    legendRefLabel = "All runs"
    ratioYAxisTitle = "Ratio: run / all runs"
    outputFilename = os.path.join(outputDirCells, "hCellEnergyTallRatio" + fileFormat)
    plotSpectra(hCellEnergyTall, hCellEnergyTallRef, 0, 0, nEvents, nEventsRef, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, "width")

    outputFilename = os.path.join(outputDirCells, "hCellEnergyTselRatio" + fileFormat)
    plotSpectra(hCellEnergyTsel, hCellEnergyTselRef, 0, 0, nEvents, nEventsRef, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, "width")


########################################################################################################
# Plot charged jet histograms    #######################################################################
########################################################################################################
def plotChargedJetQA(ispp, isPtHard, chargedJetList, outputDir, chargedJetListRef, nEvents, nEventsRef, fileFormat):
  
  # Create subdirectory for Jets
  outputDirJets = outputDir + "Jets/"
  if not os.path.exists(outputDirJets):
    os.makedirs(outputDirJets)
  
  chargedJetTHnSparse = chargedJetList.FindObject("fHistJetObservables")
  # (Centrality, eta, phi, pT, pTcorr, pT leading particle)

  if chargedJetListRef:
    chargedJetTHnSparseRef = chargedJetListRef.FindObject("fHistJetObservables")

  ROOT.gStyle.SetOptTitle(1)

  if not ispp:
    # Plot charged jet rho vs. centrality
    hChargedJetRhoVsCent = chargedJetList.FindObject("fHistRhoVsCent")
    hChargedJetRhoVsCent.SetTitle("Rho vs. Centrality, Charged Jets")
    outputFilename = os.path.join(outputDirJets, "hChargedJetRhoVsCent" + fileFormat)
    plotHist(hChargedJetRhoVsCent, outputFilename, "colz", False, True)

  # Plot charged jet eta-phi, for jet pT > threshold

  # there are ceil(250/3)=84 jet pt bins
  # (5,84) means (~12 GeV < jet pT < 250 GeV)
  # (11,84) means (~30 GeV < jet pT < 250 GeV)
  minJetPtBin = 5
  maxJetPtBin = 84
  if ispp:
    chargedJetTHnSparse.GetAxis(2).SetRange(minJetPtBin, maxJetPtBin)
  else:
    minJetPtBin = 11
    chargedJetTHnSparse.GetAxis(3).SetRange(minJetPtBin, maxJetPtBin)

  if ispp:
    hChargedJetEtaPhi = chargedJetTHnSparse.Projection(1,0)
  else:
    hChargedJetEtaPhi = chargedJetTHnSparse.Projection(2,1)
  hChargedJetEtaPhi.SetName("ChargedJetEtaPhi")
  hChargedJetEtaPhi.SetTitle("Charged Jet Occupancy, p_{T,jet} > " + str((minJetPtBin-1)*3) + " GeV")
  hChargedJetEtaPhi.GetXaxis().SetRangeUser(-0.8,0.8)
  outputFilename = os.path.join(outputDirJets, "hChargedJetEtaPhi" + fileFormat)
  plotHist(hChargedJetEtaPhi, outputFilename, "colz", False)

  # Plot ratio to reference run, if supplied
  if chargedJetListRef:
    if ispp:
      chargedJetTHnSparseRef.GetAxis(2).SetRange(minJetPtBin, maxJetPtBin)
      hChargedJetEtaPhiRef = chargedJetTHnSparseRef.Projection(1,0)
    else:
      chargedJetTHnSparseRef.GetAxis(3).SetRange(minJetPtBin, maxJetPtBin)
      hChargedJetEtaPhiRef = chargedJetTHnSparseRef.Projection(2,1)
    hChargedJetEtaPhiRef.SetName("ChargedJetEtaPhiRef")
    hChargedJetEtaPhi.Scale(1./nEvents)
    hChargedJetEtaPhiRef.Scale(1./nEventsRef)
    hChargedJetEtaPhiRatio = hChargedJetEtaPhi.Clone()
    hChargedJetEtaPhiRatio.Divide(hChargedJetEtaPhiRef)
    hChargedJetEtaPhiRatio.SetTitle("Charged Jet p_{T,jet} > " + str((minJetPtBin-1)*3) + " GeV Occupancy (per event): Current Run / All Runs")
    outputFilename = os.path.join(outputDirJets, "hChargedJetEtaPhiRatio" + fileFormat)
    plotHist(hChargedJetEtaPhiRatio, outputFilename, "colz", False, True)
    if ispp:
      chargedJetTHnSparseRef.GetAxis(2).SetRange(1, maxJetPtBin)
    else:
      chargedJetTHnSparseRef.GetAxis(3).SetRange(1, maxJetPtBin)
  if ispp:
    chargedJetTHnSparse.GetAxis(2).SetRange(1, maxJetPtBin)
  else:
    chargedJetTHnSparse.GetAxis(3).SetRange(1, maxJetPtBin)

  # Plot charged jet pT
  if ispp:
    hChargedJetPt = chargedJetTHnSparse.Projection(2)
  else:
    hChargedJetPt = chargedJetTHnSparse.Projection(3)
  hChargedJetPt.SetName("hChargedJetPt")
  
  hChargedJetPtRef = ""
  if chargedJetListRef:
    if ispp:
      hChargedJetPtRef = chargedJetTHnSparseRef.Projection(2)
    else:
      hChargedJetPtRef = chargedJetTHnSparseRef.Projection(3)
    hChargedJetPtRef.SetName("hChargedJetPt")

  outputFilename = os.path.join(outputDirJets, "hChargedJetPt" + fileFormat)
  xRangeMax = 250
  yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [GeV^{-1}]"
  legendTitle = "Charged jets"
  legendRunLabel = "Current run"
  legendRefLabel = "All runs"
  ratioYAxisTitle = "Ratio: run / all runs"
  
  plotSpectra(hChargedJetPt, hChargedJetPtRef, 0, 0, nEvents, nEventsRef, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename)

  # Plot charged jet pT leading particle vs. jet pT
  if ispp:
    hChargedJetPtLeadjetPt = chargedJetTHnSparse.Projection(3,2)
  else:
    hChargedJetPtLeadjetPt = chargedJetTHnSparse.Projection(5,3)
  hChargedJetPtLeadjetPt.SetName("fHistChJetObservables_proj_pt_leadpt")
  hChargedJetPtLeadjetPt.SetTitle("Leading pT vs. Jet pT, Charged Jets")
  outputFilename = os.path.join(outputDirJets, "hChargedJetPtLeadjetPt" + fileFormat)
  if isPtHard:
    yMin= hChargedJetPt.GetBinContent(hChargedJetPt.FindBin(200)) #find entry in bin at 200 GeV to get the right y-Axis scale
    yMax= hChargedJetPt.GetBinContent(hChargedJetPt.GetMaximumBin()) #find entry in bin at 200 GeV to get the right y-Axis scale
    hChargedJetPt.GetYaxis().SetRangeUser(yMin,yMax*1.1)
    plotHist(hChargedJetPtLeadjetPt, outputFilename, "colz", "", True)
  else:
    plotHist(hChargedJetPtLeadjetPt, outputFilename, "colz", "", True)
  ROOT.gStyle.SetOptTitle(0)

  # Plot charged jet pT, background-subtracted
  if not ispp:
    hChargedJetPtCorr = chargedJetTHnSparse.Projection(4)
    hChargedJetPtCorr.SetName("hChargedJetPtCorr")
  
    hChargedJetPtCorrRef = ""
    if chargedJetListRef:
      hChargedJetPtCorrRef = chargedJetTHnSparseRef.Projection(4)
      hChargedJetPtCorrRef.SetName("hChargedJetPtCorr")

    outputFilename = os.path.join(outputDirJets, "hChargedJetPtCorr" + fileFormat)
    xRangeMax = 150
    yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [GeV^{-1}]"
    legendTitle = "Charged jets, background subtracted"
    legendRunLabel = "Current run"
    legendRefLabel = "All runs"
    ratioYAxisTitle = "Ratio: run / all runs"
    
    plotSpectra(hChargedJetPtCorr, hChargedJetPtCorrRef, 0, 0, nEvents, nEventsRef, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename)

  # Plot charged jet pT, background-subtracted, by centrality
    chargedJetTHnSparse.GetAxis(0).SetRange(1, 1)
    hChargedJetPtCorr010 = chargedJetTHnSparse.Projection(4)
    hChargedJetPtCorr010.SetName("hChargedJetPtCorr010")

    chargedJetTHnSparse.GetAxis(0).SetRange(2, 2)
    hChargedJetPtCorr1030 = chargedJetTHnSparse.Projection(4)
    hChargedJetPtCorr1030.SetName("hChargedJetPtCorr1030")

    chargedJetTHnSparse.GetAxis(0).SetRange(3, 3)
    hChargedJetPtCorr3050 = chargedJetTHnSparse.Projection(4)
    hChargedJetPtCorr3050.SetName("hChargedJetPtCorr3050")

    chargedJetTHnSparse.GetAxis(0).SetRange(4, 4)
    hChargedJetPtCorr5090 = chargedJetTHnSparse.Projection(4)
    hChargedJetPtCorr5090.SetName("hChargedJetPtCorr5090")

    outputFilename = os.path.join(outputDirJets, "hChargedJetPtCorrCentral" + fileFormat)
    xRangeMax = 150
    yAxisTitle = "#frac{1}{N_{evts}N_{coll}}#frac{dN}{dp_{T}} [GeV^{-1}]"
    legendTitle = "Charged jets, background subtracted"
    legendRunLabel = "0-10%"
    legendRefLabel = "50-90%"
    DCalLegendLabel = "10-30%"
    PHOSLegendLabel = "30-50%"
    ratioYAxisTitle = "R_{CP}"

    # Scale by Ncoll, to compare different centralities
    # Values taken from https://twiki.cern.ch/twiki/bin/view/ALICE/CentralityCodeSnippets
    Ncoll010 = 1636.
    Ncoll1030 = 801.
    Ncoll3050 = 264.
    Ncoll5090 = 38.1
    Ncoll090 = 435.3

    hChargedJetPtCorr010.Scale(4.) # Scale by number of events in 0-10% relative to 50-90%
    hChargedJetPtCorr1030.Scale(Ncoll010/Ncoll1030 * 2.)
    hChargedJetPtCorr3050.Scale(Ncoll010/Ncoll3050 * 2.)
    hChargedJetPtCorr5090.Scale(Ncoll010/Ncoll5090)

    plotSpectra(hChargedJetPtCorr010, hChargedJetPtCorr5090, 0, 0, nEvents, nEvents, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, "", DCalLegendLabel, PHOSLegendLabel)
    chargedJetTHnSparse.GetAxis(0).SetRange(1,4)

########################################################################################################
# Plot full jet histograms     ##############################################################################
########################################################################################################
def plotFullJetQA(ispp, isPtHard, isRun2, includePhos, fullJetList, outputDir, fullJetListRef, nEvents, nEventsRef, fileFormat):
  
  # Create subdirectory for Jets
  outputDirJets = outputDir + "Jets/"
  if not os.path.exists(outputDirJets):
    os.makedirs(outputDirJets)

  fullJetTHnSparse = fullJetList.FindObject("fHistJetObservables")
  # (Centrality, eta, phi, pT, pTcorr, pT leading particle)

  if fullJetListRef:
    fullJetTHnSparseRef = fullJetListRef.FindObject("fHistJetObservables")
  
  ROOT.gStyle.SetOptTitle(1)

  if not ispp:
    # Plot full jet rho vs. centrality
    hFullJetRhoVsCent = fullJetList.FindObject("fHistRhoVsCent")
    hFullJetRhoVsCent.SetTitle("Rho vs. Centrality, Full Jets")
    outputFilename = os.path.join(outputDirJets, "hFullJetRhoVsCent" + fileFormat)
    plotHist(hFullJetRhoVsCent, outputFilename, "colz", False, True)

  # Plot Neutral Energy Fraction
  hFullJetNEF = fullJetList.FindObject("hNEFVsPtEMC")
  if not ispp:
    if hFullJetNEF:
      hFullJetNEF = hNEFVsPtCorrVsCentEMCal.Project3D("zy")
    else:
      print("hFullJetNEF not saved for PbPb in this version")
  hFullJetNEF = hFullJetNEF.ProjectionY()
  hFullJetNEFDCal = fullJetList.FindObject("hNEFVsPtDCal")
  hFullJetNEFDCal = hFullJetNEFDCal.ProjectionY()
  hFullJetNEF.SetTitle("NEF vs. p_{T,jet}, Full Jets")
  outputFilename = os.path.join(outputDirJets, "hFullJetNEF" + fileFormat)
  # plotHist(hFullJetNEF, outputFilename, "colz", True, False)
  plotNEFSpectra(hFullJetNEF,hFullJetNEFDCal, 0,nEvents, ispp, 1, "1/N_{Evt} dN/dNEF", "EMCal", outputFilename,"", "DCal")
  

  if ispp:
    # Plot Delta HadCorr vs pT
    hFullJetDeltaHcorr = fullJetList.FindObject("hDeltaEHadCorr")
    hFullJetDeltaHcorr.GetXaxis().SetRangeUser(0., 150.)
    hFullJetDeltaHcorr.SetTitle("#Delta E vs. p_{T,jet}, Full Jets")
    #outputFilename = os.path.join(outputDirJets, "hFullJetDeltaHcorr" + fileFormat)
    #plotHist(hFullJetDeltaHcorr, outputFilename, "colz", False, True)
    hFullJetDeltaHcorr.SetTitle("<#DeltaE> vs. p_{T,jet}, Full Jets")
    hDeltaEHadCorrProf = hFullJetDeltaHcorr.ProfileX()
    hDeltaEHadCorrProf.GetYaxis().SetRangeUser(0.08, 15.)
    hDeltaEHadCorrProf.SetLineColor(1)
    hDeltaEHadCorrProf.SetMarkerStyle(20)
    hDeltaEHadCorrProf.GetYaxis().SetTitleOffset(1.2)
    hDeltaEHadCorrProf.GetYaxis().SetTitle("< #sum#it{E}_{nonlincorr} - #it{E}_{hadcorr} >")
    outputFilename = os.path.join(outputDirJets, "hDeltaEHadCorrProf" + fileFormat)
    plotHist(hDeltaEHadCorrProf, outputFilename, "E", True, False)
  else:
    print("hFullJetDeltaHcorr not saved for PbPb yet") #need to project the TH3 down to 2D
  

  # Plot full jet eta-phi, for jet pT > threshold
  # there are ceil(250/3)=84 jet pt bins
  # (5,84) means (~12 GeV < jet pT < 250 GeV)
  # (11,84) means (~30 GeV < jet pT < 250 GeV)
  minJetPtBin = 5
  maxJetPtBin = 84
  if ispp:
    fullJetTHnSparse.GetAxis(2).SetRange(minJetPtBin, maxJetPtBin)
  else:
    minJetPtBin = 11
    fullJetTHnSparse.GetAxis(3).SetRange(minJetPtBin, maxJetPtBin)

  # Plot full jet eta-phi
  if ispp:
    hFullJetEtaPhi = fullJetTHnSparse.Projection(1,0)
  else:
    hFullJetEtaPhi = fullJetTHnSparse.Projection(2,1)
  hFullJetEtaPhi.SetName("FullJetEtaPhi")
  hFullJetEtaPhi.SetTitle("Full Jet Occupancy, p_{T,jet} > " + str((minJetPtBin-1)*3) + " GeV")
  outputFilename = os.path.join(outputDirJets, "hFullJetEtaPhi" + fileFormat)
  hFullJetEtaPhi.GetXaxis().SetRangeUser(-0.8,0.8)
  hFullJetEtaPhi.GetYaxis().SetRangeUser(1.2,5.8)
  plotHist(hFullJetEtaPhi, outputFilename, "colz", False)

  # Plot ratio to reference run, if supplied
  if fullJetListRef:
    if ispp:
      fullJetTHnSparseRef.GetAxis(2).SetRange(minJetPtBin, maxJetPtBin)
      hFullJetEtaPhiRef = fullJetTHnSparseRef.Projection(1,0)
    else:
      fullJetTHnSparseRef.GetAxis(3).SetRange(minJetPtBin, maxJetPtBin)
      hFullJetEtaPhiRef = fullJetTHnSparseRef.Projection(2,1)
    hFullJetEtaPhiRef.SetName("FullJetEtaPhiRef")
    hFullJetEtaPhi.Scale(1./nEvents)
    hFullJetEtaPhiRef.Scale(1./nEventsRef)
    hFullJetEtaPhiRatio = hFullJetEtaPhi.Clone()
    hFullJetEtaPhiRatio.Divide(hFullJetEtaPhiRef)
    hFullJetEtaPhiRatio.SetTitle("Full Jet p_{T,jet} > " + str((minJetPtBin-1)*3) + " GeV Occupancy (per event): Current Run / All Runs")
    outputFilename = os.path.join(outputDirJets, "hFullJetEtaPhiRatio" + fileFormat)
    plotHist(hFullJetEtaPhiRatio, outputFilename, "colz", False)
    if ispp:
      fullJetTHnSparseRef.GetAxis(2).SetRange(1, maxJetPtBin)
    else:
      fullJetTHnSparseRef.GetAxis(3).SetRange(1, maxJetPtBin)
  if ispp:
    fullJetTHnSparse.GetAxis(2).SetRange(1, maxJetPtBin)
  else:
    fullJetTHnSparse.GetAxis(3).SetRange(1, maxJetPtBin)


  ROOT.gStyle.SetOptTitle(0)


  # Plot full jet pT
  if ispp:
    hFullJetPt = fullJetTHnSparse.Projection(2)
  else:
    hFullJetPt = fullJetTHnSparse.Projection(3)
  hFullJetPt.SetName("hFullJetPt")
  
  hFullJetPtRef = ""
  if fullJetListRef:
    if ispp:
      hFullJetPtRef = fullJetTHnSparseRef.Projection(2)
    else:
      hFullJetPtRef = fullJetTHnSparseRef.Projection(3)
    hFullJetPtRef.SetName("hFullJetPt")

  outputFilename = os.path.join(outputDirJets, "hFullJetPt" + fileFormat)
  xRangeMax = 250
  yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [GeV^{-1}]"
  legendTitle = "Full jets"
  legendRunLabel = "Current run"
  legendRefLabel = "All runs"
  ratioYAxisTitle = "Ratio: run / all runs"
  
  plotSpectra(hFullJetPt, hFullJetPtRef, 0, 0, nEvents, nEventsRef, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename)

  # Plot full jet pT leading particle vs. jet pT
  if ispp:
    hFullJetPtLeadjetPt = fullJetTHnSparse.Projection(3,2)
  else:
    hFullJetPtLeadjetPt = fullJetTHnSparse.Projection(5,3)
  hFullJetPtLeadjetPt.SetName("fHistFuJetObservables_proj_pt_leadpt")
  hFullJetPtLeadjetPt.SetTitle("Leading pT vs. Jet pT, Full Jets")
  outputFilename = os.path.join(outputDirJets, "hFullJetPtLeadjetPt" + fileFormat)
  if ispp:
    hFullJetPtLeadjetPt.GetXaxis().SetRangeUser(0,200)
    hFullJetPtLeadjetPt.GetYaxis().SetRangeUser(0,100)
  if isPtHard:
    yMin  = hFullJetPt.GetBinContent(hFullJetPt.FindBin(200)) #find entry in bin at 200 GeV to get the right y-Axis scale
    maxBin= hFullJetPt.GetBinContent(hFullJetPt.GetMaximumBin())
    hFullJetPt.SetMinimum(yMin);
    hFullJetPt.SetMaximum(maxBin*1.1);
    plotHist(hFullJetPtLeadjetPt, outputFilename, "colz", "", True)
  else:
    plotHist(hFullJetPtLeadjetPt, outputFilename, "colz", "", True)

  # Plot full jet pT, background subtracted
  if not ispp:
    hFullJetPtCorr = fullJetTHnSparse.Projection(4)
    hFullJetPtCorr.SetName("hFullJetPtCorr")
    
    hFullJetPtCorrRef = ""
    if fullJetListRef:
      hFullJetPtCorrRef = fullJetTHnSparseRef.Projection(4)
      hFullJetPtCorrRef.SetName("hFullJetPtCorrRef")

    outputFilename = os.path.join(outputDirJets, "hFullJetPtCorr" + fileFormat)
    xRangeMax = 150
    yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [GeV^{-1}]"
    legendTitle = "Full jets, background subtracted"
    legendRunLabel = "Current run"
    legendRefLabel = "All runs"
    ratioYAxisTitle = "Ratio: run / all runs"
    
    plotSpectra(hFullJetPtCorr, hFullJetPtCorrRef, 0, 0, nEvents, nEventsRef, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename)
  
  # Plot full jet pT, background-subtracted, by centrality
    fullJetTHnSparse.GetAxis(0).SetRange(1, 1)
    hFullJetPtCorr010 = fullJetTHnSparse.Projection(4)
    hFullJetPtCorr010.SetName("hFullJetPtCorr010")
    
    fullJetTHnSparse.GetAxis(0).SetRange(2, 2)
    hFullJetPtCorr1030 = fullJetTHnSparse.Projection(4)
    hFullJetPtCorr1030.SetName("hFullJetPtCorr1030")
    
    fullJetTHnSparse.GetAxis(0).SetRange(3, 3)
    hFullJetPtCorr3050 = fullJetTHnSparse.Projection(4)
    hFullJetPtCorr3050.SetName("hFullJetPtCorr3050")
    
    fullJetTHnSparse.GetAxis(0).SetRange(4, 4)
    hFullJetPtCorr5090 = fullJetTHnSparse.Projection(4)
    hFullJetPtCorr5090.SetName("hFullJetPtCorr5090")
    
    outputFilename = os.path.join(outputDirJets, "hFullJetPtCorrCentral" + fileFormat)
    xRangeMax = 150
    yAxisTitle = "#propto#frac{1}{N_{evts}N_{coll}}#frac{dN}{dp_{T}} [GeV^{-1}]"
    legendTitle = "Full jets, background subtracted"
    legendRunLabel = "0-10%"
    legendRefLabel = "50-90%"
    DCalLegendLabel = "10-30%"
    PHOSLegendLabel = "30-50%"
    ratioYAxisTitle = "R_{CP}"
    
    # Scale by Ncoll, to compare different centralities
    # Values taken from https://twiki.cern.ch/twiki/bin/view/ALICE/CentralityCodeSnippets
    Ncoll010 = 1636.
    Ncoll1030 = 801.
    Ncoll3050 = 264.
    Ncoll5090 = 38.1
    Ncoll090 = 435.3
    
    hFullJetPtCorr010.Scale(4.) # Scale by number of events in 0-10% relative to 50-90%
    hFullJetPtCorr1030.Scale(Ncoll010/Ncoll1030 * 2.)
    hFullJetPtCorr3050.Scale(Ncoll010/Ncoll3050 * 2.)
    hFullJetPtCorr5090.Scale(Ncoll010/Ncoll5090)

    plotSpectra(hFullJetPtCorr010, hFullJetPtCorr5090, 0, 0, nEvents, nEvents, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, "", DCalLegendLabel, PHOSLegendLabel)
    fullJetTHnSparse.GetAxis(0).SetRange(1,4)

  # Plot full jet pT spectra separately for EMCal, DCal, PHOS jets
  if isRun2 and includePhos:
    #EMCal jets -- divide from DCal/PHOS by phi cut
    if ispp:
      phiDivideBin = fullJetTHnSparse.GetAxis(1).FindBin(4.)
      fullJetTHnSparse.GetAxis(1).SetRange(0, phiDivideBin)
      hFullJetEMCalEtaPhiPt = fullJetTHnSparse.Projection(0,1,2, "o") # "o" keep the original axis range
    else:
      phiDivideBin = fullJetTHnSparse.GetAxis(2).FindBin(4.)
      fullJetTHnSparse.GetAxis(2).SetRange(0, phiDivideBin)
      hFullJetEMCalEtaPhiPt = fullJetTHnSparse.Projection(1,2,3, "o")
      hFullJetEMCalEtaPhiPtCorr = fullJetTHnSparse.Projection(1,2,4, "o")
      hFullJetEMCalEtaPhiPtCorr.SetName("FullJetEMCalEtaPhiPtCorr");
    hFullJetEMCalEtaPhiPt.SetName("FullJetEMCalEtaPhiPt");

    hFullJetEMCalEtaPhi = hFullJetEMCalEtaPhiPt.Project3D("yx")
    outputFilename = os.path.join(outputDirJets, "hFullJetEtaPhiEMCal" + fileFormat)
    plotHist(hFullJetEMCalEtaPhi, outputFilename, "colz")

    hFullJetEMCalPt = hFullJetEMCalEtaPhiPt.Project3D("z")
    if not ispp:
      hFullJetEMCalPtCorr = hFullJetEMCalEtaPhiPtCorr.Project3D("z")

    # DCal jets -- divide from EMCal by phi cut, and divide from PHOS by |eta| > 0.22 (no fiducial cut on inner eta)
    if ispp:
      etaMinDCalBinNeg = fullJetTHnSparse.GetAxis(0).FindBin(-0.22)
      etaMinDCalBinPos = fullJetTHnSparse.GetAxis(0).FindBin(0.22)
      fullJetTHnSparse.GetAxis(1).SetRange(phiDivideBin, 101)
      fullJetTHnSparse.GetAxis(0).SetRange(1, etaMinDCalBinNeg)
      hFullJetDCalEtaPhiPtNeg = fullJetTHnSparse.Projection(0,1,2, "o")
    else:
      etaMinDCalBinNeg = fullJetTHnSparse.GetAxis(1).FindBin(-0.22)
      etaMinDCalBinPos = fullJetTHnSparse.GetAxis(1).FindBin(0.22)
      fullJetTHnSparse.GetAxis(2).SetRange(phiDivideBin, 101)
      fullJetTHnSparse.GetAxis(1).SetRange(1, etaMinDCalBinNeg)
      hFullJetDCalEtaPhiPtNeg = fullJetTHnSparse.Projection(1,2,3, "o")
      hFullJetDCalEtaPhiPtCorrNeg = fullJetTHnSparse.Projection(1,2,4, "o")
      hFullJetDCalEtaPhiPtCorrNeg.SetName("FullJetDCalEtaPhiPtCorrNeg");
    hFullJetDCalEtaPhiPtNeg.SetName("FullJetDCalEtaPhiPtNeg");

    if ispp:
      fullJetTHnSparse.GetAxis(1).SetRange(phiDivideBin, 101)
      fullJetTHnSparse.GetAxis(0).SetRange(etaMinDCalBinPos, 70)
      hFullJetDCalEtaPhiPtPos = fullJetTHnSparse.Projection(0,1,2, "o")
    else:
      fullJetTHnSparse.GetAxis(2).SetRange(phiDivideBin, 101)
      fullJetTHnSparse.GetAxis(1).SetRange(etaMinDCalBinPos, 70)
      hFullJetDCalEtaPhiPtPos = fullJetTHnSparse.Projection(1,2,3, "o")
      hFullJetDCalEtaPhiPtCorrPos = fullJetTHnSparse.Projection(1,2,4, "o")
      hFullJetDCalEtaPhiPtCorrPos.SetName("FullJetDCalEtaPhiPtCorrPos");
    hFullJetDCalEtaPhiPtPos.SetName("FullJetDCalEtaPhiPtPos");

    # Add the TH3s
    hFullJetDCalEtaPhiPt = hFullJetDCalEtaPhiPtNeg.Clone()
    hFullJetDCalEtaPhiPt.Add(hFullJetDCalEtaPhiPtPos)
    if not ispp:
      hFullJetDCalEtaPhiPtCorr = hFullJetDCalEtaPhiPtCorrNeg.Clone()
      hFullJetDCalEtaPhiPtCorr.Add(hFullJetDCalEtaPhiPtCorrPos)

    # Project to TH2 for eta-phi, and TH1 of pT
    hFullJetDCalEtaPhi = hFullJetDCalEtaPhiPt.Project3D("yx")
    outputFilename = os.path.join(outputDirJets, "hFullJetEtaPhiDCal" + fileFormat)
    plotHist(hFullJetDCalEtaPhi, outputFilename, "colz")

    hFullJetDCalPt = hFullJetDCalEtaPhiPt.Project3D("z")
    if not ispp:
      hFullJetDCalPtCorr = hFullJetDCalEtaPhiPtCorr.Project3D("z")

    # Gap jets -- divide from EMCal by phi cut, and divide from PHOS by |eta| > 0.13 and DCal by |eta| < 0.22
    if ispp:
      etaMinPHOSBin = fullJetTHnSparse.GetAxis(0).FindBin(-0.13)
      etaMaxPHOSBin = fullJetTHnSparse.GetAxis(0).FindBin(0.13)
      fullJetTHnSparse.GetAxis(1).SetRange(phiDivideBin, 101)
      fullJetTHnSparse.GetAxis(0).SetRange(etaMinDCalBinNeg, etaMinPHOSBin)
      hFullJetGapEtaPhiPtNeg = fullJetTHnSparse.Projection(0,1,2, "o")
    else:
      etaMinPHOSBin = fullJetTHnSparse.GetAxis(1).FindBin(-0.13)
      etaMaxPHOSBin = fullJetTHnSparse.GetAxis(1).FindBin(0.13)
      fullJetTHnSparse.GetAxis(2).SetRange(phiDivideBin, 101)
      fullJetTHnSparse.GetAxis(1).SetRange(etaMinDCalBinNeg, etaMinPHOSBin)
      hFullJetGapEtaPhiPtNeg = fullJetTHnSparse.Projection(1,2,3, "o")
      hFullJetGapEtaPhiPtCorrNeg = fullJetTHnSparse.Projection(1,2,4, "o")
      hFullJetGapEtaPhiPtCorrNeg.SetName("FullJetGapEtaPhiPtCorrNeg");
    hFullJetGapEtaPhiPtNeg.SetName("FullJetGapEtaPhiPtNeg");

    if ispp:
      fullJetTHnSparse.GetAxis(1).SetRange(phiDivideBin, 101)
      fullJetTHnSparse.GetAxis(0).SetRange(etaMaxPHOSBin, etaMinDCalBinPos)
      hFullJetGapEtaPhiPtPos = fullJetTHnSparse.Projection(0,1,2, "o")
    else:
      fullJetTHnSparse.GetAxis(2).SetRange(phiDivideBin, 101)
      fullJetTHnSparse.GetAxis(1).SetRange(etaMaxPHOSBin, etaMinDCalBinPos)
      hFullJetGapEtaPhiPtPos = fullJetTHnSparse.Projection(1,2,3, "o")
      hFullJetGapEtaPhiPtCorrPos = fullJetTHnSparse.Projection(1,2,4, "o")
      hFullJetGapEtaPhiPtCorrPos.SetName("FullJetGapEtaPhiPtCorrPos");
    hFullJetGapEtaPhiPtPos.SetName("FullJetGapEtaPhiPtPos");

    
    # Add the TH3s
    hFullJetGapEtaPhiPt = hFullJetGapEtaPhiPtNeg.Clone()
    hFullJetGapEtaPhiPt.Add(hFullJetGapEtaPhiPtPos)
    if not ispp:
      hFullJetGapEtaPhiPtCorr = hFullJetGapEtaPhiPtCorrNeg.Clone()
      hFullJetGapEtaPhiPtCorr.Add(hFullJetGapEtaPhiPtCorrPos)
    
    # Project to TH2 for eta-phi, and TH1 of pT
    hFullJetGapEtaPhi = hFullJetGapEtaPhiPt.Project3D("yx")
    outputFilename = os.path.join(outputDirJets, "hFullJetEtaPhiGap" + fileFormat)
    plotHist(hFullJetGapEtaPhi, outputFilename, "colz")
    
    hFullJetGapPt = hFullJetGapEtaPhiPt.Project3D("z")
    if not ispp:
      hFullJetGapPtCorr = hFullJetGapEtaPhiPtCorr.Project3D("z")

    # PHOS jets -- divide from EMCal by phi cut, and divide from DCal by eta < 0.13 (no fiducial cut on inner eta)
    #              fiducial cut on DCal (kDCALfid) ensures that remaining region is only PHOS
    if ispp:
      fullJetTHnSparse.GetAxis(1).SetRange(phiDivideBin, 101)
      fullJetTHnSparse.GetAxis(0).SetRange(etaMinPHOSBin, etaMaxPHOSBin)
      hFullJetPHOSEtaPhiPt = fullJetTHnSparse.Projection(0,1,2, "o")
    else:
      fullJetTHnSparse.GetAxis(2).SetRange(phiDivideBin, 101)
      fullJetTHnSparse.GetAxis(1).SetRange(etaMinPHOSBin, etaMaxPHOSBin)
      hFullJetPHOSEtaPhiPt = fullJetTHnSparse.Projection(1,2,3, "o")
      hFullJetPHOSEtaPhiPtCorr = fullJetTHnSparse.Projection(1,2,4, "o")
      hFullJetPHOSEtaPhiPtCorr.SetName("FullJetPHOSEtaPhiPtCorr");
    hFullJetPHOSEtaPhiPt.SetName("FullJetPHOSEtaPhiPt");
    hFullJetPHOSEtaPhi = hFullJetPHOSEtaPhiPt.Project3D("yx")
    outputFilename = os.path.join(outputDirJets, "hFullJetEtaPhiPHOS" + fileFormat)
    plotHist(hFullJetPHOSEtaPhi, outputFilename, "colz")

    hFullJetPtRef = ""
    if fullJetListRef:
      if ispp:
        hFullJetPtRef = fullJetTHnSparseRef.Projection(2)
      else:
        hFullJetPtRef = fullJetTHnSparseRef.Projection(3)
        hFullJetPtCorrRef = fullJetTHnSparseRef.Projection(4)
        hFullJetPtCorrRef.SetName("hFullJetPtCorr")
      hFullJetPtRef.SetName("hFullJetPt")

    hFullJetPHOSPt = hFullJetPHOSEtaPhiPt.Project3D("z")
    if not ispp:
      hFullJetPHOSPtCorr = hFullJetPHOSEtaPhiPtCorr.Project3D("z")

    # Now plot the EMCal/DCal/PHOS jet pT spectra and their ratio to the reference
    outputFilename = os.path.join(outputDirJets, "hFullJetPtCalo" + fileFormat)
    xRangeMax = 250
    yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [GeV^{-1}]"
    legendTitle = "Full jets"
    legendRunLabel = "EMCal jets"
    if fullJetListRef:
      legendRefLabel = "All full jets"
      ratioYAxisTitle = "Ratio to all"
      DCalLegendLabel = "DCal jets"
      PHOSLegendLabel = "PHOS jets"
      plotSpectra(hFullJetEMCalPt, hFullJetPtRef, hFullJetDCalPt, hFullJetPHOSPt, nEvents, nEventsRef, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, "", DCalLegendLabel, PHOSLegendLabel)
    else:
      legendRefLabel = "PHOS jets"
      ratioYAxisTitle = "Ratio to PHOS"
      h2LegendLabel = "DCal jets"
      plotSpectra(hFullJetEMCalPt, hFullJetPHOSPt, hFullJetDCalPt, "", nEvents, nEvents, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, "", h2LegendLabel)

    # And plot the background subtracted EMCal/DCal/PHOS jet pT spectra and their ratio to the reference
    if not ispp:
      outputFilename = os.path.join(outputDirJets, "hFullJetPtCorrCalo" + fileFormat)
      xRangeMax = 250
      yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [GeV^{-1}]"
      legendTitle = "Full jets, background subtracted"
      legendRunLabel = "EMCal jets"
      if fullJetListRef:
        legendRefLabel = "All full jets"
        ratioYAxisTitle = "Ratio to all"
        DCalLegendLabel = "DCal jets"
        PHOSLegendLabel = "PHOS jets"
        plotSpectra(hFullJetEMCalPtCorr, hFullJetPtCorrRef, hFullJetDCalPtCorr, hFullJetPHOSPtCorr, nEvents, nEventsRef, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, "", DCalLegendLabel, PHOSLegendLabel)
      else:
        legendRefLabel = "PHOS jets"
        ratioYAxisTitle = "Ratio to PHOS"
        h2LegendLabel = "DCal jets"
        plotSpectra(hFullJetEMCalPtCorr, hFullJetPHOSPtCorr, hFullJetDCalPtCorr, "", nEvents, nEvents, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, "", h2LegendLabel)

########################################################################################################
# Plot event histograms   ##############################################################################
########################################################################################################
def plotEventQA(ispp, isRun2, includePhos, qaList, outputDir, fileFormat):
  
  histNEvent = qaList.FindObject("fHistEventCount")
  nEvents = histNEvent.GetBinContent(1)
  #print("N events: %d" % nEvents)
  
  # Create subdirectory for EventQA
  outputDirEventQA = outputDir + "EventQA/"
  if not os.path.exists(outputDirEventQA):
    os.makedirs(outputDirEventQA)
  
  eventQATHnSparse = qaList.FindObject("eventQA")
  # (Centrality, N tracks, pT leading track, N clusters, leading cluster E)

  if ispp:
    # N tracks
    hEventNtracks = eventQATHnSparse.Projection(0)
    outputFilename = os.path.join(outputDirEventQA, "hEventNtracks" + fileFormat)
    plotHist(hEventNtracks, outputFilename, "hist E")
    
    # N clusters
    hEventNclusters = eventQATHnSparse.Projection(2)
    outputFilename = os.path.join(outputDirEventQA, "hEventNclusters" + fileFormat)
    plotHist(hEventNclusters, outputFilename, "hist E")
  else:
    # N tracks vs. Centrality
    hEventNtracksCentrality = eventQATHnSparse.Projection(1,0)
    outputFilename = os.path.join(outputDirEventQA, "hEventNtracksCentrality" + fileFormat)
    plotHist(hEventNtracksCentrality, outputFilename, "colz", False, True)

    # N clusters vs. Centrality
    hEventNclustersCentrality = eventQATHnSparse.Projection(3,0)
    outputFilename = os.path.join(outputDirEventQA, "hEventNclustersCentrality" + fileFormat)
    plotHist(hEventNclustersCentrality, outputFilename, "colz", False, True)

  if ispp:
    # Plot leading cluster energy
    hEventEmcalLeadClusE = eventQATHnSparse.Projection(3)
    outputFilename = os.path.join(outputDirEventQA, "hEventLeadClusE" + fileFormat)
    plotHist(hEventEmcalLeadClusE, outputFilename, "hist E", True)
  else:
    # Plot leading cluster energy vs. Centrality
    hEventLeadClusECentrality = eventQATHnSparse.Projection(4,0)
    outputFilename = os.path.join(outputDirEventQA, "hEventLeadClusECentrality" + fileFormat)
    plotHist(hEventLeadClusECentrality, outputFilename, "colz", False, True)

  # Event rejection reasons
  EventCutList = qaList.FindObject("EventCutOutput")

  hEventReject = EventCutList.FindObject("fCutStats")
  hEventReject.GetYaxis().SetTitle("N events accepted")
  outputFilename = os.path.join(outputDirEventQA, "hEventReject" + fileFormat)
  textNEvents = ROOT.TLatex()
  textNEvents.SetNDC()
  textNEvents.DrawLatex(0.65,0.87,"#it{N}_{events} = %d" % nEvents)

  plotHist(hEventReject, outputFilename, "hist", False)

########################################################################################################
# Plot Pt-hard histograms   ##############################################################################
########################################################################################################
def plotPtHard(f, qaList, nEvents, qaListRef, nEventsRef, outputDir, fileFormat):
  
  # Note: errors have not been propagated correctly for Pt-hard histos, so we do not plot them.
  
  # Create subdirectory for PtHard
  outputDirPtHard = outputDir + "PtHard/"
  if not os.path.exists(outputDirPtHard):
    os.makedirs(outputDirPtHard)
  
  ROOT.gStyle.SetOptTitle(1)

  hNEvents = f.Get("hNEventsAcc")
  outputFilename = os.path.join(outputDirPtHard, "hPtHardNEvents" + fileFormat)
  plotHist(hNEvents, outputFilename, "hist")
  
  hXSecPerEvent = f.Get("hXSecPerEvent")
  if hXSecPerEvent:
    outputFilename = os.path.join(outputDirPtHard, "hPtHardXSecPerEvent" + fileFormat)
    plotHist(hXSecPerEvent, outputFilename, "hist", True)
  
  hNTrialsPerEvent = f.Get("hNTrialsPerEvent")
  if hNTrialsPerEvent:
    outputFilename = os.path.join(outputDirPtHard, "hPtHardNTrialsPerEvent" + fileFormat)
    plotHist(hNTrialsPerEvent, outputFilename, "hist")
  
  hScaleFactor = f.Get("hScaleFactor")
  if hScaleFactor:
    outputFilename = os.path.join(outputDirPtHard, "hPtHardScaleFactor" + fileFormat)
    plotHist(hScaleFactor, outputFilename, "hist", True)

  hPtHard = qaList.FindObject("hPtHard")
  outputFilename = os.path.join(outputDirPtHard, "hPtHard" + fileFormat)
  plotHist(hPtHard, outputFilename, "hist", True)

  #if a reference is provided
  if qaListRef:
  
    hPtHardRef = qaListRef.FindObject("hPtHard")
  
    outputFilename = os.path.join(outputDirPtHard, "hPtHard_Ratio" + fileFormat)
    xRangeMax = 100
    yAxisTitle = "#frac{1}{N_{evts}}#frac{dN}{dp_{T}} [GeV^{-1}]"
    legendTitle = "pT Hard production"
    legendRunLabel = "hPtHard this run"
    legendRefLabel = "hPtHard all runs"
    ratioYAxisTitle = "Ratio: run / all runs"
    hPtHardRef.SetLineColor(1)

    ispp=1
    if nEventsRef!=0:
      plotSpectra(hPtHard, hPtHardRef,0x0, 0x0, nEvents, nEventsRef, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, "1", "2", "3")

########################################################################################################
# Plot basic histogram    ##############################################################################
########################################################################################################
def plotHist(h, outputFilename, drawOptions = "", setLogy = False, setLogz = False):
 
  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd()
  if setLogy:
    c.SetLogy()
  if setLogz:
    c.SetLogz()

  h.Draw(drawOptions)
  c.SaveAs(outputFilename)
  c.Close()

########################################################################################################
# Plot spectra (and ratio, if reference file suppled)   ################################################
########################################################################################################
def plotSpectra(h, hRef, h2, h3, nEvents, nEventsRef, ispp, xRangeMax, yAxisTitle, legendTitle, legendRunLabel, legendRefLabel, ratioYAxisTitle, outputFilename, scalingOptions = "", h2LegendLabel = "", h3LegendLabel = ""):

  h.SetLineColor(1)
  h.SetLineWidth(2)
  h.SetLineStyle(1)
  
  h.Scale(1./nEvents, scalingOptions)

  h.GetYaxis().SetTitle(yAxisTitle)
  h.GetYaxis().SetTitleSize(0.06)
  h.GetXaxis().SetRangeUser(0,xRangeMax)
  h.GetYaxis().SetRangeUser(2e-9,2e3)
  if ispp:
    h.GetYaxis().SetRangeUser(2e-11,20)
  h.GetYaxis().SetLabelFont(43)
  h.GetYaxis().SetLabelSize(20)

  if h2:
    h2.SetLineColor(2)
    h2.SetLineWidth(2)
    h2.SetLineStyle(1)
    h2.Scale(1./nEvents, scalingOptions)
    h2.GetYaxis().SetTitle(yAxisTitle)
    h2.GetYaxis().SetTitleSize(0.06)
    h2.GetXaxis().SetRangeUser(0,xRangeMax)
    h2.GetYaxis().SetRangeUser(2e-9,2e3)
    if ispp:
      h2.GetYaxis().SetRangeUser(2e-11,20)
    h2.GetYaxis().SetLabelFont(43)
    h2.GetYaxis().SetLabelSize(20)
    h2.GetXaxis().SetTitleOffset(1.4)
  if h3:
    h3.SetLineStyle(1)
    h3.SetLineColor(4)
    h3.SetLineWidth(2)
    h3.Scale(1./nEvents, scalingOptions)

  if not hRef:
    c = ROOT.TCanvas("c","c: pT",600,450)
    c.cd()

    ROOT.gPad.SetLeftMargin(0.16)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetLogy()
  
    if h2 and h3:
      h2.Draw("hist E")
      h3.Draw("hist E same")
      h.Draw("hist E same")
    elif h2:
      h2.Draw("hist E")
      h.Draw("hist E same")
    else:
      h.Draw("hist E")

  if hRef:
    c = ROOT.TCanvas("c","c: pT",800,850)
    c.cd()
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)
    pad1.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.05)
    pad1.SetTopMargin(0.05)
    pad1.SetLogy()
    pad1.Draw()
    pad1.cd()

    if h2 and h3:
      h2.Draw("hist")
      h3.Draw("hist same")
      h.Draw("hist same")
    elif h2:
      h2.Draw("hist E")
      h.Draw("hist E same")
    else:
      h.Draw("hist E")

    hRef.SetLineColor(8)
    if h2 and not h3: # hack to keep color scheme consistent in cluster spectra ratio
      hRef.SetLineColor(4)
    hRef.SetMarkerColor(1)
    hRef.SetLineStyle(1)

    hRef.Scale(1./nEventsRef, scalingOptions)
    hRef.Draw("hist E same")
    
    c.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.15)
    pad2.SetRightMargin(0.05)
    pad2.Draw()
    pad2.cd()
    
    hRatio = h.Clone()
    hRatio.Divide(hRef)
    
    hRatio.SetMarkerStyle(20)
    hRatio.SetMarkerSize(0.5)
    hRatio.SetMarkerColor(1)

    if h2:
      hRatio2 = h2.Clone()
      hRatio2.Divide(hRef)
      hRatio2.SetMarkerStyle(21)
      hRatio2.SetMarkerColor(2)
      
      hRatio2.GetYaxis().SetTitle(ratioYAxisTitle)
      hRatio2.GetYaxis().SetTitleSize(20)
      hRatio2.GetYaxis().SetTitleFont(43)
      hRatio2.GetYaxis().SetTitleOffset(2.2)
      hRatio2.GetYaxis().SetLabelFont(43)
      hRatio2.GetYaxis().SetLabelSize(20)
      hRatio2.GetYaxis().SetNdivisions(505)
      hRatio2.GetYaxis().SetRangeUser(0,2.2)
      if ratioYAxisTitle in "Ratio to all":
        hRatio2.GetYaxis().SetRangeUser(0,1.2)
      hRatio2.GetXaxis().SetRangeUser(0,xRangeMax)
      hRatio2.GetXaxis().SetTitleSize(30)
      hRatio2.GetXaxis().SetTitleFont(43)
      hRatio2.GetXaxis().SetTitleOffset(4.)
      hRatio2.GetXaxis().SetLabelFont(43)
      hRatio2.GetXaxis().SetLabelSize(20)
    if h3:
      hRatio3 = h3.Clone()
      hRatio3.Divide(hRef)
      hRatio3.SetMarkerStyle(21)
      hRatio3.SetMarkerColor(4)
    if h2 and h3:
      hRatio2.Draw("P E")
      hRatio3.Draw("P E same")
      hRatio.Draw("P E same")
    elif h2:
      hRatio2.GetYaxis().SetRangeUser(0,25)
      hRatio2.Draw("P E")
      hRatio.Draw("P E same")
    if not h2 and not h3:
      hRatio.GetYaxis().SetTitle(ratioYAxisTitle)
      hRatio.GetYaxis().SetTitleSize(20)
      hRatio.GetYaxis().SetTitleFont(43)
      hRatio.GetYaxis().SetTitleOffset(2.2)
      hRatio.GetYaxis().SetLabelFont(43)
      hRatio.GetYaxis().SetLabelSize(20)
      hRatio.GetYaxis().SetNdivisions(505)
      hRatio.GetYaxis().SetRangeUser(0,2.2)

      hRatio.GetXaxis().SetRangeUser(0,xRangeMax)
      hRatio.GetXaxis().SetTitleSize(30)
      hRatio.GetXaxis().SetTitleFont(43)
      hRatio.GetXaxis().SetTitleOffset(4.)
      hRatio.GetXaxis().SetLabelFont(43)
      hRatio.GetXaxis().SetLabelSize(20)

      hRatio.Draw("P E")

    pad1.cd()

  if nEvents > 1:
    textNEvents = ROOT.TLatex()
    textNEvents.SetNDC()
    textNEvents.DrawLatex(0.55,0.6,"#it{N}_{events} = %d" % nEvents)

  leg2 = ROOT.TLegend(0.3,0.7,0.88,0.93,legendTitle)
  leg2.SetFillColor(10)
  leg2.SetBorderSize(0)
  leg2.SetFillStyle(0)
  leg2.SetTextSize(0.04)
  leg2.AddEntry(h, legendRunLabel, "l")
  if h2:
    leg2.AddEntry(h2, h2LegendLabel, "l")
  if h3:
    leg2.AddEntry(h3, h3LegendLabel, "l")
  if hRef:
    leg2.AddEntry(hRef, legendRefLabel, "l")
  leg2.Draw("same")
  
  c.SaveAs(outputFilename)
  c.Close()

########################################################################################################
# Plot spectra and ratio (h2,h3 will be divided by h)   ################################################
########################################################################################################
def plotSpectraCent(h, h2, h3, nEvents, ispp, outputFilename, xRangeMax, yAxisTitle, ratioYAxisTitle, legendTitle, h1legendLabel, h2legendLabel, h3legendLabel = "", scalingOptions = "", yRatioMax = 32):
  
  h.SetLineColor(4)
  if not h3:
    h.SetLineColor(2)
  h.SetLineWidth(2)
  h.SetLineStyle(1)
  
  h.Scale(1./nEvents, scalingOptions)

  h.GetYaxis().SetTitle(yAxisTitle)
  h.GetYaxis().SetTitleSize(0.06)
  h.GetXaxis().SetRangeUser(0,xRangeMax)
  h.GetYaxis().SetRangeUser(2e-9,2e3)
  if ispp:
    h.GetYaxis().SetRangeUser(2e-11,20)
  h.GetYaxis().SetLabelFont(43)
  h.GetYaxis().SetLabelSize(20)
  
  h2.SetLineColor(1)
  h2.SetLineWidth(2)
  h2.SetLineStyle(1)
  h2.Scale(1./nEvents, scalingOptions)
  h2.GetYaxis().SetTitle(yAxisTitle)
  h2.GetYaxis().SetTitleSize(0.06)
  h2.GetXaxis().SetRangeUser(0,xRangeMax)
  h2.GetYaxis().SetRangeUser(2e-9,2e3)
  if ispp:
    h2.GetYaxis().SetRangeUser(2e-11,20)
  h2.GetYaxis().SetLabelFont(43)
  h2.GetYaxis().SetLabelSize(20)
  if h3:
    h3.SetLineStyle(1)
    h3.SetLineColor(2)
    h3.SetLineWidth(2)
    h3.Scale(1./nEvents, scalingOptions)

  c = ROOT.TCanvas("c","c: pT",800,850)
  c.cd()
  pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
  pad1.SetBottomMargin(0)
  pad1.SetLeftMargin(0.15)
  pad1.SetRightMargin(0.05)
  pad1.SetTopMargin(0.05)
  pad1.SetLogy()
  pad1.Draw()
  pad1.cd()
  
  if h3:
    h2.Draw("hist E")
    h3.Draw("hist E same")
    h.Draw("hist E same")
  else:
    h2.Draw("hist E")
    h.Draw("hist E same")

  c.cd()
  pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
  pad2.SetTopMargin(0)
  pad2.SetBottomMargin(0.35)
  pad2.SetLeftMargin(0.15)
  pad2.SetRightMargin(0.05)
  pad2.Draw()
  pad2.cd()

  hRatio = h2.Clone()
  hRatio.Divide(h)

  hRatio.SetMarkerStyle(21)
  hRatio.SetMarkerColor(1)

  if h3:
    hRatio2 = h3.Clone()
    hRatio2.Divide(h)
    hRatio2.SetMarkerStyle(21)
    hRatio2.SetMarkerColor(2)
    
    hRatio2.GetYaxis().SetTitle(ratioYAxisTitle)
    hRatio2.GetYaxis().SetTitleSize(20)
    hRatio2.GetYaxis().SetTitleFont(43)
    hRatio2.GetYaxis().SetTitleOffset(2.2)
    hRatio2.GetYaxis().SetLabelFont(43)
    hRatio2.GetYaxis().SetLabelSize(20)
    hRatio2.GetYaxis().SetNdivisions(505)
    hRatio2.GetYaxis().SetRangeUser(0,yRatioMax)
    hRatio2.GetXaxis().SetRangeUser(0,xRangeMax)
    hRatio2.GetXaxis().SetTitleSize(30)
    hRatio2.GetXaxis().SetTitleFont(43)
    hRatio2.GetXaxis().SetTitleOffset(4.)
    hRatio2.GetXaxis().SetLabelFont(43)
    hRatio2.GetXaxis().SetLabelSize(20)
    
    hRatio2.Draw("P E")
    hRatio.Draw("P E same")
  else:
    hRatio.GetYaxis().SetTitle(ratioYAxisTitle)
    hRatio.GetYaxis().SetTitleSize(20)
    hRatio.GetYaxis().SetTitleFont(43)
    hRatio.GetYaxis().SetTitleOffset(2.2)
    hRatio.GetYaxis().SetLabelFont(43)
    hRatio.GetYaxis().SetLabelSize(20)
    hRatio.GetYaxis().SetNdivisions(505)
    hRatio.GetYaxis().SetRangeUser(0,yRatioMax)
    
    hRatio.GetXaxis().SetRangeUser(0,xRangeMax)
    hRatio.GetXaxis().SetTitleSize(30)
    hRatio.GetXaxis().SetTitleFont(43)
    hRatio.GetXaxis().SetTitleOffset(4.)
    hRatio.GetXaxis().SetLabelFont(43)
    hRatio.GetXaxis().SetLabelSize(20)
    
    hRatio.Draw("P E")

  pad1.cd()
  leg = ROOT.TLegend(0.3,0.7,0.88,0.93,legendTitle)
  leg.SetFillColor(10)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.04)
  leg.AddEntry(h2, h2legendLabel, "l")
  if h3:
    leg.AddEntry(h3, h3legendLabel, "l")
  leg.AddEntry(h, h1legendLabel, "l")
  leg.Draw("same")
  
  c.SaveAs(outputFilename)
  c.Close()
########################################################################################################
# Plot spectra (and ratio, if reference file suppled)   ################################################
########################################################################################################
def plotNEFSpectra(h, h2, h3, nEvents, ispp, xRangeMax, yAxisTitle, h1legendLabel, outputFilename, scalingOptions = "", h2legendLabel = "", h3legendLabel = ""):
  
  h.SetLineColor(1)
  h.SetLineWidth(2)
  h.SetLineStyle(1)
  
  h.Scale(1./nEvents, scalingOptions)
  if ispp:
    h.GetYaxis().SetRangeUser(0.0000005, 0.05)
  h.GetYaxis().SetTitle(yAxisTitle)
  h.GetYaxis().SetTitleSize(0.06)
  h.GetXaxis().SetRangeUser(0,xRangeMax)
  h.GetYaxis().SetLabelFont(43)
  h.GetYaxis().SetLabelSize(20)

  if h2:
    h2.SetLineColor(2)
    h2.SetLineWidth(2)
    h2.SetLineStyle(1)
    h2.Scale(1./nEvents, scalingOptions)
    h2.GetYaxis().SetTitle(yAxisTitle)
    h2.GetYaxis().SetTitleSize(0.06)
    h2.GetXaxis().SetRangeUser(0,xRangeMax)
    #h2.GetYaxis().SetRangeUser(2e-9,2e3)
    if ispp:
      h2.GetYaxis().SetRangeUser(5e-7,0.05)
      h2.GetYaxis().SetLabelFont(43)
      h2.GetYaxis().SetLabelSize(20)
      h2.GetXaxis().SetTitleOffset(1.4)

  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd().SetLogy()

  ROOT.gPad.SetLeftMargin(0.16)
  ROOT.gPad.SetRightMargin(0.05)
  ROOT.gPad.SetBottomMargin(0.14)
  ROOT.gPad.SetTopMargin(0.05)
  
  if h3:
    h2.Draw("hist E")
    h3.Draw("hist E same")
    h.Draw("hist E same")
  else:
    h2.Draw("hist E")
    h.Draw("hist E same")

  leg = ROOT.TLegend(0.3,0.7,0.88,0.93)
  leg.SetFillColor(10)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.04)
  leg.AddEntry(h, h1legendLabel, "l")
  if h3:
    leg.AddEntry(h3, h3legendLabel, "l")
  leg.AddEntry(h2, h2legendLabel, "l")
  leg.Draw("same")
  
  c.SaveAs(outputFilename)
  c.Close()

#########################################################################################
# Function to iterate recursively through an object to set Sumw2 on all TH1/TH2/THnSparse
#########################################################################################
def SetSumw2(obj):
  if obj.InheritsFrom(ROOT.TProfile.Class()):
    pass
    #print("Sumw2 not called for TProfile %s" % obj.GetName())
  elif obj.InheritsFrom(ROOT.TH2.Class()):
    obj.Sumw2()
    #print("Sumw2 called on TH2 %s" % obj.GetName())
  elif obj.InheritsFrom(ROOT.TH1.Class()):
    obj.Sumw2()
    #print("Sumw2 called on TH1 %s" % obj.GetName())
  elif obj.InheritsFrom(ROOT.THnSparse.Class()):
    obj.Sumw2()
    #print("Sumw2 called on THnSparse %s" % obj.GetName())
  else:
    #print("Not a histogram!")
    #print obj.GetName()
    for subobj in obj:
      SetSumw2(subobj)

#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description="Compare histograms to test the new EMCal corrections framework")
  parser.add_argument("-f", "--inputFile", action="store",
                      type=str, metavar="inputFile",
                      default="AnalysisResults.root",
                      help="Path of AnalysisResults.root file")
  parser.add_argument("-o", "--outputDir", action="store",
                      type=str, metavar="outputDir",
                      default="./outputQA/",
                      help="Output directory for QA plots to be written to")
  parser.add_argument("-r", "--referenceFile", action="store",
                      type=str, metavar="referenceFile",
                      default="",
                      help="Reference root file for the inputFile histos to be compared to (when doing run-by-run QA)")
  parser.add_argument("-i", "--imageFormat", action="store",
                      type=str, metavar="imageFormat",
                      default=".pdf",
                      help="Image format to save plots in, e.g. \".pdf\" or \".png\"")

  # Parse the arguments
  args = parser.parse_args()
  
  print("Configuring...")
  print("inputFile: \"{0}\"".format(args.inputFile))
  print("ouputDir: \"{0}\"".format(args.outputDir))
  print("referenceFile: \"{0}\"".format(args.referenceFile))
  print("imageFormat: \"{0}\"".format(args.imageFormat))

  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFile):
    print("File \"{0}\" does not exist! Exiting!".format(args.inputFile))
    sys.exit(0)
  
  plotPWGJEQA(inputFile = args.inputFile, outputDir = args.outputDir, referenceFile = args.referenceFile, fileFormat = args.imageFormat)
