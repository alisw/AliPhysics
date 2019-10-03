#! /usr/bin/env python

import sys
from getopt import getopt,GetoptError
from copy import deepcopy
from ROOT import TCanvas,TFile,TH1F,TLegend,TPaveText
from ROOT import kRed,kBlack,kBlue,gROOT
from Helper import FileReaderException,HistNotFoundException,NormaliseBinWidth

gObjects = list()

def ReadHistograms(filename, trigger, pileup):
        """
        Read histogram with input spectrum from file.
        """
        infile = TFile.Open(filename)
        if not infile or infile.IsZombie():
                raise FileReaderException(filename)
        myresultlist = infile.Get("results")
        histlist = myresultlist.FindObject("List of histograms of container PtEMCalTriggerHistograms")

        resultlist = dict()
        histos = {"EventCounter":"hEvents%s" %(trigger),\
                        "SpecWithCuts":"hPt%s_%s_stdcut" %(trigger, pileup),\
                        "SpecNoCuts":"hPt%s_%s_nocut" %(trigger, pileup)}
        for key,hname in histos.iteritems():
                hist = histlist.FindObject(hname)
                if not hist:
                        raise HistNotFoundException(hname)
                else:
                        hist.SetDirectory(gROOT)
                        resultlist[key] = hist
        return resultlist

def MakePlot(specWithCuts, specWithoutCuts, options):
        """
        Creating the final plot. The plot will be a two panel plot:
        1st panel shows the comparison of the spectra. 
        2nd panel shows the ratio of the normalised spectra with and 
            without pileup rejection (as without/with)
        """
        plot = TCanvas("ComparisonPileup", "Comparison of raw spectra with/without pileup rejection", 1000, 500)
        plot.Divide(2,1)

        specpad = plot.cd(1)
        specpad.SetGrid(False, False)
        specpad.SetLogx(True)
        specpad.SetLogy(True)

        axispec = TH1F("axispec","; p_{t} (GeV/c); 1/N_{ev} dN/dp_{t} ((GeV/c)^{-1})", 1000, 0., 100.)
        axispec.SetStats(False)
        axispec.GetYaxis().SetRangeUser(1e-10,100)
        axispec.Draw("axis")
        gObjects.append(axispec)

        specWithCuts.SetMarkerColor(kRed)
        specWithCuts.SetLineColor(kRed)
        specWithCuts.SetMarkerStyle(24)
        specWithoutCuts.SetMarkerColor(kBlue)
        specWithoutCuts.SetLineColor(kBlue)
        specWithoutCuts.SetMarkerStyle(25)
        specWithCuts.Draw("epsame")
        specWithoutCuts.Draw("epsame")
        gObjects.append(specWithCuts)
        gObjects.append(specWithoutCuts)

        leg = TLegend(0.1,0.15, 0.45, 0.3)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.045)
        leg.AddEntry(specWithCuts, "With standard cuts", "lep")
        leg.AddEntry(specWithoutCuts, "Without standard cuts", "lep");
        leg.Draw()
        gObjects.append(leg)

        label = TPaveText(0.5, 0.75, 0.89, 0.89, "NDC")
        label.SetBorderSize(0)
        label.SetFillStyle(0)
        label.SetTextFont(42)
        label.SetTextAlign(12)
        label.AddText("Trigger:          %s" %(options["Trigger"]))
        label.AddText("Pileup rejection: %s" %(options["pileupreject"]))
        label.Draw()
        gObjects.append(label)

        ratiopad = plot.cd(2)
        ratiopad.SetGrid(False, False)
        ratiopad.SetLogx(True)

        axirat = TH1F("axirat", "; p_{t} (GeV/c); Without/with standard cuts", 1000, 0., 100.)
        axirat.SetStats(False)
        axirat.GetYaxis().SetRangeUser(0., 100)
        axirat.Draw("axis")
        gObjects.append(axirat)

        ratiohist = deepcopy(specWithoutCuts)
        ratiohist.SetName("ratioPileup")
        ratiohist.Divide(ratiohist, specWithCuts, 1., 1., "b")
        ratiohist.SetMarkerColor(kBlack)
        ratiohist.SetLineColor(kBlack)
        ratiohist.SetMarkerStyle(20)
        ratiohist.Draw("epsame")
        gObjects.append(ratiohist)

        plot.cd()
        gObjects.append(plot)

        return plot

def CompareTrackCuts(filename, arguments):
        try:
                opt,arg = getopt(arguments,"e:mnw", ["emcal=","minbias","nopr","withpr"])
        except GetoptError as e:
                print "Error, wrong argument(s)"
                print str(e)
                sys.exit(1)

        emctriggers = ["EMCJHigh", "EMCJLow", "EMCGHigh", "EMCGLow"]
        trigger = "MinBias"
        pileup = "wpr"
        for o,a in opt:
                if o in ("-e", "--emcal"):
                        mytrg = int(a)
                        if not mytrg in range(1, 5):
                                print "Trigger out of range"
                                sys.exit(1)
                        else:
                                trigger = emctriggers[mytrg-1]
                elif o in ("-m", "--minbias"):
                        trigger = "MinBias"
                elif o in ("-n", "--nopr"):
                        pileup = "nopr"
                elif o in ("-w", "--withpr"):
                        pileup = "wpr"

        histos = None
        try:
                histos = ReadHistograms(filename, trigger, pileup)
        except FileReaderException as f:
                print str(f)
                sys.exit(1)
        except HistNotFoundException as e:
                print str(e)
                sys.exit(1)
        
        # Normalise by binwidth and the number of events
        # Normalise spectrum with pileup rejection by number of event after pileup rejection
        # and spectrum without pileup rejection by events without pileup rejection
        specWithCuts = histos["SpecWithCuts"].ProjectionY("SpecWithPileup")
        specWithCuts.Sumw2()
        specWithoutCuts = histos["SpecNoCuts"].ProjectionY("SpecWithoutCuts")
        specWithoutCuts.Sumw2()
        NormaliseBinWidth(specWithCuts)
        NormaliseBinWidth(specWithoutCuts)
        eventbin = 2
        if pileup == "nopr":
                eventbin = 1
        specWithCuts.Scale(1./histos["EventCounter"].GetBinContent(eventbin))
        specWithoutCuts.Scale(1./histos["EventCounter"].GetBinContent(eventbin))
        options = {"Trigger" : trigger, "pileupreject" : "ON"}
        if pileup == "nopr":
                options["pileupreject"] = "OFF"
        plot = MakePlot(specWithCuts, specWithoutCuts, options)

def runComparison(filename, arguments):
        arglist = arguments.split(" ")
        CompareTrackCuts(filename, arglist)

def main():
        if len(sys.argv) < 2:
                print "At least the filename required as argument"
                sys.exit(1)
        filename = sys.argv[len(sys.argv)-1]
        arguments = sys.argv[1:len(sys.argv)-2]
        ComparePileupRejection(filename, arguments)

if __name__ == "__main__":
        main()
