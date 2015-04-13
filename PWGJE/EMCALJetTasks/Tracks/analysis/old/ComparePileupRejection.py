#! /usr/bin/env python
import sys
from getopt import getopt,GetoptError
from copy import deepcopy
from ROOT import TCanvas,TFile,TH1F,TLegend,TPaveText
from ROOT import kRed,kBlack,kBlue,gROOT
from Helper import FileReaderException,HistNotFoundException,NormaliseBinWidth

gObjects = list()

def ReadHistograms(filename, trigger, cuts):
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
                        "SpecWithPileup":"hPt%s_nopr_%s" %(trigger, cuts),\
                        "SpecNoPileup":"hPt%s_wpr_%s" %(trigger, cuts)}
        for key,hname in histos.iteritems():
                hist = histlist.FindObject(hname)
                if not hist:
                        raise HistNotFoundException(hname)
                else:
                        hist.SetDirectory(gROOT)
                        resultlist[key] = hist
        return resultlist

def MakePlot(specWithPileup, specWithoutPileup, options):
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

        specWithPileup.SetMarkerColor(kRed)
        specWithPileup.SetLineColor(kRed)
        specWithPileup.SetMarkerStyle(24)
        specWithoutPileup.SetMarkerColor(kBlue)
        specWithoutPileup.SetLineColor(kBlue)
        specWithoutPileup.SetMarkerStyle(25)
        specWithPileup.Draw("epsame")
        specWithoutPileup.Draw("epsame")
        gObjects.append(specWithPileup)
        gObjects.append(specWithoutPileup)

        leg = TLegend(0.1,0.15, 0.45, 0.3)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.045)
        leg.AddEntry(specWithPileup, "Without pileup rejection", "lep")
        leg.AddEntry(specWithoutPileup, "With pileup rejection", "lep");
        leg.Draw()
        gObjects.append(leg)

        label = TPaveText(0.5, 0.75, 0.89, 0.89, "NDC")
        label.SetBorderSize(0)
        label.SetFillStyle(0)
        label.SetTextFont(42)
        label.SetTextAlign(12)
        label.AddText("Trigger: %s" %(options["Trigger"]))
        label.AddText("Cuts:    %s" %(options["Cuts"]))
        label.Draw()
        gObjects.append(label)

        ratiopad = plot.cd(2)
        ratiopad.SetGrid(False, False)
        ratiopad.SetLogx(True)

        axirat = TH1F("axirat", "; p_{t} (GeV/c); Without/with pileup rejection", 1000, 0., 100.)
        axirat.SetStats(False)
        axirat.GetYaxis().SetRangeUser(0.8, 1.2)
        axirat.Draw("axis")
        gObjects.append(axirat)

        ratiohist = deepcopy(specWithPileup)
        ratiohist.SetName("ratioPileup")
        ratiohist.Divide(ratiohist, specWithoutPileup, 1., 1., "b")
        ratiohist.SetMarkerColor(kBlack)
        ratiohist.SetLineColor(kBlack)
        ratiohist.SetMarkerStyle(20)
        ratiohist.Draw("epsame")
        gObjects.append(ratiohist)

        plot.cd()
        gObjects.append(plot)
        return plot
        
def ComparePileupRejection(filename, arguments):
        try:
                opt,arg = getopt(arguments,"e:mns", ["emcal=","minbias","nocut","stdcut"])
        except GetoptError as e:
                print "Error, wrong argument(s)"
                print str(e)
                sys.exit(1)

        emctriggers = ["EMCJHigh", "EMCJLow", "EMCGHigh", "EMCGLow"]
        trigger = "MinBias"
        cuts = "stdcut"
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
                elif o in ("-n", "--nocut"):
                        cuts = "nocut"
                elif o in ("-s", "--stdcut"):
                        cuts = "stdcut"

        histos = None
        try:
                histos = ReadHistograms(filename, trigger, cuts)
        except FileReaderException as f:
                print str(f)
                sys.exit(1)
        except HistNotFoundException as e:
                print str(e)
                sys.exit(1)
        
        # Normalise by binwidth and the number of events
        # Normalise spectrum with pileup rejection by number of event after pileup rejection
        # and spectrum without pileup rejection by events without pileup rejection
        specWithPileup = histos["SpecWithPileup"].ProjectionY("SpecWithPileup")
        specWithPileup.Sumw2()
        specWithoutPileup = histos["SpecNoPileup"].ProjectionY("SpecWithoutPileup")
        specWithoutPileup.Sumw2()
        NormaliseBinWidth(specWithPileup)
        NormaliseBinWidth(specWithoutPileup)
        specWithPileup.Scale(1./histos["EventCounter"].GetBinContent(1))
        specWithoutPileup.Scale(1./histos["EventCounter"].GetBinContent(2))
        options = {"Trigger" : trigger, "Cuts" : "ON"}
        if cuts == "nocut":
                options["Cuts"] = "OFF"
        plot = MakePlot(specWithPileup, specWithoutPileup, options)

def run(filename, arguments):
        arglist = arguments.split(" ")
        ComparePileupRejection(filename, arglist)

def main():
        if len(sys.argv) < 2:
                print "At least the filename required as argument"
                sys.exit(1)
        filename = sys.argv[len(sys.argv)-1]
        arguments = sys.argv[1:len(sys.argv)-2]
        ComparePileupRejection(filename, arguments)

if __name__ == "__main__":
        main()
