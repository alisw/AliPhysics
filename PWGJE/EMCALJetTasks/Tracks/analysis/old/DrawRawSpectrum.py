#! /usr/bin/env python
# 
# Script plotting the raw spectrum for different trigger classes and event selections
# Can run in batch mode or in interactive mode. For the interactive mode load the script
# as module and execute the run function with the inputfile name and the command line 
# options.
#
#    Author: Markus Fasel
#

import os, sys, getopt
from ROOT import TCanvas,TFile,TH1F,TPaveText,gROOT,kRed
from Helper import NormaliseBinWidth,FileReaderException,HistNotFoundException

gObjects = list()

def ReadFromFile(filename, options):
        """
        Read histogram with input spectrum from file.
        """
        infile = TFile.Open(filename)
        if not infile or infile.IsZombie():
                raise FileReaderException(filename)
        myresultlist = infile.Get("results")
        histlist = myresultlist.FindObject("List of histograms of container PtEMCalTriggerHistograms")

        resultlist = dict()
        hEC = histlist.FindObject("hEvents%s" %(options["Trigger"]))
        if not hEC:
                infile.Close()
                raise HistNotFoundException(histnameEC)
        else:
                hEC.SetDirectory(gROOT)
                resultlist["EventCounter"] = hEC
        hSpec = histlist.FindObject("hPt%s_%s_%s" %(options["Trigger"], options["Eventselection"], options["Trackselection"]))
        if not hSpec:
                hSpec.SetDirectory(gROOT)
                infile.Close()
                raise HistNotFoundException(histnameSpec)
        else:
                resultlist["Spectrum"] = hSpec
        return resultlist

def MakePlot(spectrum, options):
        """
        Create final plot of the spectrum.
        """
        plot = TCanvas("crawspectrum", "Raw charged hadron spectrum")
        plot.cd()
        plot.SetGrid(False, False)
        plot.SetLogx()
        plot.SetLogy()

        axes = TH1F("axis", "; p_{T} (GeV/c); 1/N_{ev} dN/dp_{t} ((GeV/c)^{-1})", 1000, 0., 100.)
        axes.SetStats(False)
        axes.GetYaxis().SetRangeUser(1e-9, 100)
        axes.Draw("axis")
        gObjects.append(axes)

        spectrum.SetMarkerColor(kRed)
        spectrum.SetLineColor(kRed)
        spectrum.SetMarkerStyle(24)
        spectrum.Draw("epsame")
        gObjects.append(spectrum)

        label = TPaveText(0.65, 0.7, 0.89, 0.89, "NDC")
        label.SetBorderSize(0)
        label.SetFillStyle(0)
        label.SetTextAlign(12)
        cutstring = "ON"
        if options["Trackselection"] == "nocut":
                cutstring = "OFF"
        prstring = "ON"
        if options["Eventselection"] == "nopr":
                prstring = "OFF"
        label.AddText("Trigger:          %s" %(options["Trigger"]))
        label.AddText("Standard cuts:    %s" %(cutstring))
        label.AddText("Pileup rejection: %s" %(prstring))
        label.Draw()
        gObjects.append(label)

        gObjects.append(plot)
        return plot

def Usage():
        """
        Print help text.
        """
        print " Usage: "
        print "   "
        print "   DrawRawSpectrum.py [OPTIONS] filename"
        print "   "
        print " Options:"
        print "   -e/--emcal [1:4]:         Use EMCal triggered events. The number is the trigger type"
        print "                  1:         EMCal Jet High Threshold"
        print "                  2:         EMCal Jet Low Threshold"
        print "                  3:         EMCal Gamma High Threshold"
        print "                  4:         EMCal Gamma Low Threshold"
        print "   -h/--help:                Print help text"
        print "   -m/--minbias:             Use min. bias events"
        print "   -n/--nocuts:              No cuts applied"
        print "   -o/--outputfile name:     Write spectrum to output file"
        print "   -p/--plotname name:       Write plot to output file"
        print "   -s/--standardcuts:        Apply standard cuts"
        print "   -r/--pilupreject:         Pileup rejection applied"
        print "   -w/--withpilup:           No pileup rejection applied"

def launch(filename, arglist):
        try:
                opt, arg = getopt.getopt(arglist, "e:mno:p:rsw", ["emcal=","help","minbias","nocuts","outputfile=","plotname=","standardcuts","pileupreject","withpileup"])
        except getopt.GetoptError as e:
                print "Invalid argument(s)"
                print str(e)
                sys.exit(1)

        options = dict()
        # Set defaults
        options["Trigger"]            = "MinBias"
        options["Trackselection"]     = "stdcut"
        options["Eventselection"]     = "wpr"
        options["plotfile"]           = ""
        options["outputfile"]         = ""

        for o,a in opt:
                if o in ("--e", "--emcal"):
                        emcalTrigger = int(a)
                        if not emcalTrigger in range(1,5):
                                print "EMCal Trigger invalid"
                                Usage()
                                sys.exit(1)
                        else:
                                if emcalTrigger == 1:
                                        options["Trigger"] = "EMCJHigh"
                                elif emcalTrigger == 2:
                                        options["Trigger"] = "EMCJLow"
                                elif emcalTrigger == 3:
                                        options["Trigger"] = "EMCGHigh"
                                elif emcalTrigger == 4:
                                        options["Trigger"] = "EMCGLow"
                elif o in ("-h", "--help"):
                        Usage()
                        sys.exit(1)
                elif o in ("-m", "--minBias") :
                        options["Trigger"] = "MinBias"
                elif o in ("-n", "--nocuts") :
                        options["Trackselection"] = "nocut"
                elif o in ("-o", "--outputfile"):
                        options["outputfile"] = str(a)
                elif o in ("-p", "--outputfile"):
                        options["plotfile"] = str(a)
                elif o in ("-s", "--standardcuts"):
                        options["Trackselection"] = "stdcut"
                elif o in ("-r", "-pileupreject"):
                        options["Eventselection"] = "wpr"
                elif o in ("-w", "-withpileup"):
                        options["Eventselection"] = "nopr"

        try:
                speclist = ReadFromFile(filename, options);
        except FileReaderException as e:
                print str(e)
                sys.exit(1)
        except HistNotFoundException as e:
                print str(e)
                sys.exit(1)

        nevents = 0
        if options["Eventselection"] == "nopr":
                nevents = speclist["EventCounter"].GetBinContent(1)
        else:
                nevents = speclist["EventCounter"].GetBinContent(2)

        spectrum = speclist["Spectrum"].ProjectionY("rawspectrum")
        spectrum.Sumw2()
        spectrum.Scale(1./nevents)
        NormaliseBinWidth(spectrum)

        plot = MakePlot(spectrum, options)
        if len(options["plotfile"]):
                plot.SaveAs(options["plotfile"])
        if len(options["outputfile"]):
                outputfile = TFile(options["outputfile"], "RECREATE")
                if outputfile:
                        outputfile.cd()
                        spectrum.Write(spectrum.GetName())
                        outputfile.Close()

def run(filename,argstring):
        """
        Run plotting in interactive mode.
        """
        arglist = argstring.split(" ")
        launch(filename, arglist)

def main():
        """
        Run plotting in batch mode.
        """
        if len(sys.argv) < 2:
                Usage()
                sys.exit(1)

        filename = sys.argv[len(sys.argv) - 1]
        
        launch(filename, sys.argv[1:])

if __name__ == "__main__":
        main()
