#! /usr/bin/env python
import sys
from copy import deepcopy
from getopt import getopt, GetoptError
from ROOT import TCanvas,TFile,TH1F,TLegend,TPaveText
from ROOT import gROOT,kRed,kGreen,kBlack,kBlue,kOrange
from Helper import FileReaderException,HistNotFoundException,NormaliseBinWidth
from Graphics import Style

gObjects = list()


def ReadHistograms(filename, trigger, cuts, pr):
        """
        Read histogram with input spectrum from file.
        """
        infile = TFile.Open(filename)
        if not infile or infile.IsZombie():
                raise FileReaderException(filename)
        myresultlist = infile.Get("results")
        histlist = myresultlist.FindObject("List of histograms of container PtEMCalTriggerHistograms")

        resultlist = dict()
        histos = {"EventCounter" : "hEvents%s" %(trigger),\
                        "Spec" :"hPt%s_%s_%s" %(trigger, pr, cuts)}
        for key,hname in histos.iteritems():
                hist = histlist.FindObject(hname)
                if not hist:
                        raise HistNotFoundException(hname)
                else:
                        hist.SetDirectory(gROOT)
                        resultlist[key] = hist
        return resultlist

def SetStyle(spectrum, style):
        """
        Set the style of the graphics (line and marker)
        """
        spectrum.SetMarkerColor(style.GetColor())
        spectrum.SetMarkerStyle(style.GetMarker())
        spectrum.SetLineColor(style.GetColor())

def MakePlot(spectra, triggers, options):
        """
        Make 2-panel comparison plot (left spectra for the different trigger classes, right the ratio
        to min bias)
        """
        plot = TCanvas("comparisonTriggers", "Comparison of the raw spectra for different triggers", 1000, 500)
        plot.Divide(2,1)

        ratios = dict()

        plotstyles = {"MinBias":Style(kBlack,20), "EMCJHigh":Style(kRed,24), "EMCJLow":Style(kBlue,25), \
                        "EMCGHigh":Style(kGreen+2, 26), "EMCGLow":Style(kOrange+2, 27) }

        specpad = plot.cd(1)
        specpad.SetGrid(False, False)
        specpad.SetLogx(True)
        specpad.SetLogy(True)
        axispec = TH1F("axispec", "; p_{T} (GeV/c); 1/N_{ev} dN/dp_{t} ((GeV/c)^{-1})", 1000, 0., 100.)
        axispec.SetStats(False)
        axispec.GetYaxis().SetRangeUser(1e-10, 100)
        axispec.Draw("axis")
        gObjects.append(axispec)

        leg = TLegend(0.15, 0.15, 0.35, 0.35)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        gObjects.append(leg)

        mbspectrum = spectra["MinBias"]
        for trigger,spectrum in spectra.iteritems():
                SetStyle(spectrum, plotstyles[trigger])
                spectrum.Draw("epsame")
                gObjects.append(spectrum)
                leg.AddEntry(spectrum, trigger, "lep")
                if trigger != "MinBias":
                        ratios[trigger] = deepcopy(spectrum)
                        ratios[trigger].Divide(mbspectrum)
                        SetStyle(ratios[trigger], plotstyles[trigger])
        leg.Draw()
        label = TPaveText(0.35, 0.75, 0.89, 0.89, "NDC")
        label.SetBorderSize(0)
        label.SetFillStyle(0)
        label.SetTextFont(42)
        label.SetTextAlign(12)
        label.AddText("Cuts:        %s" %(options["cuts"]))
        label.AddText("Pileup rej.: %s" %(options["pileup"]))
        label.Draw()
        gObjects.append(label)

        ratiopad = plot.cd(2)
        ratiopad.SetLogx(True)
        ratiopad.SetGrid(False, False)
        axiratio = TH1F("axiratio", "; p_{t} (GeV/c); Ratio to Min. Bias", 1000, 0.,100.)
        axiratio.GetYaxis().SetRangeUser(0., 1000.)
        axiratio.SetStats(False)
        axiratio.Draw("axis")
        gObjects.append(axiratio)

        for ratio in ratios.values():
                ratio.Draw("epsame")
                gObjects.append(ratio)

        plot.cd()
        gObjects.append(plot)
        return plot

def NormaliseSpectra(histos, pileuprejection):
        """
        Normalise spectra by the number of events (from the event counter) 
        and the width of each bin
        """
        spectra = dict()
        eventbin = 1
        if pileuprejection == True:
                eventbin = 2
        for trgclass,hlist in histos.iteritems():
                spectra[trgclass] = hlist["Spec"].ProjectionY("spectrum%s" %(trgclass))
                spectra[trgclass].Sumw2()
                spectra[trgclass].Scale(1./hlist["EventCounter"].GetBinContent(eventbin))
                NormaliseBinWidth(spectra[trgclass])
        return spectra

def CompareTriggers(fileMB, fileTrigger, args):
        try:
                opt,arg = getopt(args, "nspr", ["nocut","stdcut","with-pileup","without-piluep"])
        except GetoptError as e:
                print "Invalid option"
                print str(e)
                sys.exit(1)
        cuts = "stdcut"
        pileup = "wpr"

        for o,a in opt:
                if o in ("-n", "--nocut"):
                        cuts = "nocut"
                elif o in ("-s", "--stdcut"):
                        cuts = "stdcut"
                elif o in ("-p", "--with-pileup"):
                        pileup = "nopr"
                elif o in ("-r" "--without-piluep"):
                        pileup = "wpr"

        histos = dict()
        EMCALTriggers = ["EMCJHigh", "EMCJLow", "EMCGHigh", "EMCGLow"]
        try:
                histos["MinBias"] = ReadHistograms(fileMB, "MinBias", cuts, pileup)
                for trg in EMCALTriggers:
                        histos[trg] = ReadHistograms(fileTrigger, trg, cuts, pileup)
        except FileReaderException as e:
                print str(e)
                sys.exit(1)
        except HistNotFoundException as e:
                print str(e)
                sys.exit(1)
        spectra = NormaliseSpectra(histos, pileup)
        options = {"cuts":cuts, "pileup":"OFF"}
        if pileup == True:
                options["pileup"] = "ON"
        MakePlot(spectra, EMCALTriggers, options)

def runComparison(fileMB, fileTrigger, argstring):
        args = argstring.split(" ")
        CompareTriggers(fileMB, fileTrigger, args)

def main():
        if len(sys.argv) < 3:
                print "At least 2 arguments (file with MB spectrum, file with triggered spectra) required"
                sys.exit(1)
        fileMB = sys.argv[len(sys.argv)-2]
        fileTrigger = sys.argv[len(sys.argv)-1]
        arglist = sys.argv[1, len(sys.argv)-3]
        CompareTriggers(fileMB, fileTrigger, arglist)

if __name__ == "__main__":
        main()
