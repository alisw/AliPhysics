#! /usr/bin/env python
import sys
from copy import deepcopy
from getopt import getopt,GetoptError
from ROOT import TCanvas,TFile,TH1F,TLegend,TPaveText
from ROOT import gROOT,kRed,kBlue,kBlack
from base.Helper import FileReaderException,HistNotFoundException,ReadHistList
from base.Graphics import Style
from base.struct.DataContainers import DataContainer

gObjects = list()

class ComparisonPlot:
    """
    Comparison plot of spectra from 2 periods
    """
    def __init__(self, spectraa, spectrab, plotname = None):
        """
        Constructor, basic initialisation for the plot
        """
        # Data
        self.__spectraA = spectraa
        self.__spectraB = spectrab
        self.__ratio = self.MakeRatios(self.__spectraA, self.__spectraB)
        self.__options = None
        
        # Graphics objects
        self.__plotname = "comparisonPlot"
        if plotname:
            self.__plotname = plotname
        self.__plotTitle = "Comparison of different periods"
        self.__comparisonPlot = None
        self.__axisspec = {}
        self.__axisrat = {}
        self.__legend = None
        self.__label = {}
        
        # Titles and ranges (default, changable) 
        tmpkey = self.__spectraA.keys()[0]
        self.__xtitle = self.__spectraA[tmpkey].GetHistogram().GetXaxis().GetTitle()
        self.__ytitle = self.__spectraA[tmpkey].GetHistogram().GetYaxis().GetTitle()
        self.__xrange = [self.__spectraA[tmpkey].GetHistogram().GetXaxis().GetXmin(), self.__spectraA[tmpkey].GetHistogram().GetXaxis().GetXmax()]
        self.__yrange = [min(self.__spectraA[tmpkey].GetHistogram().GetMinimum(), self.__spectraB[tmpkey].GetHistogram().GetMinimum()), \
                          max(self.__spectraA[tmpkey].GetHistogram().GetMaximum(), self.__spectraB[tmpkey].GetHistogram().GetMaximum())]
        self.__ratioyrange = [self.__ratio[tmpkey].GetRatioHist().GetMinimum(), self.__ratio[tmpkey].GetRatioHist().GetMaximum()]
        
    def MakeRatios(self, specnum, specden):
        result = {}
        for k in specnum.keys():
            result[k] = Ratio(specnum[k], specden[k])
        return result
        
    def SetXTitle(self, title):
        """
        Set title of the x-coordinate
        """
        self.__xtitle = title
        
    def SetYTitle(self, title):
        """
        Set title of the y-coordinate
        """
        self.__ytitle = title
        
    def SetXRange(self, xmin, xmax):
        """ 
        Set the plotting range in x-direction
        """
        self.__xrange[0] = xmin
        self.__xrange[1] = xmax
        
    def SetYRange(self, ymin, ymax):
        """
        Set the plotting range of the spectrum panel in y-direction
        """
        self.__yrange[0] = ymin
        self.__yrange[1] = ymax
        
    def SetYRangeRatio(self, ymin, ymax):
        """
        Set the plotting range of the spectrum panel in y-direction
        """
        self.__ratioyrange[0] = ymin
        self.__ratioyrange[1] = ymax
        
    def SetPlotname(self, plotname):
        """
        Name of the comparison plot
        """
        self.__plotname = plotname
        
    def SetPlottitle(self, plottitle):
        """
        Title of the comparison plot
        """
        self.__plottitle = plottitle
        
    def SaveAs(self, filenamebase):
        """
        Save plot as image file
        """
        types = ["eps", "pdf", "jpeg", "gif", "png"]
        for t in types:
            self.__comparisonPlot.SaveAs("%s.%s" %(filenamebase, t))

            
    def MakePlot(self):
        """
        Produce comparison plot
        Left panel: 2 spectra
        Right panel: ratio of the spectra
        """
        self.__comparisonPlot = TCanvas(self.__plotname, self.__plotTitle, 600, 300 *len(self.__spectraA))
        self.__comparisonPlot.Divide(2,len(self.__spectraA))
        row = 0
        for trg in self.__spectraA.keys():
            print "Doing trigger %s" %(trg)
            self.__DrawTrigger(trg, row)
            row = row + 1
        self.__comparisonPlot.cd()
        
    def __DrawTrigger(self, trg, row):
        
        padspec = self.__comparisonPlot.cd(row * 2 + 1)
        padspec.SetGrid(False, False)
        padspec.SetLogy()
        padspec.SetLogx()
        drawlegend = False
        if not self.__legend:
            self.__legend = self.__CreateLegend(0.15, 0.15, 0.35, 0.25)
            drawlegend = True
        self.__axisspec[trg] = TH1F("axisspec%s" %(self.__plotname), ";%s;%s" %(self.__xtitle, self.__ytitle), 1000, self.__xrange[0], self.__xrange[1])
        self.__axisspec[trg].SetStats(False)
        self.__axisspec[trg].GetYaxis().SetRangeUser(self.__yrange[0], self.__yrange[1])
        self.__axisspec[trg].Draw("axis")
        self.__spectraA[trg].Draw()
        self.__spectraB[trg].Draw()
        if drawlegend:
            self.__AddToLegend(self.__spectraA[trg])
            self.__AddToLegend(self.__spectraB[trg])
            self.__legend.Draw()
        
        self.__label[trg] = self.__MakeLabel(0.5,0.7,0.89,0.85,trg)
        self.__label[trg].Draw()
        
        padratio = self.__comparisonPlot.cd(row * 2 + 2)
        padratio.SetGrid(False, False)
        padratio.SetLogx()
        self.__axisrat[trg] = TH1F("axirat%s" %(self.__plotname), ";%s;%s" %(self.__xtitle, self.__ratio[trg].GetRatioTitle()), 1000, self.__xrange[0], self.__xrange[1])
        self.__axisrat[trg].SetStats(False)
        self.__axisrat[trg].GetYaxis().SetRangeUser(self.__ratioyrange[0], self.__ratioyrange[1])
        self.__axisrat[trg].Draw("axis")
        self.__ratio[trg].Draw()
    
    def __CreateLegend(self, xmin, ymin, xmax, ymax):
        """
        Create a new legend within range defined by coordinates
        """
        leg = TLegend(xmin, ymin, xmax, ymax)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        return leg
    
    def __AddToLegend(self, spectrum):
        """
        Add spectrum object to legend
        """
        self.__legend.AddEntry(spectrum.GetHistogram(), spectrum.GetTitle(), "lep")
        
    def __MakeLabel(self, xmin, xmax, ymin, ymax, trigger):
        """
        Add label with trigger, cuts and pileup
        """
        lab = TPaveText(xmin, xmax, ymin, ymax, "NDC")
        lab.SetBorderSize(0)
        lab.SetFillStyle(0)
        lab.SetTextFont(42)
        lab.SetTextAlign(12)
        lab.AddText("Trigger: %s" %(trigger))
        return lab

    
class Spectrum:
    """
    Helper class combining data of a graphics with a title and style information
    """
    def __init__(self, histogram, title, style = None):
        self.__histogram = histogram
        self.__title = title
        if style:
            self.SetStyle(style)
        else:
            self.SetStyle(Style(kBlack, 20))
        
    def GetHistogram(self):
        """
        Return Data of the spectrum
        """
        return self.__histogram
    
    def GetTitle(self):
        """
        Return title of the spectrum
        """
        return self.__title
    
    def SetHistogram(self, histogram):
        """
        Change the data of the spectrum
        """
        self.__histogram = histogram
        
    def SetTitle(self, title):
        """
        Change title of the spectrum
        """
        self.__title = title
    
    def SetStyle(self, style):
        """
        change style of the histogram
        """
        self.__histogram.SetMarkerColor(style.GetColor())
        self.__histogram.SetLineColor(style.GetColor())
        self.__histogram.SetMarkerStyle(style.GetMarker())
    
    def Draw(self):
        """
        Draw spectrum to a given Canvas
        """
        self.__histogram.Draw("epsame")

class Ratio:
    """
    Helper class creating ratio plots
    """
    def __init__(self, specnum, specden, useBinomial = False, style = None):
        self.__ratio = Spectrum(self.__calculateRatio(specnum, specden, useBinomial), "%s/%s" %(specnum.GetTitle(), specden.GetTitle()), style)
    
    def GetRatio(self):
        """
        Get the ratio result
        """
        return self.__ratio
    
    def GetRatioHist(self):
        """
        Return ratio underlying data
        """
        return self.__ratio.GetHistogram()
    
    def GetRatioTitle(self):
        """
        Return title of the ratio
        """
        return self.__ratio.GetTitle()
    
    def SetRatioTitle(self, title):
        """
        Change title of the ratio
        """
        self.__ratio.SetTitle(title)
    
    def SetStyle(self, style):
        """
        Change plotting style of the ratio hist
        """
        self.__ratio.SetStyle(style) 
        
    def Draw(self):
        """
        Draw ratio on a given canvas
        """
        self.__ratio.Draw()
        
    def __calculateRatio(self, specnum, specden, useBinomial = False):
        """ 
        Perform ratio calculation
        """
        result = deepcopy(specnum.GetHistogram())
        optstring = ""
        if useBinomial:
            optstring = "b"
        result.Divide(result, specden.GetHistogram(), 1., 1., optstring)
        return result
    
    
def MakeNormalisedSpectrum(inputdata, name):
    """
    Normalise spectrum by the number of events and by the bin width
    """
    inputdata.SetVertexRange(-10., 10.)
    inputdata.SetPileupRejection(True)
    inputdata.SelectTrackCuts(1)
    return inputdata.MakeProjection(0, "ptSpectrum%s" %(name), "p_{t} (GeV/c)", "1/N_{event} 1/(#Delta p_{t}) dN/dp_{t} ((GeV/c)^{-2})")

def ReadSpectra(filename, triggers):
    """
    Read the spectra for different trigger classes from the root file
    Returns a dictionary of triggers - spectrum container
    """
    hlist = ReadHistList(filename, "PtEMCalTriggerTask")
    result = {}
    for trg in triggers:
        result[trg] = DataContainer(eventHist = hlist.FindObject("hEventHist%s" %(trg)), trackHist = hlist.FindObject("hTrackHist%s" %(trg)))
    return result

def ComparePeriods(filea, fileb, nnum, nden, useClass = None):
    triggers = None
    plotname = None
    plottitle = None
    if not useClass:
        triggers = ["MinBias", "EMCJHigh", "EMCJLow", "EMCGHigh", "EMCGLow"]
    else:
        triggers = [useClass]
        plotname = "comparisonPlot%s"  %(useClass)
        plottitle = "Comparison of different periods for %s events" %(useClass)
    dataA = ReadSpectra(filea, triggers)
    dataB = ReadSpectra(fileb, triggers)
    
    spectraA = {}
    spectraB = {}
    for trg in triggers:
        spectraA[trg] = Spectrum(MakeNormalisedSpectrum(dataA[trg], trg), nnum, Style(kBlue, 24))
        spectraB[trg] = Spectrum(MakeNormalisedSpectrum(dataB[trg], trg), nden, Style(kRed, 25))

    resultplot = ComparisonPlot(spectraA, spectraB, plotname)
    if plottitle:
        resultplot.SetPlottitle(plottitle)
    resultplot.SetXTitle("p_{t} (GeV/c)")
    resultplot.SetYTitle("1/N_{events} 1/#delta p_{T} dN/dp_{t} ((GeV/c)^{-2})")
    resultplot.SetXRange(1., 100.)
    resultplot.SetYRange(1e-7, 100)
    resultplot.SetYRangeRatio(0., 3)
    resultplot.MakePlot()
    gObjects.append(resultplot)
    return resultplot

def runComparison(filea, fileb, argstring):
    """
    Run program in interactive mode (with ipython)
    """
    arglist = argstring.split(" ")
    ComparePeriods(filea, fileb, arglist)

def main():
    """
    Run program in batch mode
    """
    if len(sys.argv < 3):
        print "At least 2 arguments required"
    pass

if __name__ == "__main__":
    main()