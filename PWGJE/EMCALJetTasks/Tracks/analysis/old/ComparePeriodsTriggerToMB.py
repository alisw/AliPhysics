#! /usr/bin/env python

from ROOT import TCanvas, TGraphErrors, TLegend, TPaveText
from ROOT import kBlack, kBlue, kRed
from Helper import Frame, ReadHistList
from Graphics import Style
from SpectrumContainer import DataContainer
from copy import deepcopy

class PeriodComparisonPlot:
    
    def __init__(self):
        self.__comparisons = []
        
        self.__canvas = None
        self.__frames = {}
        self.__legend = None
        
    def AddComparison(self, comp):
        self.__comparisons.append(comp)
        
    def SetPlotRange(self, min ,max):
        for comp in self.__comparisons:
            comp.SetPlotRange(min, max)
        
    def Draw(self):
        self.__canvas = TCanvas("comparison%s" %(self.__comparisons[0].GetTriggerName()), "Comparison of different periods for trigger %s" %(self.__comparisons[0].GetTriggerName()), 1000, 600)
        self.__canvas.Divide(2,1)
        
        self.__legend = TLegend(0.15, 0.15, 0.45, 0.45)
        self.__legend.SetBorderSize(0)
        self.__legend.SetFillStyle(0)
        self.__legend.SetTextFont(42)
        
        specpad = self.__canvas.cd(1)
        specpad.SetGrid(False,False)
        specpad.SetLogx(True)
        specpad.SetLogy(True)
        self.__frames["Spectra"] = Frame("axisSpec%s" %(self.__comparisons[0].GetTriggerName()), 0, 100, 1e-10, 100)
        self.__frames["Spectra"].SetXtitle("p_{t} (GeV/c)")
        self.__frames["Spectra"].SetYtitle("1/N_{event} 1/(#Delta p_{t}) dN/dp_{t} ((GeV/c)^{-2})")
        self.__frames["Spectra"].Draw()
        
        self.__comparisons[0].DrawMinBiasSpectrum()
        self.__comparisons[0].AddMBtoLegend(self.__legend)
        for comp in sorted(self.__comparisons):
            comp.DrawTriggeredSpectrum()
            comp.AddTriggeredSpectrumToLegend(self.__legend)
        self.__legend.Draw()
        self.__label = self.__comparisons[0].CreateLabel(0.5, 0.75, 0.89, 0.85)
        self.__label.Draw()
        
        rpad = self.__canvas.cd(2)
        rpad.SetGrid(False, False)
        self.__frames["Ratios"] = Frame("axisRatio%s" %(self.__comparisons[0].GetTriggerName()), 0, 100, 0, 2000)
        self.__frames["Ratios"].SetXtitle("p_{t} (GeV/c)")
        self.__frames["Ratios"].SetYtitle("%s / Min. Bias" %(self.__comparisons[0].GetTriggerName()))
        self.__frames["Ratios"].Draw()
        for comp in sorted(self.__comparisons):
            comp.DrawRatioTriggeredMinBias()
        
        self.__canvas.cd()

    def SaveAs(self, filenamebase):
        """
        Save plot as image file
        """
        types = ["eps", "pdf", "jpeg", "gif", "png"]
        for t in types:
            self.__canvas.SaveAs("%s.%s" %(filenamebase, t))

class TriggerComparison:
    
    def __init__(self, trgspec, mbspec, triggername, dataname):
        self.__triggeredspectrum = trgspec
        self.__minbiasspectrum = mbspec
        self.__ratiospectra = self.__triggeredspectrum.MakeRatio(self.__minbiasspectrum)
        self.__ratiospectra.SetStyle(self.__triggeredspectrum.GetStyle())
        self.__triggername = triggername
        self.__dataname = dataname
        
    def __cmp__(self, other):
        othername = other.GetDataName()
        if self.__dataname == othername:
            return 0
        elif self.__dataname < othername:
            return -1
        else:
             return 1
        
    def SetPlotRange(self, min, max):
        self.__triggeredspectrum.SetPlotRange(min, max)
        self.__minbiasspectrum.SetPlotRange(min, max)
        self.__ratiospectra.SetPlotRange(min, max)
        
    def GetTriggerName(self):
        return self.__triggername
    
    def GetDataName(self):
        return self.__dataname

    def DrawTriggeredSpectrum(self):
        self.__triggeredspectrum.Draw()
        
    def DrawMinBiasSpectrum(self):
        self.__minbiasspectrum.Draw()
        
    def DrawRatioTriggeredMinBias(self):
        self.__ratiospectra.Draw()
        
    def AddMBtoLegend(self, leg):
        self.__minbiasspectrum.AddToLegend(leg, "MinBias")
        
    def AddTriggeredSpectrumToLegend(self, leg):
        self.__triggeredspectrum.AddToLegend(leg, self.__dataname)
        
    def CreateLabel(self, xmin, ymin, xmax, ymax):
        label = TPaveText(xmin, ymin, xmax, ymax, "NDC")
        label.SetBorderSize(0)
        label.SetFillStyle(0)
        label.SetTextFont(42)
        label.AddText("Trigger: %s" %(self.__triggername))
        return label
        
class GraphicsObject:
    
    def __init__(self, data, name):
        self._data = data
        self._graphics = None
        self._style = Style(kBlack, 20)
        self._plotrange = {"Min":None, "Max":None}
        self._name = name
        
    def SetPlotRange(self, min, max):
        self._plotrange[min] = min
        self._plotrange[max] = max
        
    def SetStyle(self, style):
        self._style = style
        
    def SetName(self, name):
        self._name = name
    
    def GetData(self):
        return self._data
        
    def GetGraphics(self):
        return self._graphics
    
    def GetStyle(self):
        return self._style
    
    def Draw(self):
        if not self._graphics:
            self._graphics = TGraphErrors()
            np = 0
            for bin in range(1, self._data.GetXaxis().GetNbins()+1):
                if self._plotrange["Min"] and self._data.GetXaxis().GetBinLowEdge(bin) < self._plotrange["Min"]:
                    continue
                if self._plotrange["Max"] and self._data.GetXaxis().GetBinUpEdge(bin) > self._plotrange["Max"]:
                    break
                self._graphics.SetPoint(np, self._data.GetXaxis().GetBinCenter(bin), self._data.GetBinContent(bin))
                self._graphics.SetPointError(np, self._data.GetXaxis().GetBinWidth(bin)/2., self._data.GetBinError(bin))
                np = np + 1
        self._graphics.SetMarkerColor(self._style.GetColor())
        self._graphics.SetLineColor(self._style.GetColor())
        self._graphics.SetMarkerStyle(self._style.GetMarker())
        self._graphics.Draw("epsame")  
        
    def AddToLegend(self, legend, title = None):
        if self._graphics:
            tit = self._name
            if title:
                tit = title
            legend.AddEntry(self._graphics, tit, "lep")  
    
class Spectrum(GraphicsObject):
    
    def __init__(self, data, name):
        GraphicsObject.__init__(self, data, name)

    def MakeRatio(self, denominator):
        result = deepcopy(self._data)
        result.Divide(denominator.GetData())
        ratio = Ratio(result)
        if self._plotrange["Min"] or self._plotrange["Max"]:
            ratio.SetPlotRange(self._plotrange["Min"], self._plotrange["Max"])
        return ratio
        
class Ratio(GraphicsObject):
    def __init__(self, data, name = None):
        GraphicsObject.__init__(self, data, name)
        
def ReadSpectra(filename, trigger):
    """
    Read the spectra for different trigger classes from the root file
    Returns a dictionary of triggers - spectrum container
    """
    hlist = ReadHistList(filename, "PtEMCalTriggerTask")
    return DataContainer(eventHist = hlist.FindObject("hEventHist%s" %(trigger)), trackHist = hlist.FindObject("hTrackHist%s" %(trigger)))

def MakeNormalisedSpectrum(inputdata, name):
    """
    Normalise spectrum by the number of events and by the bin width
    """
    inputdata.SetVertexRange(-10., 10.)
    inputdata.SetPileupRejection(True)
    inputdata.SelectTrackCuts(1)
    return inputdata.MakeProjection(0, "ptSpectrum%s" %(name), "p_{t} (GeV/c)", "1/N_{event} 1/(#Delta p_{t}) dN/dp_{t} ((GeV/c)^{-2})")

def ComparePeriods(filea, fileb, filemb, namea, nameb, trigger):
    triggers = {}
    dataA = ReadSpectra(filea, trigger)
    dataB = ReadSpectra(fileb, trigger)
    dataMB = ReadSpectra(filemb, "MinBias")
    
    specA = Spectrum(MakeNormalisedSpectrum(dataA, namea), namea)
    specA.SetStyle(Style(kBlue, 24))
    specB = Spectrum(MakeNormalisedSpectrum(dataB, nameb), nameb)
    specB.SetStyle(Style(kRed, 25))
    specMB = Spectrum(MakeNormalisedSpectrum(dataMB, "MinBias"), "MinBias")
    specMB.SetStyle(Style(kBlack, 25))
    
    plot = PeriodComparisonPlot()
    plot.AddComparison(TriggerComparison(specA, specMB, trigger, namea))
    plot.AddComparison(TriggerComparison(specB, specMB, trigger, nameb))
    plot.SetPlotRange(2., 100.)
    plot.Draw()
    return plot
    
