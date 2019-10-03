#**************************************************************************
#* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
#*                                                                        *
#* Author: The ALICE Off-line Project.                                    *
#* Contributors are mentioned in the code where appropriate.              *
#*                                                                        *
#* Permission to use, copy, modify and distribute this software and its   *
#* documentation strictly for non-commercial purposes is hereby granted   *
#* without fee, provided that the above copyright notice appears in all   *
#* copies and that both the copyright notice and this permission notice   *
#* appear in the supporting documentation. The authors make no claims     *
#* about the suitability of this software for any purpose. It is          *
#* provided "as is" without express or implied warranty.                  *
#**************************************************************************
from ROOT import kRed, kBlue, kBlack, kGreen
from copy import deepcopy
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.Helper import HistToGraph
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import Style, Frame, TwoPanelPlot, GraphicsObject
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.FileHandler import LegoTrainFileReader
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.SpectrumFitter import MinBiasFitter

class ComparisonPlot(TwoPanelPlot):
    def __init__(self, comparison):
        TwoPanelPlot.__init__(self)
        self.__comparison =  comparison
        self.__tagText = None
        
    def SetTag(self, tagText):
        self.__tagText = tagText
        
    def Create(self):
        self._CreateCanvas("comparison", "Comparison Data/Fit")

        spad = self.__canvas.cd(1)
        spad.GetPad().SetGrid(False, False)
        spad.GetPad().SetLogy(True)
        spad.GetPad().SetLogx(True)
        sframe = Frame("specframe", 0, 100, 1e-10, 100)
        sframe.SetXtitle("p_{t} (GeV/c)")
        sframe.SetYtitle("1/N_{event} 1/(#Delta p_{t}) dN/dp_{t} ((GeV/c)^{-2})")
        spad.DrawFrame(sframe)
        self.__comparison.DrawSpectra(spad)
        spad.CreateLegend()

        if self.__tagText:
            spad.DrawLabel(0.55, 0.8, 0.89, 0.88, self.__tagText)
        
        rpad = self._GetPad(2)
        rpad.GetPad().SetGrid(False, False)
        rframe = Frame("rframe", 0, 100, 0, 2)
        rframe.SetXtitle("p_{t} (GeV/c)")
        rframe.SetYtitle("Data/Patam")
        rpad.DrawFrame(rframe)
        self.__comparison.DrawRatio(rpad)
        self.__canvas.cd()
        
class MultipleFitPlot(TwoPanelPlot):
    
    def __init__(self):
        self.__data = {}
        self.__ratios = {}
        self.__reference = None
        
    def SetData(self, mycomp, fitmin, isReference):
        self.__data[fitmin] = mycomp
        if isReference:
            self.__reference = fitmin
            
    def RatiosToReference(self):
        for param in self.__data.keys():
            if param == self.__reference:
                continue
            self.__CalculateRatio(self.__data[param].GetRawParameterisation(), self.__data[self.__reference].GetRawParameterisation(), param)
    
    def Create(self):
        self._CreateCanvas("comparisonFitRange", "Comparison of the fit ranges")
        
        specpad = self.__canvas.cd(1)
        specpad.SetGrid(False, False)
        specpad.SetLogx(True)
        specpad.SetLogy(True)
        
        sframe = Frame("specframe", 0, 100, 1e-10, 100)
        sframe.SetXtitle("p_{t} (GeV/c)")
        sframe.SetYtitle("1/N_{event} 1/(#Delta p_{t}) dN/dp_{t} ((GeV/c)^{-2})")
        specpad.DraFrame(sframe)
        
        for param in sorted(self.__data.keys()):
            self.__data[param].DrawBinnedParameterisation(specpad, True, "%.1f GeV/c - 50 GeV/c" %(param))
        specpad.CreateLegend()     
        
        self.RatiosToReference()
        
        rpad = self._GetPad(2)
        rpad.GetPad().SetGrid(False, False)
        rframe = Frame("rframe", 0, 100, 0.5, 1.5)
        rframe.SetXtitle("p_{t} (GeV/c)")
        rframe.SetYtitle("Ratio to %.1f GeV/c - 50 GeV/c" %(self.__reference))
        rpad.DrawFrame(rframe)
        for ratio in sorted(self.__ratios.keys()):
            rpad.DrawGraphicsObject(self.__MakeRatioGraph(ratio, self.__data[ratio].GetStyle("Param"), self.__data[ratio].GetXrange()), False, "")
        self.__canvas.cd()
        
    def __MakeRatioGraph(self, param, style, myrange = None):
        xmin = None
        xmax = None
        if range:
            xmin = myrange["min"]
            xmax = myrange["max"]
        return GraphicsObject(HistToGraph(self.__ratios[param], xmin, xmax) , style)
          
    def __CalculateRatio(self, num, den, tag):
        self.__ratios[tag] = deepcopy(num)
        self.__ratios[tag].Divide(den)

class DataFitComparison:
    
    def __init__(self, data, fitmin):
        self.__data = data
        self.__mbfitter = MinBiasFitter(data, fitmin)
        self.__parameterised = self.__CreateBinnedParameterisation()
        self.__ratio = self.__CreateRatioDataParam()
        
        self.__xrange = {"min":None, "max":None}
        self.__graphs = {"Data" : None, "Param" : None, "Ratio" : None }
        self.__styles = {"Data" : Style(kRed, 24), "Param" : Style(kBlue, 25), "Ratio" : Style(kBlack, 20) }
        
    def SetStyle(self, what, style):
        if type in self.__styles.keys():
            self.__styles[what] = style
            
    def GetStyle(self, what):
        if not type in self.__styles.keys():
            return None
        return self.__styles[what]
            
    def SetRange(self, xmin = None, xmax = None):
        self.__xrange["min"] = xmin
        self.__xrange["max"] = xmax
        
    def GetXrange(self):
        return self.__xrange
    
    def __CreateBinnedParameterisation(self):
        print "Called"
        parameterised = deepcopy(self.__data)
        for mybin in range(1, parameterised.GetXaxis().GetNbins()+1):
            parameterised.SetBinContent(mybin, \
                            self.__mbfitter.CalculateBinned(parameterised.GetXaxis().GetBinLowEdge(mybin), parameterised.GetXaxis().GetBinUpEdge(mybin)))
            parameterised.SetBinError(mybin, 0)
        return parameterised
    
    def GetRawSpectrum(self):
        return self.__data
    
    def GetRawParameterisation(self):
        return self.__parameterised
            
    def DrawSpectra(self, pad, addToLegend = True):
        self.DrawSpectrum(pad)
        self.DrawBinnedParameterisation(pad)
        
    def DrawSpectrum(self, pad, addToLegend = True):
        if not self.__graphs["Data"]:
            self.__graphs["Data"] = self.__ConvertToGraph(self.__data)
        pad.DrawGraphicsObject(GraphicsObject(self.__graphs["Data"], self.__styles["Data"], addToLegend, "Data"))
    
    def DrawBinnedParameterisation(self, pad, addToLegend = True, legendText = "Param"):
        if not self.__graphs["Param"]:
            self.__graphs["Param"] = self.__ConvertToGraph(self.__parameterised)
        pad.DrawGraphicsObject(GraphicsObject(self.__graphs["Param"], self.__styles["Param"]), addToLegend, legendText)

    def DrawRatio(self, pad, addToLegend = False):
        self.__graphs["Ratio"] = self.__ConvertToGraph(self.__ratio)
        pad.DrawGraphicsObject(GraphicsObject(self.__DrawStyle(self.__graphs["Ratio"], self.__styles["Ratio"]), addToLegend, "Ratio Data/Param"))
                
    def __CreateRatioDataParam(self):
        ratio = deepcopy(self.__data)
        ratio.Divide(self.__parameterised)
        return ratio
    
    def __ConvertToGraph(self, hist):
        return HistToGraph(hist, self.__xrange["min"], self.__xrange["max"])


def MakeNormalisedSpectrum(inputdata, name):
    """
    Normalise spectrum by the number of events and by the bin width
    """
    inputdata.SetVertexRange(-10., 10.)
    inputdata.SetPileupRejection(True)
    inputdata.SelectTrackCuts(1)
    return inputdata.MakeProjection(0, "ptSpectrum%s" %(name), "p_{t} (GeV/c)", "1/N_{event} 1/(#Delta p_{t}) dN/dp_{t} ((GeV/c)^{-2})")

def CompareDataFit(filename, fitmin, fitmax):
    reader = LegoTrainFileReader(filename)
    mbspectrum = MakeNormalisedSpectrum(reader.ReadFile().GetData("MinBias").FindTrackContainer("tracksAll"), "MinBias")
    comparison = DataFitComparison(mbspectrum, fitmin)
    comparison.SetRange(2., 100.)
    comparisonPlot = ComparisonPlot(comparison)
    comparisonPlot.SetTag("Fit range: %.1f - %.1f Gev/c" %(fitmin, fitmax))
    comparisonPlot.Create()
    return comparisonPlot

def CheckFitRanges(filename):
    reader = LegoTrainFileReader(filename)
    mbspectrum = MakeNormalisedSpectrum(reader.ReadFile().GetData("MinBias").FindTrackContainer("tracksAll"), "MinBias")

    plot = MultipleFitPlot()
    styles = {10:Style(kBlue,24),15:Style(kBlack,25),20:Style(kRed,26),25:Style(kGreen,27)}
    for imin in range(10, 30, 5):
        comparison = DataFitComparison(mbspectrum, imin)
        comparison.SetRange(2., 100.)
        comparison.SetStyle("Param",styles[imin])
        isRef = False
        if imin == 15:
            isRef = True
        plot.SetData(comparison,imin,isRef)
    plot.Create()
    return plot