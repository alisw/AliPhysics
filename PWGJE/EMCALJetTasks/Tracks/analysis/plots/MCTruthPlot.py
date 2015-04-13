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
"""
Comparison plot for spectra in different pt-hat bins

@author: Markus Fasel
"""

from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import SinglePanelPlot, Style, GraphicsObject, Frame
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.ComparisonData import ComparisonPlot, ComparisonData, ComparisonObject
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.SpectraSum import SpectraSum
from ROOT import kBlack

class MCSpectrumPtHatBin(ComparisonObject):
    """
    Entry class of a spectrum for a given pt-hat bin
    """
    
    def __init__(self, pthatbin, spectrum, style = None):
        """
        Constructor
        """
        ComparisonObject.__init__(self, spectrum, style)
        self.__pthatbin = pthatbin
        
    def GetLegendTitle(self):
        return "Pt-hat bin %d" %(self.__pthatbin)
    
    def GetObjectName(self):
        return "SpectrumPtHat%d" %(self.__pthatbin)     
    
class MCSpectrumContainer(ComparisonData):
    """
    Container class for spectra in different pt-hat bins
    """
    
    def __init__(self):
        """
        Constructor, initialising list of bins
        """
        ComparisonData.__init__(self)
    
    def AddPtHatBin(self, pthatbin, spectrum, style = None):
        """
        Add new pt-hat bin to the container
        """
        self.AddEntry(MCSpectrumPtHatBin(pthatbin, spectrum, style))
        
    def GetSpectraSum(self):
        """
        sum up the spectra in different pt-hard bins
        """
        summer = SpectraSum()
        for pthatbin in self.GetEntries():
            summer.AddSpectrum(pthatbin.GetData())
        return summer.GetSummedSpectrum()
        
    def DrawObjects(self, pad, addtolegend = True):
        """
        Draw all spectra inside the container into a given pad
        """
        ComparisonData.DrawObjects(self, pad, addtolegend)
        # draw also sum of the different bins
        pad.DrawGraphicsObject(GraphicsObject(self.GetSpectraSum(), Style(kBlack, 20)), addtolegend, "Sum")
    
class WeightedPtSpectrumFrame(Frame):
    
    def __init__(self):
        Frame.__init__(self, "sframe", 0., 100., 1e-20, 1e-5)
        self.SetXtitle("p_{t} (GeV/c)")
        self.SetYtitle("d#sigma/dp_{t} (mb/(GeV/c))" )
    
class WeightedEnergySpectrumFrame(Frame):
    
    def __init__(self):
        Frame.__init__(self, "eframe", 0., 100., 1e-20, 1e-5)
        self.SetXtitle("E (GeV)")
        self.SetYtitle("d#sigma/dE (mb/GeV)" )
    
class MCSpectrumPlot(ComparisonPlot):
    """
    Comparison plot of spectra for different pt-hat bins
    """

    def __init__(self, plottype = "tracks"):
        """
        Constructor
        """
        ComparisonPlot.__init__(self)
        self._comparisonContainer = MCSpectrumContainer()
        self._canvasname = ""
        self._canvastitle = ""
        self.__labeltext = "MC-true spectrum"
        if plottype == "tracks":
            self.SetFrame(WeightedPtSpectrumFrame())
        else:
            self.SetFrame(WeightedEnergySpectrumFrame())
        
    def SetLabelText(self, text):
        """
        Change text of the label
        """
        self.__labeltext = text
        
    def AddMCSpectrum(self, pthatbin, spectrum, style = None):
        """
        Add new spectrum in pt-hat bin to the plot
        """
        self._comparisonContainer.AddPtHatBin(pthatbin, spectrum, style)
        
    def Create(self):
        """
        Create the plot
        """
        self.SetPadAttributes(True, True, False, False)
        self.SetLegendAttributes(0.7, 0.5, 0.89, 0.89)
        self._Create(self._canvasname, self._canvastitle)
        pad = self._GetFramedPad()
        pad.DrawLabel(0.15, 0.15, 0.45, 0.21, self.__labeltext)
        
class MCTrueSpectrumPlot(MCSpectrumPlot):
    
    def __init__(self, plottype = "tracks"):
        MCSpectrumPlot.__init__(self, plottype)
        self._canvasname = "MCtruthPlot"
        self._canvastitle = "Plot of MC-true spectra"
        self.SetLabelText("MC-true spectrum")
    
class MCRecSpectrumPlot(MCSpectrumPlot):
    
    def __init__(self, triggername, plottype = "tracks"):
        MCSpectrumPlot.__init__(self)
        self._canvasname = "MCrecPlot%s" %(triggername)
        self._canvastitle = "Plot of MC-reconstructed spectra for trigger %s" %(triggername)
        self.SetLabelText("MC-reconstructed spectrum for trigger %s" %(triggername))
    
class MCWeightPlot(SinglePanelPlot):
    """
    Class for the plot of the weights for different pt-hard bins
    """
    
    def __init__(self, weights):
        """
        Constructor
        """
        SinglePanelPlot.__init__(self)
        self.__points = weights
        
    def Create(self):
        """
        Creator function for the plot
        """
        self._OpenCanvas("weightplot", "Monte-Carlo weights")
        pad = self._GetFramedPad()
        pad.GetPad().SetLogy()
        frame = Frame("wframe", 0., 11., 1e-12, 1e-5)
        frame.SetXtitle("p_{t,hard} bin")
        frame.SetYtitle("weight factor")
        pad.DrawFrame(frame)
        pad.DrawGraphicsObject(GraphicsObject(self.__points.GetWeightingCurve(), Style(kBlack, 20)), False, "weights")