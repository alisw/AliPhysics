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
from PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import TwoPanelPlot, Style, GraphicsObject, Frame
from copy import deepcopy
from ROOT import kBlack

class Spectrum:
    """
    Graphical representation of the spectrum
    """
    
    def __init__(self, data, title, style = None):
        """
        Constructor
        """
        self.__data = data
        self.__title = title
        self.__style = None
        if style:
            self.__style = style
        else:
            self.__style = Style(kBlack, 20)
            
    def GetData(self):
        """
        Return Data of the spectrum
        """
        return self.__data
    
    def GetTitle(self):
        """
        Return title of the spectrum
        """
        return self.__title
    
    def SetData(self, histogram):
        """
        Change the data of the spectrum
        """
        self.__data = histogram
        
    def SetTitle(self, title):
        """
        Change title of the spectrum
        """
        self.__title = title
        
    def GetGraphics(self):
        """
        Return graphics object for plotting
        """
        return GraphicsObject(self.__data, self.__style)
    
    def RatioToOther(self, other, useBinomial = False):
        """
        Calculate ratio to other spectrum
        """
        result = deepcopy(self.__data)
        optstring = ""
        if useBinomial:
            optstring = "b"
        result.Divide(result, other.GetData(), 1., 1., optstring)
        return Spectrum(result, "", self.__style)


class PeriodComparisonData:
    """
    Rawspectrum data
    """
    
    class ComparisonObject:
        """
        Entry within the rawsepctrum data
        """
        
        def __init__(self, name, data, style, isReference):
            """
            Constructor
            """
            self.__data = Spectrum(data, name, style)
            self.__isReference = isReference
            
        def GetTitle(self):
            return self.__data.GetTitle()
        
        def GetData(self):
            return self.__data.GetData()
        
        def GetSpectrum(self):
            return self.__data
        
        def GetStyle(self):
            return self.__data.GetStyle()
        
        def IsReference(self):
            return self.__isReference
            
        def __cmp__(self, other):
            othername = other.GetTitle()
            if self.GetTitle() < othername:
                return -1
            if self.GetTitle() > othername:
                return +1
            return 0
            
    def __init__(self):
        """
        Constructor
        """
        self.__spectra = []
        self.__ratios = []
        
    def AddSpectrum(self, name, data, style, isReference):
        """
        Add raw sepctrum to the data collection
        """
        self.__spectra.append(self.ComparisonObject(name, data, style, isReference))
        
    def DrawSpectraInPad(self,pad):
        """
        Draw all raw spectra in the given pad
        """
        for spectrum in sorted(self.__spectra):
            pad.DrawGraphicsObject(spectrum.GetSpectrum().GetGraphics(), True, spectrum.GetTitle())
            
    def DrawRatiosInPad(self, pad):
        """
        Draw all ratios to the reference spectrum in the given pad
        """
        self.__CalculateRatios()
        for ratio in sorted(self.__ratios):
            pad.DrawGraphicsObject(ratio.GetGraphics(), True, ratio.GetTitle())
            
    def __FindReference(self):
        """
        Find the reference entry for the ratio calculation
        """
        reference = None
        for data in self.__spectra:
            if data.IsReference():
                reference = data
                break
        return reference
    
    def GetReferenceTitle(self):
        """
        GetTitle of the reference
        """
        reference = self.__FindReference()
        if not reference:
            return None
        return reference.GetTitle()
    
    def __CalculateRatios(self):
        """
        Calculate all ratios to the reference
        """
        reference = self.__FindReference()
        if not reference:
            return None
        for spec in self.__spectra:
            if spec == reference:
                continue
            ratio = spec.GetSpectrum().RatioToOther(reference.GetSpectrum())
            ratio.SetTitle(spec.GetTitle())
            self.__ratios.append(ratio)
                

class RawspectrumPeriodComparisonPlot(TwoPanelPlot):
    """
    Plot for the comparison of raw spectra in different periods
    """

    def __init__(self):
        """
        Constructor
        """
        TwoPanelPlot.__init__(self)
        self.__model = PeriodComparisonData()
        self.__label = None
        
    def AddRawSpectrum(self, title, data, style, isReference):
        """
        Add further raw spectrum to the plot
        """
        self.__model.AddSpectrum(title, data, style, isReference)
        
    def AddLabel(self, text):
        """
        Label plot
        """
        self.__label = text
        
    def Create(self):
        """
        Create the figure
        """
        self._CreateCanvas("Period comparison", "Comparison of different periods")
        
        spad = self._OpenPad(1)
        spad.GetPad().SetGrid(False, False)
        spad.GetPad().SetLogx(True)
        spad.GetPad().SetLogy(True)
        sframe = Frame("sframe", 1, 100, 1e-7, 100)
        sframe.SetXtitle("p_{t} (GeV/c)")
        sframe.SetYtitle("1/N_{events} 1/#delta p_{T} dN/dp_{t} ((GeV/c)^{-2})")
        spad.DrawFrame(sframe)
        self.__model.DrawSpectraInPad(spad)
        spad.CreateLegend(0.55, 0.75, 0.89, 0.89)
        if self.__label:
            spad.DrawLabel(0.15, 0.15, 0.35, 0.2, self.__label)
            
        rpad = self._OpenPad(2)
        rpad.GetPad().SetGrid(False, False)
        rpad.GetPad().SetLogx(True)
        rframe = Frame("rframe", 1, 100, 0, 3)
        rframe.SetXtitle("p_{t} (GeV/c)")
        rframe.SetYtitle("Ratio to %s" %(self.__model.GetReferenceTitle()))
        rpad.DrawFrame(rframe)
        self.__model.DrawRatiosInPad(spad)
        self._canvas.cd()
