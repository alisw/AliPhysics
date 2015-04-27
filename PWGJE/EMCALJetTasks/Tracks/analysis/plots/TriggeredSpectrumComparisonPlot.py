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
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.ComparisonData import ComparisonData,ComparisonObject,ComparisonPlot
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import Frame

class PtSpectrumFrame(Frame):
    
    def __init__(self, name):
        Frame.__init__(self, "ptframe%s" %(name), 0., 100., 1e-10, 100)
        self.SetXtitle("p_{t} (GeV/c)")
        self.SetYtitle("1/N_{event} dN/dp_{t} ((GeV/c)^{-1})")
    
class EnergySpectrumFrame(Frame):
    
    def __init__(self, name):
        Frame.__init__(self, "ptframe%s" %(name), 0., 100., 1e-10, 100)
        self.SetXtitle("p_{t} (GeV/c)")
        self.SetYtitle("1/N_{event} dN/dp_{t} ((GeV/c)^{-1})")
    
class SpectraComparisonObject(ComparisonObject):
    
    def __init__(self, trigger, data, style):
        ComparisonObject.__init__(self, data, style)
        self.__triggername =  trigger
        
    def GetLegendTitle(self):
        return self.__triggername
    
    def GetObjectName(self):
        return "Rawspectrum%s" %(self.__triggername)

class TriggeredSpectrumComparisonPlot(ComparisonPlot):
    """
    Comparing raw spectra of different classes
    """

    def __init__(self, frame, name = "spectrumcomparison"):
        """
        Constructor
        """
        ComparisonPlot.__init__(self)
        self._comparisonContainer = ComparisonData()
        self.SetFrame(frame)
        self.SetLegendAttributes(0.5, 0.65, 0.89, 0.89)
        self.SetPadAttributes(True, True, False, False)
        self.__name = name
        
    def AddSpectrum(self, trigger, spectrum, style):
        self._comparisonContainer.AddEntry(SpectraComparisonObject(trigger, spectrum, style))
        
    def Create(self):
        self._Create("canvas%s" %(self.__name), "Spectrum Comparison %s" %(self.__name))
        
class PtTriggeredSpectrumComparisonPlot(TriggeredSpectrumComparisonPlot):
    
    def __init__(self,name):
        TriggeredSpectrumComparisonPlot.__init__(self, PtSpectrumFrame(name), name)
        
class EnergyTriggeredSpectrumComparisonPlot(TriggeredSpectrumComparisonPlot):
    
    def __init__(self, name):
        TriggeredSpectrumComparisonPlot.__init__(self, EnergySpectrumFrame(name), name)
        