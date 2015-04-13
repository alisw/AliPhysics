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
Comparison plot of trigger efficiencies in MC in different pt-hat bins including underlying data structure

@author: Markus Fasel
"""

from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import SinglePanelPlot, GraphicsObject, Frame, Style
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.ComparisonData import ComparisonData, ComparisonObject, ComparisonPlot
from ROOT import TFile,kBlack
    
class TriggerEfficiencyClassPtHat(ComparisonObject):
    
    def __init__(self, pthatbin, triggerdata, style):
        ComparisonObject.__init__(self, triggerdata, style)
        self.__pthatbin = pthatbin
        
    def GetLegendTitle(self):
        return "p_{t}-hat bin %d" %(self.__pthatbin)
    
    def GetObjectName(self):
        return "pthat%d" %(self.__pthatbin)
    
class TriggerEfficiencyClassTriggerType(ComparisonObject):
    
    def __init__(self, triggername, triggerdata, style):
        ComparisonObject.__init__(self, triggerdata, style)
        self.__triggername =  triggername
        
    def GetLegendTitle(self):
        return self.__triggername
    
    def GetObjectName(self):
        return self.__triggername

class TriggerEfficiencyContainer(ComparisonData):
    """
    Underlying data structure for the comparison plot
    """
    
    def __init__(self):
        """
        Initialise container
        """
        ComparisonData.__init__(self)
    
    def AddEfficiency(self, trclasstype, key, efficiencyCurve, style):
        """
        Add new trigger ifno
        """
        triggerdata = None
        if trclasstype == "pthat":
            triggerdata = TriggerEfficiencyClassPtHat(key, efficiencyCurve, style)
        elif trclasstype == "triggertype":
            triggerdata = TriggerEfficiencyClassTriggerType(key, efficiencyCurve, style)
        self.AddEntry(triggerdata)
                        
class TriggerEfficiencyFrame(Frame):
    """
    Frame class for trigger efficiency plots
    """
    
    def __init__(self, name):
        """
        Constructor
        """
        Frame.__init__(self, name, 0., 100., 0., 1.)
        self.SetXtitle("p_{t} (GeV/c)")
        self.SetYtitle("Trigger efficiency")
        
class TriggerEfficiencyPlotMC(ComparisonPlot):
    """
    Comparison plot of trigger efficiencies in different pt-hat bins
    """

    def __init__(self):
        """
        Constructor
        """
        ComparisonPlot.__init__(self)
        self._comparisonContainer = TriggerEfficiencyContainer()
        self.SetFrame(TriggerEfficiencyFrame("tframe"))
        self.SetLegendAttributes(0.65, 0.15, 0.89, 0.5)
        self.__triggername = ""
        
    def SetTriggerName(self, trname):
        """
        Set triggername for the label
        """
        self.__triggername = trname
    
    def AddEfficiency(self, pthatbin, efficiency, style):
        """
        Add new efficiency container to the data structure
        """
        self._comparisonContainer.AddEfficiency("pthat", pthatbin, efficiency, style)
        
    def Create(self):
        """
        Create the plot
        """
        self._Create("triggerEfficiencyMC", "MC trigger efficiency plot")
        if len(self.__triggername):
            pad = self._GetFramedPad()
            pad.DrawLabel(0.15, 0.8, 0.5, 0.85, self.__triggername)
        
        
class TriggerEfficiencyPlotClasses(ComparisonPlot):
    """
    Plot comparing the trigger efficiency of different trigger types
    """
    
    def __init__(self):
        """
        Constructor
        """
        ComparisonPlot.__init__(self)
        self._comparisonContainer = TriggerEfficiencyContainer()
        self.SetFrame(TriggerEfficiencyFrame("tframe"))
        self.SetLegendAttributes(0.65, 0.15, 0.89, 0.5)
    
    def AddTriggerEfficiency(self, triggername, efficiency, style):
        """
        Add trigger class to the comparison data
        """
        self._comparisonContainer.AddEfficiency("triggertype", triggername, efficiency, style)
        
    def Create(self):
        self._Create("triggerclasses", "Trigger efficiencies")
            
        
class TriggerEfficiencySumPlot(SinglePanelPlot):
    """
    Plot the summed trigger efficiency from different pt-hard bins
    """
    
    def __init__(self, triggername, triggerefficiency):
        """
        Constructor
        """
        SinglePanelPlot.__init__(self)
        self.__triggername = triggername
        self.__triggereff = triggerefficiency
        
    def Create(self):
        """
        Create the plot
        """
        self._OpenCanvas("trgEffSumm", "Summed trigger efficiency")
        pad = self._GetFramedPad()
        pad.DrawFrame(TriggerEfficiencyFrame("tframe"))
        pad.DrawGraphicsObject(GraphicsObject(self.__triggereff.GetEfficiencyCurve(), Style(kBlack, 20)), False, "Trigger Eff")
        pad.DrawLabel(0.5, 0.2, 0.89, 0.25, "Trigger: %s" %(self.__triggername))