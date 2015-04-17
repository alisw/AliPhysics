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
from PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import GraphicsObject,SinglePanelPlot
from ROOT import TFile

class ComparisonObject(object):
    """
    Base entry type for object inside comparison data
    """
    
    def __init__(self, data, style):
        self.__data = data
        self.__style = style
        
    def GetData(self):
        return self.__data
        
    def GetGraphicsObject(self):
        return GraphicsObject(self.__data, self.__style)
    
    def GetRootPrimitive(self):
        self.__data.SetName(self.GetObjectName())
        return self.__data
    
    def Draw(self, pad, addToLegend = True):
        pad.DrawGraphicsObject(self.GetGraphicsObject(), addToLegend, self.GetLegendTitle())
    
    def GetLegendTitle(self):
        """
        To be implemented in inheriting classes
        """
        return ""
    
    def GetObjectName(self):
        """
        To be implemented in inheriting classes
        """
        return ""
    
class ComparisonData(object):
    """
    General comparison data collection
    """

    def __init__(self):
        """
        Constructor
        """
        self.__entries = []
        
    def GetEntries(self):
        return self.__entries
        
    def AddEntry(self, entry):
        self.__entries.append(entry)
        
    def DrawObjects(self, pad, addToLegend = True):
        for entry in self.__entries:
            entry.Draw(pad, addToLegend)
            
    def GetListOfRootObjects(self):
        """
        Get a list of root-primitive trigger efficiencies
        """
        rootprimitives = []
        for entry in self.__entries:
            rootprimitives.append(entry.GetRootPrimitive())
        return rootprimitives
    
class ComparisonPlot(SinglePanelPlot):
    """
    General comparison plot type
    """
        
    def __init__(self):
        """
        Constructor
        """
        SinglePanelPlot.__init__(self)
        self.__frame = None
        self._comparisonContainer = None    # be specified in inheriting classes
        self.__legendAttributes = None
        self.__padattributes = {"logx":False, "logy":False, "gridx":False, "gridy":False}
        
    def SetFrame(self, frame):
        self.__frame = frame
        
    def SetLegendAttributes(self, xmin, ymin, xmax, ymax):
        self.__legendAttributes = {"xmin":xmin, "xmax":xmax, "ymin":ymin, "ymax":ymax}
        
    def SetPadAttributes(self, logx, logy, gridx, gridy):
        self.__padattributes["logx"] = logx
        self.__padattributes["logy"] = logy
        self.__padattributes["gridx"] = gridx
        self.__padattributes["gridy"] = gridy
        
    def _Create(self, canvasname, canvastitle):
        """
        Make the plot
        """
        self._OpenCanvas(canvasname, canvastitle)
        pad = self._GetFramedPad()
        if self.__padattributes["logx"]:
            pad.GetPad().SetLogx()
        if self.__padattributes["logy"]:
            pad.GetPad().SetLogy()
        pad.DrawFrame(self.__frame)
        doLegend = False
        if self.__legendAttributes:
            doLegend = True
        self._comparisonContainer.DrawObjects(pad, doLegend)
        if doLegend:
            pad.CreateLegend(self.__legendAttributes["xmin"], self.__legendAttributes["ymin"], self.__legendAttributes["xmax"], self.__legendAttributes["ymax"])
        
    def WriteData(self, rootfilename):
        """
        Write out trigger efficiency curves to a root file
        """
        outputfile = TFile(rootfilename, "RECREATE")
        for rootprim in self._comparisonContainer.GetListOfRootObjects():
            rootprim.Write()
        outputfile.Close()