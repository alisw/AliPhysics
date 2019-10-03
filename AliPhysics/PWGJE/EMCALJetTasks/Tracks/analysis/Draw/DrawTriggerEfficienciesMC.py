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
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.MonteCarloFileHandler import MonteCarloFileHandler
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import Style
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.FileHandler import ResultDataBuilder
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.plots.TriggerEfficiencyPlotMC import TriggerEfficiencyPlotMC,TriggerEfficiencySumPlot,TriggerEfficiencyPlotClasses
from PWG.PWGJE.EMCALJetTasks.Tracks.analysis.base.TriggerEfficiency import TriggerEfficiency
from os import getcwd
from ROOT import kRed, kBlue, kBlack, kGreen, kMagenta, kViolet, kOrange, kTeal, kYellow, kGray, kCyan

class TriggerEffDrawerMC(object):
    """
    Controller class drawing the trigger efficiency plot
    """

    def __init__(self):
        """ 
        Constructor initialising basic variables
        """
        self.__filehandler = MonteCarloFileHandler(True)
        self.__plotter = TriggerEfficiencyPlotMC()
        self.__nbins = 1
        
    def SetNumberOfPtHatBins(self, nbins):
        """
        Set the number of pthat bins
        """
        self.__nbins = nbins
        
    def SetBaseDirectory(self, path):
        """
        Change the base directory, and read in files for all pt-hat bins
        """
        for i in range(1, self.__nbins + 1):
            self.__filehandler.AddFile("%s/%02d/AnalysisResults.root" %(path, i), i)

    def CreatePlot(self, trigger):
        """
        Create the plot
        """
        styles = [Style(kBlack, 24), Style(kBlue, 25), Style(kRed, 26), Style(kGreen, 27), Style(kMagenta, 28), Style(kOrange, 29), \
                  Style(kTeal, 30), Style(kViolet, 31), Style(kGray, 32), Style(kYellow + 2, 33), Style(kCyan+3, 34), Style(kRed-9, 35)]
        triggerCalculator = None
        for i in range(1, self.__nbins + 1):
            collection = self.__filehandler.GetCollection().GetData(i)
            tcname = "tracksWithClusters"
            #tcname = "tracksAll"
            triggerCalculator = TriggerEfficiency(trigger, collection.GetData("MinBias").FindTrackContainer(tcname), collection.GetData(trigger).FindTrackContainer(tcname))
            self.__plotter.AddEfficiency(i, triggerCalculator.GetEfficiencyCurve(), styles[i-1])
        self.__plotter.Create()
        return self.__plotter
    
class TriggerEffSumDrawer(object):
    """
    Controller class drawing the trigger efficiency plot
    """

    def __init__(self):
        """ 
        Constructor initialising basic variables
        """
        self.__filehandler = MonteCarloFileHandler(True)
        self.__plotter = None
        self.__nbins = 1
        self.__summedData = None
        
    def SetNumberOfPtHatBins(self, nbins):
        """
        Set the number of pthat bins
        """
        self.__nbins = nbins
        
    def SetBaseDirectory(self, path):
        """
        Change the base directory, and read in files for all pt-hat bins
        """
        for i in range(1, self.__nbins + 1):
            self.__filehandler.AddFile("%s/%02d/AnalysisResults.root" %(path, i), i)
            
    def WriteSummedData(self, filename):
        self.__summedData.Write(filename)

    def CreatePlot(self, trigger):
        """
        Create the plot
        """
        self.__summedData = self.__filehandler.GetCollection().SumWeightedData()
        tcname = "tracksWithClusters"
        self.__plotter = TriggerEfficiencySumPlot(trigger, TriggerEfficiency(trigger, self.__summedData.GetData("MinBias").FindTrackContainer(tcname), self.__summedData.GetData(trigger).FindTrackContainer(tcname)))
        self.__plotter.Create()
        return self.__plotter
        
def DrawTriggerEfficiency(trigger, basedir = None):
    drawer = TriggerEffDrawerMC()
    drawer.SetNumberOfPtHatBins(9)
    if not basedir:
        basedir = getcwd()
    print "Using results from directory %s" %(basedir)
    drawer.SetBaseDirectory(basedir)
    return drawer.CreatePlot(trigger)

def DrawTriggerEfficiencySummed(trigger, basedir = None):
    drawer = TriggerEffSumDrawer()
    drawer.SetNumberOfPtHatBins(4)
    if not basedir:
        basedir = getcwd()
    print "Using results from directory %s" %(basedir)
    drawer.SetBaseDirectory(basedir)
    plot = drawer.CreatePlot(trigger)
    drawer.WriteSummedData("SumPtHatBins.root")
    return plot

def DrawMergedTriggerEfficiency(inputfile, datatype="tracks",  sourcetype="tracksWithClusters"):
    reader = ResultDataBuilder("resultfile", inputfile)
    data = reader.GetResults()
    plotter = TriggerEfficiencyPlotClasses()
    styles = {"EMCJHigh":Style(kBlue, 24), "EMCJLow":Style(kRed,25), "EMCGHigh":Style(kGreen,26), "EMCGLow":Style(kOrange,27)}
    for trigger in ["EMCJHigh", "EMCJLow", "EMCGHigh", "EMCGLow"]:
        dataminbias = None
        datatriggered = None
        if datatype == "clusters":
            dataminbias = data.GetData("MinBias").FindClusterContainer(sourcetype)
            datatriggered = data.GetData(trigger).FindClusterContainer(sourcetype)
        elif datatype == "tracks":
            dataminbias = data.GetData("MinBias").FindTrackContainer(sourcetype)
            datatriggered = data.GetData(trigger).FindTrackContainer(sourcetype)
        plotter.AddTriggerEfficiency(trigger, TriggerEfficiency(trigger, dataminbias, datatriggered).GetEfficiencyCurve(), styles[trigger])
    plotter.Create()
    return plotter
