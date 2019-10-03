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
Handler module for Monte-Carlo output from the ALICE Lego Trains. Lego train output
can be min. bias events or productions in pt-hat bins

@author: Markus Fasel
"""

from PWGJE.EMCALJetTasks.Tracks.analysis.base.WeightHandler import WeightHandler
from PWGJE.EMCALJetTasks.Tracks.analysis.base.FileHandler import LegoTrainFileReader
from PWGJE.EMCALJetTasks.Tracks.analysis.base.SpectraSum import SpectraSum

class MonteCarloDataCollection(object):
    """
    Collection of Monte-Carlo based outputs
    """

    def __init__(self, isPtHat = False):
        """
        Constructor
        """
        self.__weighthandler  = None
        if isPtHat:
            self.__weighthandler = WeightHandler()
        self.__data = {"All":None}
        
    def AddData(self, results, pthatbin = -1,  weightdata = None):
        """
        Add new data (with or without pthat bins)
        """ 
        if pthatbin >= 0:
            self.__data[pthatbin] = results
            self.__weighthandler.AddPtHatBin(pthatbin, weightdata["crosssection"], weightdata["trials"])
        else:
            self.__data["All"] = results
            
    def GetData(self, pthatbin = -1):
        """
        Access to data (if necessary in a given pt-hat bin
        """
        if pthatbin >= 0:
            return self.__data[pthatbin]
        return self.__data["All"]
    
    def GetWeigthHandler(self):
        """
        Access to the weight handler
        """
        return self.__weighthandler
    
    def SumWeightedData(self):
        """
        Sum weighted containers from the different pthat bins
        """
        if not self.__weighthandler:
            print "No weight handler"
            return None
        summer = SpectraSum()
        for pthatbin in self.__data.keys():
            if pthatbin == "All":
                continue
            self.__weighthandler.ReweightSpectrum(pthatbin, self.__data[pthatbin])
            summer.AddSpectrum(self.__data[pthatbin])
        return summer.GetSummedSpectrum()
            
class MonteCarloFileHandler(object):
    """
    Class handling the reading of one file or a set of MonteCarlo files
    """
    
    def __init__(self, hasPtHardBins = False):
        """
        Constructor
        """
        self.__datacollection = MonteCarloDataCollection(hasPtHardBins)
        self.__histlist = ""
        
    def GetCollection(self):
        """
        Access to the file collection
        """
        return self.__datacollection
    
    def AddFile(self, filename, pthatbin = -1, isNew = True):
        """
        Handle new file
        """
        reader = LegoTrainFileReader(filename, isMC = True, isNew = isNew)
        if pthatbin >= 0:
            reader.SetReadWeights() 
        self.__datacollection.AddData(reader.ReadFile(), pthatbin, reader.GetWeightHistograms())
        
class MonteCarloFileMerger(object):
    """
    Class merging Monte-Carlo files in pt-hat bins, weighted by the cross section
    """
    
    def __init__(self):
        """
        Constructor
        """
        self.__reader = MonteCarloFileHandler(True)
        
    def AddFile(self, filename, pthatbin):
        """
        Add next file
        """
        self.__reader.AddFile(filename, pthatbin)
    
    def MergeAndWrite(self, outputfile):
        summed = self.__reader.GetCollection().SumWeightedData()
        summed.Write(outputfile)
        
def MergePtHardBins(outputfile, basedir, firstbin, lastbin):
    """
    Merge files from different pt-hard bins, weighted by the cross section, into one file
    """
    merger = MonteCarloFileMerger()
    for pthardbin in range(firstbin, lastbin+1):
        merger.AddFile("%s/%02d/AnalysisResults.root" %(basedir, pthardbin), pthardbin)
    merger.MergeAndWrite(outputfile)
