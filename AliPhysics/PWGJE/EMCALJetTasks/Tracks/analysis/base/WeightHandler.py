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
from PWGJE.EMCALJetTasks.Tracks.analysis.base.DataCollection import DataCollection, Datapoint

class PtHatBin:
    """
    Data entry of the pt-hat bin
    """
    
    def __init__(self, binId, weighthist, trialsHist):
        """
        Construct pt-hat bin. Automatically calculates the weight
        from the cross section and the number of trials
        """
        self.__binId = binId
        self.Set(weighthist, trialsHist)
        
    def Set(self, crossechist, ntrialshist):
        """
        Set the weight using a given cross section and number of trials hist
        """
        self.__weight = crossechist.GetBinContent(self.__binId + 1) /  ntrialshist.GetBinContent(self.__binId + 1)
        
    def ReweightSpectrum(self, spectrum):
        """
        Rescale spectrum by the weight of the given pt-hat bin
        """
        spectrum.Scale(self.__weight)
        
    def GetWeight(self):
        """
        Get value of the weight for the given bin
        """
        return self.__weight
        
    def GetBinID(self):
        """
        Get the bin id
        """
        return self.__binId
    
    def __eq__(self, other):
        if type(other) is int:
            return self.__binId == int(other)
        if type(other) is PtHatBin:
            return self.__binId == other.GetBinID()
        
    def __lt__(self, other):
        if type(other) is int:
            return self.__binId < int(other)
        if type(other) is PtHatBin:
            return self.__binId < other.GetBinID()

    def __le__(self, other):
        if type(other) is int:
            return self.__binId <= int(other)
        if type(other) is PtHatBin:
            return self.__binId <= other.GetBinID()

    def __gt__(self, other):
        if type(other) is int:
            return self.__binId > int(other)
        if type(other) is PtHatBin:
            return self.__binId > other.GetBinID()

    def __ge__(self, other):
        if type(other) is int:
            return self.__binId >= int(other)
        if type(other) is PtHatBin:
            return self.__binId >=other.GetBinID()
        
    def __str__(self):
        return "Pt-hat bin %d: weight: %e" %(self.__binId, self.__weight)
    
    def Print(self):
        """
        Print bin and weight factor
        """
        print str(self)

class WeightHandler:
    """
    Handler class for cross section weights in case the simulation was done in pt hat bins
    """
    
    def __init__(self):
        '''
        Constructor
        '''
        self.__pthatbins = []
        
    def AddPtHatBin(self, binID, crosssechist, ntrialshist):
        """
        Insert new pt-hat bin or replace weights
        """ 
        if not binID in self.__pthatbins:
            self.__pthatbins.append(PtHatBin(binID, crosssechist, ntrialshist))
        else:
            self.__pthatbins[self.__pthatbins.index(binID)].Set(crosssechist, ntrialshist)
            
    def GetWeight(self, pthardbin):
        return self.__pthatbins[self.__pthatbins.index(pthardbin)].GetWeight()
        
    def ReweightSpectrum(self, binId, spectrum):
        """ 
        Reweight spectum from a given pt-hat bin
        """
        if binId in self.__pthatbins:
            self.__pthatbins[self.__pthatbins.index(binId)].ReweightSpectrum(spectrum)
            
    def GetWeightingCurve(self):
        """
        Build graph from the different weights
        """
        result = DataCollection("weightlist")
        for mybin in self.__pthatbins:
            result.AddDataPoint(Datapoint(mybin.GetBinID(), mybin.GetWeight(), 0.5))
        return result.MakeLimitCurve(None, direction="central")
    
    def _FindPtHardBin(self, binID):
        ''' 
        Find pt-hard bin by bin id
        '''
        result = None
        for entry in self.__pthatbins:
            if entry.GetBinID() == binID:
                result = Entry
                break
        return result
            
    def Print(self):
        print "Weighting factors: "
        print "==============================="
        for entry in self.__pthatbins:
            print str(entry)