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
from ROOT import TH1D
from numpy import array

class Histogram(object):
    
    class Bin(object):
        
        def __init__(self, xmin, xmax, value, error):
            self.__xmin = xmin
            self.__xmax = xmax
            self.__value = value
            self.__error = error
      
        def GetXmin(self):
            return self.__xmin
    
        def GetXmax(self):
            return self.__xmax
    
        def GetBinCenter(self):
            return (self.__xmax + self.__xmin)/2.
    
        def GetValue(self):
            return self.__value
    
        def GetError(self):
            return self.__error
    
        def __eq__(self, other):
            return self.__xmin == other.__xmin and self.__xmax == other.__xmax
    
        def __lt__(self, other):
            return self.__xmax < other.__xmin
    
        def __le__(self, other):
            return self.__xmax <= other.__xmax
    
        def __gt__(self, other):
            return self.__xmin > other.__xmax
    
        def __ge__(self, other):
            return self.__xmin >= other.__xmin
    
    def __init__(self):
        self.__bins = []
    
    def AddBin(self, xmin, xmax, val, err):
        self.__bins.append(self.Bin(xmin, xmax, val, err))
        
    def SetData(self, datacoll, errorsource):
        graph = datacoll.MakeErrorGraphForSource(errorsource)
        for i in range(0, graph.GetN()):
            x = graph.GetX()[i]
            dx = graph.GetEXlow()[i]
            y = graph.GetY()[i]
            dy = graph.GetEYlow()[i]
            print "Doing bin from %f to %f" %(x-dx, x+dx)
            self.AddBin(x-dx, x+dx, y, dy)
    
    def __GetBinning(self):
        limits = []
        counter = 0
        for ib in sorted(self.__bins):
            if counter == 0:
                limits.append(ib.GetXmin())
            limits.append(ib.GetXmax())
            counter += 1
        return limits
  
    def MakeROOTHist(self, name, title):
        limits = self.__GetBinning()
        result = TH1D(name, title, len(limits)-1, array(limits))
        for mybin in self.__bins:
            ib = result.GetXaxis().FindBin(mybin.GetBinCenter())
            result.SetBinContent(ib, mybin.GetValue())
            result.SetBinError(ib, mybin.GetError())
        return result