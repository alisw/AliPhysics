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
Utility classes for datapoints and converters to ROOT graphics objects

@author: Markus Fasel
"""

from ROOT import TGraph, TGraphAsymmErrors
import math

class Datapoint:
    """
    Representation of data for a single point
    """
    
    def __init__(self, x, y, dx):
        """
        Constructor
        """
        self.__x = x
        self.__y = y
        self.__dx = dx
        self.__upperErrors = {}
        self.__lowerErrors = {}
        
    def AddErrorSource(self, name, lower, upper):
        """
        Add error to the datapoint
        """
        self.__upperErrors[name] = upper
        self.__lowerErrors[name] = lower
        
    def GetX(self):
        """
        Access to x-value of the data point
        """
        return self.__x
        
    def GetY(self):
        """
        Access to y-value of the data point
        """
        return self.__y
     
    def GetLowerErrorForSource(self, name):
        """
        Access to lower error for a given error source
        """
        if name == "total":
            return self.GetTotalLowerError()
        if not name in self.__lowerErrors.keys():
            return 0
        return self.__lowerErrors[name]
    
    def GetUpperErrorForSource(self, name):
        """
        Access to upper error for a given error source
        """
        if name == "total":
            return self.GetTotalUpperError()
        if not name in self.__upperErrors.keys():
            return 0
        return self.__upperErrors[name]
    
    def GetUpperLimitForSource(self, source):
        """
        Access to upper limit under a given error source
        """
        return self.__y + self.GetUpperErrorForSource(source)
    
    def GetLowerLimitForSource(self, source):
        """
        Access to lower limit under a given error source
        """
        return self.__y - self.GetLowerErrorForSource(source)

    def GetRelativeLowerError(self, source):
        return self.GetLowerErrorForSource(source)/self.__y
    
    def GetRelativeUpperError(self, source):
        return self.GetUpperErrorForSource(source)/self.__y
             
    def GetDX(self):
        """
        Access to uncertainty in x direction
        """
        return self.__dx
        
    def GetTotalLowerError(self):
        """
        calculate total error as quadratic sum of the single components
        """
        sumofsquares = 0
        for error in self.__lowerErrors.values():
            sumofsquares += math.pow(error, 2)
        return math.sqrt(sumofsquares)

    def GetTotalUpperError(self):
        """
        calculate total error as quadratic sum of the single components
        """
        sumofsquares = 0
        for error in self.__lowerErrors.values():
            sumofsquares += math.pow(error, 2)
        return math.sqrt(sumofsquares)
     
    def __eq__(self, other):
        return self.__x == other.__x
        
    def __lt__(self, other):
        return self.__x < other.__x
        
    def __le__(self, other):
        return self.__x <= other.__x
        
    def __gt__(self, other):
        return self.__x > other.__x
        
    def __ge__(self, other):
        return self.__x >= other.__x
        
    def __ne__(self, other):
        return self.__x != other.x
    
    def __str__(self):
        """
        Create string representation of the point
        """
        result = "%f +- %f GeV/c: %e" %(self.__x, self.__dx, self.__y)
        for source in self.__lowerErrors.keys():
            result += " + %e - %e (%s)" %(self.__upperErrors[source], self.__lowerErrors[source], source)
        if len(self.__lowerErrors):
            result += " [+ %e -%e (total)]" %(self.GetTotalUpperError(), self.GetTotalLowerError())
        return result

    def Print(self):
        """
        Print point representation
        """
        print str(self)

class DataCollection:
    """
    Collection of data points
    """

    def __init__(self, name):
        '''
        Constructor
        '''
        self.__name = name
        self._pointlist = []
        
    def AddDataPoint(self, point):
        self._pointlist.append(point)
        
    def AddDataXY(self, x, y):
        self._pointlist.append(Datapoint(x, y, 0.))
        
    def AddDataWithErrors(self, x, y, dx, dy):
        point = Datapoint(x,y,dx),
        point.AddErrorSource("error", dy, dy)
        self._pointlist.append(point)
        
    def GetPointList(self):
        return self._pointlist
    
    def MakeErrorGraphForSource(self, source):
        result = TGraphAsymmErrors()
        counter = 0 
        for point in sorted(self._pointlist):
            result.SetPoint(counter, point.GetX(), point.GetY())
            result.SetPointError(counter, point.GetDX(), point.GetDX(), point.GetLowerErrorForSource(source), point.GetUpperErrorForSource(source))
            counter += 1
        return result
    
    def MakeLimitCurve(self, source, direction):
        result = TGraph()
        counter = 0
        for point in self._pointlist:
            error = 0
            if direction == "upper":
                error = point.GetUpperErrorForSource(source)
            elif direction == "lower":
                error = -1. * point.GetLowerErrorForSource(source)
            elif direction == "cental":
                error = 0
            result.SetPoint(counter, point.GetX(), point.GetY() + error)
            counter += 1
        return result

    def Print(self):
        print "Data collection :s" %(self.__name)
        print "====================================================================="
        for point in sorted(self._pointlist):
            point.Print()
