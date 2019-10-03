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
Wrapper class around a ROOT THnSparse adding the functionality of cut and project in one step,
making the handling of THnSparses a lot more handy for users.

@author: Markus
"""

from copy import copy,deepcopy
from numpy import array as nparray

class AxisFormat(object):
    '''
    Definition of the axis format of a THnSparse
    '''
    
    def __init__(self, formatname):
        '''
        Constructor
        '''
        self._axes = {}
        self.__formatname = ""
        
    def GetAxes(self):
        '''
        Get the list of axes defined
        '''
        return self._axes
    
    def FindAxis(self, axisname):
        '''
        Find axis by axis name. Returns the dimension of the axis.
        '''
        result = -1
        if axisname in self._axes.keys():
            result = self._axes[axisname]
        return result
    
    def _Deepcopy(self, other, memo):
        '''
        Performing deep copy
        '''
        self._axes = deepcopy(other.GetAxes(), memo)
       
    def _Copy(self, other):
        '''
        Performing shallow copy
        '''
        self._axes = copy(other.GetAxes()) 
        
    def GetAxisName(self, dim):
        '''
        Get the name of the axis by dimension
        '''
        if not dim in self._axes.values():
            return ""
        result = ""
        for k,v in self._axes.iteritems():
            if v == dim:
                result = k 
                break
        return result
    
    def Print(self):
        for axis, dimension in self._axes.iteritems():
            print "Axis %s with dimension %d" %(axis, dimension)
        

class THnSparseCut(object):
    '''
    Cut class used in the THnSparse wrapper
    '''
    
    def __init__(self, axisname, minv, maxv):
        '''
        Constructor
        '''
        self.__axisname = axisname
        self.__minimum = minv
        self.__maximum = maxv
        
    def GetCutname(self):
        '''
        Get axis name
        '''
        return self.__axisname
        
    def GetMinimum(self):
        '''
        Get the minimum of the range
        '''
        return self.__minimum
    
    def GetMaximum(self):
        '''
        Get the maximum of the range
        '''
        return self.__maximum
    
    def SetMinimum(self, minv):
        '''
        Set the minimum of the range
        '''
        self.__minimum = minv
        
    def SetMaximum(self, maxv):
        '''
        Set the maximum of the range
        '''
        self.__maximum = maxv

class THnSparseWrapper(object):
    '''
    Wrapper class around THnSparse applying cuts on axes and performing projections
    '''
    
    def __init__(self, rootthnsparse):
        '''
        Constructor
        '''
        self._rootthnsparse = rootthnsparse
        self._axisdefinition = None
        self._cutlist = []
        
    def __deepcopy__(self, memo):
        '''
        Deep copy constructor
        '''
        result = THnSparseWrapper(deepcopy(self._rootthnsparse))
        result.CopyCuts(self._cutlist, True)
        return result
        
    def __copy__(self):
        '''
        Shallow copy constructor
        '''
        result = THnSparseWrapper(copy(self._rootthnsparse))
        result.CopyCuts(self._cutlist, False)
        return result
    
    def CopyCuts(self, reference, isDeep):
        '''
        Copy cuts into this object from a reference object
        '''
        for cut in reference.GetListOfCuts():
            newcut = None
            if isDeep:
                newcut = deepcopy(cut)
            else:
                newcut = copy(cut)
            self._cutlist.append(newcut)
        
    def GetListOfCuts(self):
        '''
        Access list of cuts
        '''
        return self._cutlist

    def GetHistogram(self):
        '''
        Access to underlying root histogram
        '''
        return self._rootthnsparse
    
    def GetHistogramName(self):
        '''
        Get the name of the underlying histogram
        '''
        return self._rootthnsparse.GetName()
    
    def Add(self, otherwrapper):
        self._rootthnsparse.Add(otherwrapper.GetHistogram())
    
    def Scale(self, scalefactor):
        self._rootthnsparse.Scale(scalefactor)
        
    def GetAxisDefinition(self):
        return self._axisdefinition
        
    def ApplyCut(self, axisname, minv, maxv):
        '''
        Apply cut on a given variable, defined by its axis name
        minv and maxv define the range. If either of them is None, the range is 
        open on one side.
        '''
        if not self._axisdefinition or self._axisdefinition.FindAxis(axisname) < 0:
            print "No axis definition or axis name (%s) not found" %(axisname)
            if self._axisdefinition:
                print "Known axes:"
                self._axisdefinition.Print()
            return
        existing = self.__FindCut(axisname)
        if not existing:
            self._cutlist.append(THnSparseCut(axisname, minv, maxv))
        else:
            existing.SetMinimum(minv)
            existing.SetMaximum(maxv)
            
    def RemoveCut(self, axisname):
        '''
        Remove cut again from the list
        '''
        for entry in self._cutlist:
            if entry.GetCutname() == axisname:
                self._cutlist.remove(entry)
        
    def ResetAxis(self, axisname):
        '''
        Reset axis range
        '''
        if not self._axisdefinition or self._axisdefinition.FindAxis(axisname) < 0:
            return
        myaxis = self._rootthnsparse.GetAxis(self._axisdefinition.FindAxis(axisname))
        myaxis.SetRange(0, myaxis.GetNbins()+1)

    def Projection1D(self, histname, axisname):
        '''
        Make projection, applying cuts defined before, and releasing the cuts afterwards.
        Projects to 1D with the axisname as dimension
        '''
        if not self._axisdefinition or self._axisdefinition.FindAxis(axisname) < 0:
            print "No axis definition or axis %s not found" %(axisname)
            return None
        self._PrepareProjection()
        result = self._rootthnsparse.Projection(self._axisdefinition.FindAxis(axisname))
        result.SetName(histname)
        self._CleanumProjection()
        return result
    
    def Projection2D(self, histname, axisdictionary):
        '''
        Make projection, applying cuts defined before, and releasing the cuts afterwards.
        Projects to 2D with the content in the axis dictionary as dimensions
        Dictionary works in the way name -> dimension, starting with 0
        '''
        if not self._axisdefinition:
            return None
        hasfound = True
        for axisname in axisdictionary.keys():
            if self._axisdefinition.FindAxis(axisname):
                hasfound = False
                break
        if not hasfound:
            return None
        self._PrepareProjection()
        xdim = None
        ydim = None
        for k,v in axisdictionary.iteritems():
            if v == 1:
                ydim = self._axisdefinition.FindAxis(k)
            else:
                xdim = self._axisdefinition.FindAxis(k)
        result = self._rootthnsparse.Projection(ydim, xdim)
        result.SetName(histname)
        self._CleanumProjection()
        return result
    
    def ProjectionND(self, histname, axisdictionary):
        '''
        Make projection, applying cuts defined before, and releasing the cuts afterwards.
        Projects to 2D with the content in the axis dictionary as dimensions
        Dictionary works in the way name -> dimension, starting with 0
        '''
        if not self._axisdefinition:
            return None
        hasfound = True
        for axisname in axisdictionary.keys():
            if self._axisdefinition.FindAxis(axisname):
                hasfound = False
                break
        if not hasfound:
            return None
        self._PrepareProjection()
        axismap = {}
        for k,v in axisdictionary.iteritems():
            axismap[v] = k
        axislist = []
        for mydim in sorted(axismap.keys()):
            axislist.append(self._axisdefinition.FindAxis(axismap[mydim]))
        result = self._rootthnsparse.Projection(len(axislist), nparray(axislist))
        result.SetName(histname)
        self._CleanumProjection()
        return result
    
    def _PrepareProjection(self):
        '''
        Apply all requested cuts before the projection
        '''
        for entry in self._cutlist:
            myaxis = self._rootthnsparse.GetAxis(self._axisdefinition.FindAxis(entry.GetCutname()))
            minv = 0 if not entry.GetMinimum() else myaxis.FindBin(entry.GetMinimum())
            maxv = myaxis.GetNbins()+1 if not entry.GetMaximum() else myaxis.FindBin(entry.GetMaximum())
            myaxis.SetRange(minv, maxv)
            
    def _CleanumProjection(self):
        '''
        Reset all possible axis cuts
        Does not remove a cut again from the list, but only releases the THnSparse
        '''
        for entry in self._cutlist:
            self.ResetAxis(entry.GetCutname())

    def __FindCut(self, cutname):
        '''
        Find cut in list by the axis name
        '''
        if not len(self._cutlist):
            return None
        result = None
        for entry in self._cutlist:
            if entry.GetCutname() == cutname:
                result = entry
                break
        return result