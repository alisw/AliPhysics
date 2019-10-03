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
from ROOT import TFile,TIter,TObject,gDirectory,gROOT
from copy import copy, deepcopy
from PWGJE.EMCALJetTasks.Tracks.analysis.base.DataSet import DataSet
from PWGJE.EMCALJetTasks.Tracks.analysis.base.MergeException import MergeException
from PWGJE.EMCALJetTasks.Tracks.analysis.base.struct.ParticleTHnSparse import ParticleTHnSparse

class ResultData(object):
    """
    Container for the result data
    Keeps the data containers for different trigger classes
    Access is implemented by the trigger class name
    """
    
    class DataException(Exception):
        """
        Exception handling for the result data
        """
        
        def __init__(self, containerName, trigger):
            """
            Initialise exception with the container name which raised the exception
            and a trigger which is not found
            """
            self.__containerName = containerName
            self.__trigger = trigger
            
        def __str__(self):
            """
            Create error message
            """
            return "Data for trigger class %s not found in container %s" %(self.__trigger, self.__containerName)
        
        def SetTrigger(self, trigger):
            """
            Set trigger name
            """
            self.__trigger = trigger
            
        def GetTrigger(self):
            """
            Access trigger name
            """
            return self.__trigger
        
        def SetContainerName(self, containerName):
            """
            Set container name
            """
            self.__containerName = containerName
            
        def GetContainerName(self):
            """
            Access container name
            """
            return self.__containerName
        
        Trigger = property(GetTrigger, SetTrigger)
        ContainerName = property(GetContainerName, SetContainerName)
    
    def __init__(self, name):
        """
        Initialise the result container
        """
        self.__name = name
        self.__data = {}
        self.__mctruth = None
        
        # for iterator
        self.__currentIndex = 0
        
    def __copy__(self):
        """
        Shallow copy constructor
        """
        print "Simple copy called from %s" %(self.__class__)
        newobject = ResultData(self.Name)
        if self.__mctruth:
            newobject.MCTruth = copy(self.__mctruth)
        for trigger in self.__data.keys():
            newobject.SetData(trigger, copy(self.__data[trigger]))
        return newobject
     
    def __deepcopy__(self, memo):
        """
        Deep copy constructor
        """
        print "deep copy called from %s" %(self.__class__)
        newobject = ResultData(self.Name)
        if self.__mctruth:
            newobject.MCTruth = deepcopy(self.__mctruth, memo)
        for trigger in self.__data.keys():
            newobject.SetData(trigger, deepcopy(self.__data[trigger], memo))
        return newobject
        
    def SetName(self, name):
        """
        Change name of the result container
        """
        self.__name = name
        
    def SetData(self, trigger, data):
        """
        Add a new trigger class with the corresponding data to the result
        container
        """
        self.__data[trigger] = data
        
    def GetName(self):
        """
        Access name of the trigger class
        """
        return self.__name
    
    def SetMCTruth(self, MCtrueSpectrum):
        """
        Set the MC truth to the result data
        """
        self.__mctruth = MCtrueSpectrum
        
    def GetMCTruth(self):
        """
        Access MC truth from the result data
        """
        return self.__mctruth
    
    def GetData(self, trigger):
        """
        Find data for a given trigger class by its name
        Raising a data exception if the trigger class is not avialable
        """
        if not self.__data.has_key(trigger):
            raise self.DataException(self.__name, trigger)
        return self.__data[trigger]
    
    # Properies
    Name = property(GetName, SetName)
    MCTruth = property(GetMCTruth, SetMCTruth)
    
    def Add(self, other):
        """
        Add MCTruth and datasets from other data to this MCTruth and the corresponding data set
        """
        if not isinstance(other, ResultData):
            raise MergeException("Type incompatibility: this(ResultData), other(%s)" %(str(other.__class__)))
        nfailure =0
        for trigger in self.GetListOfTriggers():
            if other.HasTrigger(trigger):
                self.__data[trigger].Add(other.GetData(trigger))
            else:
                nfailure += 1
        if self.__mctruth:
            self.__mctruth.Add(other.MCTruth)
        if nfailure > 0:
            raise MergeException("Unmerged histograms in this data")
    
    def Scale(self, scalefactor):
        """
        Scale all datasets and the MC truth by the scalefactor
        """
        for triggerdata in self.__data.values():
            triggerdata.Scale(scalefactor)
        # Scale also the MC truth
        self.__mctruth.Scale(scalefactor)
        
    def Write(self, rootfilename):
        """
        Write Structure to file
        """
        writer = TFile(rootfilename, "Recreate")
        writer.cd()
        for triggername, triggerdata in self.__data.iteritems():
            rootprim = triggerdata.GetRootPrimitive(triggername)
            rootprim.Write(triggername, TObject.kSingleKey)
        if self.__mctruth:
            self.__mctruth.GetRootPrimitive().Write("MCTruth", TObject.kSingleKey)
        writer.Close() 
        
    @staticmethod
    def BuildFromRootFile(filename, name):
        result = ResultData(name)
        inputfile = TFile.Open(filename)
        gROOT.cd()
        keyIter = TIter(inputfile.GetListOfKeys())
        key = keyIter.Next()
        while key:
            if key.GetName() == "MCTruth":
                result.SetMCTruth(ParticleTHnSparse(key.ReadObj()))
            else:
                result.SetData(key.GetName(), DataSet.BuildFromRootPrimitive(key.ReadObj()))
            key = keyIter.Next()
        inputfile.Close()
        print "Results successfully reconstructed from file %s" %(filename)
        #result.Print()
        return result
    
    def Print(self):
        """
        Print status of the result data sets
        """
        print "Content of result data:"
        for trg, data in self.__data.items():
            print "Trigger %s" %(trg)
            data.Print()
        
    def GetListOfTriggers(self):
        """
        Provide a list of trigger classes which are currently stored
        in the container
        """
        return self.__data.keys()
    
    def HasTrigger(self, triggername):
        """
        Check if the container has a given trigger type
        """
        return triggername in self.__data.keys()
    
    def __iter__(self):
        """
        Initialise the iterator
        """
        self.__curentIndex = 0
        return self
    
    def next(self):
        """ 
        Iterate over trigger classes
        The iterator always returns a tuple of triggerclass and data
        """
        if self.__currentIndex < len(self.__data.keys()):
            result = self.__data.keys()[self.__currentIndex], self.__data[self.__data.keys()[self.__currentIndex]]
            self.__currentIndex += 1
            return result
        else:
            raise StopIteration
    
    def __len__(self):
        """
        Return the amount of trigger classes
        """
        return len(self.__data)
    
    def __getItem__(self, key):
        """
        Access trigger item by the [] operator
        """
        return self.GetData(key)
    
    def __contains__(self, item):
        return item in self.__data.keys()
    