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

class HepDataReader(object):
    '''
    Class reading hepdata input
    '''


    def __init__(self, filename, dataname):
        '''
        Constructor
        '''
        self.__result = DataCollection(dataname)
        self.ReadFile(filename)
        
    def GetData(self):
        return self.__result        

    def ReadFile(self, filename):
        inputfile = open(filename)
        for line in inputfile:
            data = self.__ProcessDatapoint(line.replace("\n",""))
            if data:
                self.__result.AddDataPoint(data)
        inputfile.close()
        print "Successfully read in %d points" %(len(self.__result.GetPointList()))
        
    def __ProcessDatapoint(self, line):
        line = line.replace("E","e")
        tokens = line.split("\t")
        values = self.__RemoveEmpty(tokens)
        print values
        if not self.TestDigit(values[0]):
            print "%s is not a digit" %(values[0])
            return None
        result = Datapoint(float(values[0]), float(values[3]), float(values[0])-float(values[1]))
        result.AddErrorSource("stat", float(values[5]), float(values[4]))
        result.AddErrorSource("sys", float(values[7]), float(values[6]))
        result.Print()
        return result
    
    def TestDigit(self, value):
        try:
            test = float(value)
        except ValueError:
            return False
        return True
    
    def __RemoveEmpty(self, inputlist):
        output = []
        for entry in inputlist:
            if len(entry):
                output.append(entry)
        return output