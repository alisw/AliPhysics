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
from PWGJE.EMCALJetTasks.Tracks.analysis.base.FileHandler import LegoTrainFileReader, ResultStructureReader

class ResultDataBuilder(object):
    """
    General data structure builder interfacing different input file formats
    """
    
    def __init__(self, filetype, filename):
        """
        Constructor
        """
        reader = None
        if filetype == "lego":
            reader = LegoTrainFileReader(filename)
        elif filetype == "resultfile":
            reader = ResultStructureReader(filename)
        self.__results = None
        if reader:
            self.__results = reader.ReadFile()
            
    def GetResults(self):
        """
        Access to data
        """
        return self.__results
