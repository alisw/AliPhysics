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
from PWGJE.EMCALJetTasks.Tracks.analysis.base.FileHandler import LegoTrainFileReader

def Run():
    periods = ["LHC13c", "LHC13d", "LHC13e", "LHC13f", "merged"]
    for period in periods:
        GetNumberOfEventsForPeriod(period)
        
def GetNumberOfEventsForPeriod(period):
    handler = LegoTrainFileReader("%s/AnalysisResults.root" %(period))
    data = handler.ReadFile()
    print "Period %s" %(period)
    for trg in data.GetListOfTriggers():
        PrintEventCountForTrigger(data.GetData(trg).FindTrackContainer("tracksAll"), trg)
        
def PrintEventCountForTrigger(data, trigger):
    data.SetVertexRange(-10, 10)
    data.SetPileupRejection(True)
    print "%s: %d" %(trigger, data.GetEventCount())
    
if __name__ == "__main__":
    Run()