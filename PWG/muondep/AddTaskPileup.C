#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliCounterCollection.h"
#include "AliAnalysisTaskPileup.h"
#endif

AliAnalysisTaskPileup* AddTaskPileup(TString ocdbPath="alien://folder=/alice/data/2010/OCDB"){

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   
   AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskPileup", "No analysis manager to connect to.");
      return NULL;
   }   

   // This task requires an ESD or AOD output handler.
   // Check this using the analysis manager.
   //===============================================================================
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if ( ! type.Contains("ESD") && ! type.Contains("AOD") ) {
     ::Error("AddTaskPileup", "Pileup task needs the manager to have an ESD input handler.");
     return NULL;
   }

   TString baseOutName = "pileupAnalysis.root";
   TString outputfile = mgr->GetCommonFileName();
   if ( ! outputfile.IsNull() ) outputfile += ":PWG3_muon_Pileup";
   else outputfile = baseOutName;

   AliAnalysisDataContainer* coutput1 = mgr->CreateContainer("PileupCounter",AliCounterCollection::Class(),AliAnalysisManager::kOutputContainer,outputfile);
   AliAnalysisDataContainer* coutput2 = mgr->CreateContainer("PileupCorrections",TObjArray::Class(),AliAnalysisManager::kParamContainer,outputfile);


   // Create the task, add it to the manager and configure it.
   //===========================================================================   

   AliAnalysisTaskPileup* pileupTask = new AliAnalysisTaskPileup("PileupTask");
   pileupTask->SetDefaultStorage(ocdbPath);
   mgr->AddTask(pileupTask);
   
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   
   mgr->ConnectInput  (pileupTask,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (pileupTask,  1, coutput1);
   mgr->ConnectOutput (pileupTask,  2, coutput2);

   return pileupTask;
}   
