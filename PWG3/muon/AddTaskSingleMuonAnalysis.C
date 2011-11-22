#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TList.h"

#include "AliLog.h"
#include "AliVEventHandler.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliCFContainer.h"

#include "AliAnalysisTaskSingleMu.h"
#endif

AliAnalysisTaskSingleMu* AddTaskSingleMuonAnalysis(Bool_t applyPhysicsSelection=kTRUE, Int_t fillNtupleScaleDown=0, Bool_t keepAll=kFALSE, Bool_t separateSpecialOut=kFALSE){

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddtaskSingleMuonAnalysis", "No analysis manager to connect to.");
      return NULL;
   }   

   // This task requires an ESD or AOD output handler.
   // Check this using the analysis manager.
   //===============================================================================
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (!type.Contains("ESD") && !type.Contains("AOD")) {
     ::Error("AddtaskSingleMuonAnalysis", "SingleMuon task needs the manager to have an ESD or AOD input handler.");
     return NULL;
   }

   TString currName = "";
   TString outputfile = mgr->GetCommonFileName();
   if ( ! outputfile.IsNull() ) {
     currName = ( applyPhysicsSelection ) ? ":PWG3_muon_SingleMu" : ":PWG3_muon_SingleMu_NoPS";
     outputfile += currName;
   }
   else outputfile = "singleMuonAnalysis.root";

   currName = ( applyPhysicsSelection ) ? "SingleMuonContainer" : "SingleMuonContainerNoPS";
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(currName.Data(),AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile);

   currName = ( applyPhysicsSelection ) ? "SingleMuon" : "SingleMuonNoPS";
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(currName.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);

   currName = ( applyPhysicsSelection ) ? "SingleMuonQA" : "SingleMuonQANoPS";
   AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(currName.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);

   currName = ( applyPhysicsSelection ) ? "SingleMuonMC" : "SingleMuonMCNoPS";
   AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(currName.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);

   AliAnalysisDataContainer *coutput5 = 0x0;
   if ( fillNtupleScaleDown > 0 ) {
     currName = ( applyPhysicsSelection ) ? "SingleMuonNtuple" : "SingleMuonNtupleNoPS";
     TString specialOutFilename = separateSpecialOut ? currName + ".root" : outputfile;
     coutput5 = mgr->CreateContainer(currName.Data(),TTree::Class(),AliAnalysisManager::kOutputContainer,specialOutFilename);
     coutput5->SetSpecialOutput();
   }


   // Create the task, add it to the manager and configure it.
   //===========================================================================
   TString taskName =  ( applyPhysicsSelection ) ? "SingleMuonAnalysisTask" : "SingleMuonAnalysisTaskNoPS";
   AliAnalysisTaskSingleMu *singleMuonAnalysisTask = new AliAnalysisTaskSingleMu(taskName.Data(), fillNtupleScaleDown, keepAll);
   mgr->AddTask(singleMuonAnalysisTask);
   if ( applyPhysicsSelection ) 
     singleMuonAnalysisTask->SelectCollisionCandidates(AliVEvent::kAny);
   
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   
   mgr->ConnectInput  (singleMuonAnalysisTask,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (singleMuonAnalysisTask,  1, coutput1);
   mgr->ConnectOutput (singleMuonAnalysisTask,  2, coutput2);
   mgr->ConnectOutput (singleMuonAnalysisTask,  3, coutput3);
   mgr->ConnectOutput (singleMuonAnalysisTask,  4, coutput4);

   if ( coutput5 )
     mgr->ConnectOutput (singleMuonAnalysisTask,  5, coutput5);

   return singleMuonAnalysisTask;
}
