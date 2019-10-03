#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TString.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskTrigEff.h"
#include "AliAnalysisDataContainer.h"
#include "AliVEventHandler.h"

#endif


AliAnalysisTaskSE * AddTaskTrigEff(const char * outname, Bool_t isMC, Int_t ntrk) {
  // Adds my task
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskTrigEff", "No analysis manager to connect to.");
    return NULL;
  }  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskTrigEff", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
  // Configure analysis
  //===========================================================================

  TString loadTaskStr = "AliAnalysisTaskTrigEff.cxx+";
  if (gProof != NULL) {
    gProof->Load( loadTaskStr.Data() );
  }
  else {
    gROOT->LoadMacro( loadTaskStr.Data() );
  }

  AliAnalysisTaskTrigEff * task = new AliAnalysisTaskTrigEff("TaskTrigEff");

  task->SetIsMC(isMC);
  task->SetNtrkCut(2);
  mgr->AddTask(task);

  // connect input/output
  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("taskout",
							    TList::Class(),
							    AliAnalysisManager::kOutputContainer,
							    outname);

  mgr->ConnectInput(task, 0, cinput0);
  mgr->ConnectOutput(task, 1, coutput1);

  return task;

}
