#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TString.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskNanoAODESEFilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliVEventHandler.h"
#include "AliESEHelpers.h"

#endif


AliAnalysisTaskSE * AddTaskNanoAODFilter(Int_t iMC, Bool_t savecuts = 0) {
  // Adds my task
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskNanoAODESEFilter", "No analysis manager to connect to.");
    return NULL;
  }  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskNanoAODESEFilter", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
  // Configure analysis
  //===========================================================================

  //  gROOT->LoadMacro("AliAnalysisTaskNanoAODESEFilter.cxx+");

  AliAnalysisTaskNanoAODFilter * task = new AliAnalysisTaskNanoAODFilter("TaskEseFilter", savecuts);
  
  mgr->AddTask(task);
  task->SetMCMode(iMC);

  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
  //  mgr->ConnectOutput(task, 1, mgr->GetCommonOutputContainer());
  
  if(savecuts) {
    AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer("evtcuts", AliAnalysisCuts::Class(),  AliAnalysisManager::kOutputContainer,"cuts.root");
    AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer("trkcuts", AliAnalysisCuts::Class(),  AliAnalysisManager::kOutputContainer,"cuts.root");
    mgr->ConnectOutput(task, 1, coutputpt1);
    mgr->ConnectOutput(task, 2, coutputpt2);
  }

  


  return task;

}



