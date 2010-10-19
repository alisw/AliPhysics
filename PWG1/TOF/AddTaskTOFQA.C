/*  created by fbellini@cern.ch on 14/09/2010 */
/*  last modified by fbellini   on 14/09/2010 */

#include "AliAnalysisTaskTOFqa.h"
AliAnalysisTaskSE * AddTaskTOFQA() 
{
  // Task for checking TOF QA
 
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }   

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTas", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
    // Configure analysis
    //===========================================================================
  // Create the task
  AliAnalysisTaskTOFqa *task = new AliAnalysisTaskTOFqa("taskTOFqa");
  //AliLog::SetClassDebugLevel("AliAnalysisTaskTOFqa",1);
  mgr->AddTask(task);

  //--------------- set the filtering ------------
  // Barrel Tracks
  AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Standard");
  esdTrackCutsL->SetMinNClustersTPC(50); // ok 50
  esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5); // ok 3.5
  esdTrackCutsL->SetMaxCovDiagonalElements(2, 2, 0.5, 0.5, 2);//ok
  esdTrackCutsL->SetRequireTPCRefit(kTRUE);//ok (?)
  esdTrackCutsL->SetMaxDCAToVertexXY(3.0); // ok
  esdTrackCutsL->SetMaxDCAToVertexZ(3.0); // ok
  esdTrackCutsL->SetRequireSigmaToVertex(kTRUE); //ok ?
  esdTrackCutsL->SetAcceptKinkDaughters(kFALSE); // ok
  esdTrackCutsL->SetMaxNsigmaToVertex(4.0);
  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
  trackFilter->AddCuts(esdTrackCutsL);
  task->SetTrackFilter(trackFilter);
   
  
  // Create containers for input/output
  AliAnalysisDataContainer *cInputTOFqa = mgr->CreateContainer("cInputTOFqa",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *cOutGeneralTOFqa = mgr->CreateContainer("cOutGeneralTOFqa",TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName()));
   AliAnalysisDataContainer *cOutExpertsTOFqa = mgr->CreateContainer("cOutExpertsTOFqa",TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:TOF_Performance",mgr->GetCommonFileName()));
  // Attach i/o
  mgr->ConnectInput(task, 0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cOutGeneralTOFqa);
  mgr->ConnectOutput(task, 2, cOutExpertsTOFqa);
  
  return task;
}
