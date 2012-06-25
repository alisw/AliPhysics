// $Id: AddTaskEMCALPi0V2.C 56081 2012-05-01 08:57:08Z loizides $

AliAnalysisTaskPi0V2 *AddTaskEMCALPi0V2()
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALPi0V2", "No analysis manager to connect to.");
    return NULL;
  }  

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEMCALPi0V2", "This task requires an input event handler");
    return NULL;
  }

  AliEPSelectionTask *eventplaneTask = new AliEPSelectionTask("EventplaneSelection");
  eventplaneTask->SelectCollisionCandidates(AliVEvent::kSemiCentral | AliVEvent::kCentral | AliVEvent::kAnyINT);

  eventplaneTask->SetTrackType("TPC");
  eventplaneTask->SetUsePtWeight();
  eventplaneTask->SetUsePhiWeight();
  eventplaneTask->SetSaveTrackContribution();

  AliESDtrackCuts* epTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "Standard");
  epTrackCuts->SetPtRange(0.1, 4);
  eventplaneTask->SetPersonalESDtrackCuts(epTrackCuts);

  mgr->AddTask(eventplaneTask);

  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("EPStat",TList::Class(), AliAnalysisManager::kOutputContainer,"EPoutput.root");
  mgr->ConnectInput(eventplaneTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(eventplaneTask,1,coutput1);
  
  // Create the task and configure it.
  //===========================================================================
  AliAnalysisTaskPi0V2* ana = new  AliAnalysisTaskPi0V2("Pi0v2Task");
  ana->SelectCollisionCandidates( AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  
  mgr->AddTask(ana);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histosEMCALP0v2", 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1 );
   
  return ana;
}
