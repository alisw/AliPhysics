AliAnalysisTask *AddTaskHFEElecV2()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHFEElecV2", "No analysis manager found.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFEElecV2", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type=="AOD"){
    ::Error("AddTaskHFEElecV2", "The tasks exits because AODs are in input");
    return NULL;
  }
  Bool_t MCthere=kFALSE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }else{
    MCthere=kTRUE;
  }
  
  
  //Event plane task
  AliEPSelectionTask *eventplaneTask = new AliEPSelectionTask("EventplaneSelection");
  eventplaneTask->SelectCollisionCandidates(AliVEvent::kSemiCentral | AliVEvent::kEMCEGA);

  eventplaneTask->SetTrackType("TPC");
  eventplaneTask->SetUsePtWeight();
  eventplaneTask->SetUsePhiWeight();
  eventplaneTask->SetSaveTrackContribution();
  
  AliESDtrackCuts* epTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "Standard");
  epTrackCuts->SetPtRange(0.1, 4);
  eventplaneTask->SetPersonalESDtrackCuts(epTrackCuts);

  mgr->AddTask(eventplaneTask);

  TString containerName3 = mgr->GetCommonFileName();
  containerName3 += ":PWGHF_hfeCalEventPlane";
  
  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("EPStat",TList::Class(), AliAnalysisManager::kOutputContainer,containerName3.Data());
  mgr->ConnectInput(eventplaneTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(eventplaneTask,1,coutput1);

  //analysis task 
//   gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/AliAnalysisTaskElecV2.cxx++g");
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFEElecV2.C");

  AliAnalysisTaskElecV2 *taskMB = ConfigHFEElecV2(MCthere);
  AliAnalysisTaskElecV2 *taskTR = ConfigHFEElecV2(MCthere);
 
  mgr->AddTask(taskMB);
  mgr->AddTask(taskTR);
  
  // Semi-central trigger
  taskMB->SelectCollisionCandidates(AliVEvent::kSemiCentral);
  
  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWGHF_hfeCalSemiCentralV2";
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histMB", TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
  mgr->ConnectInput(taskMB, 0, cinput);
  mgr->ConnectOutput(taskMB, 1, coutput1);
  
  //L1 gamma trigger
  taskTR->SelectCollisionCandidates(AliVEvent::kEMCEGA);
  
  TString containerName2 = mgr->GetCommonFileName();
  containerName2 += ":PWGHF_hfeCalL1GammaV2";
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histTR", TList::Class(),AliAnalysisManager::kOutputContainer, containerName2.Data());
  mgr->ConnectInput(taskTR, 0, cinput);
  mgr->ConnectOutput(taskTR, 1, coutput1);
  
  if(MCthere){
    
    AliAnalysisTaskElecV2 *taskMC = ConfigHFEElecV2(MCthere);
    mgr->AddTask(taskMC);
    
    taskMC->SelectCollisionCandidates(AliVEvent::kMB);
    
    TString containerName3 = mgr->GetCommonFileName();
    containerName3 += ":PWGHF_hfeCalMCV2";
    
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histMC", TList::Class(),AliAnalysisManager::kOutputContainer, containerName3.Data());
    mgr->ConnectInput(taskMC, 0, cinput);
    mgr->ConnectOutput(taskMC, 1, coutput1);
  }
  
  
  return NULL;
}
