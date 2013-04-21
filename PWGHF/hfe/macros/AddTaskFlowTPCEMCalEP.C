AliAnalysisTask *AddTaskFlowTPCEMCalEP()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFlowTPCEMCalEP", "No analysis manager found.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFlowTPCEMCalEP", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type=="AOD"){
    ::Error("AddTaskFlowTPCEMCalEP", "The tasks exits because AODs are in input");
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
  eventplaneTask->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kSemiCentral | AliVEvent::kCentral | AliVEvent::kEMCEGA | AliVEvent::kEMCEJE);

  eventplaneTask->SetTrackType("TPC");
  eventplaneTask->SetUsePtWeight();
  eventplaneTask->SetUsePhiWeight();
  eventplaneTask->SetSaveTrackContribution();

  mgr->AddTask(eventplaneTask);

  TString containerName3 = mgr->GetCommonFileName();
  containerName3 += ":PWGHF_hfeCalEventPlane";
  
  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("EPStat",TList::Class(), AliAnalysisManager::kOutputContainer,containerName3.Data());
  mgr->ConnectInput(eventplaneTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(eventplaneTask,1,coutput1);

  //analysis task 
//   gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/AliAnalysisTaskFlowTPCEMCalEP.cxx++g");
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFE_FLOW_TPCEMCal_EP.C");

  AliAnalysisTaskFlowTPCEMCalEP *taskMB = ConfigHFE_FLOW_TPCEMCal_EP(MCthere);
  AliAnalysisTaskFlowTPCEMCalEP *taskTR = ConfigHFE_FLOW_TPCEMCal_EP(MCthere);
 
  mgr->AddTask(taskMB);
  mgr->AddTask(taskTR);
  
  // Central trigger
  taskMB->SelectCollisionCandidates(AliVEvent::kSemiCentral);
  
  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWGHF_hfeCalCentralV2";
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histMB", TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
  mgr->ConnectInput(taskMB, 0, cinput);
  mgr->ConnectOutput(taskMB, 1, coutput1);
  
  //L1 gamma and jet trigger
  taskTR->SelectCollisionCandidates(AliVEvent::kEMCEGA | AliVEvent::kEMCEJE);
  
  TString containerName2 = mgr->GetCommonFileName();
  containerName2 += ":PWGHF_hfeCalL1GammaV2";
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histTR", TList::Class(),AliAnalysisManager::kOutputContainer, containerName2.Data());
  mgr->ConnectInput(taskTR, 0, cinput);
  mgr->ConnectOutput(taskTR, 1, coutput1);
  
  if(MCthere){
    
    AliAnalysisTaskFlowTPCEMCalEP *taskMC = ConfigHFE_FLOW_TPCEMCal_EP(MCthere);
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
