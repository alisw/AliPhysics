AliAnalysisTask *AddTaskHFECal()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHFEECal", "No analysis manager found.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFECal", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type=="AOD"){
    ::Error("AddTaskHFECal", "The tasks exits because AODs are in input");
    return NULL;
  }
  Bool_t MCthere=kFALSE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }else{
    MCthere=kTRUE;
  }
  
  
  //analysis task 
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFECal.C");

  AliAnalysisTaskHFECal *hfetaskCent = ConfigHFECal(MCthere);
  AliAnalysisTaskHFECal *hfetaskTrig= ConfigHFECal(MCthere);
 
  mgr->AddTask(hfetaskCent);
  mgr->AddTask(hfetaskTrig);
  
  // Semi-central trigger
  hfetaskCent->SelectCollisionCandidates(AliVEvent::kCentral);
  
  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWGHF_hfeCalCentral";
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("HFE_Results_EMCalCentral", TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
  mgr->ConnectInput(hfetaskCent, 0, cinput);
  mgr->ConnectOutput(hfetaskCent, 1, coutput1);
  
  //L1 gamma trigger
  hfetaskTrig->SelectCollisionCandidates(AliVEvent::kEMCEGA);
  
  TString containerName2 = mgr->GetCommonFileName();
  containerName2 += ":PWGHF_hfeCalTrigEGA";
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("HFE_Results_EMCalTrigEGA", TList::Class(),AliAnalysisManager::kOutputContainer, containerName2.Data());
  mgr->ConnectInput(hfetaskTrig, 0, cinput);
  mgr->ConnectOutput(hfetaskTrig, 1, coutput1);
  
  return NULL;
}
