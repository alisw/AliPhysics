AliAnalysisTask *AddTaskHFEemcQA(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHFE", "No analysis manager found.");
    return NULL;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFE", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
 Bool_t MCthere=kFALSE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }else{
    MCthere=kTRUE;
  }


  gROOT->LoadMacro("AliAnalysisTaskHFEemcQA.cxx++g");
  
  // +++ EMCal MB
  AliAnalysisTaskHFEemcQA *hfecalqa = new AliAnalysisTaskHFEemcQA("emcqa");
  mgr->AddTask(hfecalqa);
  hfecalqa->SelectCollisionCandidates(AliVEvent::kINT8);
  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWGHF_hfeHFEemcQAINT8";
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("HFEemcQAINT8", TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
  mgr->ConnectInput(hfecalqa, 0, cinput);
  mgr->ConnectOutput(hfecalqa, 1, coutput1); 
  /*
  mgr->ConnectInput  (hfecalqa,0,cinput);
  hfecalqa->ConnectOutput(1, mgr->CreateContainer("HFE_Results_EMCALINT8", TList::Class(),
						 AliAnalysisManager::kOutputContainer, containerName.Data()));
  hfecalqa->ConnectOutput(2, mgr->CreateContainer("HFE_QA_EMCALINT8", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  */

  // EMCal EGA
  AliAnalysisTaskHFEemcQA *hfecalqaTrig0 = new AliAnalysisTaskHFEemcQA("emcqa");
  mgr->AddTask(hfecalqaTrig0);
  hfecalqaTrig0->SelectCollisionCandidates(AliVEvent::kEMCEGA);
  TString containerName1 = mgr->GetCommonFileName();
  containerName1 += ":PWGHF_hfeHFEemcQATrigGA";
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("HFEemcQATrigGA", TList::Class(),AliAnalysisManager::kOutputContainer, containerName1.Data());
  mgr->ConnectInput(hfecalqaTrig0, 0, cinput);
  mgr->ConnectOutput(hfecalqaTrig0, 1, coutput1); 

 
  // EMCal EJE
  AliAnalysisTaskHFEemcQA *hfecalqaTrig1 = new AliAnalysisTaskHFEemcQA("emcqa");
  mgr->AddTask(hfecalqaTrig1);
  hfecalqaTrig1->SelectCollisionCandidates(AliVEvent::kEMCEJE);
  TString containerName2 = mgr->GetCommonFileName();
  containerName1 += ":PWGHF_hfeHFEemcQATrigJE";
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("HFEemcQATrigJE", TList::Class(),AliAnalysisManager::kOutputContainer, containerName2.Data());
  mgr->ConnectInput(hfecalqaTrig1, 0, cinput);
  mgr->ConnectOutput(hfecalqaTrig1, 1, coutput1); 

  //return hfecalqa;
  return NULL;
}
