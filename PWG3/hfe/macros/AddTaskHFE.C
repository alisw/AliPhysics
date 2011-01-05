AliAnalysisTask *AddTaskHFE(){
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
  if (type=="AOD"){
    ::Error("AddTaskHFE", "The tasks exits because AODs are in input");
    return NULL;
  }
  Bool_t MCthere=kTRUE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }

  //============= Set Task Name ===================
  //TString taskName=("AliAnalysisTaskHFE.cxx+");
  //===============================================
  
  gROOT->LoadMacro(Form("%s/PWG3/hfe/macros/ConfigHFEstandard.C", gSystem->Getenv("ALICE_ROOT")));
  //gROOT->LoadMacro(Form("%s/PWG3/hfe/macros/ConfigHFEtrd.C", gSystem->Getenv("ALICE_ROOT")));

  AliAnalysisTaskHFE *hfetask = ConfigHFEstandard(MCthere);
  mgr->AddTask(hfetask);

  //----------------------
  //create data containers
  //----------------------
 
  //find input container
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWG3_hfe";
  
  hfetask->ConnectOutput(1, mgr->CreateContainer("HFE_Results", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  hfetask->ConnectOutput(2, mgr->CreateContainer("HFE_QA", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  mgr->ConnectInput  (hfetask,  0, cinput );
  
/*
  AliAnalysisTaskHFE *trdtask = ConfigHFEtrd(MCthere);

  //----------------------
  //create data containers
  //----------------------
 
  trdtask->ConnectOutput(1, mgr->CreateContainer("HFE_Results", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  trdtask->ConnectOutput(2, mgr->CreateContainer("HFE_QA", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  mgr->ConnectInput  (trdtask,  0, cinput );
*/

  return hfetask;
}
