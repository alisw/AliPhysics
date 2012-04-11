AliAnalysisTask *AddTaskHFECalPbPb(int TPCclust){
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
  Bool_t MCthere=kFALSE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }else{
    MCthere=kTRUE;
  }
  cout<<"AddTaskHFE - MC config is: "<<MCthere<<endl;

  //============= Set Task Name ===================
  //TString taskName=("AliAnalysisTaskHFE.cxx+");
  //===============================================
  
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFECalstandard_PbPb.C");

  //<--- task1 for EMCal trigger
  AliAnalysisTaskHFE *hfetask = ConfigHFECalstandard_PbPb(MCthere,TPCclust);
  //RequestMemory(hfetask, 250*1024);
  hfetask->SelectCollisionCandidates(AliVEvent::kEMCEGA);
  mgr->AddTask(hfetask);

  //find input container
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWGHF_hfeCalPbPbEGA";
  
  hfetask->ConnectOutput(1, mgr->CreateContainer("HFE_Results_EMCALEGA", TList::Class(),
						 AliAnalysisManager::kOutputContainer, containerName.Data()));
  hfetask->ConnectOutput(2, mgr->CreateContainer("HFE_QA_EMCALEGA", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));

  mgr->ConnectInput  (hfetask,  0, cinput );
  

  //<--- task2 for central trigger

  AliAnalysisTaskHFE *hfetask2 = ConfigHFECalstandard_PbPb(MCthere,TPCclust);
  hfetask2->SelectCollisionCandidates(AliVEvent::kCentral);
  mgr->AddTask(hfetask2);

  //find input container
  TString containerName2 = mgr->GetCommonFileName();
  containerName2 += ":PWGHF_hfeCalPbPbCent";
  
  hfetask2->ConnectOutput(1, mgr->CreateContainer("HFE_Results_EMCALCent", TList::Class(),
						 AliAnalysisManager::kOutputContainer, containerName2.Data()));
  hfetask2->ConnectOutput(2, mgr->CreateContainer("HFE_QA_EMCALCent", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName2.Data()));

  mgr->ConnectInput  (hfetask2,  0, cinput );
  

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

  //return hfetask;
  return NULL;
}
