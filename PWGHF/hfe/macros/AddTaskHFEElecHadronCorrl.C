AliAnalysisTask AddTaskHFEElecHadronCorrl()
{
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHFEElecHadronCorrl", "No analysis manager found.");
    return 0;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFEElecHadronCorrl", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type=="AOD"){
    ::Error("AddTaskHFEElecHadronCorrl", "The tasks exits because AODs are in input");
    return NULL;
  }
  Bool_t MCthere=kTRUE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }

  //analysis task 
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/AliAnalysisTaskElecHadronCorrel.cxx++g");
  //gROOT->LoadMacro("ConfigHFEemcalMod.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFEElecHadronCorl.C");
  AliAnalysisTaskElecHadronCorrel *taskHFE = ConfigHFEElecHadronCorl(MCthere);

  // output list of histos
  //TString foutputName = "ElecHadronCorrelAna.root";
  //AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("ccontainer0",TList::Class(),AliAnalysisManager::kOutputContainer,foutputName.Data());

    TString containerName = mgr->GetCommonFileName();
	containerName += ":PWGHF_hfeCalPbPbeh";
	AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("ccontainer0",TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data());
	
  mgr->ConnectInput(taskHFE,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskHFE,1,coutput3);

  mgr->AddTask(taskHFE);

  return taskHFE;
}
