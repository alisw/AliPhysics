AliAnalysisTask AddTaskHFEElecHadronCorrlPbPb()
{
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHFEElecHadronCorrlPbPb", "No analysis manager found.");
    return 0;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFEElecHadronCorrlPbPb", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
/*  if (type=="AOD"){
    ::Error("AddTaskHFEElecHadronCorrlPbPb", "The tasks exits because AODs are in input");
    return NULL;
  }
  */
  Bool_t MCthere=kTRUE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }

  //analysis task 
  gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/configs/PbPb/ConfigHFEElecHadronCorrelPbPb.C");
  AliAnalysisTaskElecHadronCorrel *taskHFEeh = ConfigHFEElecHadronCorrelPbPb(MCthere);

  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWGHF_hfeCalPbPbeh";
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("ccontainer0",TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data());

  mgr->ConnectInput(taskHFEeh,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskHFEeh,1,coutput3);

  mgr->AddTask(taskHFEeh);

  return taskHFEeh;
}
