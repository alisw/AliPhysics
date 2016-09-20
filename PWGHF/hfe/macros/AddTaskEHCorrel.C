AliAnalysisTask *AddTaskEHCorrel(Double_t centMin=0, Double_t centMax=20)
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

  AliAnalysisTaskEHCorrel *taskHFEeh = new AliAnalysisTaskEHCorrel("eh");
  taskHFEeh->SelectCollisionCandidates(AliVEvent::kINT7);
  taskHFEeh->SetCentralitySelection(centMin,centMax);
    
  TString containerName = mgr->GetCommonFileName();
  TString SubcontainerName = "EH_PbPb_INT7";
  //SubcontainerName += ContNameExt;
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(SubcontainerName,TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data());

  mgr->ConnectInput(taskHFEeh,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskHFEeh,1,coutput3);

  mgr->AddTask(taskHFEeh);
    
    
    // EMCal EGA EG1
    AliAnalysisTaskEHCorrel *taskHFEehGA01 = new AliAnalysisTaskEHCorrel("eh");
    taskHFEehGA01->SelectCollisionCandidates(AliVEvent::kEMCEGA);
    taskHFEehGA01->SetEMCalTriggerEG1(kTRUE);
    taskHFEehGA01->SetCentralitySelection(centMin,centMax);
    
    TString containerName01 = mgr->GetCommonFileName();
   // containerName01 += "EH_PbPb_GA1";
   // TString SubcontainerName01 = Form("HFEemcQATrigGAEG1_%s",calib);
     TString SubcontainerName01 = "EH_PbPb_GA1";
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName01, TList::Class(),AliAnalysisManager::kOutputContainer, containerName01.Data());
    
    mgr->ConnectInput(taskHFEehGA01, 0, cinput);
    mgr->ConnectOutput(taskHFEehGA01, 1, coutput1);
    mgr->AddTask(taskHFEehGA01);
    
  return taskHFEeh;
}
