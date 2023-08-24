////////////////////////////////////////////////////////////
//                                                        //
//  Task for Non-HFE reconstruction efficiency            //
//  Non-Photonic Electron identified with Invariant mass  //
//                                                        //
//  Author: Deepa Thomas (University of Texas at Austin)  //
//          Vivek Kumar Singh (VECC)                      //
////////////////////////////////////////////////////////////

AliAnalysisTask *AddTaskIPResolBeautyppCal(const TString ContNameExt = "", Bool_t fSwitchRIP=kTRUE)
{
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AliAnalysisTaskIPResolBeautyppCal", "No analysis manager found.");
    return 0;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEHEffiCalc", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  Bool_t MCthere=kTRUE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }
    
  TString name("ehEffiTask");
  AliAnalysisTaskIPResolBeautyppCal *task = new AliAnalysisTaskIPResolBeautyppCal(name);
  task->SelectCollisionCandidates(AliVEvent::kINT7);
  task->SwitchRecalImpPar(fSwitchRIP);
    
  TString containerName = mgr->GetCommonFileName();
  TString SubcontainerName = ContNameExt;
  SubcontainerName += "IP_test";
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(SubcontainerName,TList::Class(),AliAnalysisManager::kOutputContainer,containerName.Data());

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput3);

  return task;
}
