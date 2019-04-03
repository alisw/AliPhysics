AliAnalysisTaskPHOSEmbeddedDiffObjectCreator* AddTaskPHOSEmbeddedDiffObjectCreator(
    const char* name     = "PHOSEmbeddedDiffObjectCreator",
    const TString parname = "Pi0",
    const UInt_t trigger = AliVEvent::kINT7|AliVEvent::kPHI7,
    const Bool_t usePHOSTender = kTRUE,
    const Bool_t isMC = kFALSE,
    const Double_t BunchSpace = 100.,
    const Bool_t excludeM4     = kTRUE,
    const TString period       = "LHC15o"
    )
{
  //Add a task AliAnalysisTaskPHOSObjectCreator to the analysis train
  //Author: Daiki Sekihata
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSEmbeddedDiffObjectCreator", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSEmbeddedDiffObjectCreator", "This task requires an input event handler");
    return NULL;
  }

  TString taskname = Form("%s_%s_BS%dns",name,parname.Data(),(Int_t)BunchSpace);

  AliAnalysisTaskPHOSEmbeddedDiffObjectCreator* task = new AliAnalysisTaskPHOSEmbeddedDiffObjectCreator(taskname);
  task->SelectCollisionCandidates(trigger);
  task->SetTenderFlag(usePHOSTender);//always switch OFF tender because tender does not support embedded clusters. do it by myself
  task->SetMCFlag(isMC);
  task->SetBunchSpace(BunchSpace);//in unit of ns
  task->SetEmbeddedParticle(parname);
  task->SetEmbedding(kTRUE);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  const TString listname = Form("hist_%s",taskname.Data());
  TString outputFile = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("%s",listname.Data()), THashList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",outputFile.Data(),"PWGGA_PHOSTasks_PHOSRun2"));
  mgr->ConnectOutput(task, 1, coutput1);

  return task;
}

