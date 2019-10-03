AliAnalysisTaskGammaFlow* AddTaskPHOSGammaFlow (const char* name = "PHOSGammaFlow",
					    Int_t harmonics = 2,
					    UInt_t offlineTriggerMask = AliVEvent::kCentral)
{
  //Add a task AliAnalysisTaskGammaFlow to the analysis train
  //Author: Dmitri Peressounko
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSGammaFlow", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSGammaFlow", "This task requires an input event handler");
    return NULL;
  }

  if(harmonics!=2 && harmonics!=3){
    ::Error("AddTaskPHOSGammaFlow", Form("Only harmonics 2 and 3 allowed, you use %d",harmonics));
    return NULL;    
  }
  
  AliAnalysisTaskGammaFlow* task = new AliAnalysisTaskGammaFlow(Form("%sTask", name));
  task->SetDistCut(kFALSE) ;
  task->SetHarmonics(harmonics) ;  
  task->SelectCollisionCandidates(offlineTriggerMask);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 
  TString trigger;
  if(offlineTriggerMask==AliVEvent::kCentral) trigger = "Central" ;
  else if(offlineTriggerMask==AliVEvent::kSemiCentral) trigger = "SemiCentral" ;
  else trigger = "Other" ;
  TString cname(Form("%s%sv%d", name,trigger.Data(),harmonics));
  TString pname(Form("%s:%s", AliAnalysisManager::GetCommonFileName(), name));

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectInput(task , 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
