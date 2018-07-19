AliAnalysisPHOSResonances* AddTaskPHOSResonances(UInt_t offlineTriggerMask = AliVEvent::kAnyINT)
{
  //Add a task AliAnalysisTaskPHOSResonanses to the analysis train
  //Author: Dmitri Peresunko

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AliAnalysisPHOSResonances", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AliAnalysisPHOSResonances", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisPHOSResonances* task = new AliAnalysisPHOSResonances("PHOSResonanses");
					   
  task->SelectCollisionCandidates(offlineTriggerMask);
   
  //{gg,gpp,gpm,gppDCA,gpmDCA,gpi0,gpi0Merg,pi0pi0,pi0pi0M,pi0pp,pi0pm,pi0ppDCA,pi0pmDCA,gLam, pi0Lam,ee,gpippim,pi0pippim} ;
  int ch[18]={1,0,0,0,0,      1,   1,       1,     1,      1,    1,    1,       1,       1,    1,     0,  1,     1 } ;
  
  task->SetListOfChannels(ch,18) ;
 
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  TString cname(Form("PHOSResonancesCoutput1"));
  TString pname(Form("%s:PHOSResonanses", AliAnalysisManager::GetCommonFileName()));
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(cname.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, pname.Data());
  mgr->ConnectOutput(task, 1, coutput1);
  
  return task;
}
