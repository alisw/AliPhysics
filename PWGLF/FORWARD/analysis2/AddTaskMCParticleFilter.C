/** 
 * Add a task to storing MC particles 
 * 
 * 
 * @return Task 
 */
AliAnalysisTask*
AddTaskMCParticleFilter() 
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) Fatal("", "No analysis manager to connect to.");

  AliAnalysisTaskMCParticleFilter* mctask = 
    new AliAnalysisTaskMCParticleFilter("mcfilter");
  mgr->AddTask(mctask);
  AliAnalysisDataContainer* histOut = 
    mgr->CreateContainer("mcfilter", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 "mcParticles.root");
  mgr->ConnectInput(mctask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(mctask, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(mctask, 1, histOut);

  return mctask;
}
/*
 * EOF
 */
