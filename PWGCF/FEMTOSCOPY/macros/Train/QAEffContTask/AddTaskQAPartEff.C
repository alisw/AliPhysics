AliAnalysisTaskParticleEfficiency *AddTaskQAPartEff(const char* outfilename="AnalysisResults.root")
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhiCorrelations", "No analysis manager to connect to.");
    return NULL;
  }  

  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFemto", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  cout << "Found " <<type << " event handler" << endl;



  // Create the task and configure it.
  //===========================================================================
  AliAnalysisTaskParticleEfficiency* ana = new  AliAnalysisTaskParticleEfficiency("MyTask");

  if (!outfilename)
    outfilename = AliAnalysisManager::GetCommonFileName();

  //AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("MyList", TList::Class(),AliAnalysisManager::kOutputContainer,outfilename);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("MyList", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));

 mgr->AddTask(ana);
  

  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1 );
 
  return ana;

}
