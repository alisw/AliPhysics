AliProtonFeedDownAnalysisTask* AddTaskProtonFeedDownAnalysis(const char *analysisType="Hybrid",const char *pidMode="Bayesian")
{
  // Creates a proton analysis task and adds it to the analysis manager.
  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskProtons", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskProtons", "This task requires an input event handler");
    return NULL;
  }   
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/SPECTRA/macros/configProtonFeedDownAnalysis.C");
  AliProtonFeedDownAnalysis *pa = 0;
  if (type=="ESD") pa = GetProtonFeedDownAnalysisObject("ESD", analysisType, pidMode);
  else if (type=="AOD") pa = GetProtonFeedDownAnalysisObject("AOD", analysisType, pidMode);
  else return NULL;

  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliProtonFeedDownAnalysisTask *taskproton = new  AliProtonFeedDownAnalysisTask("TaskProtons");
  mgr->AddTask(taskproton);
  taskproton->SetAnalysisObject(pa);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWG2BaryonRatio_outputAliProtonFeedDownAnalysisTask.root";
  AliAnalysisDataContainer *cout_proton = mgr->CreateContainer("protonFeedDown", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(taskproton, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskproton, 0, cout_proton);
  
  // Return task pointer at the end
  return taskproton;
}
