//
// Macro to create the full analysis manager for Resonances
//
Bool_t AddAnalysisTaskRsnNew(const char *configMacro = "ConfigTaskRsn.C")
{
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // initialize task with all available slots, even if not all of them will be used:
  AliRsnAnalysisSE *task = new AliRsnAnalysisSE("RsnAnalysis", 1);
  task->SetZeroEventPercentWarning(100.0);

  // load and execute configuration macro
  // if we are in a PROOF environment, it must be loaded by gProof
  // otherwise if will be loaded by gROOT
  if (TProofMgr::GetListOfManagers()->GetEntries())
  {
    if (dynamic_cast<TProofLite *> gProof)
      gProof->Exec(Form(".L %s", macro.Data()));
    else
      gProof->Load(configMacro);
  }
  else
  {
    //gROOT->ProcessLine(Form(".x %s", configMacro));
    gROOT->LoadMacro(configMacro);
  }

  // configure the task
  RsnConfigTask(task);

  // add the task to manager
  mgr->AddTask(task);

  // connect input container according to source choice
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // create paths for the output in the common file
  Char_t commonPath[500];
  sprintf(commonPath, "%s:PWG2RSN", AliAnalysisManager::GetCommonFileName());

  // create containers for output
  AliAnalysisDataContainer *outputInfo = mgr->CreateContainer("RsnInfo", TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
  AliAnalysisDataContainer *outputHist = mgr->CreateContainer("RsnHist", TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
  mgr->ConnectOutput(task, 1, outputInfo);
  mgr->ConnectOutput(task, 2, outputHist);
}
