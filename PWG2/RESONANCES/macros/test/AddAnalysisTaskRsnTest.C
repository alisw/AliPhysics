//
// This macro serves to add the RSN analysis task to the steering macro.
//
// Inputs:
//   - dataLabel   = a string with informations about the type of data
//                   which could be needed to be ported to the config macro
//                   to set up some cuts
//   - configMacro = macro which configures the analysis; it has *ALWAYS*
//                   defined inside a function named 'RsnConfigTask()',
//                   whatever the name of the macro itself, whose first two
//                   arguments must have to be the task and the 'dataLabel' argument.
// 
Bool_t AddAnalysisTaskRsnTest
(const char *dataLabel, const char *configMacro = "ConfigTaskRsnTest.C")
{
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // initialize task with all available slots, even if not all of them will be used:
  AliRsnAnalysisSE *task = new AliRsnAnalysisSE("RsnAnalysis");
  task->SetZeroEventPercentWarning(100.0);
  //task->SelectCollisionCandidates();
  //task->SetMCOnly(kTRUE);

  // load and execute configuration macro
  gROOT->LoadMacro(configMacro);
  if (!RsnConfigTask(task, dataLabel)) return kFALSE;
  /*
  //                                  perf    its     expPID
  if (!RsnConfigTask(task, dataLabel, kTRUE , kTRUE , kFALSE)) return kFALSE;
  if (!RsnConfigTask(task, dataLabel, kTRUE , kFALSE, kFALSE)) return kFALSE;
  if (!RsnConfigTask(task, dataLabel, kFALSE, kTRUE , kFALSE)) return kFALSE;
  if (!RsnConfigTask(task, dataLabel, kFALSE, kFALSE, kFALSE)) return kFALSE;
  if (!RsnConfigTask(task, dataLabel, kFALSE, kTRUE , kTRUE )) return kFALSE;
  if (!RsnConfigTask(task, dataLabel, kFALSE, kFALSE, kTRUE )) return kFALSE;
  */

  // add the task to manager
  mgr->AddTask(task);

  // connect input container according to source choice
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // create paths for the output in the common file
  Char_t commonPath[500];
  sprintf(commonPath, "%s", AliAnalysisManager::GetCommonFileName());

  // create containers for output
  AliAnalysisDataContainer *outputInfo = mgr->CreateContainer("RsnInfo", TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
  AliAnalysisDataContainer *outputHist = mgr->CreateContainer("RsnHist", TList::Class(), AliAnalysisManager::kOutputContainer, commonPath);
  mgr->ConnectOutput(task, 1, outputInfo);
  mgr->ConnectOutput(task, 2, outputHist);

  return kTRUE;
}
