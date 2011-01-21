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
Bool_t AddRsnAnalysisPhi
(
  const char *options     = "", 
  const char *taskName    = "RsnAnalysis",
  const char *configMacro = "RsnConfigPhi",
  const char *configPath  = "$(ALICE_INSTALL)/PWG2/RESONANCES/macros/train/LHC2010-pp7TeV"
)
{
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    
  // create the task and connect with physics selection
  AliRsnAnalysisSE *task = new AliRsnAnalysisSE(taskName);
  task->SetZeroEventPercentWarning(100.0);
  task->SelectCollisionCandidates();

  // add the task to manager
  mgr->AddTask(task);
  
  // execute the related config with settings for adding and not adding ITS-SA
  gROOT->ProcessLine(Form(".x %s/%s(%s,%s,%s)", configPath, configMacro, taskName, options                , configPath));
  gROOT->ProcessLine(Form(".x %s/%s(%s,%s,%s)", configPath, configMacro, taskName, Form("its+%s", options), configPath));
  
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
