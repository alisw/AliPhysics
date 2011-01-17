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
Bool_t AddRsnAnalysis
(
  const char *options,
  const char *configs = "RsnConfigNoSA.C RsnConfigSA.C",
  const char *path    = "$(ALICE_ROOT)/PWG2/RESONANCES/macros/train/LHC2010-7TeV-phi"
)
{
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    
  // create the task and connect with physics selection
  AliRsnAnalysisSE *task = new AliRsnAnalysisSE("RsnAnalysis");
  task->SetZeroEventPercentWarning(100.0);
  task->SelectCollisionCandidates();

  // add the task to manager
  mgr->AddTask(task);

  // load and execute all required configuration macroes in the string (arg #2)
  TString    sList   = configs;
  TObjArray *list    = sList.Tokenize(" ");
  Int_t      nConfig = list->GetEntries();
  Int_t      iConfig = 0;
  for (iConfig = 0; iConfig < nConfig; iConfig++)
  {
    TObjString *ostr = (TObjString*)list->At(iConfig);
    
    // the config macro is assumed to be stored in the path in argument #3
    // and to have three arguments: task name, a free string of options and the path where it is stored
    // --> all of them is a string, and then it must be passed with the quote marks
    const char *macro     = ostr->GetString().Data();
    const char *argName   = Form("\"%s\"", task->GetName());
    const char *argOption = Form("\"%s\"", options);
    const char *argPath   = Form("\"%s\"", path);
    gROOT->ProcessLine(Form(".x %s/%s(%s,%s,%s)", path, macro, argName, argOption, argPath));
  }

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
