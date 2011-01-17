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
Bool_t AddAnalysisTaskRsn
(
  const char *options,
  const char *configs = "RsnConfig2010PhiFcnRealisticNoSA.C \
                         RsnConfig2010PhiFcnRealisticWithSA.C \
                         RsnConfig2010PhiFcnRealisticDipNoSA.C \
                         RsnConfig2010PhiFcnRealisticDipWithSA.C"
)
{
  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // interpret config string
  TString strDataLabel(options);
  Bool_t isSim = strDataLabel.Contains("sim");
  Bool_t isMC  = strDataLabel.Contains("MC");

  // initialize task with all available slots, even if not all of them will be used:
  AliRsnAnalysisSE *task = new AliRsnAnalysisSE("RsnAnalysis");
  task->SetZeroEventPercentWarning(100.0);
  task->SelectCollisionCandidates();
  if (isMC) task->SetMCOnly(kTRUE);

  // if not MC kinematics, set cuts for events : primary vertex range and type
  if (!isMC)
  {
    AliRsnCutPrimaryVertex *cutVertex   = new AliRsnCutPrimaryVertex("cutVertex", 10.0, 0, kFALSE);
    task->GetEventCuts()->AddCut(cutVertex);
    task->GetEventCuts()->SetCutScheme("cutVertex");
  }

  // add the task to manager
  mgr->AddTask(task);

  // load and execute configuration macroes
  TString    sList(configs);
  TObjArray *list = sList.Tokenize(" ");
  Int_t nConfig = list->GetEntries();
  Int_t iConfig = 0;
  for (iConfig = 0; iConfig < nConfig; iConfig++)
  {
    TObjString *ostr = (TObjString*)list->At(iConfig);
    cout << "***** Processing config macro '" << ostr->GetString().Data() << endl;
    gROOT->ProcessLine(Form(".x %s(\"%s\",\"%s\")", ostr->GetString().Data(), task->GetName(), options));
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
