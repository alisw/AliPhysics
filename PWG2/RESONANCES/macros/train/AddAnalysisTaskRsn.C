//
// Macro to create the full analysis manager for Resonances
//
Bool_t AddAnalysisTaskRsn
(
  Bool_t sourceESD  = kTRUE // if true, the source of data is ESD, otherwise is AOD from filter task
)
{
  AliLog::SetClassDebugLevel("AliMCEvent", AliLog::kWarning);
  //AliLog::SetClassDebugLevel("AliRsnCutStd", AliLog::kDebug+3);
  //AliLog::SetClassDebugLevel("AliRsnCutBetheBloch", AliLog::kDebug+3);

  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // initialize task with 4 slots:
  AliRsnAnalysisSE *task = new AliRsnAnalysisSE("RsnAnalysis", 2);
  
  // set prior probabilities for PID
  task->SetPriorProbability(AliPID::kElectron, 0.02);
  task->SetPriorProbability(AliPID::kMuon,     0.02);
  task->SetPriorProbability(AliPID::kPion,     0.83);
  task->SetPriorProbability(AliPID::kKaon,     0.07);
  task->SetPriorProbability(AliPID::kProton,   0.06);
  task->DumpPriors();

  // load config macros
  gROOT->LoadMacro("RsnConfig.C");

  // initialize analysis manager with pairs from config
  AliRsnAnalysisManager *anaMgr = 0x0;

  // manager #0: phi
  anaMgr = task->GetAnalysisManager(0);
  anaMgr->Add(RsnConfig("PHI_NOPID", 333, 1019.3, AliPID::kKaon, AliPID::kKaon));
  anaMgr->Add(RsnConfig("PHI_BB"   , 333, 1019.3, AliPID::kKaon, AliPID::kKaon));
  anaMgr->Add(RsnConfig("PHI_PID"  , 333, 1019.3, AliPID::kKaon, AliPID::kKaon));
  // manager #1: kstar
  anaMgr = task->GetAnalysisManager(1);
  anaMgr->Add(RsnConfig("KSTAR_NOPID", 313, 896.0, AliPID::kPion, AliPID::kKaon));
  anaMgr->Add(RsnConfig("KSTAR_BB"   , 313, 896.0, AliPID::kPion, AliPID::kKaon));
  anaMgr->Add(RsnConfig("KSTAR_PID"  , 313, 896.0, AliPID::kPion, AliPID::kKaon));

  // setup cuts for events (good primary vertex)
  AliRsnCutPrimaryVertex *cutVertex = new AliRsnCutPrimaryVertex("cutVertex", 3);
  AliRsnCutSet *cutSetEvent = new AliRsnCutSet("eventCuts");
  cutSetEvent->AddCut(cutVertex);
  cutSetEvent->SetCutScheme("cutVertex");
  task->SetEventCuts(cutSetEvent);

  // add the task to manager
  mgr->AddTask(task);

  // connect input container according to source choice
  if (sourceESD) mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  else mgr->ConnectInput(task, 0, mgr->GetCommonOutputContainer());

  // create paths for the output in the common file
  Char_t infoPath[500], phiPath[500], kstarPath[500];
  sprintf(infoPath , "%s:PWG2RSNINFO" , AliAnalysisManager::GetCommonFileName());
  sprintf(phiPath  , "%s:PWG2RSNPHI"  , AliAnalysisManager::GetCommonFileName());
  sprintf(kstarPath, "%s:PWG2RSNKSTAR", AliAnalysisManager::GetCommonFileName());

  // create containers for output
  AliAnalysisDataContainer *outputInfo = mgr->CreateContainer("RsnInfo", TList::Class(), AliAnalysisManager::kOutputContainer, infoPath);
  AliAnalysisDataContainer *outputRsn[2];
  outputRsn[0] = mgr->CreateContainer("PHI"  , TList::Class(), AliAnalysisManager::kOutputContainer, phiPath);
  outputRsn[1] = mgr->CreateContainer("KSTAR", TList::Class(), AliAnalysisManager::kOutputContainer, kstarPath);

  mgr->ConnectOutput(task, 1, outputInfo);
  mgr->ConnectOutput(task, 2, outputRsn[0]);
  mgr->ConnectOutput(task, 3, outputRsn[1]);
}
