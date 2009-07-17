//
// Macro to create the full analysis manager for Resonances
//
Bool_t AddAnalysisTaskRsn
(
  Bool_t sourceESD  = kTRUE         // if true, the source of data is ESD, otherwise is AOD from filter task
)
{
  //AliLog::SetClassDebugLevel("AliRsnCut", AliLog::kDebug+3);
  //AliLog::SetClassDebugLevel("AliRsnCutStd", AliLog::kDebug+3);
  //AliLog::SetClassDebugLevel("AliRsnCutBetheBloch", AliLog::kDebug+3);

  // retrieve analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // initialize task with 4 slots:
  AliRsnAnalysisSE *task = new AliRsnAnalysisSE("AliRsnAnalysisSE", 4);
  
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

  // manager #0: phi NO PID
  anaMgr = task->GetAnalysisManager(0);
  anaMgr->Add(RsnConfig(AliRsnPair::kNoPID, "PHI"   , 333, AliPID::kKaon, AliPID::kKaon));
  anaMgr->Add(RsnConfig(AliRsnPair::kNoPID, "PHI_BB", 333, AliPID::kKaon, AliPID::kKaon, 0.2));
  // manager #1: phi with PID (realistic & perfect)
  anaMgr = task->GetAnalysisManager(1);
  anaMgr->Add(RsnConfig(AliRsnPair::kRealisticPID, "PHI", 333, AliPID::kKaon, AliPID::kKaon));
  anaMgr->Add(RsnConfig(AliRsnPair::kPerfectPID  , "PHI"  , 333, AliPID::kKaon, AliPID::kKaon));
  // manager #2: kstar NO PID
  anaMgr = task->GetAnalysisManager(2);
  anaMgr->Add(RsnConfig(AliRsnPair::kNoPID, "KSTAR"   , 313, AliPID::kPion, AliPID::kKaon));
  anaMgr->Add(RsnConfig(AliRsnPair::kNoPID, "KSTAR_BB", 313, AliPID::kPion, AliPID::kKaon, 0.2));
  // manager #1: kstar with PID (realistic & perfect)
  anaMgr = task->GetAnalysisManager(3);
  anaMgr->Add(RsnConfig(AliRsnPair::kRealisticPID, "KSTAR", 313, AliPID::kPion, AliPID::kKaon));
  anaMgr->Add(RsnConfig(AliRsnPair::kPerfectPID  , "KSTAR"  , 313, AliPID::kPion, AliPID::kKaon));

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

  // define and connect output containers
  AliAnalysisDataContainer *outputInfo = mgr->CreateContainer("RsnInfo", TList::Class(), AliAnalysisManager::kOutputContainer, "rsn_info.root");
  AliAnalysisDataContainer *outputRsn[4];
  outputRsn[0] = mgr->CreateContainer("PHI_NOPID"  , TList::Class(), AliAnalysisManager::kOutputContainer, "phi_nopid.root");
  outputRsn[1] = mgr->CreateContainer("PHI_PID"    , TList::Class(), AliAnalysisManager::kOutputContainer, "phi_pid.root");
  outputRsn[2] = mgr->CreateContainer("KSTAR_NOPID", TList::Class(), AliAnalysisManager::kOutputContainer, "kstar_nopid.root");
  outputRsn[3] = mgr->CreateContainer("KSTAR_PID"  , TList::Class(), AliAnalysisManager::kOutputContainer, "kstar_pid.root");

  mgr->ConnectOutput(task, 1, outputInfo);
  mgr->ConnectOutput(task, 2, outputRsn[0]);
  mgr->ConnectOutput(task, 3, outputRsn[1]);
  mgr->ConnectOutput(task, 4, outputRsn[2]);
  mgr->ConnectOutput(task, 5, outputRsn[3]);
}
