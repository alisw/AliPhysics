void runRsnAnalysisSE(AliLog::EType_t type=AliLog::kInfo,Bool_t useKine = kTRUE)
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cin;
  AliAnalysisDataContainer *cout;

  cin = SetCorrectHandlersAndReturnInput(mgr,useKine);

  AliRsnAnalysisSE *task = new AliRsnAnalysisSE("AliRsnAnalysisSE");
  task->SetLogType(type,"AliRsnAnalysisManager:AliRsnPairManager:AliRsnPairManager:AliRsnPair");
  task->SetPrintInfoNumber(10);

  task->SetPriorProbability(AliRsnDaughter::kElectron, 0.02);
  task->SetPriorProbability(AliRsnDaughter::kMuon,     0.02);
  task->SetPriorProbability(AliRsnDaughter::kPion,     0.83);
  task->SetPriorProbability(AliRsnDaughter::kKaon,     0.07);
  task->SetPriorProbability(AliRsnDaughter::kProton,   0.06);
  task->DumpPriors();

  AliRsnAnalysisManager *analMgr = task->GetAnalysisManager("MyAnalysisSE");
  Int_t i=0;
  Int_t analNum=1;
  for (i=0;i<analNum;i++) {
//     analMgr->AddConfig("RsnConfigTest.C",Form("PHI%d",i),);
    analMgr->AddConfig("RsnConfig.C",Form("PHI%d",i),"RsnConfig_PHI");
    analMgr->AddConfig("RsnConfig.C",Form("KSTAR%d",i),"RsnConfig_KSTAR");
  }

  AliAnalysisDataContainer *output = mgr->CreateContainer("RSNSE", TList::Class(), AliAnalysisManager::kOutputContainer, "RSNAnalysis.root");

  mgr->AddTask(task);

  mgr->ConnectInput(task, 0, cin);
  mgr->ConnectOutput(task, 1, output);

}

AliAnalysisDataContainer *SetCorrectHandlersAndReturnInput(AliAnalysisManager *mgr,Bool_t useKine = kTRUE)
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliAnalysisDataContainer *cin;

  if (useKine) {
    AliMCEventHandler* mcInputHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!mcInputHandler) {
      Info("","Creating mcInputHandler ...");
      AliMCEventHandler* mcInputHandler = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mcInputHandler);
    }
  }

  AliESDInputHandler *esdInputHandler = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!esdInputHandler) {
    Info("","Creating esdInputHandler ...");
    esdInputHandler = new AliESDInputHandler();
    mgr->SetInputEventHandler  (esdInputHandler);
    cin = mgr->GetCommonInputContainer();
  } else {
    cin = mgr->GetCommonOutputContainer();
  }

  return cin;
}
