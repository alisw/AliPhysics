void runOnNano(const char* filePath = "AliAOD.NanoAOD.root")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("train");

  AliAODInputHandler* iH = new AliAODInputHandler();
  mgr->SetInputEventHandler(iH);

  AliAnalysisTaskSE* task = (AliAnalysisTaskSE*) gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWG/DevNanoAOD/tutorial/AddTaskSimple.C");
  task->SelectCollisionCandidates(AliVEvent::kINT7);
  
  mgr->InitAnalysis();
  mgr->PrintStatus();

  // NOTE enable this for strict mode in which you get an AliFatal also for fields in the header which are accessed but not available
  // AliNanoAODHeader::SetFatalMode();
  
  // Create chain of input files
  TChain * chain = new TChain("aodTree");
  chain->Add(filePath);

  TStopwatch watch;
  mgr->StartAnalysis("local", chain, 1000);
  watch.Print();
}
