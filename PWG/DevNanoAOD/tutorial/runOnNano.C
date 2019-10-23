void runOnNano(const char* filePath = "AliAOD.NanoAOD.root")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("train");

  AliAODInputHandler* iH = new AliAODInputHandler();
  mgr->SetInputEventHandler(iH);

  
  gInterpreter->ExecuteMacro(
        "$ALICE_PHYSICS/PWGLF/RESONANCES/PostProcessing/Xi1530/NanoCheck/"
        "AddTaskNanoCheck.C");
  /*
  AliAnalysisTaskXi1530* task_ESD =
      (AliAnalysisTaskXi1530*)gInterpreter->ExecuteMacro(
          "$ALICE_PHYSICS/PWGLF/RESONANCES/PostProcessing/Xi1530/"
          "AddTaskXi1530.C");
          */
  mgr->InitAnalysis();
  mgr->PrintStatus();

  // NOTE enable this for strict mode in which you get an AliFatal also for fields in the header which are accessed but not available
  // AliNanoAODHeader::SetFatalMode();
  
  // Create chain of input files
  TChain * chain = new TChain("aodTree");
  chain->Add(filePath);

  TStopwatch watch;
  mgr->StartAnalysis("local", chain, 100000);
  watch.Print();
    gSystem->Exec("mv AnalysisResults.root AnalysisResultsNano.root");
}
