void validateNano(TString runOn)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("train");
  
  if (runOn == "Nano") {
    AliAODInputHandler* iH = new AliAODInputHandler();
    mgr->SetInputEventHandler(iH);

    // Create chain of input files
    TChain * chain = new TChain("aodTree");
    chain->Add("AliAOD.NanoAOD.root");
  } else if (runOn == "AOD") {
    AliAODInputHandler* iH = new AliAODInputHandler();
    mgr->SetInputEventHandler(iH);

    // Multiplicity selection
    gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    
    // Create chain of input files
    TChain * chain = new TChain("aodTree");
    chain->Add("AliAOD.root");
  } else if (runOn == "ESD") {
    AliESDInputHandler* iH = new AliESDInputHandler();
    mgr->SetInputEventHandler(iH);

    // Physics selection
    gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    
    // Multiplicity selection
    gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");

    // Input files
    TChain * chain = new TChain("esdTree");
    chain->Add("AliESDs.root");
  }

  AliAnalysisTaskNanoValidator* ana = new AliAnalysisTaskNanoValidator("AliAnalysisTaskNanoValidator");
  mgr->AddTask(ana);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histograms", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", AliAnalysisManager::GetCommonFileName(), "Simple"));
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1 );
  
  mgr->InitAnalysis();
  mgr->PrintStatus();
    
  TStopwatch watch;
  mgr->StartAnalysis("local", chain, 20);
  watch.Print();
}
