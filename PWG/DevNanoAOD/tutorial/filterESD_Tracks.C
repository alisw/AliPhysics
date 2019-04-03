void filterESD_Tracks()
{
  AliAnalysisManager* mgr = new AliAnalysisManager("NanoAOD Filter", "NanoAOD filter for nanoAOD production");
    
  AliESDInputHandler* iH = new AliESDInputHandler();
  mgr->SetInputEventHandler(iH);

  // Define aod output handler
  AliAODHandler* aodOutputHandler = new AliAODHandler();
  aodOutputHandler->SetOutputFileName("AliAOD.NanoAOD.root");
  mgr->SetOutputEventHandler(aodOutputHandler);
  
  // Physics selection
  gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  
  // Multiplicity selection
  gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  
  // OCDB
  gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/AddTaskConfigOCDB.C(\"raw://\")");
  
  // ESD filter
  gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C(kFALSE, kFALSE, kFALSE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE, 1500, 3, kTRUE, kFALSE, kFALSE, kFALSE)");
  
  AliAnalysisTaskNanoAODFilter* task = (AliAnalysisTaskNanoAODFilter*) gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWG/DevNanoAOD/AddTaskNanoAODFilter.C(0, kFALSE)");
  task->AddSetter(new AliNanoAODSimpleSetter);
  
  // Event selection
  // filter bit
  task->SelectCollisionCandidates(AliVEvent::kINT7);
  
  AliAnalysisNanoAODEventCuts* evtCuts = new AliAnalysisNanoAODEventCuts;
  evtCuts->SetVertexRange(8);
  evtCuts->SetCutPileUpMV(kTRUE);

  // Track selection
  AliAnalysisNanoAODTrackCuts* trkCuts = new AliAnalysisNanoAODTrackCuts;
  trkCuts->SetBitMask((1 << 8) | (1 << 9)); // hybrid 2011
  trkCuts->SetMaxEta(0.9);
  trkCuts->SetMinPt(1.0);
  
  // Fields to store
  // event level
  // Note: vertices are kept by default
  task->SetVarListHeader("OfflineTrigger,MagField,MultSelection.RefMult08");
  // track level
  task->SetVarListTrack("pt,theta,phi");

  task->SetTrkCuts(trkCuts);
  task->SetEvtCuts(evtCuts);

  mgr->SetDebugLevel(1); // enable debug printouts
  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();
  
  // Input files
  TChain * chain = new TChain("esdTree");
  chain->Add("AliESDs.root");

  Printf("Starting Analysis....");
  mgr->StartAnalysis("local", chain, 1000);
}
