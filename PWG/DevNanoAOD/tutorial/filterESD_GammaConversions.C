void filterESD_GammaConversions()
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
  
  // PID response
  gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");

  // OCDB
  gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/AddTaskConfigOCDB.C(\"raw://\")");
  
  // V0 finder
  AliAnalysisTaskWeakDecayVertexer* v0Finder = (AliAnalysisTaskWeakDecayVertexer*) gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/AddTaskWeakDecayVertexer.C");
  v0Finder->SetUseImprovedFinding();
  
  AliCDBManager::Instance()->SetRun(285396);
  Int_t iCollision = 0;           // pp:0 , PbPb:1, pPb:2
  gInterpreter->ExecuteMacro(Form("$ALICE_PHYSICS/PWGGA/GammaConv/macros/AddTask_ConversionAODProduction.C(%d, kFALSE, \"\")", iCollision));
  mgr->RegisterExtraFile("AliAODGammaConversion.root");
  
  // ESD filter
  gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C(kFALSE, kFALSE, kFALSE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE, 1500, 3, kTRUE, kFALSE, kFALSE, kTRUE)");
  
  AliAnalysisTaskNanoAODFilter* task = (AliAnalysisTaskNanoAODFilter*) gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWG/DevNanoAOD/macros/AddTaskNanoAODFilter.C(0, kFALSE)");
  task->AddSetter(new AliNanoAODSimpleSetter);
  task->SelectCollisionCandidates(AliVEvent::kAny);
  
  // Event selection
  AliAnalysisNanoAODEventCuts* evtCuts = new AliAnalysisNanoAODEventCuts;
  // NOTE filter bit set in AliEventCuts automatically

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
  task->SetVarListTrack("pt,theta,phi,ID");

  task->SetTrkCuts(trkCuts);
  task->AddEvtCuts(evtCuts);

  // V0s
  task->SaveV0s(kTRUE, new AliAnalysisNanoAODV0Cuts);
  task->SaveCascades(kTRUE);
  
  task->SaveConversionPhotons(kTRUE);

  mgr->SetDebugLevel(1); // enable debug printouts
  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();
  
  // Input files
  TChain * chain = new TChain("esdTree");
  chain->Add("AliESDs.root");

  Printf("Starting Analysis....");
  mgr->StartAnalysis("local", chain, 100);
}
