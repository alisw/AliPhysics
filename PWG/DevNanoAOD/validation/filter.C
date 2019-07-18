void filter(TString from)
{
  AliAnalysisManager* mgr = new AliAnalysisManager("NanoAOD Filter", "NanoAOD filter for nanoAOD production");
    
  if (from == "AOD") {
    AliAODInputHandler* iH = new AliAODInputHandler();
    mgr->SetInputEventHandler(iH);

    // Define aod output handler
    AliAODHandler* aodOutputHandler = new AliAODHandler();
    aodOutputHandler->SetOutputFileName("AliAOD.NanoAOD.root");
    mgr->SetOutputEventHandler(aodOutputHandler);
    
    // Multiplicity selection
    gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    
    // Input files
    TChain * chain = new TChain("aodTree");
    chain->Add("AliAOD.root");
  } else if (from == "ESD") {
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
    
    // ESD filter
    gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C(kFALSE, kFALSE, kFALSE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE, 1500, 3, kTRUE, kFALSE, kFALSE, kFALSE)");    

    // Input files
    TChain * chain = new TChain("esdTree");
    chain->Add("AliESDs.root");
  }
  
  AliAnalysisTaskNanoAODFilter* task = (AliAnalysisTaskNanoAODFilter*) gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWG/DevNanoAOD/AddTaskNanoAODFilter.C(0, kFALSE)");
  task->AddSetter(new AliNanoAODSimpleSetter);
  
  // No event selection here to get identical events in nano and AOD
  // No track selection here to get identical events in nano and AOD

  // Fields to store
  // event level
  // Note: vertices are kept by default
  task->SetVarListHeader("BunchCrossNumber,OrbitNumber,PeriodNumber,CentrV0M,CentrTRK,CentrCL0,CentrCL1,MagField,OfflineTrigger,RunNumber,T0Spread,NumberOfESDTracks,MultSelection.RefMult08");
  // track level
  task->SetVarListTrack("pt,phi,theta,chi2perNDF,posx,posy,posz,pDCAx,pDCAy,pDCAz,posDCAx,posDCAy,posDCAz,DCA,RAtAbsorberEnd,TPCncls,ID,TPCnclsF,TPCNCrossedRows,TrackPhiOnEMCal,TrackEtaOnEMCal,TrackPtOnEMCal,ITSsignal,TPCsignal,TPCsignalTuned,TPCsignalN,TPCmomentum,TPCTgl,TOFsignal,integratedLength,TOFsignalTuned,HMPIDsignal,HMPIDoccupancy,TRDsignal,TRDChi2,TRDnSlices,TPCnclsS,FilterMap,TOFBunchCrossing,TOFchi2,TOFsignalDz,TOFsignalDx,Status");


  mgr->SetDebugLevel(1); // enable debug printouts
  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();
  
  Printf("Starting Analysis....");
  mgr->StartAnalysis("local", chain, 20);
}
