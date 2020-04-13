void filterAOD_GammaConversions()
{
  AliAnalysisManager* mgr = new AliAnalysisManager("NanoAOD Filter", "NanoAOD filter for nanoAOD production");
    
  AliAODInputHandler* iH = new AliAODInputHandler();
  mgr->SetInputEventHandler(iH);

  // Define aod output handler
  AliAODHandler* aodOutputHandler = new AliAODHandler();
  aodOutputHandler->SetOutputFileName("AliAOD.NanoAOD.root");
  mgr->SetOutputEventHandler(aodOutputHandler);
  
  // PID response
  gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");

  // Photons
  TString cutnumberPhoton = "00200078400000001240820000";
  TString cutnumberEvent = "00000000";
  Bool_t fillHistos = kTRUE;
  TString V0ReaderName = TString::Format("ConvGammaAODProduction");

  AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
  fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
  fV0ReaderV1->SetCreateAODs(kFALSE);  // AOD Output
  fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

  if (cutnumberEvent != "") {
    AliConvEventCuts *fEventCuts = new AliConvEventCuts(cutnumberEvent.Data(), cutnumberEvent.Data());
    fEventCuts->SetPreSelectionCutFlag(kTRUE);
    fEventCuts->SetV0ReaderName(V0ReaderName);
    fV0ReaderV1->SetEventCuts(fEventCuts);
    fEventCuts->SetFillCutHistograms("", kFALSE);
  }

  // Set AnalysisCut Number
  AliConversionPhotonCuts *fCuts = NULL;
  if (cutnumberPhoton != "") {
    fCuts = new AliConversionPhotonCuts(cutnumberPhoton.Data(), cutnumberPhoton.Data());
    fCuts->SetPreSelectionCutFlag(kTRUE);
    fCuts->SetIsHeavyIon(false);
    fCuts->SetV0ReaderName(V0ReaderName);
    fCuts->SetFillCutHistograms("", fillHistos);
    if (fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data()))
      fV0ReaderV1->SetConversionCuts(fCuts);
  }
  fV0ReaderV1->Init();

  // connect input V0Reader
  mgr->AddTask(fV0ReaderV1);
  mgr->ConnectInput(fV0ReaderV1, 0, mgr->GetCommonInputContainer());
  
  // Nano AOD filter task
  AliAnalysisTaskNanoAODFilter* task = (AliAnalysisTaskNanoAODFilter*) gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWG/DevNanoAOD/macros/AddTaskNanoAODFilter.C(0, kFALSE)");
  task->AddSetter(new AliNanoAODSimpleSetter);
  
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
  task->SetVarListHeader("OfflineTrigger,MagField");
  // track level
  task->SetVarListTrack("pt,theta,phi,ID");

  task->SetTrkCuts(trkCuts);
  task->AddEvtCuts(evtCuts);
  
  task->SaveConversionPhotons(kTRUE);

  mgr->SetDebugLevel(1); // enable debug printouts
  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();
  
  // Input files
  TChain * chain = new TChain("aodTree");
  chain->Add("AliAOD.root");
  
  // For conversions
  TChain* friendChain = new TChain("aodTree", "AliAODGammaConversion.root");
  friendChain->Add("AliAODGammaConversion.root");
  chain->AddFriend(friendChain);

  Printf("Starting Analysis....");
  mgr->StartAnalysis("local", chain, 1000);
}
