AliAnalysisTask *AddTask_ConversionsForNano(
    Int_t dataset = 0, Bool_t isMC = kFALSE, TString periodNameV0Reader = "",
    TString cutnumberPhoton = "00200078000000001240820000") {
  gSystem->Load("libPWGGAGammaConv");

  // get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_V0ReaderV1", "No analysis manager found.");
    return 0;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  Bool_t enableV0findingEffi = kFALSE;
  Bool_t fillHistos = kTRUE;
  Bool_t runLightOutput = kTRUE;
  AliConvEventCuts *fEventCuts = NULL;
  AliConversionPhotonCuts *fCuts = NULL;
  TString cutnumberEvent = "00000000";

  AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("ConvGammaAODProduction");

  if (periodNameV0Reader.CompareTo("") != 0)
    fV0ReaderV1->SetPeriodName(periodNameV0Reader);
  fV0ReaderV1->SetAddv0sInESDFilter(kTRUE);
  fV0ReaderV1->SetCreateAODs(kTRUE);
  fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
  fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

  if (!mgr) {
    Error("AddTask_V0ReaderV1", "No analysis manager found.");
    return NULL;
  }

  if (cutnumberEvent != "") {
    fEventCuts =
        new AliConvEventCuts(cutnumberEvent.Data(), cutnumberEvent.Data());
    fEventCuts->SetPreSelectionCutFlag(kTRUE);
    fEventCuts->SetV0ReaderName("ConvGammaAODProduction");
    if (fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())) {
      fV0ReaderV1->SetEventCuts(fEventCuts);
    }
  }

  if (cutnumberPhoton != "") {
    fCuts = new AliConversionPhotonCuts(cutnumberPhoton.Data(),
                                        cutnumberPhoton.Data());
    fCuts->SetIsHeavyIon(false);
    fCuts->SetV0ReaderName("ConvGammaAODProduction");
    if (fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())) {
      fV0ReaderV1->SetConversionCuts(fCuts);
    }
  }
  fV0ReaderV1->Init();
  AliLog::SetGlobalLogLevel(AliLog::kFatal);

  // connect input V0Reader
  mgr->AddTask(fV0ReaderV1);
  mgr->ConnectInput(fV0ReaderV1, 0, cinput);

  AliLog::SetGlobalLogLevel(AliLog::kInfo);

  //================================================
  //              data containers
  //================================================
  //            find input container
  // below the trunk version
  AliAnalysisDataContainer *coutput =
      mgr->CreateContainer("PCM onflyV0Finder container", TBits::Class(),
                           AliAnalysisManager::kExchangeContainer);

  mgr->ConnectOutput(fV0ReaderV1, 1, coutput);

  return fV0ReaderV1;
}
