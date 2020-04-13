AliAnalysisTask *AddTask_ConversionsForNano(Int_t dataset = 0, Bool_t isMC =
                                                kFALSE,
                                            TString periodNameV0Reader = "",
                                            TString cutnumberPhoton =
                                                "00200078000000001240820000") {
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
  AliConvEventCuts *fEventCutsDummy = NULL;
  AliConversionPhotonCuts *fCuts = NULL;
  AliConversionPhotonCuts *fCutsDummy = NULL;
  TString cutnumberEvent = "00000000";
  TString cutnumberEventDummy = "02132005";
  TString cutnumberPhotonDummy = "17h000i3c40104306632f21992";

  AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("ConvGammaAODProduction");
  if (periodNameV0Reader.CompareTo("") != 0)
    fV0ReaderV1->SetPeriodName(periodNameV0Reader);
  fV0ReaderV1->SetAddv0sInESDFilter(kTRUE);
  fV0ReaderV1->SetCreateAODs(kTRUE);
  fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
  fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

  if (cutnumberEvent != "") {
    fEventCuts = new AliConvEventCuts(cutnumberEvent.Data(),
                                      cutnumberEvent.Data());
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

  // This is needed because the ESD filter needs both, onthefly and offline, but those are actually not stored in the Nano...
  AliV0ReaderV1 *fV0ReaderV1Dummy = new AliV0ReaderV1(
      "ConvGammaAODProductionDummy");
  fV0ReaderV1Dummy->SetAddv0sInESDFilter(kTRUE);
  fV0ReaderV1Dummy->SetCreateAODs(kTRUE);

  if (cutnumberEventDummy != "") {
    fEventCutsDummy = new AliConvEventCuts(cutnumberEventDummy.Data(),
                                           cutnumberEventDummy.Data());
    fEventCutsDummy->SetPreSelectionCutFlag(kTRUE);
    fEventCutsDummy->SetV0ReaderName("ConvGammaAODProductionDummy");
    if (fEventCutsDummy->InitializeCutsFromCutString(
        cutnumberEventDummy.Data())) {
      fV0ReaderV1Dummy->SetEventCuts(fEventCutsDummy);
    }
  }

  if (cutnumberPhotonDummy != "") {
    fCutsDummy = new AliConversionPhotonCuts(cutnumberPhotonDummy.Data(),
                                             cutnumberPhotonDummy.Data());
    fCutsDummy->SetIsHeavyIon(false);
    fCutsDummy->SetV0ReaderName("ConvGammaAODProductionDummy");
    if (fCutsDummy->InitializeCutsFromCutString(cutnumberPhotonDummy.Data())) {
      fV0ReaderV1Dummy->SetConversionCuts(fCutsDummy);
    }
  }

  fV0ReaderV1->Init();
  fV0ReaderV1Dummy->Init();
  AliLog::SetGlobalLogLevel(AliLog::kFatal);

  // connect input V0Reader
  mgr->AddTask(fV0ReaderV1);
  mgr->AddTask(fV0ReaderV1Dummy);
  mgr->ConnectInput(fV0ReaderV1, 0, cinput);
  mgr->ConnectInput(fV0ReaderV1Dummy, 0, cinput);

  AliLog::SetGlobalLogLevel(AliLog::kInfo);

  //================================================
  //              data containers
  //================================================

  AliAnalysisDataContainer *coutput = mgr->CreateContainer(
      "PCM onflyV0Finder container", TBits::Class(),
      AliAnalysisManager::kExchangeContainer);

  mgr->ConnectOutput(fV0ReaderV1, 1, coutput);

  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(
      "PCM offlineV0Finder container", TBits::Class(),
      AliAnalysisManager::kExchangeContainer);

  mgr->ConnectOutput(fV0ReaderV1Dummy, 1, coutput2);

  return fV0ReaderV1;
}
