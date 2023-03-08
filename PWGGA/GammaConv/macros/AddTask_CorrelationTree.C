void AddTask_CorrelationTree(
    TString photonCutNumberV0Reader = "060000084001001500000000",
    TString TaskEventCutnumber = "00000003",
    TString TaskPhotonCutnumber = "090000092663743800000000",
    TString TaskCaloCutNumber = "411790109fe32230000",
    Bool_t isMC = kFALSE,
    Int_t IsHeavyIon = 0,
    Bool_t kHistograms = kTRUE,
    Double_t kTree = 1.0, // 0. / 0 / kFALSE for no, 1. / 1 / kTRUE for yes,  x > 1.0 will use only 1/x of the event statistics for the tree
    TString V0ReaderCutNumberAODBranch = "0000000060084001001500000",
    Bool_t runBasicQAWithStandardOutput = kTRUE,
    Bool_t doEtaShiftV0Reader = kFALSE,
    Bool_t enableV0findingEffi = kFALSE, // enables V0finding efficiency histograms
    TString fileNameExternalInputs = "",
    Bool_t enableElecDeDxPostCalibration = kFALSE,
    TString additionalTrainConfig = "")
{

  AliCutHandlerPCM cuts(13);
  TString fileNamedEdxPostCalib = cuts.GetSpecialFileNameFromString(fileNameExternalInputs, "FEPC:");

  TString addTaskName = "AddTask_CorrelationTree";
  TString corrTaskSetting = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF", "", addTaskName);
  if (corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  TString strTrackMatcherRunningMode  = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "TM", "", addTaskName);
  if(additionalTrainConfig.Contains("TM"))
    trackMatcherRunningMode = strTrackMatcherRunningMode.Atoi();

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    return;
  }
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler = mgr->GetInputEventHandler();

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = photonCutNumberV0Reader.Data();
  TString cutnumberEvent = "00000003";
  if (IsHeavyIon == 1)
    cutnumberEvent = "10000003";
  else if (IsHeavyIon == 2)
    cutnumberEvent = "80000003";

  //========= Check V0 Reader in  ANALYSIS manager  =====
  TString V0ReaderName = Form("V0ReaderV1_%s_%s", cutnumberEvent.Data(), cutnumberPhoton.Data());
  AliV0ReaderV1 *fV0ReaderV1 = NULL;
  if (!(AliV0ReaderV1 *)mgr->GetTask(V0ReaderName.Data()))
  {
    cout << "V0Reader: " << V0ReaderName.Data() << " not found!!" << endl;
    return;
  }
  else
  {
    cout << "V0Reader: " << V0ReaderName.Data() << " found!!" << endl;
  }

    TString caloCutPos = TaskCaloCutNumber;
    caloCutPos.Resize(1);
  TString TrackMatcherName = Form("CaloTrackMatcher_%s_%i", caloCutPos.Data(), trackMatcherRunningMode);
  if (corrTaskSetting.CompareTo(""))
  {
    TrackMatcherName = TrackMatcherName + "_" + corrTaskSetting.Data();
    cout << "Using separate track matcher for correction framework setting: " << TrackMatcherName.Data() << endl;
  }
if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi(),trackMatcherRunningMode);
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      fTrackMatcher->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }
  AliConvEventCuts *analysisEventCuts = new AliConvEventCuts();
  analysisEventCuts->SetV0ReaderName(V0ReaderName);
  analysisEventCuts->InitializeCutsFromCutString(TaskEventCutnumber.Data());
  analysisEventCuts->SetFillCutHistograms("", kFALSE);
    analysisEventCuts->SetCorrectionTaskSetting(corrTaskSetting);

  AliCaloPhotonCuts *analysisCaloCuts = new AliCaloPhotonCuts(isMC);
  analysisCaloCuts->SetV0ReaderName(V0ReaderName);
  analysisCaloCuts->InitializeCutsFromCutString(TaskCaloCutNumber.Data());
  // analysisCaloCuts->SetFillCutHistograms("",kFALSE);
  analysisCaloCuts->SetLightOutput(kTRUE);
  analysisCaloCuts->SetCaloTrackMatcherName(TrackMatcherName);
  analysisCaloCuts->SetCorrectionTaskSetting(corrTaskSetting);

  AliConversionPhotonCuts *analysisCuts = new AliConversionPhotonCuts();
  analysisCuts->SetV0ReaderName(V0ReaderName);
  if (enableElecDeDxPostCalibration)
  {
    if (isMC == 0)
    {
      if (fileNamedEdxPostCalib.CompareTo("") != 0)
      {
        analysisCuts->SetElecDeDxPostCalibrationCustomFile(fileNamedEdxPostCalib);
        cout << "Setting custom dEdx recalibration file: " << fileNamedEdxPostCalib.Data() << endl;
      }
      analysisCuts->SetDoElecDeDxPostCalibration(enableElecDeDxPostCalibration);
      cout << "Enabled TPC dEdx recalibration." << endl;
    }
    else
    {
      cout << "ERROR enableElecDeDxPostCalibration set to True even if MC file. Automatically reset to 0" << endl;
      enableElecDeDxPostCalibration = kFALSE;
      analysisCuts->SetDoElecDeDxPostCalibration(kFALSE);
    }
  }
  analysisCuts->InitializeCutsFromCutString(TaskPhotonCutnumber.Data());
  analysisCuts->SetFillCutHistograms("", kFALSE);

  AliAnalysisTaskCorrelationTree *fQA = new AliAnalysisTaskCorrelationTree(Form("%s_%s_QA", TaskEventCutnumber.Data(), TaskPhotonCutnumber.Data()));
  fQA->SetEventCuts(analysisEventCuts, IsHeavyIon);
  fQA->SetConversionCuts(analysisCuts, IsHeavyIon);
  fQA->SetCaloCuts(analysisCaloCuts, IsHeavyIon);
  fQA->SetIsMC(isMC);
  fQA->SetV0ReaderName(V0ReaderName);
  fQA->SetCorrectionTaskSetting(corrTaskSetting);
  mgr->AddTask(fQA);

  if (runBasicQAWithStandardOutput)
  {
    AliAnalysisDataContainer *coutput =
        mgr->CreateContainer(Form("GammaConv_V1QA_%s_%s", TaskEventCutnumber.Data(), TaskPhotonCutnumber.Data()), TList::Class(),
                             AliAnalysisManager::kOutputContainer, Form("%s:GammaConvV1_QA_%s_%s", AliAnalysisManager::GetCommonFileName(), TaskEventCutnumber.Data(), TaskPhotonCutnumber.Data()));
    mgr->ConnectOutput(fQA, 1, coutput);
  }
  else
  {
    AliAnalysisDataContainer *coutput =
        mgr->CreateContainer(Form("GammaConv_V1QA_%s_%s", TaskEventCutnumber.Data(), TaskPhotonCutnumber.Data()), TList::Class(),
                             AliAnalysisManager::kOutputContainer, Form("GammaConvV1_QA_%s_%s.root", TaskEventCutnumber.Data(), TaskPhotonCutnumber.Data()));
    mgr->ConnectOutput(fQA, 1, coutput);
  }
  mgr->ConnectOutput(fQA, 2, mgr->CreateContainer(Form("CorrelationTree_%s_%s", TaskEventCutnumber.Data(), TaskPhotonCutnumber.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput(fQA, 0, cinput);

  // connect containers
  return;
}
