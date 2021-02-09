void AddTask_ConvCaloTree(
  Int_t     isMC                          = 0,
  Int_t     IsHeavyIon                    = 0,
  TString   photonCutNumberV0Reader       = "00200009327000008250400000",
  TString   TaskEventCutnumber            = "00010113",
  TString   TaskEMCCutnumber              = "1111100010022700000",
  TString   TaskPHOSCutnumber             = "1111100010022700000",
  TString   TaskConversionCutnumber       = "00200009227302008254404000",
  Double_t  kTree                         = 1.0,  // 0. / 0 / kFALSE for no, 1. / 1 / kTRUE for yes,  x > 1.0 will use only 1/x of the event statistics for the tree
  TString   periodNameV0Reader            = "",
  TString   corrTaskSetting               = "",
  Int_t     enableExtMatchAndQA           = 0,
  Int_t     doSaveSurroundingTracks       = 1,        // 1: save track eta, phi and E on calo surface ; 2: Save all tracks px,py,py and eta,phi on calo
  Bool_t    doSaveMCInfo                  = 0,
  Float_t   minTrackMomentum              = 0.3,
  Bool_t    enableTriggerOverlapRej       = kTRUE,
  TString   settingMaxFacPtHard           = "3."       // maximum factor between hardest jet and ptHard generated
  ){




  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton     = photonCutNumberV0Reader.Data();
  TString cutnumberEvent      = "00000003";
  if(IsHeavyIon==1)
    cutnumberEvent = "10000003";
  else if(IsHeavyIon==2)
    cutnumberEvent = "80000003";

  //========= Check V0 Reader in  ANALYSIS manager  =====
  TString V0ReaderName        = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  AliV0ReaderV1 *fV0ReaderV1  =  NULL;
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    cout << "V0Reader: " << V0ReaderName.Data() << " not found!!"<< endl;
    return;
  } else {
    cout << "V0Reader: " << V0ReaderName.Data() << " found!!"<< endl;
  }

  TObjArray *rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_ConvCaloTree during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
  Bool_t fMinPtHardSet        = kFALSE;
  Double_t minFacPtHard       = -1;
  Bool_t fMaxPtHardSet        = kFALSE;
  Double_t maxFacPtHard       = 100;
  Bool_t fSingleMaxPtHardSet  = kFALSE;
  Double_t maxFacPtHardSingle = 100;
  for(Int_t i = 0; i<rmaxFacPtHardSetting->GetEntries() ; i++){
    TObjString* tempObjStrPtHardSetting     = (TObjString*) rmaxFacPtHardSetting->At(i);
    TString strTempSetting                  = tempObjStrPtHardSetting->GetString();
    if(strTempSetting.BeginsWith("MINPTHFAC:")){
      strTempSetting.Replace(0,10,"");
      minFacPtHard               = strTempSetting.Atof();
      cout << "running with min pT hard jet fraction of: " << minFacPtHard << endl;
      fMinPtHardSet        = kTRUE;
    } else if(strTempSetting.BeginsWith("MAXPTHFAC:")){
      strTempSetting.Replace(0,10,"");
      maxFacPtHard               = strTempSetting.Atof();
      cout << "running with max pT hard jet fraction of: " << maxFacPtHard << endl;
      fMaxPtHardSet        = kTRUE;
    } else if(strTempSetting.BeginsWith("MAXPTHFACSINGLE:")){
      strTempSetting.Replace(0,16,"");
      maxFacPtHardSingle         = strTempSetting.Atof();
      cout << "running with max single particle pT hard fraction of: " << maxFacPtHardSingle << endl;
      fSingleMaxPtHardSet        = kTRUE;
    } else if(rmaxFacPtHardSetting->GetEntries()==1 && strTempSetting.Atof()>0){
      maxFacPtHard               = strTempSetting.Atof();
      cout << "running with max pT hard jet fraction of: " << maxFacPtHard << endl;
      fMaxPtHardSet        = kTRUE;
    }
  }

  AliConvEventCuts *analysisEventCuts = new AliConvEventCuts();
  analysisEventCuts->SetV0ReaderName(V0ReaderName);
  analysisEventCuts->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
  if(fMinPtHardSet) analysisEventCuts->SetMinFacPtHard(minFacPtHard);
  if(fMaxPtHardSet) analysisEventCuts->SetMaxFacPtHard(maxFacPtHard);
  if(fSingleMaxPtHardSet) analysisEventCuts->SetMaxFacPtHardSingleParticle(maxFacPtHardSingle);
  analysisEventCuts->SetCorrectionTaskSetting(corrTaskSetting);
  if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts->SetPeriodEnum(periodNameV0Reader);
  analysisEventCuts->InitializeCutsFromCutString(TaskEventCutnumber.Data());
  analysisEventCuts->SetFillCutHistograms("",kFALSE);

  AliCaloPhotonCuts *analysisClusterCutsEMC = new AliCaloPhotonCuts();
  if(TaskEMCCutnumber.CompareTo("") != 0 ){
    analysisClusterCutsEMC->SetV0ReaderName(V0ReaderName);
    analysisClusterCutsEMC->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCutsEMC->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCutsEMC->InitializeCutsFromCutString(TaskEMCCutnumber.Data());
    analysisClusterCutsEMC->SetFillCutHistograms("");
  }

  AliCaloPhotonCuts *analysisClusterCutsPHOS = new AliCaloPhotonCuts();
  if(TaskPHOSCutnumber.CompareTo("") != 0 ){
    analysisClusterCutsPHOS->SetV0ReaderName(V0ReaderName);
    analysisClusterCutsPHOS->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCutsPHOS->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCutsPHOS->InitializeCutsFromCutString(TaskPHOSCutnumber.Data());
    analysisClusterCutsPHOS->SetFillCutHistograms("");
  }

  AliConversionPhotonCuts *analysisConversionCuts = new AliConversionPhotonCuts();
  if(TaskConversionCutnumber.CompareTo("") != 0 ){
    analysisConversionCuts->SetV0ReaderName(V0ReaderName);
    analysisConversionCuts->InitializeCutsFromCutString(TaskConversionCutnumber.Data());
    analysisConversionCuts->SetFillCutHistograms("");
  }

  AliConversionMesonCuts *analysisMesonCuts = new AliConversionMesonCuts();
  analysisMesonCuts->SetLightOutput(kTRUE);
  analysisMesonCuts->InitializeCutsFromCutString("01631031000000d0");


  AliAnalysisTaskConvCaloTree *fConvCaloTree = new AliAnalysisTaskConvCaloTree(Form("%s_%s_ConvCaloTree",TaskEventCutnumber.Data(),TaskEMCCutnumber.Data()));

  if(TaskEMCCutnumber.CompareTo("") != 0 || TaskPHOSCutnumber.CompareTo("") != 0 ) fConvCaloTree->SetSaveClusters(kTRUE);
  if(TaskConversionCutnumber.CompareTo("") != 0 )fConvCaloTree->SetSaveConversions(kTRUE);
  fConvCaloTree->SetSaveTracks(doSaveSurroundingTracks);

  fConvCaloTree->SetEventCuts(analysisEventCuts,IsHeavyIon);
  if(TaskEMCCutnumber.CompareTo("") != 0 )fConvCaloTree->SetClusterCutsEMC(analysisClusterCutsEMC);
  if(TaskPHOSCutnumber.CompareTo("") != 0 )fConvCaloTree->SetClusterCutsPHOS(analysisClusterCutsPHOS);
  if(TaskConversionCutnumber.CompareTo("") != 0 )fConvCaloTree->SetConversionCuts(analysisConversionCuts);
  fConvCaloTree->SetMesonCuts(analysisMesonCuts);
  fConvCaloTree->SetCorrectionTaskSetting(corrTaskSetting);
  fConvCaloTree->SetIsMC(isMC);
  if(isMC && doSaveMCInfo) fConvCaloTree->SetSaveMCInformation(kTRUE);
  fConvCaloTree->SetMinTrackPt(minTrackMomentum);
  fConvCaloTree->SetV0ReaderName(V0ReaderName);
  mgr->AddTask(fConvCaloTree);

  mgr->ConnectInput  (fConvCaloTree, 0,  cinput );
  mgr->ConnectOutput (fConvCaloTree, 1, mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ?  Form("GammaCaloQA_%s_%s", TaskEventCutnumber.Data(), TaskEMCCutnumber.Data()) : Form("GammaCaloQA_%s_%s_%s", TaskEventCutnumber.Data(), TaskEMCCutnumber.Data(),corrTaskSetting.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:GammaCaloQA_%s_%s", mgr->GetCommonFileName(), TaskEventCutnumber.Data(), TaskEMCCutnumber.Data())) );
  mgr->ConnectOutput (fConvCaloTree, 2, mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ?  Form("ConvCaloPhotons_%s_%s", TaskEventCutnumber.Data(), TaskEMCCutnumber.Data()) : Form("ConvCaloPhotons_%s_%s_%s", TaskEventCutnumber.Data(), TaskEMCCutnumber.Data(),corrTaskSetting.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()) );

  return;
}
