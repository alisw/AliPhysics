/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                       *
 * Author: Friederike Bock, Daniel Mühlheim                     *
 * Version 1.0                                 *
 *                                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its    *
 * documentation strictly for non-commercial purposes is hereby granted    *
 * without fee, provided that the above copyright notice appears in all    *
 * copies and that both the copyright notice and this permission notice    *
 * appear in the supporting documentation. The authors make no claims    *
 * about the suitability of this software for any purpose. It is      *
 * provided "as is" without express or implied warranty.               *
 **************************************************************************/

//***************************************************************************************
//This AddTask is supposed to set up the main task
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson.cxx) for
//PbPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_CaloMode_PbPb(
    Int_t     trainConfig                 = 1,
    Int_t     isMC                        = 0,                        //run MC
    TString   photonCutNumberV0Reader     = "",                       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
    Int_t     selectHeavyNeutralMeson     = 0,                        //run eta prime instead of omega
    Int_t     enableQAMesonTask           = 1,                        // enable QA in AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson; 0: no QA, 1: general meson QA, 2: background QA, 3: 3D histogram, 4: Dalitz plots, 5: trees, 23: enable background calculations; 
                                                                      //    combinations: 6: 1+2, 7: 1+2+3, 8: 1+2+3+5, 9: 2+3, 10: 2+3+5, 11: 1+2+3+4+5
                                                                      //    QA can't be run with light output! 
    Int_t     enableExtMatchAndQA         = 0,                        // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
    Int_t     enableTriggerMimicking      = 0,                        // enable trigger mimicking
    Bool_t    enableTriggerOverlapRej     = kFALSE,                   // enable trigger overlap rejection
    TString   fileNameExternalInputs      = "",                       // path to file for weigting input
    Bool_t    doWeighting                 = kFALSE,                   //enable Weighting
    TString   generatorName               = "HIJING",
    Double_t  tolerance                   = -1,
    TString   periodNameV0Reader          = "",                       // period Name for V0Reader
    Int_t     runLightOutput              = 0,                        // run light output option 0: no light output 1: most cut histos stiched off 2: unecessary omega hists turned off as well
    Int_t     prefilterRunFlag            = 1500,                     // flag to change the prefiltering of ESD tracks. See SetHybridTrackCutsAODFiltering() in AliPrimaryPionCuts
    Bool_t    usePtDepSelectionWindowCut  = kFALSE,                   // use pt dependent meson selection window cut
    Bool_t    enableSortingMCLabels       = kTRUE,                    // enable sorting for MC cluster labels
    TString   additionalTrainConfig       = "0"                       // additional counter for trainconfig, this has to be always the last parameter
  ) {

  AliCutHandlerPCM cuts(13);
  TString addTaskName                       = "AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_CaloMode_PbPb";

  //parse additionalTrainConfig flag
  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  TString unsmearingoutputs = "0123"; // 0: No correction, 1: One pi0 mass errer subtracted, 2: pz of pi0 corrected to fix its mass, 3: Lambda(alpha)*DeltaPi0 subtracted

  TString corrTaskSetting             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF", "", addTaskName);
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){std::cout << "ERROR during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << std::endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      if(tempStr.BeginsWith("TM")){
        TString tempType = tempStr;
        tempType.Replace(0,2,"");
        trackMatcherRunningMode = tempType.Atoi();
        cout << Form("INFO: AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_CaloMode_PbPb will use running mode '%i' for the TrackMatcher!",trackMatcherRunningMode) << endl;
      }
      if(tempStr.BeginsWith("UNSMEARING")){ // 0: No correction, 1: One pi0 mass errer subtracted, 2: pz of pi0 corrected to fix its mass, 3: Lambda(alpha)*DeltaPi0 subtracted
        TString tempType = tempStr;
        tempType.Replace(0,9,"");
        unsmearingoutputs = tempType;
        cout << "INFO: AddTask_GammaConvNeutralMesonPiPlPiMiPiZero_CaloMode_pPb will output the following minv_pT histograms:" << endl;
        if(unsmearingoutputs.Contains("0")) cout << "- Uncorrected" << endl;
        if(unsmearingoutputs.Contains("1")) cout << "- SubNDM" << endl;
        if(unsmearingoutputs.Contains("2")) cout << "- Fixpz" << endl;
        if(unsmearingoutputs.Contains("3")) cout << "- SubLambda" << endl;
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    std::cout << "INFO: AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_CaloMode_PbPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << std::endl;
  }

  Int_t isHeavyIon = 1;
  Int_t neutralPionMode = 2;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_CaloMode_PbPb_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();


  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = photonCutNumberV0Reader.Data();
  TString cutnumberEvent = "10000003";
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName        = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  AliV0ReaderV1 *fV0ReaderV1  =  NULL;
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    std::cout << "V0Reader: " << V0ReaderName.Data() << " not found!!"<< std::endl;
    return;
  } else {
    std::cout << "V0Reader: " << V0ReaderName.Data() << " found!!"<< std::endl;
  }

  TString PionCuts      = "000000200";
  //================================================
  //========= Add Pion Selector ====================
  if( !(AliPrimaryPionSelector*)mgr->GetTask("PionSelector") ){
    AliPrimaryPionSelector *fPionSelector = new AliPrimaryPionSelector("PionSelector");
    AliPrimaryPionCuts *fPionCuts=0;
    if( PionCuts!=""){
      fPionCuts= new AliPrimaryPionCuts(PionCuts.Data(),PionCuts.Data());
      fPionCuts->SetPrefilterRunFlag(prefilterRunFlag);
      if(runLightOutput>0) fPionCuts->SetLightOutput(kTRUE);
      fPionCuts->SetPeriodName(periodNameV0Reader);
      if(fPionCuts->InitializeCutsFromCutString(PionCuts.Data())){
        fPionSelector->SetPrimaryPionCuts(fPionCuts);
        fPionCuts->SetFillCutHistograms("",kTRUE);
      }
    }

    fPionSelector->Init();
    mgr->AddTask(fPionSelector);

    AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
    mgr->ConnectInput (fPionSelector,0,cinput1);
  }

  AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson *task=NULL;
  task= new AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson(Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i",neutralPionMode, trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  if(runLightOutput>1) task->SetLightOutput(kTRUE);
  task->SetTolerance(tolerance);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);
  task->SetCorrectionTaskSetting(corrTaskSetting);



  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                          ETA MESON
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if( trainConfig == 1 ) {
    cuts.AddCutHeavyMesonCalo("00000113","1111113047032230000","000010400","0103503a00000000","0403503000000000");
    cuts.AddCutHeavyMesonCalo("00000113","1111113047032230000","000010400","0103503a00000000","0403503000000000");
    cuts.AddCutHeavyMesonCalo("10130a13","411798305k0a2220000","000010400","0103503a00000000","0403503000000000");
    cuts.AddCutHeavyMesonCalo("11310a13","411798305k0b2220000","000010400","0103503a00000000","0403503000000000");
    cuts.AddCutHeavyMesonCalo("13530a13","411798305k032220000","000010400","0103503a00000000","0403503000000000");
    cuts.AddCutHeavyMesonCalo("15910a13","411798305k032220000","000010400","0103503a00000000","0403503000000000");
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                          OMEGA MESON
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  } else if( trainConfig == 100 ) {
    cuts.AddCutHeavyMesonCalo("10130a13","411798305k0a2220000","32c51070a","0103603700000000","0453503000000000"); //
    cuts.AddCutHeavyMesonCalo("11310a13","411798305k0b2220000","32c51070a","0103603700000000","0453503000000000"); //
    cuts.AddCutHeavyMesonCalo("13530a13","411798305k032220000","32c51070a","0103603700000000","0453503000000000"); //
    cuts.AddCutHeavyMesonCalo("15910a13","411798305k032220000","32c51070a","0103603700000000","0453503000000000"); //

  }else if (trainConfig == 101){ // Standard 5 TeV omega INT7 cutstring in pPb -> First look into PbPb
    cuts.AddCutHeavyMesonCalo("10910a13","411790105fe30220000","32c51079a","0000003100000000","0400503000000000"); // 0-90% MB
    cuts.AddCutHeavyMesonCalo("10130a13","411790105fe30220000","32c51079a","0000003100000000","0400503000000000"); // 0-10% Triggered
    cuts.AddCutHeavyMesonCalo("11310a13","411790105fe30220000","32c51079a","0000003100000000","0400503000000000"); // 10-30% MB
    cuts.AddCutHeavyMesonCalo("13530a13","411790105fe30220000","32c51079a","0000003100000000","0400503000000000"); // 30-50% Triggeres
    cuts.AddCutHeavyMesonCalo("15910a13","411790105fe30220000","32c51079a","0000003100000000","0400503000000000"); // 50-90% MB

  }else if (trainConfig == 102){ // Standard 5 TeV omega INT7 cutstring in pPb -> Semicentral
    cuts.AddCutHeavyMesonCalo("12410a13","411790105fe30220000","32c51079a","0000003100000000","0400503000000000"); // 20-40% 
  }else if (trainConfig == 103){ // Standard 5 TeV omega INT7 cutstring in pPb -> Peripheral
    cuts.AddCutHeavyMesonCalo("16810a13","411790105fe30220000","32c51079a","0000003100000000","0400503000000000"); // 60-80% 
 
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                          ETA PRIME MESON
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  } else if( trainConfig == 200 ) {
    cuts.AddCutHeavyMesonCalo("10130a13","411798305k0a2220000","32c510700","0103603l00000000","0453503000000000");
    cuts.AddCutHeavyMesonCalo("11310a13","411798305k0b2220000","32c510700","0103603l00000000","0453503000000000");
    cuts.AddCutHeavyMesonCalo("13530a13","411798305k032220000","32c510700","0103603l00000000","0453503000000000");
    cuts.AddCutHeavyMesonCalo("15910a13","411798305k032220000","32c510700","0103603l00000000","0453503000000000");
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                          D0 MESON
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  } else if( trainConfig == 300 ) {
    cuts.AddCutHeavyMesonCalo("10130a13","411798305k0a2220000","000010400","0103503a00000000","0103503000000000");
    cuts.AddCutHeavyMesonCalo("11310a13","411798305k0b2220000","000010400","0103503a00000000","0103503000000000");
    cuts.AddCutHeavyMesonCalo("13530a13","411798305k032220000","000010400","0103503a00000000","0103503000000000");
    cuts.AddCutHeavyMesonCalo("15910a13","411798305k032220000","000010400","0103503a00000000","0103503000000000");
  } else {
    Error(Form("GammaConvNeutralMeson_CaloMode_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    std::cout << "\n\n****************************************************" << std::endl;
    std::cout << "ERROR: No valid cuts stored in CutHandlerNeutralCalo! Returning..." << std::endl;
    std::cout << "****************************************************\n\n" << std::endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ClusterCutList  = new TList();
  TList *NeutralPionCutList = new TList();
  TList *MesonCutList = new TList();
  TList *PionCutList  = new TList();

  TList *HeaderList = new TList();
  TObjString *Header1 = new TObjString("pi0_1");
  HeaderList->Add(Header1);
  TObjString *Header3 = new TObjString("eta_2");
  HeaderList->Add(Header3);

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
  NeutralPionCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisNeutralPionCuts   = new AliConversionMesonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];
  PionCutList->SetOwner(kTRUE);
  AliPrimaryPionCuts **analysisPionCuts     = new AliPrimaryPionCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    //create AliCaloTrackMatcher instance, if there is none present
    TString caloCutPos = cuts.GetClusterCut(i);
    caloCutPos.Resize(1);
    TString TrackMatcherName = Form("CaloTrackMatcher_%s_%i",caloCutPos.Data(),trackMatcherRunningMode);
    if(corrTaskSetting.CompareTo("")){
      TrackMatcherName = TrackMatcherName+"_"+corrTaskSetting.Data();
      cout << "Using separate track matcher for correction framework setting: " << TrackMatcherName.Data() << endl;
    }
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi(),trackMatcherRunningMode);
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      fTrackMatcher->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    if(runLightOutput>0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);


    analysisClusterCuts[i] = new AliCaloPhotonCuts();
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    if(runLightOutput>0) analysisClusterCuts[i]->SetLightOutput(kTRUE);

    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if( ! analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data()) ) {
      std::cout<<"ERROR: analysisClusterCuts [" <<i<<"]"<<std::endl;
      return 0;
    } else {
      analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
      ClusterCutList->Add(analysisClusterCuts[i]);
      analysisClusterCuts[i]->SetFillCutHistograms("");
    }

    analysisNeutralPionCuts[i] = new AliConversionMesonCuts();
    if(runLightOutput>0) analysisNeutralPionCuts[i]->SetLightOutput(kTRUE);
    if( ! analysisNeutralPionCuts[i]->InitializeCutsFromCutString((cuts.GetNDMCut(i)).Data()) ) {
      std::cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<std::endl;
      return 0;
    } else {
      NeutralPionCutList->Add(analysisNeutralPionCuts[i]);
      analysisNeutralPionCuts[i]->SetFillCutHistograms("");
    }

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if(runLightOutput>0) analysisMesonCuts[i]->SetLightOutput(kTRUE);
    if( ! analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data()) ) {
      std::cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<std::endl;
      return 0;
    } else {
      MesonCutList->Add(analysisMesonCuts[i]);
      analysisMesonCuts[i]->SetFillCutHistograms("");
    }
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);

    TString cutName( Form("%s_%s_%s_%s_%s",(cuts.GetEventCut(i)).Data(), (cuts.GetClusterCut(i)).Data(),(cuts.GetPionCut(i)).Data(),(cuts.GetNDMCut(i)).Data(), (cuts.GetMesonCut(i)).Data() ) );
    analysisPionCuts[i] = new AliPrimaryPionCuts();
    analysisPionCuts[i]->SetPrefilterRunFlag(prefilterRunFlag);
    analysisPionCuts[i]->SetPeriodName(periodNameV0Reader);
    if(runLightOutput>0) analysisPionCuts[i]->SetLightOutput(kTRUE);

        if( !analysisPionCuts[i]->InitializeCutsFromCutString((cuts.GetPionCut(i)).Data())) {
      std::cout<< "ERROR:  analysisPionCuts [ " <<i<<" ] "<<std::endl;
      return 0;
    } else {
      PionCutList->Add(analysisPionCuts[i]);
      analysisPionCuts[i]->SetFillCutHistograms("",kFALSE,cutName);
    }
  }

  task->SetNDMRecoMode(neutralPionMode);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetClusterCutList(ClusterCutList);
  task->SetNeutralPionCutList(NeutralPionCutList);
  task->SetMesonCutList(MesonCutList);
  task->SetPionCutList(PionCutList);
  task->SetCorrectionTaskSetting(corrTaskSetting);

  task->SetMoveParticleAccordingToVertex(kFALSE);
  task->SetSelectedHeavyNeutralMeson(selectHeavyNeutralMeson);

  task->SetDoMesonQA(enableQAMesonTask );

  task->SetUnsmearedOutputs(unsmearingoutputs);

  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig) : Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i_%s.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig, corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
