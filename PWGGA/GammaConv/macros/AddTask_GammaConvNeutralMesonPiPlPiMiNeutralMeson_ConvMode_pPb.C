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
//pPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_ConvMode_pPb(
    Int_t     trainConfig                 = 1,
    Int_t     isMC                        = 0,                        //run MC
    TString   photonCutNumberV0Reader     = "",                       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
    Int_t     selectHeavyNeutralMeson     = 0,                        //run eta prime instead of omega
    Int_t     enableQAMesonTask           = 1,                        //enable QA in AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson
    Int_t     enableTriggerMimicking      = 0,                        // enable trigger mimicking
    Bool_t    enableTriggerOverlapRej     = kFALSE,                   // enable trigger overlap rejection
    TString   fileNameInputForWeighting   = "MCSpectraInput.root",    // path to file for weigting input
    Bool_t    doWeighting                 = kFALSE,                   //enable Weighting
    TString   generatorName               = "HIJING",
    Double_t  tolerance                   = -1,
    TString   periodNameV0Reader          = "",                       // period Name for V0Reader
    Int_t     runLightOutput              = 0,                        // run light output option 0: no light output 1: most cut histos stiched off 2: unecessary omega hists turned off as well
    Int_t     prefilterRunFlag            = 1500,                     // flag to change the prefiltering of ESD tracks. See SetHybridTrackCutsAODFiltering() in AliPrimaryPionCuts
    TString   additionalTrainConfig       = "0"                       // additional counter for trainconfig, this has to be always the last parameter
  ) {

  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      cout << "INFO: nothing to do, no definition available!" << endl;
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_ConvMode_pPb running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  Int_t isHeavyIon = 2;
  Int_t neutralPionMode = 0;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvNeutralMesonPiPlPiMiNeutralMeson_ConvMode_pPb_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = photonCutNumberV0Reader.Data();
  TString cutnumberEvent = "80000003";
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

    //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName        = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  AliV0ReaderV1 *fV0ReaderV1  =  NULL;
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    cout << "V0Reader: " << V0ReaderName.Data() << " not found!!"<< endl;
    return;
  } else {
    cout << "V0Reader: " << V0ReaderName.Data() << " found!!"<< endl;
  }

  //================================================
  //========= Add Pion Selector ====================
  TString PionCuts          = "000000200";

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
    //connect input V0Reader
    mgr->ConnectInput (fPionSelector,0,cinput1);

  }

  AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson *task=NULL;
  task= new AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson(Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i",neutralPionMode, trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  if(runLightOutput>1) task->SetLightOutput(kTRUE);
  task->SetTolerance(tolerance);
  AliCutHandlerPCM cuts(13);

  // ******************************************************************************************************
  // ++++++++++++++++++++++++++++++   N A M I N G  C O N V E N T I O N   ++++++++++++++++++++++++++++++++++
  //
  //  config = meson * 1000  + method * 100 + variation
  //
  //  meson:  eta = 0; omega = 1; etaPrime = 2
  //  method: EMC = 0; PHOS = 5;
  //
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ******************************************************************************************************


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                          ETA MESON
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   //************************************************ PCM- EDC analysis 5 TeV pPb *********************************************
  // no event mixing background
  if (trainConfig == 1){ // INT7 run1 & run2
    cuts.AddCutHeavyMesonPCM("80010113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, INT7
  } else if (trainConfig == 2){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80052113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EMC7
  } else if (trainConfig == 3){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80083113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG1
  } else if (trainConfig == 4){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80085113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG2

    // EMCal + EDC triggers
  } else if (trainConfig == 5){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("8008d113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG1
  } else if (trainConfig == 6){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("8008e113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG2

 //************************************************ PCM- PHOS analysis 5 TeV pPb ********************************************
  } else if (trainConfig == 501){ // PHOS  PHI7 run1
    cuts.AddCutHeavyMesonPCM("80062113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000");  // 0-100%

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                          OMEGA MESON
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //************************************************ PCM- EDC analysis 5 TeV pPb *********************************************

  // no event mixing background
  }else if (trainConfig == 1001){ // EMC  INT7 run1 & run2
    cuts.AddCutHeavyMesonPCM("80010113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%
  } else if (trainConfig == 1002){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80052113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EMC7
  } else if (trainConfig == 1003){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80083113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG1
  } else if (trainConfig == 1004){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80085113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG2

    // EMCal + EDC triggers
  } else if (trainConfig == 1005){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("8008d113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG1
  } else if (trainConfig == 1006){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("8008e113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000"); // 0-100%, EG2
 //************************************************ PCM- PHOS analysis 5 TeV pPb ********************************************
  } else if (trainConfig == 1501){ // PHOS  PHI7 run1
    cuts.AddCutHeavyMesonPCM("80062113","00200009f9730000dge0400000","32c51070a","0103603s00000000","0453503000000000");  // 0-100%


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                          ETA PRIME MESON
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //************************************************ PCM- EDC analysis 5 TeV pPb *********************************************

  // no event mixing background
  }else if (trainConfig == 2001){ // EMC  INT7 run1 & run2
    cuts.AddCutHeavyMesonPCM("80010113","00200009f9730000dge0400000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL
  } else if (trainConfig == 2002){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80052113","00200009f9730000dge0400000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EMC7
  } else if (trainConfig == 2003){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80083113","00200009f9730000dge0400000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EG1
  } else if (trainConfig == 2004){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("80085113","00200009f9730000dge0400000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EG2

    // EMCal + EDC triggers
  } else if (trainConfig == 2005){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("8008d113","00200009f9730000dge0400000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EG1
  } else if (trainConfig == 2006){ // EMC EMC triggers
    cuts.AddCutHeavyMesonPCM("8008e113","00200009f9730000dge0400000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL just EMC, EG2

  //************************************************ PCM- PHOS analysis 5 TeV pPb ********************************************
  } else if (trainConfig == 2501){ // PHOS  PHI7 run1
    cuts.AddCutHeavyMesonPCM("80062113","00200009f9730000dge0400000","32c510700","0103603l00000000","0453503000000000"); // 0-100% without NL
  } else {
    Error(Form("GammaConvNeutralMeson_ConvMode_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerNeutralConv! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList  = new TList();
  TList *NeutralPionCutList = new TList();
  TList *MesonCutList = new TList();
  TList *PionCutList  = new TList();

  TList *HeaderList = new TList();
  TObjString *Header1 = new TObjString("pi0_1");
  HeaderList->Add(Header1);
  TObjString *Header3 = new TObjString("eta_2");
  HeaderList->Add(Header3);

  if (periodNameV0Reader.Contains("LHC17g6a2") || periodNameV0Reader.Contains("LHC17g6a3") ){
    TObjString *HeaderPMB = new TObjString("Dpmjet_0");
    TObjString *HeaderP8J = new TObjString("Pythia8JetsGammaTrg_1");
    if (doWeightingPart==4) { // all headers
      HeaderList->Add(HeaderPMB);
      HeaderList->Add(HeaderP8J);
    } else if (doWeightingPart==5) { // only MB header
      HeaderList->Add(HeaderPMB);
    } else { // only JJ header
      HeaderList->Add(HeaderP8J);
    }
  }
  
  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  NeutralPionCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisNeutralPionCuts   = new AliConversionMesonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts   = new AliConversionMesonCuts*[numberOfCuts];
  PionCutList->SetOwner(kTRUE);
  AliPrimaryPionCuts **analysisPionCuts     = new AliPrimaryPionCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    if(runLightOutput>0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    if(runLightOutput>0) analysisCuts[i]->SetLightOutput(kTRUE);
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    if( ! analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data()) ) {
      cout<<"ERROR: analysisCuts [" <<i<<"]"<<endl;
      return 0;
    } else {
      ConvCutList->Add(analysisCuts[i]);
      analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    }

    analysisNeutralPionCuts[i] = new AliConversionMesonCuts();
    if(runLightOutput>0) analysisNeutralPionCuts[i]->SetLightOutput(kTRUE);
    if( ! analysisNeutralPionCuts[i]->InitializeCutsFromCutString((cuts.GetNDMCut(i)).Data()) ) {
      cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      NeutralPionCutList->Add(analysisNeutralPionCuts[i]);
      analysisNeutralPionCuts[i]->SetFillCutHistograms("");
    }

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if(runLightOutput>0) analysisMesonCuts[i]->SetLightOutput(kTRUE);
    if( ! analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data()) ) {
      cout<<"ERROR: analysisMesonCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      MesonCutList->Add(analysisMesonCuts[i]);
      analysisMesonCuts[i]->SetFillCutHistograms("");
    }
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);

    TString cutName( Form("%s_%s_%s_%s_%s",(cuts.GetEventCut(i)).Data(), (cuts.GetPhotonCut(i)).Data(),(cuts.GetPionCut(i)).Data(),(cuts.GetNDMCut(i)).Data(), (cuts.GetMesonCut(i)).Data() ) );
    analysisPionCuts[i] = new AliPrimaryPionCuts();
    analysisPionCuts[i]->SetPrefilterRunFlag(prefilterRunFlag);
    analysisPionCuts[i]->SetPeriodName(periodNameV0Reader);
    if(runLightOutput>0) analysisPionCuts[i]->SetLightOutput(kTRUE);
    if( !analysisPionCuts[i]->InitializeCutsFromCutString((cuts.GetPionCut(i)).Data())) {
      cout<< "ERROR:  analysisPionCuts [ " <<i<<" ] "<<endl;
      return 0;
    } else {
      PionCutList->Add(analysisPionCuts[i]);
      analysisPionCuts[i]->SetFillCutHistograms("",kFALSE,cutName);
    }
  }

  task->SetNDMRecoMode(neutralPionMode);
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(ConvCutList);
  task->SetNeutralPionCutList(NeutralPionCutList);
  task->SetMesonCutList(MesonCutList);
  task->SetPionCutList(PionCutList);

  task->SetMoveParticleAccordingToVertex(kTRUE);

  task->SetSelectedHeavyNeutralMeson(selectHeavyNeutralMeson);

  task->SetDoMesonQA(enableQAMesonTask);

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvNeutralMesonPiPlPiMiNeutralMeson_%i_%i_%i.root",selectHeavyNeutralMeson,neutralPionMode, trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
