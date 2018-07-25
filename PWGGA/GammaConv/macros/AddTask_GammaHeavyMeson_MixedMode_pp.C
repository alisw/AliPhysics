/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Sandro Wenzel                                 *
 * Version 1.0                                                            *
 *                                                                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//***************************************************************************************
//This AddTask is supposed to set up the main task
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskHeavyNeutralMesonToGG.cxx) for
//pPb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
// no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************
class CutHandlerHeavyMesonConvCalo{
  public:
    CutHandlerHeavyMesonConvCalo(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; clusterCutArray[i] = ""; mesonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString clusterCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerHeavyMesonConvCalo: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || clusterCut.Length()!=19 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerHeavyMesonConvCalo: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; clusterCutArray[nCuts]=clusterCut; mesonCutArray[nCuts]=mesonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerHeavyMesonConvCalo: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerHeavyMesonConvCalo: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerHeavyMesonConvCalo: GetClusterCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerHeavyMesonConvCalo: GetMesonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* clusterCutArray;
    TString* mesonCutArray;
};

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaHeavyMeson_MixedMode_pp(  Int_t     selectedMeson                 = 0,                    // select the corresponding meson: 0 pi0, 1 eta, 2 eta'
                                            Int_t     trainConfig                   = 1,                    // change different set of cuts
                                            Int_t     isMC                          = 0,                    // run MC
                                            Int_t     enableQAMesonTask             = 0,                    // enable QA in AliAnalysisTaskGammaConvV1
                                            Int_t     enableQAPhotonTask            = 0,                    // enable additional QA task
                                            TString   fileNameInputForWeighting     = "MCSpectraInput.root",// path to file for weigting input / modified acceptance
                                            Int_t     doWeightingPart               = 0,                    // enable Weighting
                                            TString   periodName                    = "DPMJET",             // generator Name
                                            TString   cutnumberAODBranch            = "800000006008400000001500000",  // cutnumber for AOD branch
                                            Int_t     enableExtMatchAndQA           = 0,                    // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
                                            Bool_t    isUsingTHnSparse              = kTRUE,                // enable or disable usage of THnSparses for background estimation
                                            Bool_t    enableV0findingEffi           = kFALSE,               // enables V0finding efficiency histograms
                                            Bool_t    enableTriggerMimicking        = kFALSE,               // enable trigger mimicking
                                            Bool_t    enableTriggerOverlapRej       = kFALSE,               // enable trigger overlap rejection
                                            Float_t   maxFacPtHard                  = 3,                    // maximum factor between hardest jet and ptHard generated
                                            TString   periodNameV0Reader            = "",                   // period Name for V0Reader
                                            Bool_t    doMultiplicityWeighting       = kFALSE,               // enable multiplicity weights
                                            TString   fileNameInputForMultWeighing  = "Multiplicity.root",  // file for multiplicity weights
                                            TString   periodNameAnchor              = "",                   // anchor period name for mult weighting
                                            Bool_t    enableSortingMCLabels         = kTRUE,                // enable sorting for MC cluster labels
                                            Int_t     runLightOutput                = 0,                    // switch to run light output 0 (disabled), 1 (for CutClasses), 2 (for cutClasses and task)
                                            Bool_t    doPrimaryTrackMatching        = kTRUE,                // enable basic track matching for all primary tracks to cluster
                                            TString   additionalTrainConfig         = "0"                   // additional counter for trainconfig, this has to be always the last parameter
                                          ) {

  TH1S* histoAcc = 0x0;                     // histo for modified acceptance
  TString corrTaskSetting = ""; // select which correction task setting to use
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaHeavyMeson_MixedMode_pp during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0) rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    else{
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      if(tempStr.BeginsWith("MODIFYACC")){
        cout << "INFO: AddTask_GammaHeavyMeson_MixedMode_pp activating 'MODIFYACC'" << endl;
        TString tempType = tempStr;
        tempType.Replace(0,9,"");
        cout << "INFO: connecting to alien..." << endl;
        TGrid::Connect("alien://");
        cout << "done!" << endl;
        TFile *w = TFile::Open(fileNameInputForWeighting.Data());
        if(!w){cout << "ERROR: Could not open file: " << fileNameInputForWeighting.Data() << endl;return;}
        histoAcc = (TH1S*) w->Get(tempType.Data());
        if(!histoAcc) {cout << "ERROR: Could not find histo: " << tempType.Data() << endl;return;}
        cout << "found: " << histoAcc << endl;
      }else if(tempStr.BeginsWith("CF")){
        cout << "INFO: AddTask_GammaHeavyMeson_MixedMode_pp will use custom branch from Correction Framework!" << endl;
        corrTaskSetting = tempStr;
        corrTaskSetting.Replace(0,2,"");
      }
    }
  }
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: AddTask_GammaHeavyMeson_MixedMode_pp running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t isHeavyIon    = 2;
  // meson reco mode: 0 - PCM-PCM, 1 - PCM-Calo, 2 - Calo-Calo
  Int_t mesonRecoMode = 1;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaHeavyMeson_MixedMode_pp_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  Bool_t isMCForOtherTasks = kFALSE;
  if (isMC > 0) isMCForOtherTasks = kTRUE;


  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMCForOtherTasks);
  }

  Printf("here \n");

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton   = "00000008400100001500000000";
  if (  periodNameV0Reader.CompareTo("LHC16f") == 0 || periodNameV0Reader.CompareTo("LHC17g")==0 || periodNameV0Reader.CompareTo("LHC18c")==0 ||
        periodNameV0Reader.CompareTo("LHC17d1") == 0  || periodNameV0Reader.CompareTo("LHC17d12")==0 ||
        periodNameV0Reader.CompareTo("LHC17h3")==0 || periodNameV0Reader.CompareTo("LHC17k1")==0 ||
        periodNameV0Reader.CompareTo("LHC17f8b") == 0 ||
        periodNameV0Reader.CompareTo("LHC16P1JJLowB") == 0 || periodNameV0Reader.CompareTo("LHC16P1Pyt8LowB") == 0 )
    cutnumberPhoton         = "00000088400000000100000000";

  TString cutnumberEvent    = "00000003";
  Bool_t doEtaShift         = kFALSE;
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1(V0ReaderName.Data());
    if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderV1->SetPeriodName(periodNameV0Reader);
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
    fV0ReaderV1->SetProduceV0FindingEfficiency(enableV0findingEffi);

    AliConvEventCuts *fEventCuts=NULL;
    if(cutnumberEvent!=""){
      fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
      if (periodNameV0Reader.CompareTo("") != 0) fEventCuts->SetPeriodEnum(periodNameV0Reader);
      if (runLightOutput > 0) fEventCuts->SetLightOutput(kTRUE);
      if(fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())){
        fEventCuts->DoEtaShift(doEtaShift);
        fV0ReaderV1->SetEventCuts(fEventCuts);
        fEventCuts->SetFillCutHistograms("",kTRUE);
      }
    }

    // Set AnalysisCut Number
    AliConversionPhotonCuts *fCuts=NULL;
    if(cutnumberPhoton!=""){
      fCuts= new AliConversionPhotonCuts(cutnumberPhoton.Data(),cutnumberPhoton.Data());
      fCuts->SetPreSelectionCutFlag(kTRUE);
      fCuts->SetIsHeavyIon(isHeavyIon);
      fCuts->SetV0ReaderName(V0ReaderName);
      if (runLightOutput > 0) fCuts->SetLightOutput(kTRUE);
      if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
        fV0ReaderV1->SetConversionCuts(fCuts);
        fCuts->SetFillCutHistograms("",kTRUE);
      }
    }
    if(inputHandler->IsA()==AliAODInputHandler::Class()){
    // AOD mode
      fV0ReaderV1->AliV0ReaderV1::SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
    }
    fV0ReaderV1->Init();

    AliLog::SetGlobalLogLevel(AliLog::kFatal);

    //connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);

  }

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  AliAnalysisTaskHeavyNeutralMesonToGG *task=NULL;
  task= new AliAnalysisTaskHeavyNeutralMesonToGG(Form("HeavyNeutralMesonToGG_%i_%i_%i", mesonRecoMode, selectedMeson, trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);

  task->SetMesonRecoMode(mesonRecoMode); // meson reco mode: 0 - PCM-PCM, 1 - PCM-Calo, 2 - Calo-Calo
  if (runLightOutput > 1) task->SetLightOutput(kTRUE);
  task->SetDoPrimaryTrackMatching(doPrimaryTrackMatching);

  //create cut handler
  CutHandlerHeavyMesonConvCalo cuts;

  // *****************************************************************************************************
  // pp 2.76 TeV EMC configurations, pi0/eta paper cuts -
  // *****************************************************************************************************
  if (trainConfig == 1){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1) final analysis cuts
    cuts.AddCut("00003113","00200009327000008250400000","1111121057032230000","0163103000000010");
    cuts.AddCut("00051013","00200009327000008250400000","1111121057032230000","0163103000000010");
  } else if (trainConfig == 2){  // LHC13g final analysis cuts
    cuts.AddCut("00010113","00200009327000008250400000","1111121067032230000","0163103000000010"); // INT7
    cuts.AddCut("00052013","00200009327000008250400000","1111121067032230000","0163103000000010"); // EMC7
    cuts.AddCut("00083013","00200009327000008250400000","1111121067032230000","0163103000000010"); // EMCEG1,
    cuts.AddCut("00085013","00200009327000008250400000","1111121067032230000","0163103000000010"); // EMCEG2,
  } else if (trainConfig == 3){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1) final analysis cuts
    cuts.AddCut("00003113","00200009327000008250400000","1111121057032230000","0163103b00000010");
    cuts.AddCut("00051013","00200009327000008250400000","1111121057032230000","0163103b00000010");
  } else if (trainConfig == 4){  // LHC13g final analysis cuts
    cuts.AddCut("00010113","00200009327000008250400000","1111121067032230000","0163103b00000010"); // INT7
    cuts.AddCut("00052013","00200009327000008250400000","1111121067032230000","0163103b00000010"); // EMC7
    cuts.AddCut("00083013","00200009327000008250400000","1111121067032230000","0163103b00000010"); // EMCEG1,
    cuts.AddCut("00085013","00200009327000008250400000","1111121067032230000","0163103b00000010"); // EMCEG2,

  // *****************************************************************************************************
  // pp 8 TeV EMC configurations, pi0/eta paper cuts
  // *****************************************************************************************************
  } else if (trainConfig == 100){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103000000010"); // std
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103000000010"); // std
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103000000010"); // std
  } else if (trainConfig == 101){ // EMCAL clusters 8 TeV LHC12
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103b00000010"); // std
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103b00000010"); // std
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103b00000010"); // std

  // *****************************************************************************************************
  // pp 7 TeV LHC10x EMC configurations
  // *****************************************************************************************************
  } else if (trainConfig == 200){
    cuts.AddCut("00000113","00200009327000008250400000","11111110b7032230000","0163103000000010"); // std
    cuts.AddCut("00000113","00200009327000008250400000","1111111007032230000","0163103000000010"); // std
  } else if (trainConfig == 201){
    cuts.AddCut("00000113","00200009327000008250400000","11111110b7032230000","0163103b00000010"); // std
    cuts.AddCut("00000113","00200009327000008250400000","1111111007032230000","0163103b00000010"); // std

  // *********************************************************************************************************
  // 5 TeV  pp Run2 - EMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 300){ // EMCAL clusters // -50ns, 30ns timing cut
    cuts.AddCut("00010113","00200009327000008250400000","1111100017032230000","0163103000000010");
    cuts.AddCut("00052013","00200009327000008250400000","1111100017032230000","0163103000000010");
    cuts.AddCut("00085013","00200009327000008250400000","1111100017032230000","0163103000000010");
    cuts.AddCut("00083013","00200009327000008250400000","1111100017032230000","0163103000000010");
  } else if (trainConfig == 301){ // EMCAL clusters // -50ns, 30ns timing cut
    cuts.AddCut("00010113","00200009327000008250400000","1111100067032230000","0163103000000010");
    cuts.AddCut("00052013","00200009327000008250400000","1111100067032230000","0163103000000010");
    cuts.AddCut("00085013","00200009327000008250400000","1111100067032230000","0163103000000010");
    cuts.AddCut("00083013","00200009327000008250400000","1111100067032230000","0163103000000010");
  } else if (trainConfig == 302){ // EMCAL clusters // -50ns, 30ns timing cut
    cuts.AddCut("00010113","00200009327000008250400000","1111100017032230000","0163103b00000010");
    cuts.AddCut("00052013","00200009327000008250400000","1111100017032230000","0163103b00000010");
    cuts.AddCut("00085013","00200009327000008250400000","1111100017032230000","0163103b00000010");
    cuts.AddCut("00083013","00200009327000008250400000","1111100017032230000","0163103b00000010");
  } else if (trainConfig == 303){ // EMCAL clusters // -50ns, 30ns timing cut
    cuts.AddCut("00010113","00200009327000008250400000","1111100067032230000","0163103b00000010");
    cuts.AddCut("00052013","00200009327000008250400000","1111100067032230000","0163103b00000010");
    cuts.AddCut("00085013","00200009327000008250400000","1111100067032230000","0163103b00000010");
    cuts.AddCut("00083013","00200009327000008250400000","1111100067032230000","0163103b00000010");

  // *********************************************************************************************************
  // 5 TeV  pp Run2 - DMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 350){ // DCAL clusters
    cuts.AddCut("00010113","00200009327000008250400000","3885500081041220000","0163103000000010"); //
    cuts.AddCut("00010113","00200009327000008250400000","3885500011041220000","0163103000000010"); //
  } else if (trainConfig == 351){ // DCAL clusters // -1000ns, 1000ns timing cut // no NL
    cuts.AddCut("00055113","00200009327000008250400000","3885500011041220000","0163103000000010");
    cuts.AddCut("00089113","00200009327000008250400000","3885500011041220000","0163103000000010");
    cuts.AddCut("0008b113","00200009327000008250400000","3885500011041220000","0163103000000010");
  } else if (trainConfig == 352){ // DCAL clusters // -50ns, 30ns timing cut // no NL
    cuts.AddCut("00055113","00200009327000008250400000","3885500081041220000","0163103000000010");
    cuts.AddCut("00089113","00200009327000008250400000","3885500011041220000","0163103000000010");
    cuts.AddCut("0008b113","00200009327000008250400000","3885500081041220000","0163103000000010");
  } else if (trainConfig == 353){ // DCAL clusters
    cuts.AddCut("00010113","00200009327000008250400000","3885500081041220000","0163103b00000010"); //
    cuts.AddCut("00010113","00200009327000008250400000","3885500011041220000","0163103b00000010"); //
  } else if (trainConfig == 354){ // DCAL clusters // -1000ns, 1000ns timing cut // no NL
    cuts.AddCut("00055113","00200009327000008250400000","3885500011041220000","0163103b00000010");
    cuts.AddCut("00089113","00200009327000008250400000","3885500011041220000","0163103b00000010");
    cuts.AddCut("0008b113","00200009327000008250400000","3885500011041220000","0163103b00000010");
  } else if (trainConfig == 355){ // DCAL clusters // -50ns, 30ns timing cut // no NL
    cuts.AddCut("00055113","00200009327000008250400000","3885500081041220000","0163103b00000010");
    cuts.AddCut("00089113","00200009327000008250400000","3885500011041220000","0163103b00000010");
    cuts.AddCut("0008b113","00200009327000008250400000","3885500081041220000","0163103b00000010");

  // *********************************************************************************************************
  // 13 TeV  pp Run2 - EMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 400){ // EMCAL clusters
    cuts.AddCut("00010113","00200009327000008250400000","1111100017032230000","0163103000000010");
    cuts.AddCut("00010113","00200009327000008250400000","1111100067032230000","0163103000000010");
  } else if (trainConfig == 401){ // EMCAL clusters
    cuts.AddCut("00052013","00200009327000008250400000","1111100017032230000","0163103000000010");
    cuts.AddCut("00085013","00200009327000008250400000","1111100017032230000","0163103000000010");
    cuts.AddCut("00083013","00200009327000008250400000","1111100017032230000","0163103000000010");
  } else if (trainConfig == 402){ // EMCAL clusters
    cuts.AddCut("00052013","00200009327000008250400000","1111100067032230000","0163103000000010");
    cuts.AddCut("00085013","00200009327000008250400000","1111100067032230000","0163103000000010");
    cuts.AddCut("00083013","00200009327000008250400000","1111100067032230000","0163103000000010");
  } else if (trainConfig == 403){ // EMCAL clusters, PCM-EMC NL // -50ns, 30ns timing cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103000000010");
    cuts.AddCut("00052013","00200009327000008250400000","1111111067032230000","0163103000000010");
    cuts.AddCut("00085013","00200009327000008250400000","1111111067032230000","0163103000000010");
    cuts.AddCut("00083013","00200009327000008250400000","1111111067032230000","0163103000000010");
  } else if (trainConfig == 404){ // EMCAL clusters
    cuts.AddCut("00010113","00200009327000008250400000","1111100017032230000","0163103b00000010");
    cuts.AddCut("00010113","00200009327000008250400000","1111100067032230000","0163103b00000010");
  } else if (trainConfig == 405){ // EMCAL clusters
    cuts.AddCut("00052013","00200009327000008250400000","1111100017032230000","0163103b00000010");
    cuts.AddCut("00085013","00200009327000008250400000","1111100017032230000","0163103b00000010");
    cuts.AddCut("00083013","00200009327000008250400000","1111100017032230000","0163103b00000010");
  } else if (trainConfig == 406){ // EMCAL clusters
    cuts.AddCut("00052013","00200009327000008250400000","1111100067032230000","0163103b00000010");
    cuts.AddCut("00085013","00200009327000008250400000","1111100067032230000","0163103b00000010");
    cuts.AddCut("00083013","00200009327000008250400000","1111100067032230000","0163103b00000010");

  // *********************************************************************************************************
  // 13 TeV  DMC configurations
  // *********************************************************************************************************
  } else if (trainConfig == 450){ // DCAL clusters
    cuts.AddCut("00010113","00200009327000008250400000","3885500081041220000","0163103000000010"); //
    cuts.AddCut("00010113","00200009327000008250400000","3885500011041220000","0163103000000010"); //
  } else if (trainConfig == 451){ // DCAL clusters // -1000ns, 1000ns timing cut // no NL
    cuts.AddCut("00055113","00200009327000008250400000","3885500011041220000","0163103000000010");
    cuts.AddCut("00089113","00200009327000008250400000","3885500011041220000","0163103000000010");
    cuts.AddCut("0008b113","00200009327000008250400000","3885500011041220000","0163103000000010");
  } else if (trainConfig == 452){ // DCAL clusters // -50ns, 30ns timing cut // no NL
    cuts.AddCut("00055113","00200009327000008250400000","3885500081041220000","0163103000000010");
    cuts.AddCut("00089113","00200009327000008250400000","3885500011041220000","0163103000000010");
    cuts.AddCut("0008b113","00200009327000008250400000","3885500081041220000","0163103000000010");
  } else if (trainConfig == 453){ // DCAL clusters
    cuts.AddCut("00010113","00200009327000008250400000","3885500081041220000","0163103b00000010"); //
    cuts.AddCut("00010113","00200009327000008250400000","3885500011041220000","0163103b00000010"); //
  } else if (trainConfig == 454){ // DCAL clusters // -1000ns, 1000ns timing cut // no NL
    cuts.AddCut("00055113","00200009327000008250400000","3885500011041220000","0163103b00000010");
    cuts.AddCut("00089113","00200009327000008250400000","3885500011041220000","0163103b00000010");
    cuts.AddCut("0008b113","00200009327000008250400000","3885500011041220000","0163103b00000010");
  } else if (trainConfig == 455){ // DCAL clusters // -50ns, 30ns timing cut // no NL
    cuts.AddCut("00055113","00200009327000008250400000","3885500081041220000","0163103b00000010");
    cuts.AddCut("00089113","00200009327000008250400000","3885500011041220000","0163103b00000010");
    cuts.AddCut("0008b113","00200009327000008250400000","3885500081041220000","0163103b00000010");

  // *****************************************************************************************************
  // 2.76 TeV pp Run1 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 500) { //PHOS clusters
    cuts.AddCut("00003113","00200009327000008250400000","2444400041033200000","0163103000000010");
    cuts.AddCut("00061113","00200009327000008250400000","2444400041033200000","0163103000000010");
    // LHC13g & LHC12x
  } else if (trainConfig == 501) { //PHOS clusters
    cuts.AddCut("00010113","00200009327000008250400000","2444400041033200000","0163103000000010");
    cuts.AddCut("00062113","00200009327000008250400000","2444400041033200000","0163103000000010");
  } else if (trainConfig == 502) { //PHOS clusters
    cuts.AddCut("00003113","00200009327000008250400000","2444400041033200000","0163103b00000010");
    cuts.AddCut("00061113","00200009327000008250400000","2444400041033200000","0163103b00000010");
    // LHC13g & LHC12x
  } else if (trainConfig == 503) { //PHOS clusters
    cuts.AddCut("00010113","00200009327000008250400000","2444400041033200000","0163103b00000010");
    cuts.AddCut("00062113","00200009327000008250400000","2444400041033200000","0163103b00000010");

  // *****************************************************************************************************
  // 8 TeV pp Run1 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 600){ // PHOS clusters
    cuts.AddCut("00010113","00200009327000008250400000","2444400072023200000","0163103000000010"); // 600 MeV cluster min energy, t<|30|ns
    cuts.AddCut("00062113","00200009327000008250400000","2444400072023200000","0163103000000010"); // 600 MeV cluster min energy, t<|30|ns
  } else if (trainConfig == 601){ // PHOS clusters
    cuts.AddCut("00010113","00200009327000008250400000","2444400072023200000","0163103b00000010"); // 600 MeV cluster min energy, t<|30|ns
    cuts.AddCut("00062113","00200009327000008250400000","2444400072023200000","0163103b00000010"); // 600 MeV cluster min energy, t<|30|ns

  // *****************************************************************************************************
  // 7 TeV pp Run1 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 700){
    cuts.AddCut("00000113","00200009327000008250400000","2444400000013300000","0163103000000010"); // QA
  } else if (trainConfig == 701){
    cuts.AddCut("00000113","00200009327000008250400000","2444400040013300000","0163103000000010"); // 100ns timing cut, no track matching
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","0163103000000010"); // 100ns timing cut
  } else if (trainConfig == 702){ // train config for bad channels and NonLin Variation
    cuts.AddCut("00000113","00200009327000008250400000","2444411000013300000","0163103000000010"); // own constant ConvCalo NonLin
    cuts.AddCut("00000113","00200009327000008250400000","2444421000013300000","0163103000000010"); // ex
  } else if (trainConfig == 703){
    cuts.AddCut("00000113","00200009327000008250400000","2444400000013300000","0163103b00000010"); // QA
  } else if (trainConfig == 704){
    cuts.AddCut("00000113","00200009327000008250400000","2444400040013300000","0163103b00000010"); // 100ns timing cut, no track matching
    cuts.AddCut("00000113","00200009327000008250400000","2444400043013300000","0163103b00000010"); // 100ns timing cut
  } else if (trainConfig == 705){ // train config for bad channels and NonLin Variation
    cuts.AddCut("00000113","00200009327000008250400000","2444411000013300000","0163103b00000010"); // own constant ConvCalo NonLin
    cuts.AddCut("00000113","00200009327000008250400000","2444421000013300000","0163103b00000010"); // ex

  // *****************************************************************************************************
  // 5 TeV pp Run2 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 800){ // PHOS clusters with larger acceptance
    cuts.AddCut("00010113","00200009327000008250400000","2446600040013300000","0163103000000010"); // INT7
    cuts.AddCut("00062113","00200009327000008250400000","2446600040013300000","0163103000000010"); // PHI7
  } else if (trainConfig == 801){ // PHOS clusters with larger acceptance
    cuts.AddCut("00010113","00200009327000008250400000","2446600040013300000","0163103b00000010"); // INT7
    cuts.AddCut("00062113","00200009327000008250400000","2446600040013300000","0163103b00000010"); // PHI7

  // *****************************************************************************************************
  // 13 TeV pp Run2 - PHOS configurations
  // *****************************************************************************************************
  } else if (trainConfig == 900){ // PHOS clusters with larger acceptance
    cuts.AddCut("00010113","00200009327000008250400000","2446600000013300000","0163103000000010"); // QA
    cuts.AddCut("00010113","00200009327000008250400000","2446600040013300000","0163103000000010"); // QA, 100ns timing
    cuts.AddCut("00010113","00200009327000008250400000","2446600043013300000","0163103000000010"); // QA, 100ns timing, TM on with default EMC params
  } else if (trainConfig == 901){ // PHI7
    cuts.AddCut("00062113","00200009327000008250400000","2446600000013300000","0163103000000010"); // QA
    cuts.AddCut("00062113","00200009327000008250400000","2446600040013300000","0163103000000010"); // QA, 100ns timing
    cuts.AddCut("00062113","00200009327000008250400000","2446600043013300000","0163103000000010"); // QA, 100ns timing, TM on with default EMC params
  } else if (trainConfig == 902){ // PHOS clusters with larger acceptance
    cuts.AddCut("00010113","00200009327000008250400000","2446600000013300000","0163103b00000010"); // QA
    cuts.AddCut("00010113","00200009327000008250400000","2446600040013300000","0163103b00000010"); // QA, 100ns timing
    cuts.AddCut("00010113","00200009327000008250400000","2446600043013300000","0163103b00000010"); // QA, 100ns timing, TM on with default EMC params
  } else if (trainConfig == 903){ // PHI7
    cuts.AddCut("00062113","00200009327000008250400000","2446600000013300000","0163103b00000010"); // QA
    cuts.AddCut("00062113","00200009327000008250400000","2446600040013300000","0163103b00000010"); // QA, 100ns timing
    cuts.AddCut("00062113","00200009327000008250400000","2446600043013300000","0163103b00000010"); // QA, 100ns timing, TM on with default EMC params

  } else {
    Error(Form("HeavyNeutralMesonToGG_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerHeavyMesonConvCalo! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  TList *EventCutList   = new TList();
  TList *ConvCutList    = new TList();
  TList *ClusterCutList = new TList();
  TList *MesonCutList   = new TList();

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *HeaderList     = new TList();
  if (doWeightingPart==1) {
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
  }
  if (doWeightingPart==2){
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }
  if (doWeightingPart==3) {
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts      = new AliConversionPhotonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    //create AliCaloTrackMatcher instance, if there is none present
    TString caloCutPos = cuts.GetClusterCut(i);
    caloCutPos.Resize(1);
    TString TrackMatcherName = Form("CaloTrackMatcher_%s",caloCutPos.Data());
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi());
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      fTrackMatcher->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i] = new AliConvEventCuts();

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    TString triggerString       = cuts.GetEventCut(i);
    triggerString               = triggerString(3,2);

    dataInputMultHisto          = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto            = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());

    if (doMultiplicityWeighting){
      cout << "enabling mult weighting" << endl;
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameInputForMultWeighing, dataInputMultHisto, mcInputMultHisto );
    }

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    if (runLightOutput > 0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);


    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    if (runLightOutput > 0) analysisCuts[i]->SetLightOutput(kTRUE);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i] = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if (runLightOutput > 0) analysisClusterCuts[i]->SetLightOutput(kTRUE);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if (runLightOutput > 0) analysisMesonCuts[i]->SetLightOutput(kTRUE);
    analysisMesonCuts[i]->SetRunningMode(2);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetMesonType(selectedMeson);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetDoClusterQA(1);  //Attention new switch small for Cluster QA
  task->SetUseTHnSparse(isUsingTHnSparse);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("HeavyNeutralMesonToGG_%i_%i_%i", mesonRecoMode, selectedMeson, trainConfig)
                                                        :  Form("HeavyNeutralMesonToGG_%i_%i_%i_%s", mesonRecoMode, selectedMeson, trainConfig, corrTaskSetting.Data()), TList::Class(),
                        AliAnalysisManager::kOutputContainer,Form("HeavyNeutralMesonToGG_%i_%i.root",mesonRecoMode,trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
