/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Daniel Mühlheim, Marvin Hemmer                *
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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskOmegaToPiZeroGamma.cxx) for
//pp together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
//no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************
class CutHandlerPiZeroGamma{
  public:
    CutHandlerPiZeroGamma(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts]; neutralPionCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; clusterCutArray[i] = ""; neutralPionCutArray[i] = ""; mesonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString clusterCut, TString neutralPionCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerPiZeroGamma: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || clusterCut.Length()!=19 || neutralPionCut.Length()!=16 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerPiZeroGamma: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; clusterCutArray[nCuts]=clusterCut; neutralPionCutArray[nCuts]=neutralPionCut; mesonCutArray[nCuts]=mesonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerPiZeroGamma: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerPiZeroGamma: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerPiZeroGamma: GetClusterCut wrong index i" << endl;return "";}}
    TString GetNeutralPionCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return neutralPionCutArray[i]; else {cout << "ERROR in CutHandlerPiZeroGamma: GetNeutralPionCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerPiZeroGamma: GetMesonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* clusterCutArray;
    TString* neutralPionCutArray;
    TString* mesonCutArray;
};

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_OmegaToPiZeroGamma_pp(
                                Int_t     trainConfig                   = 1,                      // change different set of cuts
                                Int_t     isMC                          = 0,                      // run MC
                                TString   photonCutNumberV0Reader       = "",                     // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
                                Int_t     enableQAMesonTask             = 1,                      // enable QA in AliAnalysisTaskGammaConvV1
                                Int_t     enableQAPhotonTask            = 1,                      // enable additional QA task
                                Int_t     DoPiZeroGammaAngleCut         = kFALSE,                 // flag for enabling cut on pi0-gamma angle
                                Double_t  lowerFactor                   = 0.75,                   // scale factor for lower limit in pi0-gamma angle cut
                                Double_t  upperFactor                   = 2.5,                    // scale factor for upper limit in pi0-gamma angle cut
                                TString   cutnumberAODBranch            = "000000006008400001001500000",
                                Int_t     enableExtMatchAndQA           = 0,                      // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
                                Bool_t    enableV0findingEffi           = kFALSE,                 // enables V0finding efficiency histograms
                                Int_t     enableTriggerMimicking        = 0,                      // enable trigger mimicking
                                Bool_t    enableTriggerOverlapRej       = kFALSE,                 // enable trigger overlap rejection
                                TString   settingMaxFacPtHard           = "3.",                   // maximum factor between hardest jet and ptHard generated
                                TString   periodNameV0Reader            = "",                     // period Name for V0 Reader
                                Bool_t    doMultiplicityWeighting       = kFALSE,                 // enable multiplicity weights
                                TString   fileNameInputForMultWeighing  = "Multiplicity.root",    // file for multiplicity weights
                                TString   periodNameAnchor              = "",                     // anchor period name for mult weighting
                                Bool_t    enableSortingMCLabels         = kTRUE,                  // enable sorting for MC cluster labels
                                Bool_t    runLightOutput                = kFALSE,                 // switch to run light output (only essential histograms for afterburner)
                                Bool_t    doSmear                       = kFALSE,                 // switches to run user defined smearing
                                Double_t  bremSmear                     = 1.,
                                Double_t  smearPar                      = 0.,                     // conv photon smearing params
                                Double_t  smearParConst                 = 0.,                     // conv photon smearing params
                                Bool_t    usePtDepSelectionWindowCut    = kFALSE,                 // use pt dependent meson selection window cut
                                Int_t     PhotonSelectionMode           = 0,                      // 0 none, 1 normal (Cal-Cal, PCM-PCM), 2 strict (normal + Cal-PCM)
                                Bool_t    useDalitzCut                  = kFALSE,                 // flag to cut on Dalitz distribution or not
                                TString   additionalTrainConfig         = "0"                     // additional counter for trainconfig, this has to be always the last parameter
              ) {


  // CHANGED // TODO From AddTask_GammaHeavyMeson_CaloMode_pp
  AliCutHandlerPCM cutsPCM(13);


  TString addTaskName                 = "AddTask_OmegaToPiZeroGamma_pp";
  TString sAdditionalTrainConfig      = cutsPCM.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  TString corrTaskSetting             = cutsPCM.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF", "", addTaskName);
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
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
        cout << Form("INFO: AddTask_OmegaToPiZeroGamma_pp will use running mode '%i' for the TrackMatcher!",trackMatcherRunningMode) << endl;
      }
    }
  }
  // TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  // if (sAdditionalTrainConfig.Atoi() > 0){
  //   trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
  //   cout << "INFO: AddTask_OmegaToPiZeroGamma_pp running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  // }

  TObjArray *rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_OmegaToPiZeroGamma_pp during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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


  Int_t isHeavyIon = 0;

  // fReconMethod = first digit of trainconfig
  Int_t ReconMethod = trainConfig/1000;
  if(ReconMethod == 6) ReconMethod = 2;                                         // use EDC-EDC-EDC for 6000 trainconfigs

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_OmegaToPiZeroGamma_%i_%i", trainConfig, isMC), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  Bool_t isMCForOtherSettings = 0;
  if (isMC > 0) isMCForOtherSettings = 1;

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = photonCutNumberV0Reader.Data();
  TString cutnumberEvent = "00000003";
  Bool_t doEtaShift = kFALSE;
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

    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
    }

    AliConvEventCuts *fEventCuts=NULL;
    if(cutnumberEvent!=""){
      fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
      fEventCuts->SetV0ReaderName(V0ReaderName);
      fEventCuts->SetLightOutput(runLightOutput);
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
      fCuts->SetLightOutput(runLightOutput);
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
  AliAnalysisTaskOmegaToPiZeroGamma *task=NULL;
  task= new AliAnalysisTaskOmegaToPiZeroGamma(Form("OmegaToPiZeroGamma_%i_%i", trainConfig, isMC));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetV0ReaderName(V0ReaderName);
  task->SetReconMethod(ReconMethod);
  task->SetDoPiZeroGammaAngleCut(DoPiZeroGammaAngleCut);
  task->SetlowerFactor(lowerFactor);
  task->SetupperFactor(upperFactor);
  task->SetTrackMatcherRunningMode(trackMatcherRunningMode);
  task->SetPhotonSelectionMode(PhotonSelectionMode);

  //create cut handler
  CutHandlerPiZeroGamma cuts;

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"

  // ************************************* EMCAL cuts ****************************************************

  // 7 TeV
  // cuts for ReconMethod==0 PCM-Cal-Cal
  if(trainConfig == 1){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // pion mass (0.08,0.145), last string is not used as of now

  } else if(trainConfig == 2){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 3){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 4){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 51){ // std 8TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 52){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 53){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 54){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
  } else if( trainConfig == 60) {
    //std MB 13TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00085113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00083113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
  } else if( trainConfig == 61) {
    //only MB 13TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  // cuts for ReconMethod==1 PCM-Cal-PCM
  } else if(trainConfig == 1001){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // pion mass (0.08,0.145)

  } else if(trainConfig == 1002){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 1003){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 1004){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 1051){ // std 8TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 1052){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 1053){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 1054){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
  } else if(trainConfig == 1060) {
    //std MB 13TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00085113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00083113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
  } else if(trainConfig == 1061) {
    //only MB 13TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
  } else if(trainConfig == 1062) {
    //only MB 13TeV EMcal + Dcal
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0");
  } else if(trainConfig == 1063) {
    //EG2 13TeV EMcal + Dcal
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0");
  } else if(trainConfig == 1064) {
    //EG1 13TeV EMcal + Dcal
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0");

  // cuts for ReconMethod==2 Cal-Cal-Cal
  } else if(trainConfig == 2001){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // pion mass (0.08,0.145)

  } else if(trainConfig == 2002){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 2003){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 2004){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 2051){ // std 8TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 2052){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 2053){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 2054){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 2055){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111100067032230000","0163103100000010","0163103000000010"); // same as standard but without non lin
    cuts.AddCut("00052113","00200009327000008250400000","1111100067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111100067032230000","0163103100000010","0163103000000010");


    /////////////////////////////////
    // 13 TeV
    /////////////////////////////////

    // EMCal + DCal
  } else if( trainConfig == 2060) {
    // MB 13TeV EMcal + DCal
    cuts.AddCut("00010113","00200009327000008250400000","411791106f032230000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00052113","00200009327000008250400000","411791106f032230000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00085113","00200009327000008250400000","411791106f032230000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00083113","00200009327000008250400000","411791106f032230000","01631031000000d0","01631030000000d0");
  } else if( trainConfig == 2061) {
    // MB 13TeV only EMCal
    cuts.AddCut("00010113","00200009327000008250400000","1111100067032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2062) {
    // MB std NL 12 EMCal + DCal
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2063) {
    // MB 13TeV EMCal + DCal, Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 2064) {
    // MB 13TeV EMCal + DCal Background Variation (Swapping Method by Joshua)
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation (omega)
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints (omega)
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints (omega)
  } else if( trainConfig == 2065) {
    // MB 13TeV EMCal + DCal, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631071000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631051000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631041000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631081000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 2066) {
    // MB 13TeV EMCal + DCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 2067) {
    // MB 13TeV EMCal + DCal Background Variation (Swapping Method by Joshua), AP-like cut
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, no AP like cut
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031010000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 1 sigma
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031020000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 2 sigma
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031030000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 3 sigma
  } else if( trainConfig == 2068) {
    // MB 13TeV EMCal + DCal Background Variation (Swapping around Pi0, Method by Joshua)
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints (omega)
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0y631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation (Pi0)
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0z631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints (Pi0)
  } else if( trainConfig == 2071) {
    // EG2 13TeV EMcal
    cuts.AddCut("0008e113","00200009327000008250400000","1111100067032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2072) {
    // EG2 std NL 12 EMCal + DCal
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2073) {
    // EG2 13TeV EMcal + Dcal, Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 2074) {
    // EG2 13TeV EMcal + Dcal, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 2075) {
    // EG2 13TeV EMCal + DCal, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631071000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631051000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631041000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631081000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 2076) {
    // EG2 13TeV EMCal + DCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 2077) {
    // MB 13TeV EMCal + DCal Background Variation (Swapping Method by Joshua), AP-like cut
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, no AP like cut
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031010000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 1 sigma
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031020000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 2 sigma
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031030000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 3 sigma
  } else if( trainConfig == 2078) {
    // EG2 13TeV EMcal + Dcal, Background Variation (Swapping around Pi0, Method by Joshua)
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints (omega)
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0y631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation (Pi0)
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0z631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints (Pi0)
  } else if( trainConfig == 2081) {
    // EG1 13TeV EMcal
    cuts.AddCut("0008d113","00200009327000008250400000","1111100067032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2082) {
    // EG1 13TeV EMcal + Dcal
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2083) {
    // EG1 13TeV EMcal + Dcal, Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 2084) {
    // EG1 13TeV EMcal + Dcal, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 2085) {
    // EG1 13TeV EMCal + DCal, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631071000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631051000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631041000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631081000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 2086) {
    // EG1 13TeV EMCal + DCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 2087) {
    // EG1 13TeV EMCal + DCal Background Variation (Swapping Method by Joshua), AP-like cut
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, no AP like cut
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031010000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 1 sigma
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031020000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 2 sigma
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031030000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 3 sigma
  } else if( trainConfig == 2088) {
    // EG2 13TeV EMcal + Dcal, Background Variation (Swapping around Pi0, Method by Joshua)
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints (omega)
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0y631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation (Pi0)
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0z631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints (Pi0)

  } else if( trainConfig == 2160) {
    // MB 13TeV PHOS
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00052113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00085113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00083113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631030000000d0");
  } else if( trainConfig == 2162) {
    // MB std NL 12 PHOS
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2163) {
    // MB 13TeV PHOS, Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 2164) {
    // MB 13TeV PHOS Background Variation (Swapping Method by Joshua)
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 2165) {
    // MB 13TeV PHOS, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631071000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631051000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631041000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631081000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 2166) {
    // MB 13TeV PHOS, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 2167) {
    // MB 13TeV EMCal + DCal Background Variation (Swapping Method by Joshua), AP-like cut
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, no AP like cut
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031010000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 1 sigma
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031020000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 2 sigma
    cuts.AddCut("00010113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031030000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 3 sigma
  } else if( trainConfig == 2172) {
    // EG2 std NL 12 PHOS
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2173) {
    // EG2 13TeV PHOS, Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 2174) {
    // EG2 13TeV PHOS, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 2175) {
    // EG2 13TeV PHOS, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631071000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631051000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631041000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631081000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 2176) {
    // EG2 13TeV PHOS, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 2177) {
    // EG1 13TeV PHOS Background Variation (Swapping Method by Joshua), AP-like cut
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, no AP like cut
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031010000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 1 sigma
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031020000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 2 sigma
    cuts.AddCut("0008e113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031030000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 3 sigma
  } else if( trainConfig == 2182) {
    // EG1 13TeV PHOS
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2183) {
    // EG1 13TeV PHOS, Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 2184) {
    // EG1 13TeV PHOS, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 2185) {
    // EG1 13TeV PHOS, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631071000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631051000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631041000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631081000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 2186) {
    // EG1 13TeV PHOS, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 2187) {
    // EG1 13TeV PHOS Background Variation (Swapping Method by Joshua), AP-like cut
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, no AP like cut
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031010000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 1 sigma
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031020000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 2 sigma
    cuts.AddCut("0008d113","00200009327000008250400000","24466190sa01cc00000","01631031000000d0","0x631031030000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 3 sigma

  } else if( trainConfig == 2262) {
    // MB std NL 12 EMCal
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2263) {
    // MB 13TeV EMCal, Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 2264) {
    // MB 13TeV EMCal Background Variation (Swapping Method by Joshua)
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 2265) {
    // MB 13TeV EMCal, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","01631071000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","01631051000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","01631041000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","01631081000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 2266) {
    // MB 13TeV EMCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 2272) {
    // EG2 std NL 12 EMCal
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2273) {
    // EG2 13TeV EMcal, Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 2274) {
    // EG2 13TeV EMcal, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 2275) {
    // EG2 13TeV EMCal, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","01631071000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","01631051000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","01631041000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","01631081000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 2276) {
    // EG2 13TeV EMCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 2282) {
    // EG1 13TeV EMcal
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2283) {
    // EG1 13TeV EMcal, Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 2284) {
    // EG1 13TeV EMcal, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 2285) {
    // EG1 13TeV EMCal, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","01631071000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","01631051000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","01631041000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","01631081000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 2286) {
    // EG1 13TeV EMCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","111791206f032230000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega

  } else if( trainConfig == 2360) {
    // MB 13TeV DCal
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00052113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00085113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00083113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631030000000d0");
  } else if( trainConfig == 2362) {
    // MB std NL 12 DCal
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2363) {
    // MB 13TeV DCal, Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 2364) {
    // MB 13TeV DCal Background Variation (Swapping Method by Joshua)
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 2365) {
    // MB 13TeV DCal, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","01631071000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","01631051000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","01631041000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","01631081000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 2366) {
    // MB 13TeV DCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 2372) {
    // EG2 std NL 12 PHOS
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2373) {
    // EG2 13TeV DCal, Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 2374) {
    // EG2 13TeV DCal, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 2375) {
    // EG2 13TeV DCal, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","01631071000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","01631051000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","01631041000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","01631081000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 2376) {
    // EG2 13TeV DCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 2382) {
    // EG1 13TeV DCal
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 2383) {
    // EG1 13TeV DCal, Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 2384) {
    // EG1 13TeV DCal, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 2385) {
    // EG1 13TeV DCal, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","01631071000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","01631051000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","01631041000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","01631081000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 2386) {
    // EG1 13TeV DCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","311791206f032230000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega

  } else if( trainConfig == 2463) {
    // MB 13TeV EMCal + DCal, for strong Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout


  // cuts for ReconMethod==3 Cal-Cal-PCM
  } else if(trainConfig == 3001){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // pion mass (0.08,0.145)

  } else if(trainConfig == 3002){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 3003){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 3004){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 3051){ // std 8TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 3052){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 3053){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 3054){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 3055){ // std 8TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111100067032230000","0163103100000010","0163103000000010"); // same as standard but wo non lin
    cuts.AddCut("00052113","00200009327000008250400000","1111100067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111100067032230000","0163103100000010","0163103000000010");

  /////////////////////////////////
  // 13 TeV
  /////////////////////////////////
  } else if( trainConfig == 3060) {
    // MB 13TeV PCM Photon with EMcal + DCal Pi0
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791106f032230000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00052113","0dm00009f9730000dge0404000","411791106f032230000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00085113","0dm00009f9730000dge0404000","411791106f032230000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00083113","0dm00009f9730000dge0404000","411791106f032230000","01631031000000d0","01631030000000d0");
  } else if( trainConfig == 3061) {
    // MB 13TeV PCM Photon with EMcal Pi0
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","1111100067032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 3062) {
    // MB std NL 12 PCM Photon with EMcal + DCal Pi0
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 3063) {
    // MB 13TeV PCM Photon with EMcal + DCal Pi0, Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 3064) {
    // MB 13TeV PCM Photon with EMcal + DCal Pi0, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","0r631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","0v631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","0x631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 3065) {
    // MB 13TeV PCM Photon with EMcal + DCal, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163107d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163105d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163104d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163108d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 3066) {
    // MB 13TeV PCM Photon with EMcal + DCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631071000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631051000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631041000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631081000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 3067) {
    // MB 13TeV PCM Photon with EMCal + DCal Background Variation (Swapping Method by Joshua), AP-like cut
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, no AP like cut
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","0x631031010000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 1 sigma
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","0x631031020000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 2 sigma
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","0x631031030000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 3 sigma
  } else if( trainConfig == 3071) {
    // EG2 13TeV PCM Photon with EMcal Pi0
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","1111100067032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 3072) {
    // EG2 std NL 12 PCM Photon with EMcal + DCal Pi0
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 3073) {
    // EG2 13TeV EMcal + Dcal, Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 3074) {
    // EG2 13TeV EMcal + Dcal, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","01631034000000d0","0r631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","01631034000000d0","0v631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","01631034000000d0","0x631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 3075) {
    // EG2 13TeV PCM Photon with EMcal + DCal, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","0163107d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","0163105d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","0163104d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","0163108d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 3076) {
    // EG2 13TeV PCM Photon with EMcal + DCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631071000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631051000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631041000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631081000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 3077) {
    // EG2 13TeV PCM Photon with EMCal + DCal Background Variation (Swapping Method by Joshua), AP-like cut
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, no AP like cut
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","01631031000000d0","0x631031010000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 1 sigma
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","01631031000000d0","0x631031020000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 2 sigma
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","411791206f032230000","01631031000000d0","0x631031030000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 3 sigma
  } else if( trainConfig == 3081) {
    // EG1 13TeV EMcal
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","1111100067032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 3082) {
    // EG1 13TeV EMcal + Dcal
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 3083) {
    // EG1 13TeV EMcal + Dcal, Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 3084) {
    // EG1 13TeV EMcal + Dcal, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","0r631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","0v631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","0x631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 3085) {
    // EG1 13TeV PCM Photon with EMcal + DCal, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163107d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163105d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163104d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163108d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 3086) {
    // EG1 13TeV PCM Photon with EMcal + DCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631071000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631051000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631041000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","01631081000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 3087) {
    // EG1 13TeV PCM Photonwith  EMCal + DCal Background Variation (Swapping Method by Joshua), AP-like cut
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, no AP like cut
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","0x631031010000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 1 sigma
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","0x631031020000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 2 sigma
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","411791206f032230000","0163103w000000d0","0x631031030000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 3 sigma

    // PCM(gamma) - PHOS (pi0)
  } else if( trainConfig == 3162) {
    // MB std NL 12
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 3163) {
    // MB 13TeV, Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 3164) {
    // MB 13TeV, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","0r631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","0v631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","0x631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 3165) {
    // MB 13TeV, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163107d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163105d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163104d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163108d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 3166) {
    // MB 13TeV, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631071000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631051000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631041000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("00010113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631081000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 3172) {
    // EG2 std NL 12
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 3173) {
    // EG2 13TeV, Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 3174) {
    // EG2 13TeV, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","0r631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","0v631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","0x631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 3175) {
    // EG2 13TeV, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163107d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163105d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163104d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163108d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 3176) {
    // EG2 13TeV, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631071000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631051000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631041000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("0008e113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631081000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 3182) {
    // EG1 13TeV
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 3183) {
    // EG1 13TeV, Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 3184) {
    // EG1 13TeV, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","0r631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","0v631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","0x631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 3185) {
    // EG1 13TeV, Pi0 selection plus Gamma dropout, alpha cut pi0
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163107d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for pi0
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163105d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for pi0
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163104d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for pi0
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163108d000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for pi0
  } else if( trainConfig == 3186) {
    // EG1 13TeV, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631071000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631051000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631041000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("0008d113","0dm00009f9730000dge0404000","24466190sa01cc00000","0163103w000000d0","01631081000000d0"); // 4 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega


  // cuts for ReconMethod==4 PCM-PCM-Cal
  } else if(trainConfig == 4001){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // pion mass (0.08,0.145)

  } else if(trainConfig == 4002){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 4003){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 4004){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 4051){ // std 8TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 4052){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 4053){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 4054){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 4054){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
  } else if( trainConfig == 4060) {
    //std MB 13TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00085113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00083113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
  } else if( trainConfig == 4061) {
    //only MB 13TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
  } else if( trainConfig == 4062) {
    //only MB 13TeV EMcal + DCal
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0163103000000010"); // NL 12
  } else if( trainConfig == 4063) {
    //EG2 13TeV EMcal + Dcal
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0"); // NL 12
  } else if( trainConfig == 4064) {
    //EG1 13TeV EMcal + Dcal
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0"); // NL 12

  // cuts for ReconMethod==5 PCM-PCM-PCM
  } else if(trainConfig == 5001){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // pion mass (0.08,0.145)

  } else if(trainConfig == 5002){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 5003){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 5004){ // EMCAL clusters pp 7 TeV
    cuts.AddCut("00000113","00200009327000008250400000","1111111017032230000","0163103a00000010","0163103000000010"); // same as standard, for varying angle cut

  } else if(trainConfig == 5051){ // std 8TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 5052){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 5053){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if(trainConfig == 5054){ // same as std 8TeV, for varying angle cut
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00081113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");

  } else if( trainConfig == 5060) {
    //std MB 13TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00052113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00085113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
    cuts.AddCut("00083113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
  } else if( trainConfig == 5061) {
    //only MB 13TeV
    cuts.AddCut("00010113","00200009327000008250400000","1111111067032230000","0163103100000010","0163103000000010");
  } else if( trainConfig == 5062) {
    //only MB 13TeV EMcal + DCal
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0163103000000010"); // NL 12
  } else if( trainConfig == 5063) {
    //EG2 13TeV EMcal + Dcal
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0"); // NL 12
  } else if( trainConfig == 5064) {
    //EG1 13TeV EMcal + Dcal
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0"); // NL 12

    /////////////////////////////////
    // 13 TeV 2nd batch of trainconfigs
    /////////////////////////////////
  } else if( trainConfig == 6034) {
    // EG2 13TeV EMcal + Dcal, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f030230000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f030230000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f030230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 6044) {
    // EG1 13TeV EMcal + Dcal, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f030230000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f030230000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f030230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints
    // EMCal + DCal
  } else if( trainConfig == 6060) {
    // MB 13TeV EMcal + DCal
    cuts.AddCut("00010113","00200009327000008250400000","411791106f032230000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00052113","00200009327000008250400000","411791106f032230000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00085113","00200009327000008250400000","411791106f032230000","01631031000000d0","01631030000000d0");
    cuts.AddCut("00083113","00200009327000008250400000","411791106f032230000","01631031000000d0","01631030000000d0");
  } else if( trainConfig == 6061) {
    // MB 13TeV only EMCal
    cuts.AddCut("00010113","00200009327000008250400000","1111100067032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 6062) {
    // MB std NL 12 EMCal + DCal
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // background scheme swapping method trough TGPS with constraints (omega)
  } else if( trainConfig == 6063) {
    // MB 13TeV EMCal + DCal, Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 6064) {
    // MB 13TeV EMCal + DCal Background Variation (Swapping Method by Joshua)
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0r631031000000d0"); // background scheme swapping method trough roation (omega)
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0v631031000000d0"); // background scheme swapping method trough TGPS without constraints (omega)
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // background scheme swapping method trough TGPS with constraints (omega)
  } else if( trainConfig == 6065) {
    // MB 13TeV EMCal + DCal, Pi0 selection plus Gamma dropout, pi0 sigma range
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","0163103u000000d0","0x631031000000d0"); // 1 sigma Pi0, background scheme swapping method trough TGPS with constraints
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0, background scheme swapping method trough TGPS with constraints
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","0163103v000000d0","0x631031000000d0"); // 3 sigma Pi0, background scheme swapping method trough TGPS with constraints
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","0163103w000000d0","0x631031000000d0"); // 4 sigma Pi0, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 6066) {
    // MB 13TeV EMCal + DCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 6067) {
    // MB 13TeV EMCal + DCal Background Variation (Swapping Method by Joshua), AP-like cut
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, no AP like cut
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031010000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 1 sigma
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031020000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 2 sigma
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031030000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 3 sigma
  } else if( trainConfig == 6068) {
    // MB 13TeV EMCal + DCal Background Variation (Swapping around Pi0, Method by Joshua)
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints (omega)
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0y631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation (Pi0)
    cuts.AddCut("00010113","00200009327000008250400000","411791206f032230000","01631031000000d0","0z631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints (Pi0)
  } else if( trainConfig == 6071) {
    // EG2 13TeV EMcal
    cuts.AddCut("0008e113","00200009327000008250400000","1111100067032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 6072) {
    // EG2 std NL 12 EMCal + DCal
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 6073) {
    // EG2 13TeV EMcal + Dcal, Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 6074) {
    // EG2 13TeV EMcal + Dcal, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0r631031000000d0"); // background scheme swapping method trough roation
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0v631031000000d0"); // background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 6075) {
    // EG2 13TeV EMCal + DCal, Pi0 selection plus Gamma dropout, pi0 sigma range
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","0163103u000000d0","0v631031000000d0"); // 1 sigma Pi0, background scheme swapping method trough TGPS
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0, background scheme swapping method trough TGPS
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","0163103v000000d0","0v631031000000d0"); // 3 sigma Pi0, background scheme swapping method trough TGPS
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","0163103w000000d0","0v631031000000d0"); // 4 sigma Pi0, background scheme swapping method trough TGPS
  } else if( trainConfig == 6076) {
    // EG2 13TeV EMCal + DCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 6077) {
    // MB 13TeV EMCal + DCal Background Variation (Swapping Method by Joshua), AP-like cut
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, no AP like cut
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031010000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 1 sigma
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031020000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 2 sigma
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031030000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 3 sigma
  } else if( trainConfig == 6078) {
    // EG2 13TeV EMcal + Dcal, Background Variation (Swapping around Pi0, Method by Joshua)
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints (omega)
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0y631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation (Pi0)
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0z631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints (Pi0)
  } else if( trainConfig == 6079) {
    // EG2 13TeV EMcal + Dcal, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0r631031000000d0"); // background scheme swapping method trough roation
    cuts.AddCut("0008e113","00200009327000008250400000","411791206f032230000","01631031000000d0","0s631031000000d0"); // background scheme swapping method trough roation with weighting
  } else if( trainConfig == 6081) {
    // EG1 13TeV EMcal
    cuts.AddCut("0008d113","00200009327000008250400000","1111100067032230000","01631031000000d0","01631031000000d0");
  } else if( trainConfig == 6082) {
    // EG1 13TeV EMcal + Dcal
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 6083) {
    // EG1 13TeV EMcal + Dcal, Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","0163103u000000d0","01631031000000d0"); // 1 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","0163103v000000d0","01631031000000d0"); // 3 sigma Pi0 selection plus Gamma dropout
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","0163103w000000d0","01631031000000d0"); // 4 sigma Pi0 selection plus Gamma dropout
  } else if( trainConfig == 6084) {
    // EG1 13TeV EMcal + Dcal, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0, background scheme swapping method trough roation
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0, background scheme swapping method trough TGPS without constraints
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0, background scheme swapping method trough TGPS with constraints
  } else if( trainConfig == 6085) {
    // EG1 13TeV EMCal + DCal, Pi0 selection plus Gamma dropout, pi0 sigma range
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","0163103u000000d0","0v631031000000d0"); // 1 sigma Pi0, background scheme swapping method trough TGPS
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0v631031000000d0"); // 2 sigma Pi0, background scheme swapping method trough TGPS
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","0163103v000000d0","0v631031000000d0"); // 3 sigma Pi0, background scheme swapping method trough TGPS
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","0163103w000000d0","0v631031000000d0"); // 4 sigma Pi0, background scheme swapping method trough TGPS
  } else if( trainConfig == 6086) {
    // EG1 13TeV EMCal + DCal, Pi0 selection plus Gamma dropout, alpha cut omega
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631071000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.85 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631051000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.75 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631041000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.65 alpha cut for omega
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","01631081000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, 0.0-0.60 alpha cut for omega
  } else if( trainConfig == 6087) {
    // EG1 13TeV EMCal + DCal Background Variation (Swapping Method by Joshua), AP-like cut
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, no AP like cut
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031010000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 1 sigma
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031020000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 2 sigma
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031030000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints, AP like cut 3 sigma
  } else if( trainConfig == 6088) {
    // EG2 13TeV EMcal + Dcal, Background Variation (Swapping around Pi0, Method by Joshua)
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0x631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints (omega)
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0y631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough roation (Pi0)
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0z631031000000d0"); // 2 sigma Pi0 selection plus Gamma dropout, background scheme swapping method trough TGPS with constraints (Pi0)
  } else if( trainConfig == 6089) {
    // EG1 13TeV EMcal + Dcal, Background Variation (Swapping Method by Joshua)
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0r631031000000d0"); // 2 sigma Pi0, background scheme swapping method trough roation
    cuts.AddCut("0008d113","00200009327000008250400000","411791206f032230000","01631031000000d0","0s631031000000d0"); // 2 sigma Pi0, background scheme swapping method trough roation with weighting
  } else {
    Error(Form("OmegaToPiZeroGamma_%i_%i", trainConfig, isMC), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerPiZeroGamma! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();
  TList *ClusterCutList = new TList();
  TList *neutralPionCutList = new TList();
  TList *MesonCutList = new TList();

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
  neutralPionCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisNeutralPionCuts = new AliConversionMesonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    //create AliCaloTrackMatcher instance, if there is none present
    TString caloCutPos = cuts.GetClusterCut(i);
    caloCutPos.Resize(1);
    TString TrackMatcherName = Form("CaloTrackMatcher_%s_%i",caloCutPos.Data(),trackMatcherRunningMode);
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
      AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi(),trackMatcherRunningMode);
      fTrackMatcher->SetV0ReaderName(V0ReaderName);
      mgr->AddTask(fTrackMatcher);
      mgr->ConnectInput(fTrackMatcher,0,cinput);
    }

    analysisEventCuts[i] = new AliConvEventCuts();

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    TString triggerString   = cuts.GetEventCut(i);
    triggerString           = triggerString(3,2);
    if (triggerString.CompareTo("03")==0)
      triggerString         = "00";
    if (periodNameAnchor.CompareTo("LHC13g") == 0 && triggerString.CompareTo("10")== 0 )
      triggerString         = "00";

    dataInputMultHisto      = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto        = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());

    if (doMultiplicityWeighting){
      cout << "enableling mult weighting" << endl;
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameInputForMultWeighing, dataInputMultHisto, mcInputMultHisto );
    }


    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    if(fMinPtHardSet)
      analysisEventCuts[i]->SetMinFacPtHard(minFacPtHard);
    if(fMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    if(fSingleMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHardSingleParticle(maxFacPtHardSingle);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetLightOutput(runLightOutput);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i] = new AliCaloPhotonCuts((isMC==2));
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    analysisClusterCuts[i]->SetLightOutput(runLightOutput);
    if(ReconMethod == 2) analysisClusterCuts[i]->SetIsPureCaloCut(2);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    // analysisClusterCuts[i]->SetFillCutHistograms("");

    analysisNeutralPionCuts[i] = new AliConversionMesonCuts();
    analysisNeutralPionCuts[i]->SetLightOutput(runLightOutput);
    analysisNeutralPionCuts[i]->SetRunningMode(2);
    if(usePtDepSelectionWindowCut) analysisNeutralPionCuts[i]->SetUsePtDepSelectionWindow(usePtDepSelectionWindowCut);
    if(doSmear) analysisNeutralPionCuts[i]->SetDefaultSmearing(bremSmear,smearPar,smearParConst);
    analysisNeutralPionCuts[i]->InitializeCutsFromCutString((cuts.GetNeutralPionCut(i)).Data());
    neutralPionCutList->Add(analysisNeutralPionCuts[i]);
    analysisNeutralPionCuts[i]->SetFillCutHistograms("NeutralPionCuts");

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(runLightOutput);
    analysisMesonCuts[i]->SetRunningMode(2);
    analysisMesonCuts[i]->SetEnableOmegaAPlikeCut(kTRUE);                       // enables the overloaded ToCloseV0sCut to work as an AP like cut
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("MesonCuts");
    if(doSmear) analysisMesonCuts[i]->SetDefaultSmearing(bremSmear,smearPar,smearParConst);
    if(analysisMesonCuts[i]->DoGammaSwappForBg()) analysisClusterCuts[i]->SetUseEtaPhiMapForBackCand(kTRUE);
    analysisClusterCuts[i]->SetFillCutHistograms("");
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetNeutralPionCutList(numberOfCuts,neutralPionCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetCorrectionTaskSetting(corrTaskSetting);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetLightOutput(runLightOutput);
  task->SetDalitzCut(useDalitzCut);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);

  // task->SetDoMesonAnalysis(kTRUE); // ADDED

  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}

  //connect containers
  AliAnalysisDataContainer *coutput =
    // mgr->CreateContainer(Form("OmegaToPiZeroGamma_%i_%i", trainConfig, isMC), TList::Class(),
    //           AliAnalysisManager::kOutputContainer,Form("OmegaToPiZeroGamma_%i_%i.root", trainConfig, isMC));
    mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("OmegaToPiZeroGamma_%i",trainConfig) : Form("OmegaToPiZeroGamma_%i_%s",trainConfig,corrTaskSetting.Data()), TList::Class(),
              AliAnalysisManager::kOutputContainer, Form("OmegaToPiZeroGamma_%i.root",trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
