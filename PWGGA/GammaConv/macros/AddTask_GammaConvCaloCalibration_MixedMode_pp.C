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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaConvCaloCalibration.cxx) for
//for pp together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvCaloCalibration_MixedMode_pp(
  Int_t     selectedMeson                 = 0,        // select the corresponding meson: 0 pi0, 1 eta, 2 eta'
  Int_t     trainConfig                   = 1,        // change different set of cuts
  Int_t     isMC                          = 0,        // run MC
  TString   photonCutNumberV0Reader       = "",       // 00000008400000000100000000 nom. B, 00000088400000000100000000 low B
  TString   periodNameV0Reader            = "",
  // general setting for task
  Int_t     enableQAMesonTask             = 0,        // enable QA in AliAnalysisTaskGammaConvV1
  Int_t     enableQAPhotonTask            = 0,        // enable additional QA task
  Int_t     enableExtMatchAndQA           = 0,        // disabled (0), extMatch (1), extQA_noCellQA (2), extMatch+extQA_noCellQA (3), extQA+cellQA (4), extMatch+extQA+cellQA (5)
  Int_t     enableLightOutput             = 0,        // switch to run light output (only essential histograms for afterburner)
  Bool_t    enableTHnSparse               = kFALSE,   // switch on THNsparse
  Int_t     enableTriggerMimicking        = 0,        // enable trigger mimicking
  Bool_t    enableTriggerOverlapRej       = kFALSE,   // enable trigger overlap rejection
  TString   settingMaxFacPtHard           = "3.",       // maximum factor between hardest jet and ptHard generated
  Int_t     debugLevel                    = 0,        // introducing debug levels for grid running
  // settings for weights
  // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, separate with ;
  TString   fileNameExternalInputs        = "",
  Int_t     doWeightingPart               = 0,        // enable Weighting
  TString   generatorName                 = "DPMJET", // generator Name
  Bool_t    enableMultiplicityWeighting   = kFALSE,   //
  TString   periodNameAnchor              = "",       //
  // special settings
  Bool_t    enableSortingMCLabels         = kTRUE,    // enable sorting for MC cluster labels
  Bool_t    enableTreeConvGammaShape      = kFALSE,   // enable additional tree for conversion properties for clusters
  Bool_t    doSmear                       = kFALSE,   // switches to run user defined smearing
  Double_t  bremSmear                     = 1.,
  Double_t  smearPar                      = 0.,       // conv photon smearing params
  Double_t  smearParConst                 = 0.,       // conv photon smearing params
  Bool_t    doPrimaryTrackMatching        = kTRUE,    // enable basic track matching for all primary tracks to cluster
  Int_t     isRun2                        = kTRUE,    // enables different number of SM
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
) {

  AliCutHandlerPCM cuts;

  TString fileNamePtWeights           = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights         = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");

  TString addTaskName                 = "AddTask_GammaConvCaloCalibration_MixedMode_pp";
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }
  TString corrTaskSetting             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF", "", addTaskName);
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;

  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  TString strTrackMatcherRunningMode  = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "TM", "", addTaskName);
  if(additionalTrainConfig.Contains("TM"))
    trackMatcherRunningMode = strTrackMatcherRunningMode.Atoi();

  Bool_t doTreeEOverP = kFALSE; // switch to produce EOverP tree
  TString strdoTreeEOverP             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "EPCLUSTree", "", addTaskName);
  if(strdoTreeEOverP.Atoi()==1)
    doTreeEOverP = kTRUE;

  Bool_t doTreeClusterShowerShape = kFALSE; // switch to produce EOverP tree
  TString strdoTreeClusterShowerShape             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "INVMASSCLUSTree", "", addTaskName);
  if(strdoTreeClusterShowerShape.Atoi()==1)
    doTreeClusterShowerShape = kTRUE;

  TH1S* histoAcc = 0x0;         // histo for modified acceptance
  TString strModifiedAcc              = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "MODIFYACC", "", addTaskName);
  if(strModifiedAcc.Contains("MODIFYACC")){
    TString tempType = strModifiedAcc;
    tempType.Replace(0,9,"");
    cout << "INFO: connecting to alien..." << endl;
    TGrid::Connect("alien://");
    cout << "done!" << endl;
    TFile *w = TFile::Open(fileNamePtWeights.Data());
    if(!w){cout << "ERROR: Could not open file: " << fileNamePtWeights.Data() << endl;return;}
    histoAcc = (TH1S*) w->Get(tempType.Data());
    if(!histoAcc) {cout << "ERROR: Could not find histo: " << tempType.Data() << endl;return;}
    cout << "found: " << histoAcc << endl;
  }

  Int_t localDebugFlag = 0;
  TString strLocalDebugFlag              = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "LOCALDEBUGFLAG", "", addTaskName);
  if(strLocalDebugFlag.Atoi()>0)
    localDebugFlag = strLocalDebugFlag.Atoi();

  TObjArray *rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaConvCaloCalibration_MixedMode_pp during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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
  // meson reco mode: 0 - PCM-PCM, 1 - PCM-Calo, 2 - Calo-Calo
  Int_t mesonRecoMode = 1;

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("%s_%i", addTaskName.Data(),  trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = photonCutNumberV0Reader.Data();
  TString cutnumberEvent = "00000003";
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
  //========= Add task to the ANALYSIS manager =====
  //================================================
  AliAnalysisTaskConvCaloCalibration *task=NULL;
  task= new AliAnalysisTaskConvCaloCalibration(Form("ConvCaloCalibration_%i_%i_%i", mesonRecoMode, selectedMeson, trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetCorrectionTaskSetting(corrTaskSetting);

  task->SetMesonRecoMode(mesonRecoMode); // meson reco mode: 0 - PCM-PCM, 1 - PCM-Calo, 2 - Calo-Calo
  if (enableLightOutput > 1) task->SetLightOutput(kTRUE);
  task->SetDoPrimaryTrackMatching(doPrimaryTrackMatching);

  // *****************************************************************************************************
  // 13 TeV  pp Run2 - EDC configurations
  // *****************************************************************************************************
  if (trainConfig == 1){ // EMCAL + DCal clusters 13 TeV
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117900007032220000","0163103100000010"); // no timing cut, no NL INT7
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117900067032220000","0163103100000010"); // -50ns, 30ns timing cut, no NL INT7
  } else if (trainConfig == 2){ // EMCAL + DCal clusters 13 TeV
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","41179000a7032230000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("0008e113","00200009327000008250400000","41179000a7032230000","0163103100000010"); // EG2+DG2
    cuts.AddCutPCMCalo("0008d113","00200009327000008250400000","41179000a7032230000","0163103100000010"); // EG1+DG1
    cuts.AddCutPCMCalo("0009b113","00200009327000008250400000","41179060a7032230000","0163103100000010"); // EGJ+DJ1
  } else if (trainConfig == 3){ // EMCAL + DCal clusters 5 TeV - iteration2 test
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117917067032220000","0163103100000010"); // -50ns, 30ns timing cut, no NL INT7
  } else if (trainConfig == 4){ // EMCAL + DCal clusters 13 TeV - iteration 2
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117900067032230000","0163103100000010"); // INT7 NoNL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","4117911067032230000","0163103100000010"); // INT7 NL11
  } else if (trainConfig == 5){ // EMCAL + DCal clusters 13 TeV - iteration 2
    cuts.AddCutPCMCalo("0008e113","00200009327000008250400000","4117900067032230000","0163103100000010"); // EG2  No NL
    cuts.AddCutPCMCalo("0008e113","00200009327000008250400000","4117911067032230000","0163103100000010"); // EG2  NL11
  } else if (trainConfig == 6){ // EMCAL + DCal clusters 13 TeV - iteration 2
    cuts.AddCutPCMCalo("0008d113","00200009327000008250400000","4117911067032230000","0163103100000010"); // EG1  No NL
    cuts.AddCutPCMCalo("0008d113","00200009327000008250400000","4117911067032230000","0163103100000010"); // EG1  NL11
  } else if (trainConfig == 7){ // EMCAL + DCal clusters 13 TeV - NCell study
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","411790106fe30220000","0163103100000010"); // INT7 NoNL + FT
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","411790106f030000000","0163103100000010"); // INT7 No Clus cuts TBNL + FT
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","411793706f030000000","0163103100000010"); // INT7 No Clus cuts TBNL
    cuts.AddCutPCMCalo("00010113","00200009327000008250400000","411790406f030000000","0163103100000010"); // INT7 No Clus cuts TBNL no Scale
  } else if (trainConfig == 8){ // EMCAL + DCal clusters 13 TeV - NCell study
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0163103100000010"); // no NCell cut  NL01
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe32220000","0163103100000010"); // NCell >= 2 cut NL01
  } else if (trainConfig == 9){ // EMCAL + DCal clusters 13 TeV - NCell study
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe30220000","0163103100000010"); // no NCell cut  NL21
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe32220000","0163103100000010"); // NCell >= 2 cut NL21

  // diff clusterizer
  } else if (trainConfig == 10){ // std
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe30220000","0163103100000010"); // no NCell cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe32220000","0163103100000010"); // NCell >= 2 cut
  } else if (trainConfig == 11){ // var1
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe30220000","0163103100000010"); // no NCell cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe32220000","0163103100000010"); // NCell >= 2 cut
  } else if (trainConfig == 12){ // var2
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe30220000","0163103100000010"); // no NCell cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe32220000","0163103100000010"); // NCell >= 2 cut
  } else if (trainConfig == 13){ // var3
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe30220000","0163103100000010"); // no NCell cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe32220000","0163103100000010"); // NCell >= 2 cut
  } else if (trainConfig == 14){ // var4
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe30220000","0163103100000010"); // no NCell cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109fe32220000","0163103100000010"); // NCell >= 2 cut

  // diff clusterizer without M02 and exotics cut
  } else if (trainConfig == 15){ // std
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109f020000000","0163103100000010"); // no NCell cut
  } else if (trainConfig == 16){ // var1
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109f020000000","0163103100000010"); // no NCell cut
  } else if (trainConfig == 17){ // var2
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109f020000000","0163103100000010"); // no NCell cut
  } else if (trainConfig == 18){ // var3
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109f020000000","0163103100000010"); // no NCell cut
  } else if (trainConfig == 19){ // var4
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411792109f020000000","0163103100000010"); // no NCell cut

  // different NL scales
  } else if (trainConfig == 20){ // no scale in NL
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799309f020000000","0163103100000010"); // no NCell cut
  } else if (trainConfig == 21){ // 1.5% scale in NL
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799409f020000000","0163103100000010"); // no NCell cut
  } else if (trainConfig == 22){ // 3.5% scale in NL
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799509f020000000","0163103100000010"); // no NCell cut

  // different clusterization settings and different NonLins (applied in CF)
  } else if (trainConfig == 25){ // S500A100,  min energy = 700MeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe32220000","0163103100000010"); // INT7
  } else if (trainConfig == 26){ // S500A100, min energy = 700MeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe32220000","0163103100000010"); // INT7
  } else if (trainConfig == 27){ // S500A100, min energy = 700MeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe32220000","0163103100000010"); // INT7
  } else if (trainConfig == 28){ // S100A50, min energy = 700MeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe32220000","0163103100000010"); // INT7
  } else if (trainConfig == 29){ // S100A50 min energy = 200MeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fei2220000","0163103100000010"); // INT7
  } else if (trainConfig == 30){ // S100A50 min energy = 300MeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009feg2220000","0163103100000010"); // INT7

  } else if (trainConfig == 31){ // S500A100,  min energy = 700MeV
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","411790009fe32220000","0163103100000010"); // EMC7
  } else if (trainConfig == 32){ // S500A100, min energy = 700MeV
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","411790009fe32220000","0163103100000010"); // EMC7
  } else if (trainConfig == 33){ // S500A100, min energy = 700MeV
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","411790009fe32220000","0163103100000010"); // EMC7
  } else if (trainConfig == 34){ // S100A50, min energy = 700MeV
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","411790009fe32220000","0163103100000010"); // EMC7
  } else if (trainConfig == 35){ // S100A50 min energy = 200MeV
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","411790009fei2220000","0163103100000010"); // EMC7
  } else if (trainConfig == 36){ // S100A50 min energy = 300MeV
    cuts.AddCutPCMCalo("00052113","0dm00009f9730000dge0404000","411790009feg2220000","0163103100000010"); // EMC7

    // NonLin settings with only test beam WO scale (no fine tuning)
  } else if (trainConfig == 40){ // NL 96, nominal Bfield setting
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799609fe32220000","0163103100000010"); // INT7
  } else if (trainConfig == 41){ // NL 96, nominal Bfield setting
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411799609fe32220000","0163103100000010"); // EG2
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411799609fe32220000","0163103100000010"); // EG1

    // NonLin settings with only test beam WO scale (PCMEDC fine tuning)
  } else if (trainConfig == 42){ // NL 97, nominal Bfield setting
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799709fe32220000","0163103100000010"); // INT7
  } else if (trainConfig == 43){ // NL 97, nominal Bfield setting
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411799709fe32220000","0163103100000010"); // EG2
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411799709fe32220000","0163103100000010"); // EG1

    // NonLin settings with only test beam WO scale (EDC fine tuning)
  } else if (trainConfig == 44){ // NL 98, nominal Bfield setting
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799809fe32220000","0r63103100000010"); // INT7
  } else if (trainConfig == 45){ // NL 98, nominal Bfield setting
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411799809fe32220000","0r63103100000010"); // EG2
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411799809fe32220000","0r63103100000010"); // EG1


  // Same as above (4x) settings but with rotation instead of mixing
    // NonLin settings with only test beam WO scale (no fine tuning)
  } else if (trainConfig == 50){ // NL 96, nominal Bfield setting
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799609fe32220000","0r63103100000010"); // INT7
  } else if (trainConfig == 51){ // NL 96, nominal Bfield setting
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411799609fe32220000","0r63103100000010"); // EG2
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411799609fe32220000","0r63103100000010"); // EG1

    // NonLin settings with only test beam WO scale (PCMEDC and EMC fine tuning)
  } else if (trainConfig == 52){ // NL 97, 98, nominal Bfield setting
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799709fe32220000","0r63103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799809fe32220000","0r63103100000010"); // INT7
  } else if (trainConfig == 53){ // NL 97, 98, nominal Bfield setting
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411799709fe32220000","0r63103100000010"); // EG2
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411799809fe32220000","0r63103100000010"); // EG2
  } else if (trainConfig == 54){ // NL 97, 98, nominal Bfield setting
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411799709fe32220000","0r63103100000010"); // EG1
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411799809fe32220000","0r63103100000010"); // EG1

    // NonLin settings with only test beam WO scale (interpolation between PCM-EMC and EDC fine tuning)
  } else if (trainConfig == 55){ // NL 99, nominal Bfield setting
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909fe32220000","0r63103100000010"); // INT7
  } else if (trainConfig == 56){ // NL 99, nominal Bfield setting
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411799909fe32220000","0r63103100000010"); // EG2
  } else if (trainConfig == 57){ // NL 99, nominal Bfield setting
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411799909fe32220000","0r63103100000010"); // EG1


  } else if (trainConfig == 60){ // NL 96, low Bfield setting
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","411799609fe32220000","0163103100000010"); // INT7
  } else if (trainConfig == 61){ // NL 97, 98, low Bfield setting
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","411799709fe32220000","0163103100000010"); // INT7
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","411799809fe32220000","0163103100000010"); // INT7
  } else if (trainConfig == 62){ // NL 99, low Bfield setting
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","411799909fe32220000","0163103100000010"); // INT7


    // NonLin settings with only test beam WO scale (interpolation between PCM-EMC and EDC fine tuning)
  } else if (trainConfig == 63){ // NL 99, nominal Bfield setting, without NCell cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909fe30220000","0r63103100000010"); // INT7
  } else if (trainConfig == 64){ // NL 99, nominal Bfield settingwithout NCell, M02 and exotics cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909f030000000","0r63103100000010"); // INT7


  } else if (trainConfig == 70){ // NL applied in CF, nominal Bfield setting
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe32220000","0r63103100000010"); // INT7, NonLin applied in CF
  } else if (trainConfig == 71){ // NL applied in CF, nominal Bfield setting
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790009fe32220000","0r63103100000010"); // EG2, NonLin applied in CF
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790009fe32220000","0r63103100000010"); // EG1, NonLin applied in CF
  } else if (trainConfig == 72){ // NL applied in CF, nominal Bfield setting
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe32220000","0r63103100000010"); // INT7, NonLin applied in CF
  } else if (trainConfig == 73){ // NL applied in CF, nominal Bfield setting
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790009fe32220000","0r63103100000010"); // EG2, NonLin applied in CF
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790009fe32220000","0r63103100000010"); // EG1, NonLin applied in CF
  } else if (trainConfig == 74){ // NL applied in CF, nominal Bfield setting
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe32220000","0r63103100000010"); // INT7, NonLin applied in CF
  } else if (trainConfig == 75){ // NL applied in CF, nominal Bfield setting
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790009fe32220000","0r63103100000010"); // EG2, NonLin applied in CF
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790009fe32220000","0r63103100000010"); // EG1, NonLin applied in CF
  } else if (trainConfig == 76){ // NL applied in CF, nominal Bfield, with lower min energy (400MeV)
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009feh2220000","0r63103100000010"); // INT7, NonLin applied in CF
  } else if (trainConfig == 77){ // NL applied in CF, nominal Bfield, with lower min energy (400MeV)
    cuts.AddCutPCMCalo("0008e113","0dm00009f9730000dge0404000","411790009feh2220000","0r63103100000010"); // EG2, NonLin applied in CF
    cuts.AddCutPCMCalo("0008d113","0dm00009f9730000dge0404000","411790009feh2220000","0r63103100000010"); // EG1, NonLin applied in CF

  } else if (trainConfig == 100){ //
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909fe30220000","0r63103100000010"); // INT7 with M02/exotic cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909f030000000","0r63103100000010"); // INT7 without M02/exotic cut
  } else if (trainConfig == 101){ //
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799609fe30220000","0r63103100000010"); // INT7 with M02/exotic cut, no FT
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799609f030000000","0r63103100000010"); // INT7 without M02/exotic cut, no FT

  } else if (trainConfig == 102){ // with isolated cluster cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909fe302200i0","0r63103100000010"); // INT7 with M02/exotic cut + isolated cluster
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909fe302200j0","0r63103100000010"); // INT7 with M02/exotic cut + isolated cluster
  } else if (trainConfig == 102){ // with isolated cluster cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909f0300000i0","0r63103100000010"); // INT7 without M02/exotic cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909f0300000j0","0r63103100000010"); // INT7 without M02/exotic cut

    // lowB settings
  } else if (trainConfig == 110){ // lowB
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109fe30220000","0r63103100000010"); // INT7 with M02/exotic cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790109f030000000","0r63103100000010"); // INT7 without M02/exotic cut
  } else if (trainConfig == 111){ // lowB
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","411799609fe30220000","0r63103100000010"); // INT7 with M02/exotic cut, no FT
    cuts.AddCutPCMCalo("00010113","0dm00089f9730000iih0404000","411799609f030000000","0r63103100000010"); // INT7 without M02/exotic cut, no FT


  } else if (trainConfig == 120){ // minE = 700 MeV, noNL
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe30220000","0163103100000010"); // no NCell cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe32220000","0163103100000010"); // NCell >= 2 cut
  } else if (trainConfig == 121){ // minE = 500 MeV, noNL
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe20220000","0163103100000010"); // no NCell cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe22220000","0163103100000010"); // NCell >= 2 cut
  } else if (trainConfig == 122){ // minE = 400 MeV, noNL
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe80220000","0163103100000010"); // no NCell cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe82220000","0163103100000010"); // NCell >= 2 cut
  } else if (trainConfig == 123){ // minE = 100 MeV, noNL
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe00220000","0163103100000010"); // no NCell cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fe02220000","0163103100000010"); // NCell >= 2 cut
  } else if (trainConfig == 124){ // minE = 200 MeV, noNL
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fei0220000","0163103100000010"); // no NCell cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411790009fei2220000","0163103100000010"); // NCell >= 2 cut

  // Cluster efficiency with lower cluster Emin
  } else if (trainConfig == 125){ // NL 99, nominal Bfield setting, without NCell cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909feg09v0000","0r63103100000010"); // INT7, Emin = 300 MeV
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909feg0200000","0r63103100000010"); // INT7, Emin = 300 MeV, open M02
  // same as 125 but with exotics > 6
  } else if (trainConfig == 126){ // NL 99, nominal Bfield setting, without NCell cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909flg09v0000","0r63103100000010"); // INT7, Emin = 300 MeV, exotic_E > 6
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909flg0200000","0r63103100000010"); // INT7, Emin = 300 MeV, open M02, exotic_E > 6

  } else if (trainConfig == 130){ // NL 99, nominal Bfield setting, without NCell cut + isolated cluster cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909feg09v00i0","0r63103100000010"); // INT7, Emin = 300 MeV + isolated clusters
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909feg09v00i0","0r63103100000010"); // INT7, Emin = 300 MeV + isolated clusters
  } else if (trainConfig == 131){ // NL 99, nominal Bfield setting, without NCell cut + isolated cluster cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909feg02000i0","0r63103100000010"); // INT7, Emin = 300 MeV, open M02 + isolated clusters
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909feg02000j0","0r63103100000010"); // INT7, Emin = 300 MeV, open M02 + isolated clusters
  } else if (trainConfig == 132){ // NL 99, nominal Bfield setting, without NCell cut + isolated cluster cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909flg09v00i0","0r63103100000010"); // INT7, Emin = 300 MeV + isolated clusters, exotic_E > 6
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909flg09v00i0","0r63103100000010"); // INT7, Emin = 300 MeV + isolated clusters, exotic_E > 6
  } else if (trainConfig == 133){ // NL 99, nominal Bfield setting, without NCell cut + isolated cluster cut
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909flg02000i0","0r63103100000010"); // INT7, Emin = 300 MeV, open M02 + isolated clusters, exotic_E > 6
    cuts.AddCutPCMCalo("00010113","0dm00009f9730000dge0404000","411799909flg02000j0","0r63103100000010"); // INT7, Emin = 300 MeV, open M02 + isolated clusters, exotic_E > 6

  } else {
    Error(Form("AddTask_GammaConvCaloCalibration_MixedMode_pp%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
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
    if(corrTaskSetting.CompareTo("")){
      TrackMatcherName = TrackMatcherName+"_"+corrTaskSetting.Data();
      cout << "Using separate track matcher for correction framework setting: " << TrackMatcherName.Data() << endl;
    }
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

    if (enableMultiplicityWeighting){
      cout << "enabling mult weighting" << endl;
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameMultWeights, dataInputMultHisto, mcInputMultHisto );
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
    analysisEventCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    if (enableLightOutput > 0) analysisEventCuts[i]->SetLightOutput(kTRUE);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    EventCutList->Add(analysisEventCuts[i]);


    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    if (enableLightOutput > 0) analysisCuts[i]->SetLightOutput(kTRUE);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisClusterCuts[i] = new AliCaloPhotonCuts(isMC);
    analysisClusterCuts[i]->SetHistoToModifyAcceptance(histoAcc);
    analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisClusterCuts[i]->SetCorrectionTaskSetting(corrTaskSetting);
    analysisClusterCuts[i]->SetCaloTrackMatcherName(TrackMatcherName);
    if (enableLightOutput > 0) analysisClusterCuts[i]->SetLightOutput(kTRUE);
    analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedMatchAndQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    if (enableLightOutput > 0) analysisMesonCuts[i]->SetLightOutput(kTRUE);
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
  task->SetUseTHnSparse(enableTHnSparse);
  task->SetEnableSortingOfMCClusLabels(enableSortingMCLabels);
  if(enableExtMatchAndQA > 1){ task->SetPlotHistsExtQA(kTRUE);}
  if(isRun2){ task->SetNumOfCaloModules(20); }
  if(!isRun2){ task->SetNumOfCaloModules(10); }

  //connect containers
  AliAnalysisDataContainer *coutput =
  mgr->CreateContainer(!(corrTaskSetting.CompareTo("")) ? Form("ConvCaloCalibration_%i_%i_%i", mesonRecoMode, selectedMeson, trainConfig)
                                                        :  Form("ConvCaloCalibration_%i_%i_%i_%s", mesonRecoMode, selectedMeson, trainConfig, corrTaskSetting.Data()), TList::Class(),
                        AliAnalysisManager::kOutputContainer,Form("ConvCaloCalibration_%i_%i.root",mesonRecoMode,trainConfig) );

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
