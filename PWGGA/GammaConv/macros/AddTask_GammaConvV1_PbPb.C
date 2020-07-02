/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Lucia Leardini                                *
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
//($ALIPHYSICS/PWGGA/GammaConv/AliAnalysisTaskGammaConvV1.cxx) for
//Pb-Pb together with all supporting classes
//***************************************************************************************

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvV1_PbPb(
  Int_t     trainConfig                   = 1,        // change different set of cuts
  Int_t     isMC                          = 0,        // run MC
  TString   photonCutNumberV0Reader       = "",
  TString   periodNameV0Reader            = "",
  // general setting for task
  Int_t     enableQAMesonTask             = 0,        // enable QA in AliAnalysisTaskGammaConvV1
  Int_t     enableQAPhotonTask            = 0,        // enable additional QA task
  Bool_t    enableLightOutput             = kFALSE,   // switch to run light output (only essential histograms for afterburner)
  Bool_t    enableTHnSparse               = kFALSE,   // switch on THNsparse
  Int_t     enableTriggerMimicking        = 0,        // enable trigger mimicking
  Bool_t    enableTriggerOverlapRej       = kFALSE,   // enable trigger overlap rejection
  TString   settingMaxFacPtHard           = "3.",     // maximum factor between hardest jet and ptHard generated
  Int_t     debugLevel                    = 0,        // introducing debug levels for grid running
  // settings for weights
  // FPTW:fileNamePtWeights, FMUW:fileNameMultWeights, FMAW:fileNameMatBudWeights, FEPC:fileNamedEdxPostCalib, FCEF:fileNameCentFlattening, separate with ;
  TString   fileNameExternalInputs        = "",
  Int_t     doWeightingPart               = 0,        // enable Weighting
  Bool_t    enablePtWeighting             = kFALSE,   // enable Weighting
  TString   generatorName                 = "DPMJET", // generator Name
  Bool_t    enableMultiplicityWeighting   = kFALSE,   //
  TString   periodNameAnchor              = "",       //
  Int_t     enableMatBudWeightsPi0        = 0,        // 1 = three radial bins, 2 = 10 radial bins
  Bool_t    enableElecDeDxPostCalibration = kFALSE,
  Int_t     enableFlattening              = 0,
  // special settings
  Bool_t    enableChargedPrimary          = kFALSE,
  Bool_t    enablePlotVsCentrality        = kFALSE,
  Bool_t    processAODcheckForV0s         = kFALSE,   // flag for AOD check if V0s contained in AliAODs.root and AliAODGammaConversion.root
  // subwagon config
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig + special settings
)  {


  AliCutHandlerPCM cuts;

  Int_t isHeavyIon                    = 1;
  TString fileNamePtWeights           = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights         = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameMatBudWeights       = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMAW:");
  TString fileNamedEdxPostCalib       = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FEPC:");
  TString fileNameCentFlattening      = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FCEF:");
  TString fileNameBDT                 = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FBDT:");
  Bool_t enableBDT = kFALSE;
  if(fileNameBDT.CompareTo("")!= 0){
    enableBDT = kTRUE;
    cout << "enabling BDT !!!! " << fileNameBDT.Data() << endl;
  }

  TString addTaskName                 = "AddTask_GammaConvV1_PbPb";
  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: " << addTaskName.Data() << " running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("%s_%i", addTaskName.Data(),  trainConfig), "No analysis manager found.");
    return ;
  }

  TObjArray *rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaConvV1_PbPb during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
  Bool_t fMinPtHardSet        = kFALSE;
  Double_t minFacPtHard       = -1;
  Bool_t fMaxPtHardSet        = kFALSE;
  Double_t maxFacPtHard       = 100;
  Bool_t fSingleMaxPtHardSet  = kFALSE;
  Double_t maxFacPtHardSingle = 100;
  Bool_t fJetFinderUsage      = kFALSE;
  Bool_t fUsePtHardFromFile      = kFALSE;
  Bool_t fUseAddOutlierRej      = kFALSE;
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
    } else if(strTempSetting.BeginsWith("USEJETFINDER:")){
      strTempSetting.Replace(0,13,"");
      if(strTempSetting.Atoi()==1){
        cout << "using MC jet finder for outlier removal" << endl;
        fJetFinderUsage        = kTRUE;
      }
    } else if(strTempSetting.BeginsWith("PTHFROMFILE:")){
      strTempSetting.Replace(0,12,"");
      if(strTempSetting.Atoi()==1){
        cout << "using MC jet finder for outlier removal" << endl;
        fUsePtHardFromFile        = kTRUE;
      }
    } else if(strTempSetting.BeginsWith("ADDOUTLIERREJ:")){
      strTempSetting.Replace(0,14,"");
      if(strTempSetting.Atoi()==1){
        cout << "using path based outlier removal" << endl;
        fUseAddOutlierRej        = kTRUE;
      }
    } else if(rmaxFacPtHardSetting->GetEntries()==1 && strTempSetting.Atof()>0){
      maxFacPtHard               = strTempSetting.Atof();
      cout << "running with max pT hard jet fraction of: " << maxFacPtHard << endl;
      fMaxPtHardSet        = kTRUE;
    }
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = photonCutNumberV0Reader.Data();
  TString cutnumberEvent = "10000003";


  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  TString V0ReaderName        = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  AliV0ReaderV1 *fV0ReaderV1  =  NULL;
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    cout << "ERROR: V0Reader: " << V0ReaderName.Data() << " not found!!"<< endl;
    return;
  } else {
    cout << "V0Reader: " << V0ReaderName.Data() << " found!!"<< endl;
  }

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  AliAnalysisTaskGammaConvV1 *task=NULL;
  task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(enableLightOutput);
  if(enableBDT) task->SetFileNameBDT(fileNameBDT.Data());

  //****************************************************************************************************
  // 2.76TeV Pb-Pb LHC10h
  //****************************************************************************************************
  if (trainConfig == 1){ // Standard cuts
    cuts.AddCutPCM("60100013", "04200009297002003220000000", "0152204500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "04200009297002003220000000", "0152204500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "04200009297002003220000000", "0152204500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "04200009297002003220000000", "0152204500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "04200009297002003220000000", "0152204500900000"); // 0-20%
  } else if (trainConfig == 2) { // Standard cuts
    cuts.AddCutPCM("52400013", "04200009297002003220000000", "0152204500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "04200009297002003220000000", "0152206500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "04200009297002003220000000", "0152206500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "04200009297002003220000000", "0152206500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "04200009297002003220000000", "0152206500900000"); // 40-90%
  } else if (trainConfig == 3) { // Standard cuts only added signals
    cuts.AddCutPCM("60100023", "04200009297002003220000000", "0152204500900000"); // 0-5%
    cuts.AddCutPCM("61200023", "04200009297002003220000000", "0152204500900000"); // 5-10%
    cuts.AddCutPCM("50100023", "04200009297002003220000000", "0152204500900000"); // 0-10%
    cuts.AddCutPCM("51200023", "04200009297002003220000000", "0152204500900000"); // 10-20%
    cuts.AddCutPCM("50200023", "04200009297002003220000000", "0152204500900000"); // 0-20%
  } else if (trainConfig == 4) { // Standard cuts only added signals
    cuts.AddCutPCM("52400023", "04200009297002003220000000", "0152204500900000"); // 20-40%
    cuts.AddCutPCM("54600023", "04200009297002003220000000", "0152206500900000"); // 40-60%
    cuts.AddCutPCM("56800023", "04200009297002003220000000", "0152206500900000"); // 60-80%
    cuts.AddCutPCM("54800023", "04200009297002003220000000", "0152206500900000"); // 20-40%
    cuts.AddCutPCM("54900023", "04200009297002003220000000", "0152206500900000"); // 40-90%
  } else if (trainConfig == 5){ // R-minCut 7.5 cm
    cuts.AddCutPCM("60100013", "04900009297002003220000000", "0152204500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "04900009297002003220000000", "0152204500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "04900009297002003220000000", "0152204500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "04900009297002003220000000", "0152204500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "04900009297002003220000000", "0152204500900000"); // 0-20%
  } else if (trainConfig == 6) { // R-minCut 7.5 cm
    cuts.AddCutPCM("52400013", "04900009297002003220000000", "0152204500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "04900009297002003220000000", "0152206500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "04900009297002003220000000", "0152206500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "04900009297002003220000000", "0152206500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "04900009297002003220000000", "0152206500900000"); // 40-90%
  } else if (trainConfig == 7) {// R-minCut 7.5 cm
    cuts.AddCutPCM("60100023", "04900009297002003220000000", "0152204500900000"); // 0-5%
    cuts.AddCutPCM("61200023", "04900009297002003220000000", "0152204500900000"); // 5-10%
    cuts.AddCutPCM("50100023", "04900009297002003220000000", "0152204500900000"); // 0-10%
    cuts.AddCutPCM("51200023", "04900009297002003220000000", "0152204500900000"); // 10-20%
    cuts.AddCutPCM("50200023", "04900009297002003220000000", "0152204500900000"); // 0-20%
  } else if (trainConfig == 8) { // R-minCut 7.5 cm
    cuts.AddCutPCM("52400023", "04900009297002003220000000", "0152204500900000"); // 20-40%
    cuts.AddCutPCM("54600023", "04900009297002003220000000", "0152206500900000"); // 40-60%
    cuts.AddCutPCM("56800023", "04900009297002003220000000", "0152206500900000"); // 60-80%
    cuts.AddCutPCM("54800023", "04900009297002003220000000", "0152206500900000"); // 20-40%
    cuts.AddCutPCM("54900023", "04900009297002003220000000", "0152206500900000"); // 40-90%
  } else if (trainConfig == 9){ // R-minCut 12.5 cm
    cuts.AddCutPCM("60100013", "04800009297002003220000000", "0152204500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "04800009297002003220000000", "0152204500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "04800009297002003220000000", "0152204500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "04800009297002003220000000", "0152204500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "04800009297002003220000000", "0152204500900000"); // 0-20%
  } else if (trainConfig == 10) { // R-minCut 12.5 cm
    cuts.AddCutPCM("52400013", "04800009297002003220000000", "0152204500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "04800009297002003220000000", "0152206500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "04800009297002003220000000", "0152206500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "04800009297002003220000000", "0152206500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "04800009297002003220000000", "0152206500900000"); // 40-90%
  } else if (trainConfig == 11) {// R-minCut 12.5 cm
    cuts.AddCutPCM("60100023", "04800009297002003220000000", "0152204500900000"); // 0-5%
    cuts.AddCutPCM("61200023", "04800009297002003220000000", "0152204500900000"); // 5-10%
    cuts.AddCutPCM("50100023", "04800009297002003220000000", "0152204500900000"); // 0-10%
    cuts.AddCutPCM("51200023", "04800009297002003220000000", "0152204500900000"); // 10-20%
    cuts.AddCutPCM("50200023", "04800009297002003220000000", "0152204500900000"); // 0-20%
  } else if (trainConfig == 12) { // R-minCut 12.5 cm
    cuts.AddCutPCM("52400023", "04800009297002003220000000", "0152204500900000"); // 20-40%
    cuts.AddCutPCM("54600023", "04800009297002003220000000", "0152206500900000"); // 40-60%
    cuts.AddCutPCM("56800023", "04800009297002003220000000", "0152206500900000"); // 60-80%
    cuts.AddCutPCM("54800023", "04800009297002003220000000", "0152206500900000"); // 20-40%
    cuts.AddCutPCM("54900023", "04800009297002003220000000", "0152206500900000"); // 40-90%
  }else  if (trainConfig == 13){ // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCutPCM("60100013", "03200009297002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "03200009297002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "03200009297002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "03200009297002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "03200009297002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 14) {  // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCutPCM("52400013", "03200009297002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "03200009297002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "03200009297002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "03200009297002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "03200009297002003220000000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 15) { // LHC10h standard, eta 0.65, y = 0.6  added signals
    cuts.AddCutPCM("60100023", "03200009297002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200023", "03200009297002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100023", "03200009297002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200023", "03200009297002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200023", "03200009297002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 16) { // LHC10h standard, eta 0.65, y = 0.6  added signals
    cuts.AddCutPCM("52400023", "03200009297002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600023", "03200009297002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800023", "03200009297002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800023", "03200009297002003220000000", "0152306500900000"); // 20-40%
    cuts.AddCutPCM("54900023", "03200009297002003220000000", "0152306500900000"); // 40-90%
  }else  if (trainConfig == 17){ // LHC10h standard, eta 0.65, y = 0.6, photon quality 1
    cuts.AddCutPCM("60100013", "03200009297002003220020000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "03200009297002003220020000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "03200009297002003220020000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "03200009297002003220020000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "03200009297002003220020000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 18) {  // LHC10h standard, eta 0.65, y = 0.6, photon quality 1
    cuts.AddCutPCM("52400013", "03200009297002003220020000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "03200009297002003220020000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "03200009297002003220020000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "03200009297002003220020000", "0152306500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "03200009297002003220020000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 19) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 1  added signal
    cuts.AddCutPCM("60100023", "03200009297002003220020000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200023", "03200009297002003220020000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100023", "03200009297002003220020000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200023", "03200009297002003220020000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200023", "03200009297002003220020000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 20) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 1 added signal
    cuts.AddCutPCM("52400023", "03200009297002003220020000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600023", "03200009297002003220020000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800023", "03200009297002003220020000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800023", "03200009297002003220020000", "0152306500900000"); // 20-40%
    cuts.AddCutPCM("54900023", "03200009297002003220020000", "0152306500900000"); // 40-90%
  }else  if (trainConfig == 21){ // LHC10h standard, eta 0.65, y = 0.6, photon quality 2
    cuts.AddCutPCM("60100013", "03200009297002003220030000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "03200009297002003220030000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "03200009297002003220030000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "03200009297002003220030000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "03200009297002003220030000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 22) {  // LHC10h standard, eta 0.65, y = 0.6, photon quality 2
    cuts.AddCutPCM("52400013", "03200009297002003220030000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "03200009297002003220030000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "03200009297002003220030000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "03200009297002003220030000", "0152306500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "03200009297002003220030000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 23) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 2  added signal
    cuts.AddCutPCM("60100023", "03200009297002003220030000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200023", "03200009297002003220030000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100023", "03200009297002003220030000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200023", "03200009297002003220030000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200023", "03200009297002003220030000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 24) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 2 added signal
    cuts.AddCutPCM("52400023", "03200009297002003220030000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600023", "03200009297002003220030000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800023", "03200009297002003220030000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800023", "03200009297002003220030000", "0152306500900000"); // 20-40%
    cuts.AddCutPCM("54900023", "03200009297002003220030000", "0152306500900000"); // 40-90%
  }else  if (trainConfig == 25){ // LHC10h standard, eta 0.65, y = 0.6, photon quality 3
    cuts.AddCutPCM("60100013", "03200009297002003220040000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "03200009297002003220040000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "03200009297002003220040000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "03200009297002003220040000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "03200009297002003220040000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 26) {  // LHC10h standard, eta 0.65, y = 0.6, photon quality 3
    cuts.AddCutPCM("52400013", "03200009297002003220040000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "03200009297002003220040000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "03200009297002003220040000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "03200009297002003220040000", "0152306500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "03200009297002003220040000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 27) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 3  added signal
    cuts.AddCutPCM("60100023", "03200009297002003220040000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200023", "03200009297002003220040000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100023", "03200009297002003220040000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200023", "03200009297002003220040000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200023", "03200009297002003220040000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 28) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 3 added signal
    cuts.AddCutPCM("52400023", "03200009297002003220040000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600023", "03200009297002003220040000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800023", "03200009297002003220040000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800023", "03200009297002003220040000", "0152306500900000"); // 20-40%
    cuts.AddCutPCM("54900023", "03200009297002003220040000", "0152306500900000"); // 40-90%
  }else  if (trainConfig == 29){ // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm
    cuts.AddCutPCM("60100013", "03700009297002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "03700009297002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "03700009297002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "03700009297002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "03700009297002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 30) {  // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm
    cuts.AddCutPCM("52400013", "03700009297002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "03700009297002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "03700009297002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "03700009297002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "03700009297002003220000000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 31) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm  added signals
    cuts.AddCutPCM("60100023", "03700009297002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200023", "03700009297002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100023", "03700009297002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200023", "03700009297002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200023", "03700009297002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 32) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm  added signals
    cuts.AddCutPCM("52400023", "03700009297002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600023", "03700009297002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800023", "03700009297002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800023", "03700009297002003220000000", "0152306500900000"); // 20-40%
    cuts.AddCutPCM("54900023", "03700009297002003220000000", "0152306500900000"); // 40-90%
  }else  if (trainConfig == 33){ // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1
    cuts.AddCutPCM("60100013", "03700009297002003220020000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "03700009297002003220020000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "03700009297002003220020000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "03700009297002003220020000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "03700009297002003220020000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 34) {  // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1
    cuts.AddCutPCM("52400013", "03700009297002003220020000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "03700009297002003220020000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "03700009297002003220020000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "03700009297002003220020000", "0152306500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "03700009297002003220020000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 35) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1  added signals
    cuts.AddCutPCM("60100023", "03700009297002003220020000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200023", "03700009297002003220020000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100023", "03700009297002003220020000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200023", "03700009297002003220020000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200023", "03700009297002003220020000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 36) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1  added signals
    cuts.AddCutPCM("52400023", "03700009297002003220020000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600023", "03700009297002003220020000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800023", "03700009297002003220020000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800023", "03700009297002003220020000", "0152306500900000"); // 20-40%
    cuts.AddCutPCM("54900023", "03700009297002003220020000", "0152306500900000"); // 40-90%
  }else  if (trainConfig == 37){ // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3
    cuts.AddCutPCM("60100013", "03700009297002003220040000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "03700009297002003220040000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "03700009297002003220040000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "03700009297002003220040000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "03700009297002003220040000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 38) {  // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3
    cuts.AddCutPCM("52400013", "03700009297002003220040000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "03700009297002003220040000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "03700009297002003220040000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "03700009297002003220040000", "0152306500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "03700009297002003220040000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 39) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3  added signals
    cuts.AddCutPCM("60100023", "03700009297002003220040000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200023", "03700009297002003220040000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100023", "03700009297002003220040000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200023", "03700009297002003220040000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200023", "03700009297002003220040000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 40) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3  added signals
    cuts.AddCutPCM("52400023", "03700009297002003220040000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600023", "03700009297002003220040000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800023", "03700009297002003220040000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800023", "03700009297002003220040000", "0152306500900000"); // 20-40%
    cuts.AddCutPCM("54900023", "03700009297002003220040000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 41){ // Standard cuts, eta 0.9, only to be run on data :kSemiCentral
    cuts.AddCutPCM("60160013", "00200009297002003220000000", "0152504500900000"); // 0-5%
    cuts.AddCutPCM("61260013", "00200009297002003220000000", "0152504500900000"); // 5-10%
    cuts.AddCutPCM("50160013", "00200009297002003220000000", "0152504500900000"); // 0-10%
    cuts.AddCutPCM("61260013", "00200009297002003220000000", "0152504500900000"); // 10-20%
    cuts.AddCutPCM("50260013", "00200009297002003220000000", "0152504500900000"); // 0-20%
  } else if (trainConfig == 42) { // Standard cuts, eta 0.9, only to be run on data :kSemiCentral
    cuts.AddCutPCM("52360013", "00200009297002003220000000", "0152506500900000");
    cuts.AddCutPCM("53460013", "00200009297002003220000000", "0152506500900000");
    cuts.AddCutPCM("54560013", "00200009297002003220000000", "0152506500900000");
    cuts.AddCutPCM("55660013", "00200009297002003220000000", "0152506500900000");
    cuts.AddCutPCM("56760013", "00200009297002003220000000", "0152506500900000");
  } else if (trainConfig == 43){ // Standard cuts, eta 0.9, only to be run on data
    cuts.AddCutPCM("52300013", "00200009297002003220000000", "0152504500900000");
    cuts.AddCutPCM("53400013", "00200009297002003220000000", "0152506500900000");
    cuts.AddCutPCM("54500013", "00200009297002003220000000", "0152506500900000");
    cuts.AddCutPCM("55600013", "00200009297002003220000000", "0152506500900000");
    cuts.AddCutPCM("52500013", "00200009297002003220000000", "0152506500900000");
  } else if ( trainConfig == 44){ // qt elipse cut 0.05
    cuts.AddCutPCM("60100013", "00200009297002008220000000", "0152506500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009297002008220000000", "0152506500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009297002008220000000", "0152506500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "00200009297002008220000000", "0152506500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "00200009297002008220000000", "0152506500900000"); // 0-20%
  } else if ( trainConfig == 45) { // qt elipse cut 0.05
    cuts.AddCutPCM("52400013", "00200009297002008220000000", "0152506500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "00200009297002008220000000", "0152506500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "00200009297002008220000000", "0152506500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "00200009297002008220000000", "0152506500900000"); // 40-80%
    cuts.AddCutPCM("53500013", "00200009297002008220000000", "0152506500900000"); // 30-50%
  } else if ( trainConfig == 46){ // qt elipse cut 0.05
    cuts.AddCutPCM("52300013", "00200009297002008220000000", "0152506500900000"); // 20-30%
    cuts.AddCutPCM("53400013", "00200009297002008220000000", "0152506500900000"); // 30-40%
    cuts.AddCutPCM("54500013", "00200009297002008220000000", "0152506500900000"); // 40-50%
    cuts.AddCutPCM("55600013", "00200009297002008220000000", "0152506500900000"); // 50-60%
    cuts.AddCutPCM("52500013", "00200009297002008220000000", "0152506500900000"); // 60-70%
  } else if ( trainConfig == 47){ // cos(theta_point) cut 0.85
    cuts.AddCutPCM("60100013", "00200009297002003220400000", "0152506500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009297002003220400000", "0152506500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009297002003220400000", "0152506500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "00200009297002003220400000", "0152506500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "00200009297002003220400000", "0152506500900000"); // 0-20%
  } else if ( trainConfig == 48) { // cos(theta_point) cut 0.85
    cuts.AddCutPCM("52400013", "00200009297002003220400000", "0152506500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "00200009297002003220400000", "0152506500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "00200009297002003220400000", "0152506500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "00200009297002003220400000", "0152506500900000"); // 40-80%
    cuts.AddCutPCM("53500013", "00200009297002003220400000", "0152506500900000"); // 30-50%
  } else if ( trainConfig == 49){ // cos(theta_point) cut 0.85
    cuts.AddCutPCM("52300013", "00200009297002003220400000", "0152506500900000"); // 20-30%
    cuts.AddCutPCM("53400013", "00200009297002003220400000", "0152506500900000"); // 30-40%
    cuts.AddCutPCM("54500013", "00200009297002003220400000", "0152506500900000"); // 40-50%
    cuts.AddCutPCM("55600013", "00200009297002003220400000", "0152506500900000"); // 50-60%
    cuts.AddCutPCM("52500013", "00200009297002003220400000", "0152506500900000"); // 60-70%
  } else if ( trainConfig == 50){ // psi pair 2D 0.05
    cuts.AddCutPCM("60100013", "00200009297002003260000000", "0152506500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009297002003260000000", "0152506500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009297002003260000000", "0152506500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "00200009297002003260000000", "0152506500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "00200009297002003260000000", "0152506500900000"); // 0-20%
  } else if ( trainConfig == 51) { // psi pair 2D 0.05
    cuts.AddCutPCM("52400013", "00200009297002003260000000", "0152506500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "00200009297002003260000000", "0152506500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "00200009297002003260000000", "0152506500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "00200009297002003260000000", "0152506500900000"); // 40-80%
    cuts.AddCutPCM("53500013", "00200009297002003260000000", "0152506500900000"); // 30-50%
  } else if ( trainConfig == 52){ // psi pair 2D 0.05
    cuts.AddCutPCM("52300013", "00200009297002003260000000", "0152506500900000"); // 20-30%
    cuts.AddCutPCM("53400013", "00200009297002003260000000", "0152506500900000"); // 30-40%
    cuts.AddCutPCM("54500013", "00200009297002003260000000", "0152506500900000"); // 40-50%
    cuts.AddCutPCM("55600013", "00200009297002003260000000", "0152506500900000"); // 50-60%
    cuts.AddCutPCM("52500013", "00200009297002003260000000", "0152506500900000"); // 60-70%
  } else if ( trainConfig == 53){ // psi pair 2D 0.1
    cuts.AddCutPCM("60100013", "00200009297002003250000000", "0152506500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009297002003250000000", "0152506500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009297002003250000000", "0152506500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "00200009297002003250000000", "0152506500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "00200009297002003250000000", "0152506500900000"); // 0-20%
  } else if ( trainConfig == 54) { // psi pair 2D 0.1
    cuts.AddCutPCM("52400013", "00200009297002003250000000", "0152506500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "00200009297002003250000000", "0152506500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "00200009297002003250000000", "0152506500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "00200009297002003250000000", "0152506500900000"); // 40-80%
    cuts.AddCutPCM("53500013", "00200009297002003250000000", "0152506500900000"); // 30-50%
  } else if ( trainConfig == 55){ // psi pair 2D 0.1
    cuts.AddCutPCM("52300013", "00200009297002003250000000", "0152506500900000"); // 20-30%
    cuts.AddCutPCM("53400013", "00200009297002003250000000", "0152506500900000"); // 30-40%
    cuts.AddCutPCM("54500013", "00200009297002003250000000", "0152506500900000"); // 40-50%
    cuts.AddCutPCM("55600013", "00200009297002003250000000", "0152506500900000"); // 50-60%
    cuts.AddCutPCM("52500013", "00200009297002003250000000", "0152506500900000"); // 60-70%
  } else if ( trainConfig == 56){ // cleaner cuts central classes
    cuts.AddCutPCM("60100013", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("51200013", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("50200013", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 57){ // cleaner cuts central classes - added signal
    cuts.AddCutPCM("60100023", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("51200023", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("50200023", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 58) { // cleaner cuts semicentral classes
    cuts.AddCutPCM("52400013", "00200009247602008250400000", "0152506500000000"); // 20-40%
    cuts.AddCutPCM("54600013", "00200009247602008250400000", "0152506500000000"); // 40-60%
    cuts.AddCutPCM("56800013", "00200009247602008250400000", "0152506500000000"); // 60-80%
    cuts.AddCutPCM("54800013", "00200009247602008250400000", "0152506500000000"); // 40-80%
    cuts.AddCutPCM("53500013", "00200009247602008250400000", "0152506500000000"); // 30-50%
  } else if ( trainConfig == 59) { // cleaner cuts semicentral classes - added signal
    cuts.AddCutPCM("52400023", "00200009247602008250400000", "0152506500000000"); // 20-40%
    cuts.AddCutPCM("54600023", "00200009247602008250400000", "0152506500000000"); // 40-60%
    cuts.AddCutPCM("56800023", "00200009247602008250400000", "0152506500000000"); // 60-80%
    cuts.AddCutPCM("54800023", "00200009247602008250400000", "0152506500000000"); // 40-80%
    cuts.AddCutPCM("53500023", "00200009247602008250400000", "0152506500000000"); // 30-50%
  } else if ( trainConfig == 60){ // cleaner cuts finer centralities
    cuts.AddCutPCM("52300013", "00200009247602008250400000", "0152506500000000"); // 20-30%
    cuts.AddCutPCM("53400013", "00200009247602008250400000", "0152506500000000"); // 30-40%
    cuts.AddCutPCM("54500013", "00200009247602008250400000", "0152506500000000"); // 40-50%
    cuts.AddCutPCM("55600013", "00200009247602008250400000", "0152506500000000"); // 50-60%
    cuts.AddCutPCM("52500013", "00200009247602008250400000", "0152506500000000"); // 60-70%
  } else if ( trainConfig == 61){ // cleaner cuts finer centralities - added signal
    cuts.AddCutPCM("52300023", "00200009247602008250400000", "0152506500000000"); // 20-30%
    cuts.AddCutPCM("53400023", "00200009247602008250400000", "0152506500000000"); // 30-40%
    cuts.AddCutPCM("54500023", "00200009247602008250400000", "0152506500000000"); // 40-50%
    cuts.AddCutPCM("55600023", "00200009247602008250400000", "0152506500000000"); // 50-60%
    cuts.AddCutPCM("52500023", "00200009247602008250400000", "0152506500000000"); // 60-70%
  } else if ( trainConfig == 62){ // cleaner cuts finer centralities
    cuts.AddCutPCM("62300013", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("63400013", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("64500013", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("65600013", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("66700013", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 63){ // cleaner cuts finer centralities - added signal
    cuts.AddCutPCM("62300023", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("63400023", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("64500023", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("65600023", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("66700023", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 64){ // cleaner cuts
    cuts.AddCutPCM("67800013", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("68900013", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("56700013", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("57800013", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("58900013", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 65){ // cleaner cuts added signal
    cuts.AddCutPCM("67800023", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("68900023", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("56700023", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("57800023", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("58900023", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 66){ // cleaner cuts
    // 0 -0%, 1-5%, 2-10%, 3-15%, 4-20%, 5-25%, 6-30%, 7-35%, 8-40%, 9-45%, a-50%, b-55%, c-60%, d-65%, e-70%, f-75%, g-80%, h-85%, i-90%, j-95%, k-100%
    cuts.AddCutPCM("69a00013", "00200009247602008250400000", "0152506500000000"); // 45-50%
    cuts.AddCutPCM("6ab00013", "00200009247602008250400000", "0152506500000000"); // 50-55%
    cuts.AddCutPCM("6bc00013", "00200009247602008250400000", "0152506500000000"); // 55-60%
    cuts.AddCutPCM("6cd00013", "00200009247602008250400000", "0152506500000000"); // 60-65%
    cuts.AddCutPCM("6de00013", "00200009247602008250400000", "0152506500000000"); // 65-70%
  } else if ( trainConfig == 67){ // cleaner cuts added signal
    cuts.AddCutPCM("69a00023", "00200009247602008250400000", "0152506500000000"); // 45-50%
    cuts.AddCutPCM("6ab00023", "00200009247602008250400000", "0152506500000000"); // 50-55%
    cuts.AddCutPCM("6bc00023", "00200009247602008250400000", "0152506500000000"); // 55-60%
    cuts.AddCutPCM("6cd00023", "00200009247602008250400000", "0152506500000000"); // 60-65%
    cuts.AddCutPCM("6de00023", "00200009247602008250400000", "0152506500000000"); // 65-70%
  } else if ( trainConfig == 68){ // cleaner cuts
    cuts.AddCutPCM("6ef00013", "00200009247602008250400000", "0152506500000000"); // 70-75%
    cuts.AddCutPCM("6fg00013", "00200009247602008250400000", "0152506500000000"); // 75-80%
    cuts.AddCutPCM("6gh00013", "00200009247602008250400000", "0152506500000000"); // 80-85%
    cuts.AddCutPCM("6hi00013", "00200009247602008250400000", "0152506500000000"); // 85-90%
    cuts.AddCutPCM("6ij00013", "00200009247602008250400000", "0152506500000000"); // 90-95%
  } else if ( trainConfig == 69){ // cleaner cuts added signal
    cuts.AddCutPCM("6ef00023", "00200009247602008250400000", "0152506500000000"); // 70-75%
    cuts.AddCutPCM("6fg00023", "00200009247602008250400000", "0152506500000000"); // 75-80%
    cuts.AddCutPCM("6gh00023", "00200009247602008250400000", "0152506500000000"); // 80-85%
    cuts.AddCutPCM("6hi00023", "00200009247602008250400000", "0152506500000000"); // 85-90%
    cuts.AddCutPCM("6ij00023", "00200009247602008250400000", "0152506500000000"); // 90-95%

  //****************************************************************************************************
  // 2.76TeV Pb-Pb LHC11h
  //****************************************************************************************************
    //----------- here start the syst var for LHC11h (and not only) -------------------
  } else if ( trainConfig == 70){ // min R = 35 cm
    cuts.AddCutPCM("60100013", "00700009247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00700009247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00700009247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00700009247602008850404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("52500013", "00700009247602008850404000", "0152501500000000"); // 0-20%
  } else if ( trainConfig == 71){ // min R = 35 cm added signal
    cuts.AddCutPCM("60100023", "00700009247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00700009247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00700009247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00700009247602008850404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("52500023", "00700009247602008850404000", "0152501500000000"); // 0-20%
  } else if ( trainConfig == 72){ // min R = 35 cm with phi cut
    cuts.AddCutPCM("60100013", "00716609247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00716609247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00716609247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00716609247602008850404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("52500013", "00716609247602008850404000", "0152501500000000"); // 0-20%
  } else if ( trainConfig == 73){ // min R = 35 cm with phi cut - added signal
    cuts.AddCutPCM("60100023", "00716609247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00716609247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00716609247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00716609247602008850404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("52500023", "00716609247602008850404000", "0152501500000000"); // 0-20%
  } else if ( trainConfig == 74){ // single pt 0.075 --------------------------------------
    cuts.AddCutPCM("60100013", "00200049247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200049247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200049247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200049247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200049247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 75){ // single pt 0.075 - added signal
    cuts.AddCutPCM("60100023", "00200049247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200049247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200049247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200049247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200049247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 76){ // single pt 0.075 with phi cut
    cuts.AddCutPCM("60100013", "00216649247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216649247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216649247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216649247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216649247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 77){ // single pt 0.075 with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216649247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216649247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216649247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216649247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216649247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 78){ // single pt 0.1 ----------------------------------------
    cuts.AddCutPCM("60100013", "00200019247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200019247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200019247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200019247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200019247602008850404000", "0152501500000000"); //20-50%
  } else if ( trainConfig == 79){ // single pt 0.1 - added signal
    cuts.AddCutPCM("60100023", "00200019247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200019247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200019247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200019247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200019247602008850404000", "0152501500000000"); //20-50%
  } else if ( trainConfig == 80){ // single pt 0.1 with phi cut
    cuts.AddCutPCM("60100013", "00216619247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216619247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216619247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216619247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216619247602008850404000", "0152501500000000"); //20-50%
  } else if ( trainConfig == 81){ // single pt 0.1  with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216619247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216619247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216619247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216619247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216619247602008850404000", "0152501500000000"); //20-50%
  } else if ( trainConfig == 82){ // variation TPC cls 0.7 --------------------------------
    cuts.AddCutPCM("60100013", "00200006247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200006247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200006247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200006247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200006247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 83){ // variation TPC cls 0.7 added signal
    cuts.AddCutPCM("60100023", "00200006247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200006247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200006247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200006247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200006247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 84){ // variation TPC cls 0.7 with phi cut
    cuts.AddCutPCM("60100013", "00216606247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216606247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216606247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216606247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216606247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 85){ // variation TPC cls 0.7 with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216606247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216606247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216606247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216606247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216606247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 86){ // variation TPC cls 0.35 -------------------------------
    cuts.AddCutPCM("60100013", "00200008247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200008247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200008247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200008247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200008247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 87){ // variation TPC cls 0.35 - added signal
    cuts.AddCutPCM("60100023", "00200008247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200008247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200008247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200008247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200008247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 88){ // variation TPC cls 0.35 with phi cut
    cuts.AddCutPCM("60100013", "00216608247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216608247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216608247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216608247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216608247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 89){ // variation TPC cls 0.35  with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216608247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216608247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216608247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216608247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216608247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 90){ // variation edEdx -4,5 ---------------------------------
    cuts.AddCutPCM("60100013", "00200009347602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009347602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009347602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009347602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009347602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 91){ // variation edEdx  -4,5 - added signal
    cuts.AddCutPCM("60100023", "00200009347602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009347602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009347602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009347602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009347602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 92){ // variation edEdx  -4,5 with phi cut
    cuts.AddCutPCM("60100013", "00216609347602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609347602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609347602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609347602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609347602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 93){ // variation edEdx  -4,5 with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609347602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609347602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609347602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609347602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609347602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 94){ // variation edEdx  -2.5,5 ------------------------------
    cuts.AddCutPCM("60100013", "00200009947602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009947602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009947602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009947602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009947602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 95){ // variation edEdx  -2.5,5 - added signal
    cuts.AddCutPCM("60100023", "00200009947602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009947602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009947602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009947602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009947602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 96){ // variation edEdx  -2.5,5 with phi cut
    cuts.AddCutPCM("60100013", "00216609947602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609947602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609947602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609947602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609947602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 97){ // variation edEdx  -2.5,5 with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609947602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609947602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609947602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609947602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609947602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 98){ // variation pdEdx 2.0 sigma, 1 sigma high pt -----------
    cuts.AddCutPCM("60100013", "00200009287602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009287602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009287602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009287602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009287602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 99){ // variation pdEdx 2.0 sigma, 1 sigma high pt - added signal
    cuts.AddCutPCM("60100023", "00200009287602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009287602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009287602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009287602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009287602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 100){ // variation pdEdx 2.0 sigma, 1 sigma high pt with phi cut
    cuts.AddCutPCM("60100013", "00216609287602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609287602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609287602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609287602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609287602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 101){ // variation pdEdx 2.0 sigma, 1 sigma high pt with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609287602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609287602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609287602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609287602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609287602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 102){ // variation pdEdx 3.0 sigma, no high pt ---------------
    cuts.AddCutPCM("60100013", "00200009247002008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009247002008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009247002008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009247002008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009247002008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 103){ // variation pdEdx 3.0 sigma, no high pt - added signal
    cuts.AddCutPCM("60100023", "00200009247002008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009247002008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009247002008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009247002008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009247002008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 104){ // variation pdEdx 3.0 sigma, no high pt with phi cut
    cuts.AddCutPCM("60100013", "00216609247002008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609247002008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609247002008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609247002008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609247002008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 105){ // variation pdEdx 3.0 sigma, no high pt  with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609247002008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609247002008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609247002008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609247002008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609247002008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 106){ // variation pdEdx std (3.0) sigma, but 0.3-3. ----------
    cuts.AddCutPCM("60100013", "00200009245402008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009245402008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009245402008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009245402008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009245402008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 107){ //variation pdEdx 3.0 sigma,  0.3-3. - added signal
    cuts.AddCutPCM("60100023", "00200009245402008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009245402008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009245402008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009245402008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009245402008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 108){ //variation pdEdx 3.0 sigma, 0.3-3. with phi cut
    cuts.AddCutPCM("60100013", "00216609245402008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609245402008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609245402008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609245402008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609245402008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 109){ //variation pdEdx 3.0 sigma, 0.3-3. with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609245402008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609245402008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609245402008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609245402008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609245402008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 110){ // TOF el. PID -3,5 --------------------------------------
    cuts.AddCutPCM("60100013", "00200009247603008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009247603008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009247603008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009247603008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009247603008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 111){ // TOF el. PID -3,5 - added signal
    cuts.AddCutPCM("60100023", "00200009247603008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009247603008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009247603008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009247603008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009247603008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 112){ // TOF el. PID -3,5 with phi cut
    cuts.AddCutPCM("60100013", "00216609247603008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609247603008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609247603008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609247603008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609247603008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 113){ // TOF el. PID -3,5 with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609247603008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609247603008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609247603008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609247603008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609247603008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 114){ // TOF el. PID -2,3 ---------------------------------------
    cuts.AddCutPCM("60100013", "00200009247604008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009247604008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009247604008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009247604008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009247604008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 115){ // TOF el. PID -2,3 - added signal
    cuts.AddCutPCM("60100023", "00200009247604008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009247604008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009247604008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009247604008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009247604008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 116){ // TOF el. PID -2,3 with phi cut
    cuts.AddCutPCM("60100013", "00216609247604008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609247604008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609247604008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609247604008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609247604008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 117){ // TOF el. PID -2,3  with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609247604008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609247604008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609247604008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609247604008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609247604008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 118){ // qt 0.03 2D ----------------------------------------
    cuts.AddCutPCM("60100013", "00200009247602009850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009247602009850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009247602009850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009247602009850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009247602009850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 119){ // qt 0.03 2D - added signal
    cuts.AddCutPCM("60100023", "00200009247602009850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009247602009850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009247602009850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009247602009850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009247602009850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 120){ // qt 0.03 2D with phi cut
    cuts.AddCutPCM("60100013", "00216609247602009850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609247602009850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609247602009850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609247602009850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609247602009850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 121){ // qt 0.03 2D with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609247602009850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609247602009850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609247602009850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609247602009850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609247602009850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 122){ // qt 0.06 2D ---------------------------------------
    cuts.AddCutPCM("60100013", "00200009247602002850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009247602002850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009247602002850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009247602002850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009247602002850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 123){ // qt 0.06 2D - added signal
    cuts.AddCutPCM("60100023", "00200009247602002850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009247602002850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009247602002850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009247602002850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009247602002850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 124){ // qt 0.06 2D with phi cut
    cuts.AddCutPCM("60100013", "00216609247602002850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609247602002850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609247602002850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609247602002850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609247602002850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 125){ // qt 0.06 2D with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609247602002850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609247602002850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609247602002850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609247602002850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609247602002850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 126){ // chi2  50. ----------------------------------------
    cuts.AddCutPCM("60100013", "00200009247602008150404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009247602008150404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009247602008150404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009247602008150404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009247602008150404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 127){ // chi2  50.  added signal
    cuts.AddCutPCM("60100023", "00200009247602008150404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009247602008150404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009247602008150404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009247602008150404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009247602008150404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 128){ // chi2  50. with phi cut
    cuts.AddCutPCM("60100013", "00216609247602008150404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609247602008150404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609247602008150404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609247602008150404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609247602008150404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 129){ // chi2  50. with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609247602008150404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609247602008150404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609247602008150404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609247602008150404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609247602008150404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 130){ // chi2  30. -----------------------------------------
    cuts.AddCutPCM("60100013", "00200009247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 131){ // chi2  30.  added signal
    cuts.AddCutPCM("60100023", "00200009247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 132){ // chi2  30. with phi cut
    cuts.AddCutPCM("60100013", "00216609247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 133){ // chi2  30. with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 134){ // psi pair 0.05 2D ---------------------------------
    cuts.AddCutPCM("60100013", "00200009247602008860404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009247602008860404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009247602008860404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009247602008860404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009247602008860404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 135){ // psi pair 0.05 2D - added signal
    cuts.AddCutPCM("60100023", "00200009247602008860404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009247602008860404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009247602008860404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009247602008860404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009247602008860404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 136){ // psi pair 0.05 2D with phi cut
    cuts.AddCutPCM("60100013", "00216609247602008860404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609247602008860404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609247602008860404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609247602008860404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609247602008860404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 137){ // psi pair 0.05 2D with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609247602008860404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609247602008860404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609247602008860404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609247602008860404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609247602008860404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 138){ // psi pair 0.2 2D ----------------------------------
    cuts.AddCutPCM("60100013", "00200009247602008880404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009247602008880404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009247602008880404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009247602008880404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009247602008880404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 139){ // psi pair 0.2 2D - added signal
    cuts.AddCutPCM("60100023", "00200009247602008880404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009247602008880404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009247602008880404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009247602008880404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009247602008880404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 140){ // psi pair 0.2 2D with phi cut
    cuts.AddCutPCM("60100013", "00216609247602008880404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609247602008880404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609247602008880404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609247602008880404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609247602008880404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 141){ // psi pair 0.2 2D with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609247602008880404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609247602008880404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609247602008880404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609247602008880404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609247602008880404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 142){ // cosPA -1 ----------------------------------------
    cuts.AddCutPCM("60100013", "00200009247602008850004000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009247602008850004000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009247602008850004000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009247602008850004000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009247602008850004000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 143){ // cosPA -1 - added signal
    cuts.AddCutPCM("60100023", "00200009247602008850004000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009247602008850004000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009247602008850004000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009247602008850004000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009247602008850004000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 144){ // cosPA -1 with phi cut
    cuts.AddCutPCM("60100013", "00216609247602008850004000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609247602008850004000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609247602008850004000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609247602008850004000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609247602008850004000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 145){ // cosPA -1 with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609247602008850004000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609247602008850004000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609247602008850004000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609247602008850004000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609247602008850004000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 146){ // variation alpha 0.75 ------------------------------
    cuts.AddCutPCM("60100013", "00200009247602008850404000", "0152505500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009247602008850404000", "0152505500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009247602008850404000", "0152505500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009247602008850404000", "0152505500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009247602008850404000", "0152505500000000"); // 20-50%
  } else if ( trainConfig == 147){ // variation alpha 0.75 - added signal
    cuts.AddCutPCM("60100023", "00200009247602008850404000", "0152505500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009247602008850404000", "0152505500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009247602008850404000", "0152505500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009247602008850404000", "0152505500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009247602008850404000", "0152505500000000"); // 20-50%
  } else if ( trainConfig == 148){ // variation alpha 0.75 with phi cut
    cuts.AddCutPCM("60100013", "00216609247602008850404000", "0152505500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609247602008850404000", "0152505500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609247602008850404000", "0152505500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609247602008850404000", "0152505500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609247602008850404000", "0152505500000000"); // 20-50%
  } else if ( trainConfig == 149){ // variation alpha 0.75 with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609247602008850404000", "0152505500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609247602008850404000", "0152505500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609247602008850404000", "0152505500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609247602008850404000", "0152505500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609247602008850404000", "0152505500000000"); // 20-50%
  } else if ( trainConfig == 150){ // variation alpha 1. --------------------------------
    cuts.AddCutPCM("60100013", "00200009247602008850404000", "0152503500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009247602008850404000", "0152503500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009247602008850404000", "0152503500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009247602008850404000", "0152503500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009247602008850404000", "0152503500000000"); // 20-50%
  } else if ( trainConfig == 151){ // variation alpha 1. - added signal
    cuts.AddCutPCM("60100023", "00200009247602008850404000", "0152503500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009247602008850404000", "0152503500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009247602008850404000", "0152503500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009247602008850404000", "0152503500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009247602008850404000", "0152503500000000"); // 20-50%
  } else if ( trainConfig == 152){ // variation alpha 1. with phi cut
    cuts.AddCutPCM("60100013", "00216609247602008850404000", "0152503500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609247602008850404000", "0152503500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609247602008850404000", "0152503500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609247602008850404000", "0152503500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609247602008850404000", "0152503500000000"); // 20-50%
  } else if ( trainConfig == 153){ // variation alpha 1. with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609247602008850404000", "0152503500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609247602008850404000", "0152503500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609247602008850404000", "0152503500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609247602008850404000", "0152503500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609247602008850404000", "0152503500000000"); // 20-50%
  } else if ( trainConfig == 154){ // variation with phi cut at 2.0 - 4.0 --------------
    cuts.AddCutPCM("60100013", "00215509247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00215509247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00215509247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00215509247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00215509247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 155){ // variation with phi cut at 2.0 - 4.0 - added signal
    cuts.AddCutPCM("60100023", "00215509247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00215509247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00215509247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00215509247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00215509247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 156){ // variation with phi cut at 2.4 - 3.6 --------------
    cuts.AddCutPCM("60100013", "00217709247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00217709247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00217709247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00217709247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00217709247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 157){ // variation with phi cut at 2.4 - 3.6 - added signal
    cuts.AddCutPCM("60100023", "00217709247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00217709247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00217709247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00217709247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00217709247602008850404000", "0152501500000000"); // 20-50%
    //-------------------- end LHC11h syst cut variations ---------------------------

  } else if ( trainConfig == 158){ // standard LHC11h cut selection
    cuts.AddCutPCM("60100013", "00200009247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 159){ // standard LHC11h cut selection - double rejec - added signal
    cuts.AddCutPCM("60100023", "00200009247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 160){ // standard LHC11h cut selection - double rejec with phi cut
    cuts.AddCutPCM("60100013", "00216609247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 161){ // standard LHC11h cut selection - double rejec with phi cut - added signal
    cuts.AddCutPCM("60100023", "00216609247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 162){ // RP EM V0 mult background
    cuts.AddCutPCM("60100013", "00200009247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("61200013", "00200009247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("50100013", "00200009247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("52400013", "00200009247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("52500013", "00200009247602008850404000", "0652501500000000"); //
  } else if ( trainConfig == 163){ // added signals
    cuts.AddCutPCM("60100023", "00200009247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("61200023", "00200009247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("50100023", "00200009247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("52400023", "00200009247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("52500023", "00200009247602008850404000", "0652501500000000"); //
  } else if ( trainConfig == 164){ //  with phi cut
    cuts.AddCutPCM("60100013", "00216609247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("61200013", "00216609247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("50100013", "00216609247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("52400013", "00216609247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("52500013", "00216609247602008850404000", "0652501500000000"); //
  } else if ( trainConfig == 165){ //  with phi cut - added signals
    cuts.AddCutPCM("60100023", "00216609247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("61200023", "00216609247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("50100023", "00216609247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("52400023", "00216609247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("52500023", "00216609247602008850404000", "0652501500000000"); //

  } else if ( trainConfig == 166){ // cleaner cuts photon Quality 1
    cuts.AddCutPCM("60100013", "00200009297002008250420000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009297002008250420000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009297002008250420000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009297002008250420000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("52500013", "00200009297002008250420000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 167){ // cleaner cuts added signal photon Quality 1
    cuts.AddCutPCM("60100023", "00200009297002008250420000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009297002008250420000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009297002008250420000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009297002008250420000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("52500023", "00200009297002008250420000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 168){ // cleaner cuts photon Quality 2
    cuts.AddCutPCM("60100013", "00200009297002008250430000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009297002008250430000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009297002008250430000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009297002008250430000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("52500013", "00200009297002008250430000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 169){ // cleaner cuts added signal photon Quality 2
    cuts.AddCutPCM("60100023", "00200009297002008250430000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009297002008250430000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009297002008250430000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009297002008250430000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("52500023", "00200009297002008250430000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 170){ // cleaner cuts photon Quality 3
    cuts.AddCutPCM("60100013", "00200009297002008250440000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009297002008250440000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009297002008250440000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009297002008250440000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("52500013", "00200009297002008250440000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 171){ // cleaner cuts added signal photon Quality 3
    cuts.AddCutPCM("60100023", "00200009297002008250440000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009297002008250440000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009297002008250440000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009297002008250440000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("52500023", "00200009297002008250440000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 172){ // cleaner cuts, photon Quality 1, min R = 35 cm
    cuts.AddCutPCM("60100013", "00700009297002008250420000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00700009297002008250420000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00700009297002008250420000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00700009297002008250420000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("52500013", "00700009297002008250420000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 173){ // cleaner cuts added signal, photon Quality 1, min R = 35 cm
    cuts.AddCutPCM("60100023", "00700009297002008250420000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00700009297002008250420000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00700009297002008250420000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00700009297002008250420000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("52500023", "00700009297002008250420000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 174){ // cleaner cuts, photon Quality 3, min R = 35 cm
    cuts.AddCutPCM("60100013", "00700009297002008250440000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00700009297002008250440000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00700009297002008250440000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00700009297002008250440000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("52500013", "00700009297002008250440000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 175){ // cleaner cuts added signal, photon Quality 3, min R = 35 cm
    cuts.AddCutPCM("60100023", "00700009297002008250440000", "0152506500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00700009297002008250440000", "0152506500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00700009297002008250440000", "0152506500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00700009297002008250440000", "0152506500000000"); // 10-20%
    cuts.AddCutPCM("52500023", "00700009297002008250440000", "0152506500000000"); // 0-20%

  } else if ( trainConfig == 176){ // flow cuts with eta = 0.9, y = 0.85
    cuts.AddCutPCM("60100013", "00200009297002008250400000", "0152506500000000");
    cuts.AddCutPCM("61200013", "00200009297002008250400000", "0152506500000000");
    cuts.AddCutPCM("51200013", "00200009297002008250400000", "0152506500000000");
    cuts.AddCutPCM("52300013", "00200009297002008250400000", "0152506500000000");
    cuts.AddCutPCM("53400013", "00200009297002008250400000", "0152506500000000");
    cuts.AddCutPCM("54600013", "00200009297002008250400000", "0152506500000000");
    cuts.AddCutPCM("56800013", "00200009297002008250400000", "0152506500000000");
  } else if ( trainConfig == 177){ // flow cuts with eta = 0.65, y = 0.6
    cuts.AddCutPCM("60100013", "03200009297002008250400000", "0152306500000000");
    cuts.AddCutPCM("61200013", "03200009297002008250400000", "0152306500000000");
    cuts.AddCutPCM("51200013", "03200009297002008250400000", "0152306500000000");
    cuts.AddCutPCM("52300013", "03200009297002008250400000", "0152306500000000");
    cuts.AddCutPCM("53400013", "03200009297002008250400000", "0152306500000000");
    cuts.AddCutPCM("54600013", "03200009297002008250400000", "0152306500000000");
    cuts.AddCutPCM("56800013", "03200009297002008250400000", "0152306500000000");
  } else if ( trainConfig == 178){ // flow cuts with eta = 0.6, y = 0.5
    cuts.AddCutPCM("60100013", "01200009297002008250400000", "0152406500000000");
    cuts.AddCutPCM("61200013", "01200009297002008250400000", "0152406500000000");
    cuts.AddCutPCM("51200013", "01200009297002008250400000", "0152406500000000");
    cuts.AddCutPCM("52300013", "01200009297002008250400000", "0152406500000000");
    cuts.AddCutPCM("53400013", "01200009297002008250400000", "0152406500000000");
    cuts.AddCutPCM("54600013", "01200009297002008250400000", "0152406500000000");
    cuts.AddCutPCM("56800013", "01200009297002008250400000", "0152406500000000");
  } else if ( trainConfig == 179){ // flow cuts with eta = 0.9, y = 0.85
    cuts.AddCutPCM("60100013", "00200009297002208250400000", "0152506500000000");
    cuts.AddCutPCM("61200013", "00200009297002208250400000", "0152506500000000");
    cuts.AddCutPCM("51200013", "00200009297002208250400000", "0152506500000000");
    cuts.AddCutPCM("52300013", "00200009297002208250400000", "0152506500000000");
    cuts.AddCutPCM("53400013", "00200009297002208250400000", "0152506500000000");
    cuts.AddCutPCM("54600013", "00200009297002208250400000", "0152506500000000");
    cuts.AddCutPCM("56800013", "00200009297002208250400000", "0152506500000000");
  } else if ( trainConfig == 180){ // flow cuts with eta = 0.65, y = 0.6
    cuts.AddCutPCM("60100013", "03200009297002208250400000", "0152306500000000");
    cuts.AddCutPCM("61200013", "03200009297002208250400000", "0152306500000000");
    cuts.AddCutPCM("51200013", "03200009297002208250400000", "0152306500000000");
    cuts.AddCutPCM("52300013", "03200009297002208250400000", "0152306500000000");
    cuts.AddCutPCM("53400013", "03200009297002208250400000", "0152306500000000");
    cuts.AddCutPCM("54600013", "03200009297002208250400000", "0152306500000000");
    cuts.AddCutPCM("56800013", "03200009297002208250400000", "0152306500000000");
  } else if ( trainConfig == 181){ // flow cuts with eta = 0.6, y = 0.5
    cuts.AddCutPCM("60100013", "01200009297002208250400000", "0152406500000000");
    cuts.AddCutPCM("61200013", "01200009297002208250400000", "0152406500000000");
    cuts.AddCutPCM("51200013", "01200009297002208250400000", "0152406500000000");
    cuts.AddCutPCM("52300013", "01200009297002208250400000", "0152406500000000");
    cuts.AddCutPCM("53400013", "01200009297002208250400000", "0152406500000000");
    cuts.AddCutPCM("54600013", "01200009297002208250400000", "0152406500000000");
    cuts.AddCutPCM("56800013", "01200009297002208250400000", "0152406500000000");

  } else if ( trainConfig == 182){ // standard LHC11h cut selection - no dEdx
    cuts.AddCutPCM("60100013", "00200009000002008250400000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009000002008250400000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009000002008250400000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00200009000002008250400000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00200009000002008250400000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 183){ // standard LHC11h cut selection - no dEdx add signals
    cuts.AddCutPCM("60100023", "00200009000002008250400000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00200009000002008250400000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00200009000002008250400000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00200009000002008250400000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00200009000002008250400000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 184){ // standard LHC11h cut selection - no dEdx with phi cut
    cuts.AddCutPCM("60100013", "00216609000002008250400000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00216609000002008250400000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00216609000002008250400000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00216609000002008250400000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00216609000002008250400000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 185){ // standard LHC11h cut selection - no dEdx with phi cut add signals
    cuts.AddCutPCM("60100023", "00216609000002008250400000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00216609000002008250400000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00216609000002008250400000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00216609000002008250400000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00216609000002008250400000", "0152501500000000"); // 20-50%

  } else if ( trainConfig == 186){ // dir photons
    cuts.AddCutPCM("50100013", "00200009247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 187){ //
    cuts.AddCutPCM("50100013", "00216609247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 188){ // dir photons
    cuts.AddCutPCM("52400013", "00200009247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("52500013", "00200009247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 189){ //
    cuts.AddCutPCM("52400013", "00216609247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("52500013", "00216609247002008850404000", "0152501500000000"); //

  } else if ( trainConfig == 190){ // phi variation 1
    cuts.AddCutPCM("50100013", "00215509247002008850404000", "0152501500000000"); // 2 - 4 rad
  } else if ( trainConfig == 191){ // phi variation 2
    cuts.AddCutPCM("50100013", "00217709247002008850404000", "0152501500000000"); // 2.6 - 3.4 rad
  } else if ( trainConfig == 192){ // phi variation 1
    cuts.AddCutPCM("52400013", "00215509247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("52500013", "00215509247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 193){ // phi variation 2
    cuts.AddCutPCM("52400013", "00217709247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("52500013", "00217709247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 194){ // single pt variation
    cuts.AddCutPCM("50100013", "00200049247002008850404000", "0152501500000000"); // 0.075
    cuts.AddCutPCM("50100013", "00200019247002008850404000", "0152501500000000"); // 0.1
  } else if ( trainConfig == 195){ //  single pt variation - phi
    cuts.AddCutPCM("50100013", "00216649247002008850404000", "0152501500000000"); // 0.075
    cuts.AddCutPCM("50100013", "00216619247002008850404000", "0152501500000000"); // 0.1
  } else if ( trainConfig == 196){ // single pt variation
    cuts.AddCutPCM("52400013", "00200049247002008850404000", "0152501500000000"); // 0.075
    cuts.AddCutPCM("52400013", "00200019247002008850404000", "0152501500000000"); // 0.1
    cuts.AddCutPCM("52500013", "00200049247002008850404000", "0152501500000000"); // 0.075
    cuts.AddCutPCM("52500013", "00200019247002008850404000", "0152501500000000"); // 0.1
  } else if ( trainConfig == 197){ //  single pt variation - phi
    cuts.AddCutPCM("52400013", "00216649247002008850404000", "0152501500000000"); // 0.075
    cuts.AddCutPCM("52400013", "00216619247002008850404000", "0152501500000000"); // 0.1
    cuts.AddCutPCM("52500013", "00216649247002008850404000", "0152501500000000"); // 0.075
    cuts.AddCutPCM("52500013", "00216619247002008850404000", "0152501500000000"); // 0.1
  } else if ( trainConfig == 198){ // TPC clusters
    cuts.AddCutPCM("50100013", "00200008247002008850404000", "0152501500000000"); // 0.35
    cuts.AddCutPCM("50100013", "00200006247002008850404000", "0152501500000000"); // 0.7
  } else if ( trainConfig == 199){ // TPC clusters - phi
    cuts.AddCutPCM("50100013", "00216608247002008850404000", "0152501500000000"); // 0.35
    cuts.AddCutPCM("50100013", "00216606247002008850404000", "0152501500000000"); // 0.7
  } else if ( trainConfig == 200){ // TPC clusters
    cuts.AddCutPCM("52400013", "00200008247002008850404000", "0152501500000000"); // 0.35
    cuts.AddCutPCM("52400013", "00200006247002008850404000", "0152501500000000"); // 0.7
    cuts.AddCutPCM("52500013", "00200008247002008850404000", "0152501500000000"); // 0.35
    cuts.AddCutPCM("52500013", "00200006247002008850404000", "0152501500000000"); // 0.7
  } else if ( trainConfig == 201){ // TPC clusters - phi
    cuts.AddCutPCM("52400013", "00216608247002008850404000", "0152501500000000"); // 0.35
    cuts.AddCutPCM("52400013", "00216606247002008850404000", "0152501500000000"); // 0.7
    cuts.AddCutPCM("52500013", "00216608247002008850404000", "0152501500000000"); // 0.35
    cuts.AddCutPCM("52500013", "00216606247002008850404000", "0152501500000000"); // 0.7
  } else if ( trainConfig == 202){ // edEdx
    cuts.AddCutPCM("50100013", "00200009347002008850404000", "0152501500000000"); // -4, 5
    cuts.AddCutPCM("50100013", "00200009947002008850404000", "0152501500000000"); // -2.5, 5
  } else if ( trainConfig == 203){ // edEdx - phi
    cuts.AddCutPCM("50100013", "00216609347002008850404000", "0152501500000000"); // -4, 5
    cuts.AddCutPCM("50100013", "00216609947002008850404000", "0152501500000000"); // -2.5, 5
  } else if ( trainConfig == 204){ // edEdx
    cuts.AddCutPCM("52400013", "00200009347002008850404000", "0152501500000000"); // -4, 5
    cuts.AddCutPCM("52400013", "00200009947002008850404000", "0152501500000000"); // -2.5, 5
    cuts.AddCutPCM("52500013", "00200009347002008850404000", "0152501500000000"); // -4, 5
    cuts.AddCutPCM("52500013", "00200009947002008850404000", "0152501500000000"); // -2.5, 5
  } else if ( trainConfig == 205){ // edEdx - phi
    cuts.AddCutPCM("52400013", "00216609347002008850404000", "0152501500000000"); // -4, 5
    cuts.AddCutPCM("52400013", "00216609947002008850404000", "0152501500000000"); // -2.5, 5
    cuts.AddCutPCM("52500013", "00216609347002008850404000", "0152501500000000"); // -4, 5
    cuts.AddCutPCM("52500013", "00216609947002008850404000", "0152501500000000"); // -2.5, 5
  } else if ( trainConfig == 206){ // pdEdx
    cuts.AddCutPCM("50100013", "00200009240002008850404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCutPCM("50100013", "00200009247602008850404000", "0152501500000000"); // 3&1sigma
  } else if ( trainConfig == 207){ // pdEdx - phi
    cuts.AddCutPCM("50100013", "00216609240002008850404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCutPCM("50100013", "00216609247602008850404000", "0152501500000000"); // 3&1sigma
  } else if ( trainConfig == 208){ // pdEdx
    cuts.AddCutPCM("52400013", "00200009240002008850404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCutPCM("52400013", "00200009247602008850404000", "0152501500000000"); // 3&1sigma
    cuts.AddCutPCM("52500013", "00200009240002008850404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCutPCM("52500013", "00200009247602008850404000", "0152501500000000"); // 3&1sigma
  } else if ( trainConfig == 209){ // pdEdx - phi
    cuts.AddCutPCM("52400013", "00216609240002008850404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCutPCM("52400013", "00216609247602008850404000", "0152501500000000"); // 3&1sigma
    cuts.AddCutPCM("52500013", "00216609240002008850404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCutPCM("52500013", "00216609247602008850404000", "0152501500000000"); // 3&1sigma
  } else if ( trainConfig == 210){ // pdEdx & TOF
    cuts.AddCutPCM("50100013", "00200009237002008850404000", "0152501500000000"); // 2.5sigma
    cuts.AddCutPCM("50100013", "00200009247005008850404000", "0152501500000000"); // TOF -3, 3
  } else if ( trainConfig == 211){ // pdEdx & TOF - phi
    cuts.AddCutPCM("50100013", "00216609237002008850404000", "0152501500000000"); // 2.5sigma
    cuts.AddCutPCM("50100013", "00216609247005008850404000", "0152501500000000"); // TOF -3, 3
  } else if ( trainConfig == 212){ // pdEdx & TOF
    cuts.AddCutPCM("52400013", "00200009237002008850404000", "0152501500000000"); // 2.5sigma
    cuts.AddCutPCM("52400013", "00200009247005008850404000", "0152501500000000"); // TOF -3, 3
    cuts.AddCutPCM("52500013", "00200009237002008850404000", "0152501500000000"); // 2.5sigma
    cuts.AddCutPCM("52500013", "00200009247005008850404000", "0152501500000000"); // TOF -3, 3
  } else if ( trainConfig == 213){ // pdEdx & TOF - phi
    cuts.AddCutPCM("52400013", "00216609237002008850404000", "0152501500000000"); // 2.5sigma
    cuts.AddCutPCM("52400013", "00216609247005008850404000", "0152501500000000"); // TOF -3, 3
    cuts.AddCutPCM("52500013", "00216609237002008850404000", "0152501500000000"); // 2.5sigma
    cuts.AddCutPCM("52500013", "00216609247005008850404000", "0152501500000000"); // TOF -3, 3
  } else if ( trainConfig == 214){ // qT
    cuts.AddCutPCM("50100013", "00200009247002002850404000", "0152501500000000"); // 0.06 2D
    cuts.AddCutPCM("50100013", "00200009247002009850404000", "0152501500000000"); // 0.03 2D
  } else if ( trainConfig == 215){ // qT - phi
    cuts.AddCutPCM("50100013", "00216609247002002850404000", "0152501500000000"); // 0.06 2D
    cuts.AddCutPCM("50100013", "00216609247002009850404000", "0152501500000000"); // 0.03 2D
  } else if ( trainConfig == 216){ // qT
    cuts.AddCutPCM("52400013", "00200009247002002850404000", "0152501500000000"); // 0.06 2D
    cuts.AddCutPCM("52400013", "00200009247002009850404000", "0152501500000000"); // 0.03 2D
    cuts.AddCutPCM("52500013", "00200009247002002850404000", "0152501500000000"); // 0.06 2D
    cuts.AddCutPCM("52500013", "00200009247002009850404000", "0152501500000000"); // 0.03 2D
  } else if ( trainConfig == 217){ // qT - phi
    cuts.AddCutPCM("52400013", "00216609247002002850404000", "0152501500000000"); // 0.06 2D
    cuts.AddCutPCM("52400013", "00216609247002009850404000", "0152501500000000"); // 0.03 2D
    cuts.AddCutPCM("52500013", "00216609247002002850404000", "0152501500000000"); // 0.06 2D
    cuts.AddCutPCM("52500013", "00216609247002009850404000", "0152501500000000"); // 0.03 2D
  } else if ( trainConfig == 218){ // chi2 + psi pair
    cuts.AddCutPCM("50100013", "00200009247002008750404000", "0152501500000000"); // 10+0.1 2D
    cuts.AddCutPCM("50100013", "00200009247002008a50404000", "0152501500000000"); // 25+0.1 2D
  } else if ( trainConfig == 219){ // chi2 + psi pair - phi
    cuts.AddCutPCM("50100013", "00216609247002008750404000", "0152501500000000"); // 10+0.1 2D
    cuts.AddCutPCM("50100013", "00216609247002008a50404000", "0152501500000000"); // 25+0.1 2D
  } else if ( trainConfig == 220){ // chi2 + psi pair
    cuts.AddCutPCM("52400013", "00200009247002008750404000", "0152501500000000"); // 10+0.1 2D
    cuts.AddCutPCM("52400013", "00200009247002008a50404000", "0152501500000000"); // 25+0.1 2D
    cuts.AddCutPCM("52500013", "00200009247002008750404000", "0152501500000000"); // 10+0.1 2D
    cuts.AddCutPCM("52500013", "00200009247002008a50404000", "0152501500000000"); // 25+0.1 2D
  } else if ( trainConfig == 221){ // chi2 + psi pair - phi
    cuts.AddCutPCM("52400013", "00216609247002008750404000", "0152501500000000"); // 10+0.1 2D
    cuts.AddCutPCM("52400013", "00216609247002008a50404000", "0152501500000000"); // 25+0.1 2D
    cuts.AddCutPCM("52500013", "00216609247002008750404000", "0152501500000000"); // 10+0.1 2D
    cuts.AddCutPCM("52500013", "00216609247002008a50404000", "0152501500000000"); // 25+0.1 2D
  } else if ( trainConfig == 222){ // chi2 + psi pair
    cuts.AddCutPCM("50100013", "00200009247002008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCutPCM("50100013", "00200009247002008870404000", "0152501500000000"); // 20+0.07 2D
  } else if ( trainConfig == 223){ // chi2 + psi pair - phi
    cuts.AddCutPCM("50100013", "00216609247002008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCutPCM("50100013", "00216609247002008870404000", "0152501500000000"); // 20+0.07 2D
  } else if ( trainConfig == 224){ // chi2 + psi pair
    cuts.AddCutPCM("52400013", "00200009247002008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCutPCM("52400013", "00200009247002008870404000", "0152501500000000"); // 20+0.07 2D
    cuts.AddCutPCM("52500013", "00200009247002008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCutPCM("52500013", "00200009247002008870404000", "0152501500000000"); // 20+0.07 2D
  } else if ( trainConfig == 225){ // chi2 + psi pair - phi
    cuts.AddCutPCM("52400013", "00216609247002008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCutPCM("52400013", "00216609247002008870404000", "0152501500000000"); // 20+0.07 2D
    cuts.AddCutPCM("52500013", "00216609247002008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCutPCM("52500013", "00216609247002008870404000", "0152501500000000"); // 20+0.07 2D
  } else if ( trainConfig == 226){ // asym 1D & cosPA
    cuts.AddCutPCM("50100013", "00200009247002008857404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCutPCM("50100013", "00200009247002008850004000", "0152501500000000"); // cosPA off
  } else if ( trainConfig == 227){ // asym 1D & cosPA - phi
    cuts.AddCutPCM("50100013", "00216609247002008857404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCutPCM("50100013", "00216609247002008850004000", "0152501500000000"); // cosPA off
  } else if ( trainConfig == 228){ // asym 1D & cosPA
    cuts.AddCutPCM("52400013", "00200009247002008857404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCutPCM("52400013", "00200009247002008850004000", "0152501500000000"); // cosPA off
    cuts.AddCutPCM("52500013", "00200009247002008857404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCutPCM("52500013", "00200009247002008850004000", "0152501500000000"); // cosPA off
  } else if ( trainConfig == 229){ // asym 1D & cosPA - phi
    cuts.AddCutPCM("52400013", "00216609247002008857404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCutPCM("52400013", "00216609247002008850004000", "0152501500000000"); // cosPA off
    cuts.AddCutPCM("52500013", "00216609247002008857404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCutPCM("52500013", "00216609247002008850004000", "0152501500000000"); // cosPA off
  } else if ( trainConfig == 230){ // asym 1D  & proton rejection
    cuts.AddCutPCM("50100013", "00200009247002008856404000", "0152501500000000"); // asym > 6GeV
    cuts.AddCutPCM("50100013", "00200009247082008850404000", "0152501500000000"); // asym > 8GeV
  } else if ( trainConfig == 231){ // asym 1D  & proton rejection - phi
    cuts.AddCutPCM("50100013", "00216609247002008856404000", "0152501500000000"); // asym > 6GeV
    cuts.AddCutPCM("50100013", "00216609247082008850404000", "0152501500000000"); // asym > 8GeV
  } else if ( trainConfig == 232){ // asym 1D  & proton rejection
    cuts.AddCutPCM("52400013", "00200009247002008856404000", "0152501500000000"); // asym > 6GeV
    cuts.AddCutPCM("52400013", "00200009247082008850404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCutPCM("52500013", "00200009247002008856404000", "0152501500000000"); // asym > 6GeV
    cuts.AddCutPCM("52500013", "00200009247082008850404000", "0152501500000000"); // asym > 8GeV
  } else if ( trainConfig == 233){ // asym 1D  & proton rejection - phi
    cuts.AddCutPCM("52400013", "00216609247002008856404000", "0152501500000000"); // asym > 6GeV
    cuts.AddCutPCM("52400013", "00216609247082008850404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCutPCM("52500013", "00216609247002008856404000", "0152501500000000"); // asym > 6GeV
    cuts.AddCutPCM("52500013", "00216609247082008850404000", "0152501500000000"); // asym > 8GeV

  } else if ( trainConfig == 234){ // standard LHC11h cut selection with R bins
    cuts.AddCutPCM("60100013", "00a00009247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00a00009247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00a00009247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00a00009247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00a00009247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 235){ //
    cuts.AddCutPCM("60100023", "00a00009247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00a00009247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00a00009247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00a00009247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00a00009247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 236){ //
    cuts.AddCutPCM("60100013", "00a16609247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00a16609247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00a16609247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00a16609247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00a16609247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 237){ //
    cuts.AddCutPCM("60100023", "00a16609247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00a16609247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00a16609247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00a16609247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00a16609247602008850404000", "0152501500000000"); // 20-50%

  } else if ( trainConfig == 238){ // standard LHC11h cut selection with R bins
    cuts.AddCutPCM("60100013", "00b00009247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00b00009247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00b00009247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00b00009247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00b00009247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 239){ //
    cuts.AddCutPCM("60100023", "00b00009247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00b00009247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00b00009247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00b00009247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00b00009247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 240){ //
    cuts.AddCutPCM("60100013", "00b16609247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00b16609247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00b16609247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00b16609247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00b16609247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 241){ //
    cuts.AddCutPCM("60100023", "00b16609247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00b16609247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00b16609247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00b16609247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00b16609247602008850404000", "0152501500000000"); // 20-50%

  } else if ( trainConfig == 242){ // standard LHC11h cut selection with R bins
    cuts.AddCutPCM("60100013", "00c00009247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00c00009247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00c00009247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00c00009247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00c00009247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 243){ //
    cuts.AddCutPCM("60100023", "00c00009247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00c00009247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00c00009247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00c00009247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00c00009247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 244){ //
    cuts.AddCutPCM("60100013", "00c16609247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200013", "00c16609247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100013", "00c16609247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400013", "00c16609247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500013", "00c16609247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 245){ //
    cuts.AddCutPCM("60100023", "00c16609247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCutPCM("61200023", "00c16609247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("50100023", "00c16609247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("52400023", "00c16609247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("52500023", "00c16609247602008850404000", "0152501500000000"); // 20-50%

  //****************************************************************************************************
  // 5.02TeV Pb-Pb LHC15o
  //****************************************************************************************************
  } else if (trainConfig == 246){ // LHC15o, kINT7, cent. from V0M, reject added particles
    cuts.AddCutPCM("10110013", "00200009247602008250404000", "0652501500000000"); //  0-10%
    cuts.AddCutPCM("11210013", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCutPCM("12510013", "00200009247602008250404000", "0652501500000000"); // 20-50%
    cuts.AddCutPCM("15910013", "00200009247602008250404000", "0652501500000000"); // 50-90%
    cuts.AddCutPCM("10910013", "00200009247602008250404000", "0652501500000000"); //  0-90%
  } else if (trainConfig == 247){  // LHC15o, kINT7, cent. from V0M, user defined header
    cuts.AddCutPCM("10110023", "00200009247602008250404000", "0652501500000000"); //  0-10%
    cuts.AddCutPCM("11210023", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCutPCM("12510023", "00200009247602008250404000", "0652501500000000"); // 20-50%
    cuts.AddCutPCM("15910023", "00200009247602008250404000", "0652501500000000"); // 50-90%
    cuts.AddCutPCM("10910023", "00200009247602008250404000", "0652501500000000"); //  0-90%
  } else if (trainConfig == 248){   // LHC15o, kINT7,  cent. from V0M, reject added particles
    cuts.AddCutPCM("30110013", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCutPCM("31210013", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCutPCM("32310013", "00200009247602008250404000", "0652501500000000"); // 10-15%
    cuts.AddCutPCM("33410013", "00200009247602008250404000", "0652501500000000"); // 15-20%
    cuts.AddCutPCM("34510013", "00200009247602008250404000", "0652501500000000"); // 20-25%
  } else if (trainConfig == 249){   //  LHC15o, kINT7,  cent. from V0M, user defined header
    cuts.AddCutPCM("30110023", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCutPCM("31210023", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCutPCM("32310023", "00200009247602008250404000", "0652501500000000"); // 10-15%
    cuts.AddCutPCM("33410023", "00200009247602008250404000", "0652501500000000"); // 15-20%
    cuts.AddCutPCM("34510023", "00200009247602008250404000", "0652501500000000"); // 20-25%
  } else if (trainConfig == 250){   // LHC15o, kINT7,  cent. from V0M, reject added particles
    cuts.AddCutPCM("35610013", "00200009247602008250404000", "0652501500000000"); // 25-30%
    cuts.AddCutPCM("36710013", "00200009247602008250404000", "0652501500000000"); // 30-35%
    cuts.AddCutPCM("37810013", "00200009247602008250404000", "0652501500000000"); // 35-40%
    cuts.AddCutPCM("38910013", "00200009247602008250404000", "0652501500000000"); // 40-45%
  } else if (trainConfig == 251){   // LHC15o, kINT7,  cent. from V0M, user defined header
    cuts.AddCutPCM("35610023", "00200009247602008250404000", "0652501500000000"); // 25-30%
    cuts.AddCutPCM("36710023", "00200009247602008250404000", "0652501500000000"); // 30-35%
    cuts.AddCutPCM("37810023", "00200009247602008250404000", "0652501500000000"); // 35-40%
    cuts.AddCutPCM("38910023", "00200009247602008250404000", "0652501500000000"); // 40-45%
  } else if (trainConfig == 252){   // LHC15o, kINT7,  cent. from V0M, reject added particles
    // 0 -0%, 1-5%, 2-10%, 3-15%, 4-20%, 5-25%, 6-30%, 7-35%, 8-40%, 9-45%, a-50%, b-55%, c-60%, d-65%, e-70%, f-75%, g-80%, h-85%, i-90%, j-95%, k-100%
    cuts.AddCutPCM("39a10013", "00200009247602008250404000", "0652501500000000"); // 45-50%
    cuts.AddCutPCM("3ab10013", "00200009247602008250404000", "0652501500000000"); // 50-55%
    cuts.AddCutPCM("3bc10013", "00200009247602008250404000", "0652501500000000"); // 55-60%
    cuts.AddCutPCM("3cd10013", "00200009247602008250404000", "0652501500000000"); // 60-65%
    cuts.AddCutPCM("3de10013", "00200009247602008250404000", "0652501500000000"); // 65-70%
  } else if (trainConfig == 253){   // LHC15o, kINT7,  cent. from V0M, user defined header
    cuts.AddCutPCM("39a10023", "00200009247602008250404000", "0652501500000000"); // 45-50%
    cuts.AddCutPCM("3ab10023", "00200009247602008250404000", "0652501500000000"); // 50-55%
    cuts.AddCutPCM("3bc10023", "00200009247602008250404000", "0652501500000000"); // 55-60%
    cuts.AddCutPCM("3cd10023", "00200009247602008250404000", "0652501500000000"); // 60-65%
    cuts.AddCutPCM("3de10023", "00200009247602008250404000", "0652501500000000"); // 65-70%
  } else if (trainConfig == 254){   // LHC15o, kINT7,  cent. from V0M, reject added particles
    cuts.AddCutPCM("3ef10013", "00200009247602008250404000", "0652501500000000"); // 70-75%
    cuts.AddCutPCM("3fg10013", "00200009247602008250404000", "0652501500000000"); // 75-80%
    cuts.AddCutPCM("3gh10013", "00200009247602008250404000", "0652501500000000"); // 80-85%
    cuts.AddCutPCM("3hi10013", "00200009247602008250404000", "0652501500000000"); // 85-90%
  } else if (trainConfig == 255){   // LHC15o, kINT7,  cent. from V0M, user defined header
    cuts.AddCutPCM("3ef10023", "00200009247602008250404000", "0652501500000000"); // 70-75%
    cuts.AddCutPCM("3fg10023", "00200009247602008250404000", "0652501500000000"); // 75-80%
    cuts.AddCutPCM("3gh10023", "00200009247602008250404000", "0652501500000000"); // 80-85%
    cuts.AddCutPCM("3hi10023", "00200009247602008250404000", "0652501500000000"); // 85-90%
  } else if (trainConfig == 256){ // LHC15o, kINT7,  cent. from track mult., reject added particles
    cuts.AddCutPCM("50110013", "00200009247602008250404000", "0652501500000000"); //  0-10%
    cuts.AddCutPCM("51210013", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCutPCM("52510013", "00200009247602008250404000", "0652501500000000"); // 20-50%
    cuts.AddCutPCM("55910013", "00200009247602008250404000", "0652501500000000"); // 50-90%
    cuts.AddCutPCM("50910013", "00200009247602008250404000", "0652501500000000"); //  0-90%
  } else if (trainConfig == 257){  // LHC15o, kINT7,  cent. from track mult., user defined header
    cuts.AddCutPCM("50110023", "00200009247602008250404000", "0652501500000000"); //  0-10%
    cuts.AddCutPCM("51210023", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCutPCM("52510023", "00200009247602008250404000", "0652501500000000"); // 20-50%
    cuts.AddCutPCM("55910023", "00200009247602008250404000", "0652501500000000"); // 50-90%
    cuts.AddCutPCM("50910023", "00200009247602008250404000", "0652501500000000"); //  0-90%
  } else if (trainConfig == 258){   //  LHC15o, kINT7,  cent. from track mult., reject added particles
    cuts.AddCutPCM("60110013", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCutPCM("61210013", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCutPCM("62310013", "00200009247602008250404000", "0652501500000000"); // 10-15%
    cuts.AddCutPCM("63410013", "00200009247602008250404000", "0652501500000000"); // 15-20%
    cuts.AddCutPCM("64510013", "00200009247602008250404000", "0652501500000000"); // 20-25%
  } else if (trainConfig == 259){   //  LHC15o, kINT7,  cent. from track mult., user defined header
    cuts.AddCutPCM("60110023", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCutPCM("61210023", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCutPCM("62310023", "00200009247602008250404000", "0652501500000000"); // 10-15%
    cuts.AddCutPCM("63410023", "00200009247602008250404000", "0652501500000000"); // 15-20%
    cuts.AddCutPCM("64510023", "00200009247602008250404000", "0652501500000000"); // 20-25%
  } else if (trainConfig == 260){   // LHC15o, kINT7,  cent. from track mult., reject added particles
    cuts.AddCutPCM("65610013", "00200009247602008250404000", "0652501500000000"); // 25-30%
    cuts.AddCutPCM("66710013", "00200009247602008250404000", "0652501500000000"); // 30-35%
    cuts.AddCutPCM("67810013", "00200009247602008250404000", "0652501500000000"); // 35-40%
    cuts.AddCutPCM("68910013", "00200009247602008250404000", "0652501500000000"); // 40-45%
  } else if (trainConfig == 261){   // LHC15o, kINT7,  cent. from track mult., user defined header
    cuts.AddCutPCM("65610023", "00200009247602008250404000", "0652501500000000"); // 25-30%
    cuts.AddCutPCM("66710023", "00200009247602008250404000", "0652501500000000"); // 30-35%
    cuts.AddCutPCM("67810023", "00200009247602008250404000", "0652501500000000"); // 35-40%
    cuts.AddCutPCM("68910023", "00200009247602008250404000", "0652501500000000"); // 40-45%
  } else if (trainConfig == 262){   // LHC15o, kINT7,  cent. from track mult., reject added particles
    cuts.AddCutPCM("69a10013", "00200009247602008250404000", "0652501500000000"); // 45-50%
    cuts.AddCutPCM("6ab10013", "00200009247602008250404000", "0652501500000000"); // 50-55%
    cuts.AddCutPCM("6bc10013", "00200009247602008250404000", "0652501500000000"); // 55-60%
    cuts.AddCutPCM("6cd10013", "00200009247602008250404000", "0652501500000000"); // 60-65%
    cuts.AddCutPCM("6de10013", "00200009247602008250404000", "0652501500000000"); // 65-70%
  } else if (trainConfig == 263){   // LHC15o, kINT7, cent. from track mult., user defined header
    cuts.AddCutPCM("69a10023", "00200009247602008250404000", "0652501500000000"); // 45-50%
    cuts.AddCutPCM("6ab10023", "00200009247602008250404000", "0652501500000000"); // 50-55%
    cuts.AddCutPCM("6bc10023", "00200009247602008250404000", "0652501500000000"); // 55-60%
    cuts.AddCutPCM("6cd10023", "00200009247602008250404000", "0652501500000000"); // 60-65%
    cuts.AddCutPCM("6de10023", "00200009247602008250404000", "0652501500000000"); // 65-70%
  } else if (trainConfig == 264){   // LHC15o, kINT7, cent. from track mult., reject added particles
    cuts.AddCutPCM("6ef10013", "00200009247602008250404000", "0652501500000000"); // 70-75%
    cuts.AddCutPCM("6fg10013", "00200009247602008250404000", "0652501500000000"); // 75-80%
    cuts.AddCutPCM("6gh10013", "00200009247602008250404000", "0652501500000000"); // 80-85%
    cuts.AddCutPCM("6hi10013", "00200009247602008250404000", "0652501500000000"); // 85-90%
  } else if (trainConfig == 265){   // LHC15o, kINT7, cent. from track mult., user defined header
    cuts.AddCutPCM("6ef10023", "00200009247602008250404000", "0652501500000000"); // 70-75%
    cuts.AddCutPCM("6fg10023", "00200009247602008250404000", "0652501500000000"); // 75-80%
    cuts.AddCutPCM("6gh10023", "00200009247602008250404000", "0652501500000000"); // 80-85%
    cuts.AddCutPCM("6hi10023", "00200009247602008250404000", "0652501500000000"); // 85-90%
  } else if (trainConfig == 266){ // LHC15o, kINT7,  cent. from track mult., reject added particles
    cuts.AddCutPCM("60110013", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCutPCM("61210013", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCutPCM("51210013", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCutPCM("52310013", "00200009247602008250404000", "0652501500000000"); // 20-30%
    cuts.AddCutPCM("53410013", "00200009247602008250404000", "0652501500000000"); // 30-40%
  } else if (trainConfig == 267){ // LHC15o, kINT7,  cent. from track mult., user defined header
    cuts.AddCutPCM("60110023", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCutPCM("61210023", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCutPCM("51210023", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCutPCM("52310023", "00200009247602008250404000", "0652501500000000"); // 20-30%
    cuts.AddCutPCM("53410023", "00200009247602008250404000", "0652501500000000"); // 30-40%
  } else if (trainConfig == 268){ // LHC15o, kINT7,  cent. from track mult., reject added particles
    cuts.AddCutPCM("54610013", "00200009247602008250404000", "0652501500000000"); // 40-60%
    cuts.AddCutPCM("56810013", "00200009247602008250404000", "0652501500000000"); // 60-80%
    cuts.AddCutPCM("58910013", "00200009247602008250404000", "0652501500000000"); // 80-90%
    cuts.AddCutPCM("58010013", "00200009247602008250404000", "0652501500000000"); // 80-100%
    cuts.AddCutPCM("50910013", "00200009247602008250404000", "0652501500000000"); //  0-90%
  } else if (trainConfig == 269){ // LHC15o, kINT7,  cent. from track mult., user defined header
    cuts.AddCutPCM("54610023", "00200009247602008250404000", "0652501500000000"); // 40-60%
    cuts.AddCutPCM("56810023", "00200009247602008250404000", "0652501500000000"); // 60-80%
    cuts.AddCutPCM("58910023", "00200009247602008250404000", "0652501500000000"); // 80-90%
    cuts.AddCutPCM("58010023", "00200009247602008250404000", "0652501500000000"); // 80-100%
    cuts.AddCutPCM("50910023", "00200009247602008250404000", "0652501500000000"); //  0-90%
  } else if (trainConfig == 270){ // LHC15o, kINT7,  cent. from V0M, reject added particles
    cuts.AddCutPCM("30110013", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCutPCM("31210013", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCutPCM("11210013", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCutPCM("12310013", "00200009247602008250404000", "0652501500000000"); // 20-30%
    cuts.AddCutPCM("13410013", "00200009247602008250404000", "0652501500000000"); // 30-40%
  } else if (trainConfig == 271){ // LHC15o, kINT7,  cent. from V0M, user defined header
    cuts.AddCutPCM("30110023", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCutPCM("31210023", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCutPCM("11210023", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCutPCM("12310023", "00200009247602008250404000", "0652501500000000"); // 20-30%
    cuts.AddCutPCM("13410023", "00200009247602008250404000", "0652501500000000"); // 30-40%
  } else if (trainConfig == 272){ // LHC15o, kINT7,  cent. from V0M, reject added particles
    cuts.AddCutPCM("14610013", "00200009247602008250404000", "0652501500000000"); // 40-60%
    cuts.AddCutPCM("16810013", "00200009247602008250404000", "0652501500000000"); // 60-80%
    cuts.AddCutPCM("18910013", "00200009247602008250404000", "0652501500000000"); // 80-90%
    cuts.AddCutPCM("18010013", "00200009247602008250404000", "0652501500000000"); // 80-100%
    cuts.AddCutPCM("10910013", "00200009247602008250404000", "0652501500000000"); //  0-90%
  } else if (trainConfig == 273){ // LHC15o, kINT7,  cent. from V0M, user defined header
    cuts.AddCutPCM("14610023", "00200009247602008250404000", "0652501500000000"); // 40-60%
    cuts.AddCutPCM("16810023", "00200009247602008250404000", "0652501500000000"); // 60-80%
    cuts.AddCutPCM("18910023", "00200009247602008250404000", "0652501500000000"); // 80-90%
    cuts.AddCutPCM("18010023", "00200009247602008250404000", "0652501500000000"); // 80-100%
    cuts.AddCutPCM("10910023", "00200009247602008250404000", "0652501500000000"); //  0-90%
  } else if (trainConfig == 274){ // LHC15o, kINT7,  cent. from track mult., reject added particles, with PU cut
    cuts.AddCutPCM("60110613", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCutPCM("61210613", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCutPCM("51210613", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCutPCM("52310613", "00200009247602008250404000", "0652501500000000"); // 20-30%
    cuts.AddCutPCM("53410613", "00200009247602008250404000", "0652501500000000"); // 30-40%
  } else if (trainConfig == 275){ // LHC15o, kINT7,  cent. from track mult., user defined header, with PU cut
    cuts.AddCutPCM("60110623", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCutPCM("61210623", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCutPCM("51210623", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCutPCM("52310623", "00200009247602008250404000", "0652501500000000"); // 20-30%
    cuts.AddCutPCM("53410623", "00200009247602008250404000", "0652501500000000"); // 30-40%
  } else if (trainConfig == 276){ // LHC15o, kINT7,  cent. from track mult., reject added particles, with PU cut
    cuts.AddCutPCM("54610613", "00200009247602008250404000", "0652501500000000"); // 40-60%
    cuts.AddCutPCM("56810613", "00200009247602008250404000", "0652501500000000"); // 60-80%
    cuts.AddCutPCM("58910613", "00200009247602008250404000", "0652501500000000"); // 80-90%
    cuts.AddCutPCM("58010613", "00200009247602008250404000", "0652501500000000"); // 80-100%
    cuts.AddCutPCM("50910613", "00200009247602008250404000", "0652501500000000"); //  0-90%
  } else if (trainConfig == 277){ // LHC15o, kINT7,  cent. from track mult., user defined header, with PU cut
    cuts.AddCutPCM("54610623", "00200009247602008250404000", "0652501500000000"); // 40-60%
    cuts.AddCutPCM("56810623", "00200009247602008250404000", "0652501500000000"); // 60-80%
    cuts.AddCutPCM("58910623", "00200009247602008250404000", "0652501500000000"); // 80-90%
    cuts.AddCutPCM("58010623", "00200009247602008250404000", "0652501500000000"); // 80-100%
    cuts.AddCutPCM("50910623", "00200009247602008250404000", "0652501500000000"); //  0-90%
  } else if (trainConfig == 278){ // LHC15o, kINT7,  cent. from V0M, reject added particles, with PU cut
    cuts.AddCutPCM("30110613", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCutPCM("31210613", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCutPCM("11210613", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCutPCM("12310613", "00200009247602008250404000", "0652501500000000"); // 20-30%
    cuts.AddCutPCM("13410613", "00200009247602008250404000", "0652501500000000"); // 30-40%
  } else if (trainConfig == 279){ // LHC15o, kINT7,  cent. from V0M, user defined header, with PU cut
    cuts.AddCutPCM("30110623", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCutPCM("31210623", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCutPCM("11210623", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCutPCM("12310623", "00200009247602008250404000", "0652501500000000"); // 20-30%
    cuts.AddCutPCM("13410623", "00200009247602008250404000", "0652501500000000"); // 30-40%
  } else if (trainConfig == 280){ // LHC15o, kINT7,  cent. from V0M, reject added particles, with PU cut
    cuts.AddCutPCM("14610613", "00200009247602008250404000", "0652501500000000"); // 40-60%
    cuts.AddCutPCM("16810613", "00200009247602008250404000", "0652501500000000"); // 60-80%
    cuts.AddCutPCM("18910613", "00200009247602008250404000", "0652501500000000"); // 80-90%
    cuts.AddCutPCM("18010613", "00200009247602008250404000", "0652501500000000"); // 80-100%
    cuts.AddCutPCM("10910613", "00200009247602008250404000", "0652501500000000"); //  0-90%
  } else if (trainConfig == 281){ // LHC15o, kINT7,  cent. from V0M, user defined header, with PU cut
    cuts.AddCutPCM("14610623", "00200009247602008250404000", "0652501500000000"); // 40-60%
    cuts.AddCutPCM("16810623", "00200009247602008250404000", "0652501500000000"); // 60-80%
    cuts.AddCutPCM("18910623", "00200009247602008250404000", "0652501500000000"); // 80-90%
    cuts.AddCutPCM("18010623", "00200009247602008250404000", "0652501500000000"); // 80-100%
    cuts.AddCutPCM("10910623", "00200009247602008250404000", "0652501500000000"); //  0-90%
  } else if (trainConfig == 282){ // LHC15o, kINT7, test Pileup cuts
    cuts.AddCutPCM("10910013", "00200009247602008250404000", "0652501500000000"); // no pileup rejection
    cuts.AddCutPCM("10910113", "00200009247602008250404000", "0652501500000000"); // pileup rejection using SPD (strict cut)
    cuts.AddCutPCM("10910613", "00200009247602008250404000", "0652501500000000"); // pileup rejection using SPD (strict cut) and V0+TPC
    cuts.AddCutPCM("10910a13", "00200009247602008250404000", "0652501500000000"); // pileup rejection using SPD (open cut) and V0+TPC
    cuts.AddCutPCM("10910b13", "00200009247602008250404000", "0652501500000000"); // pileup rejection using V0+TPC
  } else if (trainConfig == 283){ // LHC15o, kINT7, cent from V0M, pileup rejection using V0+TPC, reject added particles
    cuts.AddCutPCM("30110a13", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCutPCM("31210a13", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCutPCM("11210a13", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCutPCM("12310a13", "00200009247602008250404000", "0652501500000000"); // 20-30%
    cuts.AddCutPCM("13410a13", "00200009247602008250404000", "0652501500000000"); // 30-40%
  } else if (trainConfig == 284){ // LHC15o, kINT7, cent from V0M, pileup rejection using V0+TPC, user defined header
    cuts.AddCutPCM("30110a23", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCutPCM("31210a23", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCutPCM("11210a23", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCutPCM("12310a23", "00200009247602008250404000", "0652501500000000"); // 20-30%
    cuts.AddCutPCM("13410a23", "00200009247602008250404000", "0652501500000000"); // 30-40%
  } else if (trainConfig == 285){ // LHC15o, kINT7, cent from V0M, pileup rejection using V0+TPC, reject added particles
    cuts.AddCutPCM("14610a13", "00200009247602008250404000", "0652501500000000"); // 40-60%
    cuts.AddCutPCM("16810a13", "00200009247602008250404000", "0652501500000000"); // 60-80%
    cuts.AddCutPCM("18910a13", "00200009247602008250404000", "0652501500000000"); // 80-90%
    cuts.AddCutPCM("10110a13", "00200009247602008250404000", "0652501500000000"); // 0-10%
    cuts.AddCutPCM("10910a13", "00200009247602008250404000", "0652501500000000"); // 0-90%
  } else if (trainConfig == 286){ // LHC15o, kINT7, cent from V0M, pileup rejection using V0+TPC, user defined header
    cuts.AddCutPCM("14610a23", "00200009247602008250404000", "0652501500000000"); // 40-60%
    cuts.AddCutPCM("16810a23", "00200009247602008250404000", "0652501500000000"); // 60-80%
    cuts.AddCutPCM("18910a23", "00200009247602008250404000", "0652501500000000"); // 80-90%
    cuts.AddCutPCM("10110a23", "00200009247602008250404000", "0652501500000000"); // 0-10%
    cuts.AddCutPCM("10910a23", "00200009247602008250404000", "0652501500000000"); // 0-90%
  } else if (trainConfig == 287){ // for test for inv mass cut
    cuts.AddCutPCM("10910013", "00200009247602008250404000", "0652501500000000"); // analysis cuts
    cuts.AddCutPCM("10910013", "00200018357602002140000000", "0652501500000000"); // open for systematic checks
  } else if (trainConfig == 288){ // LHC15o, kINT7, cent from V0M, pileup rejection using V0+TPC, reject added particles
    cuts.AddCutPCM("30110a13", "00200009247602008250404000", "0152501500000000"); //  0-5%
    cuts.AddCutPCM("12410a13", "00200009247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("14610a13", "00200009247602008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("13410a13", "00200009247602008250404000", "0152501500000000"); // 30-40%
    cuts.AddCutPCM("16910a13", "00200009247602008250404000", "0152501500000000"); // 60-90%
  } else if (trainConfig == 289){ // LHC15o, kINT7, cent from V0M, pileup rejection using V0+TPC, reject added particles
    cuts.AddCutPCM("31210a13", "00200009247602008250404000", "0152501500000000"); //  5-10%
    cuts.AddCutPCM("12510a13", "00200009247602008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("15710a13", "00200009247602008250404000", "0152501500000000"); // 50-70%
    cuts.AddCutPCM("16810a13", "00200009247602008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("10910a13", "00200009247602008250404000", "0152501500000000"); // 0-90%
  } else if (trainConfig == 290){ // LHC15o, kINT7, cent from V0M, pileup rejection using V0+TPC, reject added particles
    cuts.AddCutPCM("10110a13", "00200009247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200009247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12310a13", "00200009247602008250404000", "0152501500000000"); // 20-30%
    cuts.AddCutPCM("17910a13", "00200009247602008250404000", "0152501500000000"); // 70-90%
    cuts.AddCutPCM("18910a13", "00200009247602008250404000", "0152501500000000"); // 80-90%
  } else if (trainConfig == 291){ // LHC15o, kINT7, cent from V0M, pileup rejection using V0+TPC, user defined header
    cuts.AddCutPCM("30110a23", "00200009247602008250404000", "0152501500000000"); //  0-5%
    cuts.AddCutPCM("12410a23", "00200009247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("14610a23", "00200009247602008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("13410a23", "00200009247602008250404000", "0152501500000000"); // 30-40%
    cuts.AddCutPCM("16910a23", "00200009247602008250404000", "0152501500000000"); // 60-90%
  } else if (trainConfig == 292){ // LHC15o, kINT7, cent from V0M, pileup rejection using V0+TPC, user defined header
    cuts.AddCutPCM("31210a23", "00200009247602008250404000", "0152501500000000"); //  5-10%
    cuts.AddCutPCM("12510a23", "00200009247602008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("15710a23", "00200009247602008250404000", "0152501500000000"); // 50-70%
    cuts.AddCutPCM("16810a23", "00200009247602008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("10910a23", "00200009247602008250404000", "0152501500000000"); // 0-90%
  } else if (trainConfig == 293){ // LHC15o, kINT7, cent from V0M, pileup rejection using V0+TPC, user defined header
    cuts.AddCutPCM("10110a23", "00200009247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200009247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12310a23", "00200009247602008250404000", "0152501500000000"); // 20-30%
    cuts.AddCutPCM("17910a23", "00200009247602008250404000", "0152501500000000"); // 70-90%
    cuts.AddCutPCM("18910a23", "00200009247602008250404000", "0152501500000000"); // 80-90%
    // LHC15o standard cuts, LHC15o, kINT7, cent from V0M, pileup rejection using V0+TPC
  } else if (trainConfig == 294){ // without added particles
    cuts.AddCutPCM("10110a13", "00200009247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200009247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200009247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200009247602008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 295){
    cuts.AddCutPCM("14610a13", "00200009247602008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200009247602008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200009247602008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200009247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a13", "00200009247602008250404000", "0152501500000000"); // 0-20%
  } else if (trainConfig == 296){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200009247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200009247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200009247602008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 297){
    cuts.AddCutPCM("14610a23", "00200009247602008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200009247602008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200009247602008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200009247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a23", "00200009247602008250404000", "0152501500000000"); // 0-20%

  } else if (trainConfig == 298){
    cuts.AddCutPCM("10110013", "00200009247602008250404000", "0152501500000000");
    cuts.AddCutPCM("11310013", "00200009247602008250404000", "0152501500000000");
    cuts.AddCutPCM("13510013", "00200009247602008250404000", "0152501500000000");
    cuts.AddCutPCM("15910013", "00200009247602008250404000", "0152501500000000");
    cuts.AddCutPCM("10910013", "00200009247602008250404000", "0152501500000000");
  } else if (trainConfig == 299){
    cuts.AddCutPCM("10110113", "00200009247602008250404000", "0152501500000000");
    cuts.AddCutPCM("11310113", "00200009247602008250404000", "0152501500000000");
    cuts.AddCutPCM("13510113", "00200009247602008250404000", "0152501500000000");
    cuts.AddCutPCM("15910113", "00200009247602008250404000", "0152501500000000");
    cuts.AddCutPCM("10910113", "00200009247602008250404000", "0152501500000000");

  } else  if (trainConfig == 300){ // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCutPCM("60100013", "03200009300002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "03200009300002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "03200009300002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "03200009300002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "03200009300002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 301) {  // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCutPCM("52400013", "03200009300002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "03200009300002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "03200009300002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "03200009300002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "03200009300002003220000000", "0152306500900000"); // 40-90%
  } else  if (trainConfig == 302){ // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCutPCM("60100013", "03200009400002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "03200009400002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "03200009400002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "03200009400002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "03200009400002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 303) {  // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCutPCM("52400013", "03200009400002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "03200009400002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "03200009400002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "03200009400002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "03200009400002003220000000", "0152306500900000"); // 40-90%
  } else  if (trainConfig == 304){ // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCutPCM("60100013", "03200009500002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "03200009500002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "03200009500002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "03200009500002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "03200009500002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 305) {  // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCutPCM("52400013", "03200009500002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "03200009500002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "03200009500002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "03200009500002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "03200009500002003220000000", "0152306500900000"); // 40-90%
  } else  if (trainConfig == 306){ // LHC10h standard, eta 0.9, y = 0.6
    cuts.AddCutPCM("60100013", "00200009300002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009300002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009300002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "00200009300002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "00200009300002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 307) {  // LHC10h standard, eta 0.9, y = 0.6
    cuts.AddCutPCM("52400013", "00200009300002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "00200009300002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "00200009300002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "00200009300002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "00200009300002003220000000", "0152306500900000"); // 40-90%
  } else  if (trainConfig == 308){ // LHC10h standard, eta 0.9, y = 0.6
    cuts.AddCutPCM("60100013", "00200009400002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009400002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009400002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "00200009400002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "00200009400002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 309) { // LHC10h standard, eta 0.9, y = 0.6
    cuts.AddCutPCM("52400013", "00200009400002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "00200009400002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "00200009400002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "00200009400002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "00200009400002003220000000", "0152306500900000"); // 40-90%
  } else  if (trainConfig == 310){ // LHC10h standard, eta 0.9, y = 0.6
    cuts.AddCutPCM("60100013", "00200009500002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCutPCM("61200013", "00200009500002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCutPCM("50100013", "00200009500002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCutPCM("51200013", "00200009500002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCutPCM("50200013", "00200009500002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 311) {  // LHC10h standard, eta 0.9, y = 0.6
    cuts.AddCutPCM("52400013", "00200009500002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCutPCM("54600013", "00200009500002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCutPCM("56800013", "00200009500002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCutPCM("54800013", "00200009500002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCutPCM("54900013", "00200009500002003220000000", "0152306500900000"); // 40-90%
  } else  if (trainConfig == 312){ // LHC10h standard direct photon flow cuts
    cuts.AddCutPCM("50200013", "00200009307000008250400000", "0152304500900000");
    cuts.AddCutPCM("52400013", "00200009307000008250400000", "0152304500900000");
    cuts.AddCutPCM("54800013", "00200009307000008250400000", "0152304500900000");

    // added sign. gamma only rejection for dir photons ------------------------
  } else if ( trainConfig == 313){ // dir photons (tighter dEdx+TOF)
    cuts.AddCutPCM("50100023", "00200009247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 314){ //
    cuts.AddCutPCM("50100023", "00216609247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 315){ // dir photons (tighter dEdx+TOF)
    cuts.AddCutPCM("52400023", "00200009247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("52500023", "00200009247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 316){ //
    cuts.AddCutPCM("52400023", "00216609247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("52500023", "00216609247002008850404000", "0152501500000000"); //

  } else if ( trainConfig == 317){ // phi variation 1
    cuts.AddCutPCM("50100023", "00215509247002008850404000", "0152501500000000"); // 2 - 4 rad
  } else if ( trainConfig == 318){ // phi variation 2
    cuts.AddCutPCM("50100023", "00217709247002008850404000", "0152501500000000"); // 2.6 - 3.4 rad
  } else if ( trainConfig == 319){ // phi variation 1
    cuts.AddCutPCM("52400023", "00215509247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("52500023", "00215509247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 320){ // phi variation 2
    cuts.AddCutPCM("52400023", "00217709247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("52500023", "00217709247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 321){ // single pt variation
    cuts.AddCutPCM("50100023", "00200049247002008850404000", "0152501500000000"); // 0.075
    cuts.AddCutPCM("50100023", "00200019247002008850404000", "0152501500000000"); // 0.1
  } else if ( trainConfig == 322){ //  single pt variation - phi
    cuts.AddCutPCM("50100023", "00216649247002008850404000", "0152501500000000"); // 0.075
    cuts.AddCutPCM("50100023", "00216619247002008850404000", "0152501500000000"); // 0.1
  } else if ( trainConfig == 323){ // single pt variation
    cuts.AddCutPCM("52400023", "00200049247002008850404000", "0152501500000000"); // 0.075
    cuts.AddCutPCM("52400023", "00200019247002008850404000", "0152501500000000"); // 0.1
    cuts.AddCutPCM("52500023", "00200049247002008850404000", "0152501500000000"); // 0.075
    cuts.AddCutPCM("52500023", "00200019247002008850404000", "0152501500000000"); // 0.1
  } else if ( trainConfig == 324){ //  single pt variation - phi
    cuts.AddCutPCM("52400023", "00216649247002008850404000", "0152501500000000"); // 0.075
    cuts.AddCutPCM("52400023", "00216619247002008850404000", "0152501500000000"); // 0.1
    cuts.AddCutPCM("52500023", "00216649247002008850404000", "0152501500000000"); // 0.075
    cuts.AddCutPCM("52500023", "00216619247002008850404000", "0152501500000000"); // 0.1
  } else if ( trainConfig == 325){ // TPC clusters
    cuts.AddCutPCM("50100023", "00200008247002008850404000", "0152501500000000"); // 0.35
    cuts.AddCutPCM("50100023", "00200006247002008850404000", "0152501500000000"); // 0.7
  } else if ( trainConfig == 326){ // TPC clusters - phi
    cuts.AddCutPCM("50100023", "00216608247002008850404000", "0152501500000000"); // 0.35
    cuts.AddCutPCM("50100023", "00216606247002008850404000", "0152501500000000"); // 0.7
  } else if ( trainConfig == 327){ // TPC clusters
    cuts.AddCutPCM("52400023", "00200008247002008850404000", "0152501500000000"); // 0.35
    cuts.AddCutPCM("52400023", "00200006247002008850404000", "0152501500000000"); // 0.7
    cuts.AddCutPCM("52500023", "00200008247002008850404000", "0152501500000000"); // 0.35
    cuts.AddCutPCM("52500023", "00200006247002008850404000", "0152501500000000"); // 0.7
  } else if ( trainConfig == 328){ // TPC clusters - phi
    cuts.AddCutPCM("52400023", "00216608247002008850404000", "0152501500000000"); // 0.35
    cuts.AddCutPCM("52400023", "00216606247002008850404000", "0152501500000000"); // 0.7
    cuts.AddCutPCM("52500023", "00216608247002008850404000", "0152501500000000"); // 0.35
    cuts.AddCutPCM("52500023", "00216606247002008850404000", "0152501500000000"); // 0.7
  } else if ( trainConfig == 329){ // edEdx
    cuts.AddCutPCM("50100023", "00200009347002008850404000", "0152501500000000"); // -4, 5
    cuts.AddCutPCM("50100023", "00200009947002008850404000", "0152501500000000"); // -2.5, 5
  } else if ( trainConfig == 330){ // edEdx - phi
    cuts.AddCutPCM("50100023", "00216609347002008850404000", "0152501500000000"); // -4, 5
    cuts.AddCutPCM("50100023", "00216609947002008850404000", "0152501500000000"); // -2.5, 5
  } else if ( trainConfig == 331){ // edEdx
    cuts.AddCutPCM("52400023", "00200009347002008850404000", "0152501500000000"); // -4, 5
    cuts.AddCutPCM("52400023", "00200009947002008850404000", "0152501500000000"); // -2.5, 5
    cuts.AddCutPCM("52500023", "00200009347002008850404000", "0152501500000000"); // -4, 5
    cuts.AddCutPCM("52500023", "00200009947002008850404000", "0152501500000000"); // -2.5, 5
  } else if ( trainConfig == 332){ // edEdx - phi
    cuts.AddCutPCM("52400023", "00216609347002008850404000", "0152501500000000"); // -4, 5
    cuts.AddCutPCM("52400023", "00216609947002008850404000", "0152501500000000"); // -2.5, 5
    cuts.AddCutPCM("52500023", "00216609347002008850404000", "0152501500000000"); // -4, 5
    cuts.AddCutPCM("52500023", "00216609947002008850404000", "0152501500000000"); // -2.5, 5
  } else if ( trainConfig == 333){ // pdEdx
    cuts.AddCutPCM("50100023", "00200009240002008850404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCutPCM("50100023", "00200009247602008850404000", "0152501500000000"); // 3&1sigma
  } else if ( trainConfig == 334){ // pdEdx - phi
    cuts.AddCutPCM("50100023", "00216609240002008850404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCutPCM("50100023", "00216609247602008850404000", "0152501500000000"); // 3&1sigma
  } else if ( trainConfig == 335){ // pdEdx
    cuts.AddCutPCM("52400023", "00200009240002008850404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCutPCM("52400023", "00200009247602008850404000", "0152501500000000"); // 3&1sigma
    cuts.AddCutPCM("52500023", "00200009240002008850404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCutPCM("52500023", "00200009247602008850404000", "0152501500000000"); // 3&1sigma
  } else if ( trainConfig == 336){ // pdEdx - phi
    cuts.AddCutPCM("52400023", "00216609240002008850404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCutPCM("52400023", "00216609247602008850404000", "0152501500000000"); // 3&1sigma
    cuts.AddCutPCM("52500023", "00216609240002008850404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCutPCM("52500023", "00216609247602008850404000", "0152501500000000"); // 3&1sigma
  } else if ( trainConfig == 337){ // pdEdx & TOF
    cuts.AddCutPCM("50100023", "00200009237002008850404000", "0152501500000000"); // 2.5sigma
    cuts.AddCutPCM("50100023", "00200009247005008850404000", "0152501500000000"); // TOF -3, 3
  } else if ( trainConfig == 338){ // pdEdx & TOF - phi
    cuts.AddCutPCM("50100023", "00216609237002008850404000", "0152501500000000"); // 2.5sigma
    cuts.AddCutPCM("50100023", "00216609247005008850404000", "0152501500000000"); // TOF -3, 3
  } else if ( trainConfig == 339){ // pdEdx & TOF
    cuts.AddCutPCM("52400023", "00200009237002008850404000", "0152501500000000"); // 2.5sigma
    cuts.AddCutPCM("52400023", "00200009247005008850404000", "0152501500000000"); // TOF -3, 3
    cuts.AddCutPCM("52500023", "00200009237002008850404000", "0152501500000000"); // 2.5sigma
    cuts.AddCutPCM("52500023", "00200009247005008850404000", "0152501500000000"); // TOF -3, 3
  } else if ( trainConfig == 340){ // pdEdx & TOF - phi
    cuts.AddCutPCM("52400023", "00216609237002008850404000", "0152501500000000"); // 2.5sigma
    cuts.AddCutPCM("52400023", "00216609247005008850404000", "0152501500000000"); // TOF -3, 3
    cuts.AddCutPCM("52500023", "00216609237002008850404000", "0152501500000000"); // 2.5sigma
    cuts.AddCutPCM("52500023", "00216609247005008850404000", "0152501500000000"); // TOF -3, 3
  } else if ( trainConfig == 341){ // qT
    cuts.AddCutPCM("50100023", "00200009247002002850404000", "0152501500000000"); // 0.06 2D
    cuts.AddCutPCM("50100023", "00200009247002009850404000", "0152501500000000"); // 0.03 2D
  } else if ( trainConfig == 342){ // qT - phi
    cuts.AddCutPCM("50100023", "00216609247002002850404000", "0152501500000000"); // 0.06 2D
    cuts.AddCutPCM("50100023", "00216609247002009850404000", "0152501500000000"); // 0.03 2D
  } else if ( trainConfig == 343){ // qT
    cuts.AddCutPCM("52400023", "00200009247002002850404000", "0152501500000000"); // 0.06 2D
    cuts.AddCutPCM("52400023", "00200009247002009850404000", "0152501500000000"); // 0.03 2D
    cuts.AddCutPCM("52500023", "00200009247002002850404000", "0152501500000000"); // 0.06 2D
    cuts.AddCutPCM("52500023", "00200009247002009850404000", "0152501500000000"); // 0.03 2D
  } else if ( trainConfig == 344){ // qT - phi
    cuts.AddCutPCM("52400023", "00216609247002002850404000", "0152501500000000"); // 0.06 2D
    cuts.AddCutPCM("52400023", "00216609247002009850404000", "0152501500000000"); // 0.03 2D
    cuts.AddCutPCM("52500023", "00216609247002002850404000", "0152501500000000"); // 0.06 2D
    cuts.AddCutPCM("52500023", "00216609247002009850404000", "0152501500000000"); // 0.03 2D
  } else if ( trainConfig == 345){ // chi2 + psi pair
    cuts.AddCutPCM("50100023", "00200009247002008750404000", "0152501500000000"); // 10+0.1 2D
    cuts.AddCutPCM("50100023", "00200009247002008a50404000", "0152501500000000"); // 25+0.1 2D
  } else if ( trainConfig == 346){ // chi2 + psi pair - phi
    cuts.AddCutPCM("50100023", "00216609247002008750404000", "0152501500000000"); // 10+0.1 2D
    cuts.AddCutPCM("50100023", "00216609247002008a50404000", "0152501500000000"); // 25+0.1 2D
  } else if ( trainConfig == 347){ // chi2 + psi pair
    cuts.AddCutPCM("52400023", "00200009247002008750404000", "0152501500000000"); // 10+0.1 2D
    cuts.AddCutPCM("52400023", "00200009247002008a50404000", "0152501500000000"); // 25+0.1 2D
    cuts.AddCutPCM("52500023", "00200009247002008750404000", "0152501500000000"); // 10+0.1 2D
    cuts.AddCutPCM("52500023", "00200009247002008a50404000", "0152501500000000"); // 25+0.1 2D
  } else if ( trainConfig == 348){ // chi2 + psi pair - phi
    cuts.AddCutPCM("52400023", "00216609247002008750404000", "0152501500000000"); // 10+0.1 2D
    cuts.AddCutPCM("52400023", "00216609247002008a50404000", "0152501500000000"); // 25+0.1 2D
    cuts.AddCutPCM("52500023", "00216609247002008750404000", "0152501500000000"); // 10+0.1 2D
    cuts.AddCutPCM("52500023", "00216609247002008a50404000", "0152501500000000"); // 25+0.1 2D
  } else if ( trainConfig == 349){ // chi2 + psi pair
    cuts.AddCutPCM("50100023", "00200009247002008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCutPCM("50100023", "00200009247002008870404000", "0152501500000000"); // 20+0.07 2D
  } else if ( trainConfig == 350){ // chi2 + psi pair - phi
    cuts.AddCutPCM("50100023", "00216609247002008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCutPCM("50100023", "00216609247002008870404000", "0152501500000000"); // 20+0.07 2D
  } else if ( trainConfig == 351){ // chi2 + psi pair
    cuts.AddCutPCM("52400023", "00200009247002008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCutPCM("52400023", "00200009247002008870404000", "0152501500000000"); // 20+0.07 2D
    cuts.AddCutPCM("52500023", "00200009247002008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCutPCM("52500023", "00200009247002008870404000", "0152501500000000"); // 20+0.07 2D
  } else if ( trainConfig == 352){ // chi2 + psi pair - phi
    cuts.AddCutPCM("52400023", "00216609247002008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCutPCM("52400023", "00216609247002008870404000", "0152501500000000"); // 20+0.07 2D
    cuts.AddCutPCM("52500023", "00216609247002008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCutPCM("52500023", "00216609247002008870404000", "0152501500000000"); // 20+0.07 2D
  } else if ( trainConfig == 353){ // asym 1D & cosPA
    cuts.AddCutPCM("50100023", "00200009247002008857404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCutPCM("50100023", "00200009247002008850004000", "0152501500000000"); // cosPA off
  } else if ( trainConfig == 354){ // asym 1D & cosPA - phi
    cuts.AddCutPCM("50100023", "00216609247002008857404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCutPCM("50100023", "00216609247002008850004000", "0152501500000000"); // cosPA off
  } else if ( trainConfig == 355){ // asym 1D & cosPA
    cuts.AddCutPCM("52400023", "00200009247002008857404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCutPCM("52400023", "00200009247002008850004000", "0152501500000000"); // cosPA off
    cuts.AddCutPCM("52500023", "00200009247002008857404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCutPCM("52500023", "00200009247002008850004000", "0152501500000000"); // cosPA off
  } else if ( trainConfig == 356){ // asym 1D & cosPA - phi
    cuts.AddCutPCM("52400023", "00216609247002008857404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCutPCM("52400023", "00216609247002008850004000", "0152501500000000"); // cosPA off
    cuts.AddCutPCM("52500023", "00216609247002008857404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCutPCM("52500023", "00216609247002008850004000", "0152501500000000"); // cosPA off

  } else if ( trainConfig == 357){ // meson STD cut with pileup removal
    cuts.AddCutPCM("60100713", "00200009247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("61200713", "00200009247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("50100713", "00200009247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("52400713", "00200009247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("52500713", "00200009247602008850404000", "0652501500000000"); //
  } else if ( trainConfig == 358){ //  with phi cut
    cuts.AddCutPCM("60100713", "00216609247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("61200713", "00216609247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("50100713", "00216609247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("52400713", "00216609247602008850404000", "0652501500000000"); //
    cuts.AddCutPCM("52500713", "00216609247602008850404000", "0652501500000000"); //

  } else if ( trainConfig == 359){ // chi2 + psi pair
    cuts.AddCutPCM("50100013", "00200009247002008250404000", "0152501500000000"); // 30+0.1 2D
  } else if ( trainConfig == 360){ // chi2 + psi pair - phi
    cuts.AddCutPCM("50100013", "00216609247002008250404000", "0152501500000000"); // 30+0.1 2D
  } else if ( trainConfig == 361){ // chi2 + psi pair
//     cuts.AddCutPCM("52400013", "00200009247002008250404000", "0152501500000000"); // 30+0.1 2D
    cuts.AddCutPCM("52500013", "00200009247002008250404000", "0152501500000000"); // 30+0.1 2D
  } else if ( trainConfig == 362){ // chi2 + psi pair - phi
    cuts.AddCutPCM("52400013", "00216609247002008250404000", "0152501500000000"); // 30+0.1 2D
    cuts.AddCutPCM("52500013", "00216609247002008250404000", "0152501500000000"); // 30+0.1 2D

  } else if ( trainConfig == 363){ // chi2 + psi pair
    cuts.AddCutPCM("50100023", "00200009247002008250404000", "0152501500000000"); // 30+0.1 2D
  } else if ( trainConfig == 364){ // chi2 + psi pair - phi
    cuts.AddCutPCM("50100023", "00216609247002008250404000", "0152501500000000"); // 30+0.1 2D
  } else if ( trainConfig == 365){ // chi2 + psi pair
    cuts.AddCutPCM("52400023", "00200009247002008250404000", "0152501500000000"); // 30+0.1 2D
    cuts.AddCutPCM("52500023", "00200009247002008250404000", "0152501500000000"); // 30+0.1 2D
  } else if ( trainConfig == 366){ // chi2 + psi pair - phi
    cuts.AddCutPCM("52400023", "00216609247002008250404000", "0152501500000000"); // 30+0.1 2D
    cuts.AddCutPCM("52500023", "00216609247002008250404000", "0152501500000000"); // 30+0.1 2D

  } else if ( trainConfig == 367){ // dir photons finer cent
    cuts.AddCutPCM("52300013", "00200009247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("53400013", "00200009247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("53500013", "00200009247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 368){ // with phi cut
    cuts.AddCutPCM("52300013", "00216609247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("53400013", "00216609247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("53500013", "00216609247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 369){ //
    cuts.AddCutPCM("52300023", "00200009247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("53400023", "00200009247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("53500023", "00200009247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 370){ //
    cuts.AddCutPCM("52300023", "00216609247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("53400023", "00216609247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("53500023", "00216609247002008850404000", "0152501500000000"); //

  } else if ( trainConfig == 371){ // dir photons finer cent
    cuts.AddCutPCM("60100013", "00200009247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("61200013", "00200009247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("51200013", "00200009247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 372){ // with phi cut
    cuts.AddCutPCM("60100013", "00216609247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("61200013", "00216609247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("51200013", "00216609247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 373){ //
    cuts.AddCutPCM("60100023", "00200009247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("61200023", "00200009247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("51200023", "00200009247002008850404000", "0152501500000000"); //
  } else if ( trainConfig == 374){ //
    cuts.AddCutPCM("60100023", "00216609247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("61200023", "00216609247002008850404000", "0152501500000000"); //
    cuts.AddCutPCM("51200023", "00216609247002008850404000", "0152501500000000"); //

  //****************************************************************************************************
  // 5.44TeV Xe-Xe LHC17n
  //****************************************************************************************************
  } else if (trainConfig == 400){
    cuts.AddCutPCM("10810013","00200089f97300008250400000","0163103100000000"); // 0-80
    cuts.AddCutPCM("10810013","00200089f9730000iih0400000","0163103100000000"); // 0-80
  } else if (trainConfig == 401){
    cuts.AddCutPCM("10210013","00200089f97300008250400000","0163103100000000"); // 0-20
    cuts.AddCutPCM("12410013","00200089f97300008250400000","0163103100000000"); // 20-40
    cuts.AddCutPCM("10410013","00200089f97300008250400000","0163103100000000"); // 0-40
    cuts.AddCutPCM("14810013","00200089f97300008250400000","0163103100000000"); // 40-80
  } else if (trainConfig == 402){ // optimized psi pair & chi2
    cuts.AddCutPCM("10210013","00200089f97300008ih0400000","0163103100000000"); // 0-20
    cuts.AddCutPCM("12410013","00200089f97300008ih0400000","0163103100000000"); // 20-40
    cuts.AddCutPCM("10410013","00200089f97300008ih0400000","0163103100000000"); // 0-40
    cuts.AddCutPCM("14810013","00200089f97300008ih0400000","0163103100000000"); // 40-80
  } else if (trainConfig == 403){ // optimized qt vs pt
    cuts.AddCutPCM("10210013","00200089f9730000i250400000","0163103100000000"); // 0-20
    cuts.AddCutPCM("12410013","00200089f9730000i250400000","0163103100000000"); // 20-40
    cuts.AddCutPCM("10410013","00200089f9730000i250400000","0163103100000000"); // 0-40
    cuts.AddCutPCM("14810013","00200089f9730000i250400000","0163103100000000"); // 40-80
  } else if (trainConfig == 404){ // optimized psi pair & chi2 & qt vs pt
    cuts.AddCutPCM("10210013","00200089f9730000iih0400000","0163103100000000"); // 0-20
    cuts.AddCutPCM("12410013","00200089f9730000iih0400000","0163103100000000"); // 20-40
    cuts.AddCutPCM("10410013","00200089f9730000iih0400000","0163103100000000"); // 0-40
    cuts.AddCutPCM("14810013","00200089f9730000iih0400000","0163103100000000"); // 40-80

  //****************************************************************************************************
  // 5.02TeV Pb-Pb LHC15o
  //****************************************************************************************************
    // LHC15o systematics,  kINT7 trigger, MC centrality from V0M, pileup rejection using V0+TPC
    // standard: Chi2 < 30 (2)
    // Chi2 < 20 (8)
  } else if (trainConfig == 500){ //-------------_-----------------------------------------
    cuts.AddCutPCM("10110a13", "00200009247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200009247602008850404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200009247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200009247602008850404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 501){
    cuts.AddCutPCM("14610a13", "00200009247602008850404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200009247602008850404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200009247602008850404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200009247602008850404000", "0152501500000000"); // 5-10%
  } else if (trainConfig == 502){ // added particles
    cuts.AddCutPCM("10110a23", "00200009247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200009247602008850404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200009247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200009247602008850404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 503){ // added particles
    cuts.AddCutPCM("14610a23", "00200009247602008850404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200009247602008850404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200009247602008850404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200009247602008850404000", "0152501500000000"); // 5-10%
    // Chi2 < 50 (1)
  } else if (trainConfig == 504){ //-------------_-----------------------------------------
    cuts.AddCutPCM("10110a13", "00200009247602008150404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200009247602008150404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200009247602008150404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200009247602008150404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 505){
    cuts.AddCutPCM("14610a13", "00200009247602008150404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200009247602008150404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200009247602008150404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200009247602008150404000", "0152501500000000"); // 5-10%
  } else if (trainConfig == 506){ // added particles
    cuts.AddCutPCM("10110a23", "00200009247602008150404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200009247602008150404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200009247602008150404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200009247602008150404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 507){ // added particles
    cuts.AddCutPCM("14610a23", "00200009247602008150404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200009247602008150404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200009247602008150404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200009247602008150404000", "0152501500000000"); // 5-10%

    // standard: PsiPair < 0.1 (5)
    // PsiPair < 0.05 (6)
  } else if (trainConfig == 508){ //--------------_---------------------------------------
    cuts.AddCutPCM("10110a13", "00200009247602008260404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200009247602008260404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200009247602008260404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200009247602008260404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 509){
    cuts.AddCutPCM("14610a13", "00200009247602008260404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200009247602008260404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200009247602008260404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200009247602008260404000", "0152501500000000"); // 5-10%
  } else if (trainConfig == 510){ // added particles
    cuts.AddCutPCM("10110a23", "00200009247602008260404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200009247602008260404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200009247602008260404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200009247602008260404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 511){ // added particles
    cuts.AddCutPCM("14610a23", "00200009247602008260404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200009247602008260404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200009247602008260404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200009247602008260404000", "0152501500000000"); // 5-10%
    // PsiPair < 0.2 (8)
  } else if (trainConfig == 512){ //--------------_---------------------------------------
    cuts.AddCutPCM("10110a13", "00200009247602008280404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200009247602008280404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200009247602008280404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200009247602008280404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 513){
    cuts.AddCutPCM("14610a13", "00200009247602008280404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200009247602008280404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200009247602008280404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200009247602008280404000", "0152501500000000"); // 5-10%
  } else if (trainConfig == 514){ // added particles
    cuts.AddCutPCM("10110a23", "00200009247602008280404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200009247602008280404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200009247602008280404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200009247602008280404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 515){ // added particles
    cuts.AddCutPCM("14610a23", "00200009247602008280404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200009247602008280404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200009247602008280404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200009247602008280404000", "0152501500000000"); // 5-10%

    // standard: electron pT > 0.05 GeV and photon pT > 0.02 GeV (0)
    // electron pT > 0.075 (4)
  } else if (trainConfig == 516){ //--_---------------------------------------------------
    cuts.AddCutPCM("10110a13", "00200049247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200049247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200049247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200049247602008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 517){
    cuts.AddCutPCM("14610a13", "00200049247602008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200049247602008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200049247602008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200049247602008250404000", "0152501500000000"); // 5-10%
  } else if (trainConfig == 518){ // added particles
    cuts.AddCutPCM("10110a23", "00200049247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200049247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200049247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200049247602008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 519){ // added particles
    cuts.AddCutPCM("14610a23", "00200049247602008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200049247602008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200049247602008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200049247602008250404000", "0152501500000000"); // 5-10%
    // electron pT > 0.1 (1)
  } else if (trainConfig == 520){ //--_---------------------------------------------------
    cuts.AddCutPCM("10110a13", "00200019247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200019247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200019247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200019247602008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 521){
    cuts.AddCutPCM("14610a13", "00200019247602008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200019247602008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200019247602008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200019247602008250404000", "0152501500000000"); // 5-10%
  } else if (trainConfig == 522){ // added particles
    cuts.AddCutPCM("10110a23", "00200019247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200019247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200019247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200019247602008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 523){ // added particles
    cuts.AddCutPCM("14610a23", "00200019247602008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200019247602008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200019247602008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200019247602008250404000", "0152501500000000"); // 5-10%
    // photon pT > 0.1 (9)
  } else if (trainConfig == 524){ //--_---------------------------------------------------
    cuts.AddCutPCM("10110a13", "00200099247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200099247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200099247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200099247602008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 525){
    cuts.AddCutPCM("14610a13", "00200099247602008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200099247602008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200099247602008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200099247602008250404000", "0152501500000000"); // 5-10%
  } else if (trainConfig == 526){ // with added particles
    cuts.AddCutPCM("10110a23", "00200099247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200099247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200099247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200099247602008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 527){ // with added particles
    cuts.AddCutPCM("14610a23", "00200099247602008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200099247602008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200099247602008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200099247602008250404000", "0152501500000000"); // 5-10%

    // standard: TPC cluster > 60% (9)
    // TPC clusters > 70 % (6)
  } else if (trainConfig == 528){ //---_--------------------------------------------------
    cuts.AddCutPCM("10110a13", "00200006247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200006247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200006247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200006247602008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 529){
    cuts.AddCutPCM("14610a13", "00200006247602008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200006247602008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200006247602008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200006247602008250404000", "0152501500000000"); // 5-10%
  } else if (trainConfig == 530){ // with added particles
    cuts.AddCutPCM("10110a23", "00200006247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200006247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200006247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200006247602008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 531){ // with added particles
    cuts.AddCutPCM("14610a23", "00200006247602008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200006247602008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200006247602008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200006247602008250404000", "0152501500000000"); // 5-10%
    // TPC clusters > 35 % (8)
  } else if (trainConfig == 532){ //---_--------------------------------------------------
    cuts.AddCutPCM("10110a13", "00200008247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200008247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200008247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200008247602008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 533){
    cuts.AddCutPCM("14610a13", "00200008247602008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200008247602008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200008247602008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200008247602008250404000", "0152501500000000"); // 5-10%
  } else if (trainConfig == 534){ // with added particles
    cuts.AddCutPCM("10110a23", "00200008247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200008247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200008247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200008247602008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 535){ // with added particles
    cuts.AddCutPCM("14610a23", "00200008247602008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200008247602008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200008247602008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200008247602008250404000", "0152501500000000"); // 5-10%

    // standard: qT < 0.05 2D with alpha < 0.95 (8)
    // qT < 0.03 2D (9)
  } else if (trainConfig == 536){ //------------_-----------------------------------------
    cuts.AddCutPCM("10110a13", "00200009247600009250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00200009247600009250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 537){
    cuts.AddCutPCM("16810a13", "00200009247600009250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 538){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600009250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00200009247600009250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 539){ // with added particles
    cuts.AddCutPCM("16810a23", "00200009247600009250404000", "0152501500000000"); // 60-80%
    // qT < 0.06 2D (2)
  } else if (trainConfig == 540){ //------------_-----------------------------------------
    cuts.AddCutPCM("10110a13", "00200009247600002250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00200009247600002250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 541){
    cuts.AddCutPCM("16810a13", "00200009247600002250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 542){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600002250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00200009247600002250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 543){ // with added particles
    cuts.AddCutPCM("16810a23", "00200009247600002250404000", "0152501500000000"); // 60-80%
    // qT < 0.05 1D i.e. alpha cut off (3)
  } else if (trainConfig == 544){ //------------_-----------------------------------------
    cuts.AddCutPCM("10110a13", "00200009247600003250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00200009247600003250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 545){
    cuts.AddCutPCM("16810a13", "00200009247600003250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 546){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600003250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00200009247600003250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 547){ // with added particles
    cuts.AddCutPCM("16810a23", "00200009247600003250404000", "0152501500000000"); // 60-80%

    // standard TPC PID electron line -3,5 (2)
    // TPC PID -2.5,4 (6)
  } else if (trainConfig == 548){ //----_-------------------------------------------------
    cuts.AddCutPCM("10110a13", "00200009647600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00200009647600008250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 549){
    cuts.AddCutPCM("16810a13", "00200009647600008250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 550){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009647600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00200009647600008250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 551){  // with added particles
    cuts.AddCutPCM("16810a23", "00200009647600008250404000", "0152501500000000"); // 60-80%
    // TPC PID -4,5 (3)
  } else if (trainConfig == 552){ //----_-------------------------------------------------
    cuts.AddCutPCM("10110a13", "00200009347600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00200009347600008250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 553){
    cuts.AddCutPCM("16810a13", "00200009347600008250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 554){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009347600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00200009347600008250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 555){  // with added particles
    cuts.AddCutPCM("16810a23", "00200009347600008250404000", "0152501500000000"); // 60-80%

    // TPC PID standard 3/1 sigma above pion line between 0.4GeV - 2.0GeV
    // pion rejection range 0.3-3GeV (454)
  } else if (trainConfig == 556){ //-----___----------------------------------------------
    cuts.AddCutPCM("10110a13", "00200009245400008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00200009245400008250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 557){
    cuts.AddCutPCM("16810a13", "00200009245400008250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 558){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009245400008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00200009245400008250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 559){ // with added particles
    cuts.AddCutPCM("16810a23", "00200009245400008250404000", "0152501500000000"); // 60-80%
    // 2 sigmas above pion line (576)
  } else if (trainConfig == 560){ //-----___----------------------------------------------
    cuts.AddCutPCM("10110a13", "00200009257600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00200009257600008250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 561){
    cuts.AddCutPCM("16810a13", "00200009257600008250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 562){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009257600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00200009257600008250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 563){ // with added particles
    cuts.AddCutPCM("16810a23", "00200009257600008250404000", "0152501500000000"); // 60-80%

    // Chi2 varied in small steps for 0-10% centrality class
  } else if (trainConfig == 564){ //-------------_-----------------------------------------
    cuts.AddCutPCM("10110a13", "00200009247600008a50404000", "0152501500000000"); // 25
    cuts.AddCutPCM("10110a13", "00200009247600008950404000", "0152501500000000"); // 15
    cuts.AddCutPCM("10110a13", "00200009247600008750404000", "0152501500000000"); // 10
    cuts.AddCutPCM("10110a13", "00200009247600008650404000", "0152501500000000"); // 5
  } else if (trainConfig == 565){
    cuts.AddCutPCM("10110a13", "00200009247600008b50404000", "0152501500000000"); // 35
    cuts.AddCutPCM("10110a13", "00200009247600008c50404000", "0152501500000000"); // 40
    cuts.AddCutPCM("10110a13", "00200009247600008d50404000", "0152501500000000"); // 45
    cuts.AddCutPCM("10110a13", "00200009247600008e50404000", "0152501500000000"); // 55
  } else if (trainConfig == 566){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600008a50404000", "0152501500000000"); // 25
    cuts.AddCutPCM("10110a23", "00200009247600008950404000", "0152501500000000"); // 15
    cuts.AddCutPCM("10110a23", "00200009247600008750404000", "0152501500000000"); // 10
    cuts.AddCutPCM("10110a23", "00200009247600008650404000", "0152501500000000"); // 5
  } else if (trainConfig == 567){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600008b50404000", "0152501500000000"); // 35
    cuts.AddCutPCM("10110a23", "00200009247600008c50404000", "0152501500000000"); // 40
    cuts.AddCutPCM("10110a23", "00200009247600008d50404000", "0152501500000000"); // 45
    cuts.AddCutPCM("10110a23", "00200009247600008e50404000", "0152501500000000"); // 55

    // new standard cut without TOF cut
  } else if (trainConfig == 568){ //_________-____________
    cuts.AddCutPCM("10110a13", "00200009247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200009247600008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200009247600008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200009247600008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 569){
    cuts.AddCutPCM("14610a13", "00200009247600008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200009247600008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200009247600008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200009247600008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a13", "00200009247600008250404000", "0152501500000000"); // 0-20%
  } else if (trainConfig == 570){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200009247600008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200009247600008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200009247600008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 571){
    cuts.AddCutPCM("14610a23", "00200009247600008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200009247600008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200009247600008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200009247600008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a23", "00200009247600008250404000", "0152501500000000"); // 0-20%
    // standard with material budget weights
 } else if (trainConfig == 572){ // without added particles
    cuts.AddCutPCM("10110a13", "00200009247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200009247600008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200009247600008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200009247600008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 573){
    cuts.AddCutPCM("14610a13", "00200009247600008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200009247600008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200009247600008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200009247600008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a13", "00200009247600008250404000", "0152501500000000"); // 0-20%
  } else if (trainConfig == 574){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200009247600008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200009247600008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200009247600008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 575){
    cuts.AddCutPCM("14610a23", "00200009247600008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200009247600008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200009247600008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200009247600008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a23", "00200009247600008250404000", "0152501500000000"); // 0-20%

    // double counting cut variations. Standard:4
    // without cut (0)
  } else if (trainConfig == 576){ //__________________-___
    cuts.AddCutPCM("10110a13", "00200009247600008250400000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00200009247600008250400000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 577){
    cuts.AddCutPCM("16810a13", "00200009247600008250400000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 578){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600008250400000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00200009247600008250400000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 579){
    cuts.AddCutPCM("16810a23", "00200009247600008250400000", "0152501500000000"); // 60-80%
    // fminV0Dist = 2 (2)
  } else if (trainConfig == 580){ //__________________-___
    cuts.AddCutPCM("10110a13", "00200009247600008250402000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00200009247600008250402000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 581){
    cuts.AddCutPCM("16810a13", "00200009247600008250402000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 582){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600008250402000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00200009247600008250402000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 583){
    cuts.AddCutPCM("16810a23", "00200009247600008250402000", "0152501500000000"); // 60-80%
    // fminV0Dist = 3 (3)
  } else if (trainConfig == 584){ //__________________-___
    cuts.AddCutPCM("10110a13", "00200009247600008250403000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00200009247600008250403000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 585){
    cuts.AddCutPCM("16810a13", "00200009247600008250403000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 586){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600008250403000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00200009247600008250403000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 587){
    cuts.AddCutPCM("16810a23", "00200009247600008250403000", "0152501500000000"); // 60-80%

    // phi cut vaiations. Standard: no cut (000)
    // accept only photons in good phi regions (399)
  } else if (trainConfig == 588){//---___________________
    cuts.AddCutPCM("10110a13", "00239909247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00239909247600008250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 589){
    cuts.AddCutPCM("16810a13", "00239909247600008250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 590){ // with added particles
    cuts.AddCutPCM("10110a23", "00239909247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00239909247600008250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 591){
    cuts.AddCutPCM("16810a23", "00239909247600008250404000", "0152501500000000"); // 60-80%
    // accept only photons in phi regions with strong distortions (4aa)
  } else if (trainConfig == 592){//---___________________
    cuts.AddCutPCM("10110a13", "0024aa09247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "0024aa09247600008250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 593){
    cuts.AddCutPCM("16810a13", "0024aa09247600008250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 594){ // with added particles
    cuts.AddCutPCM("10110a23", "0024aa09247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "0024aa09247600008250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 595){
    cuts.AddCutPCM("16810a23", "0024aa09247600008250404000", "0152501500000000"); // 60-80%

    // meson alpha cut variations. Standard: 0.5 (1)
    // 0.75 (5)
  } else if (trainConfig == 596){//_________________________________-_________
    cuts.AddCutPCM("10110a13", "00200009247600008250404000", "0152505500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00200009247600008250404000", "0152505500000000"); // 20-40%
  } else if (trainConfig == 597){
    cuts.AddCutPCM("16810a13", "00200009247600008250404000", "0152505500000000"); // 60-80%
  } else if (trainConfig == 598){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600008250404000", "0152505500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00200009247600008250404000", "0152505500000000"); // 20-40%
  } else if (trainConfig == 599){
    cuts.AddCutPCM("16810a23", "00200009247600008250404000", "0152505500000000"); // 60-80%
    // 1.0  (3)
  } else if (trainConfig == 600){//_________________________________-_________
    cuts.AddCutPCM("10110a13", "00200009247600008250404000", "0152503500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00200009247600008250404000", "0152503500000000"); // 20-40%
  } else if (trainConfig == 601){
    cuts.AddCutPCM("16810a13", "00200009247600008250404000", "0152503500000000"); // 60-80%
  } else if (trainConfig == 602){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600008250404000", "0152503500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00200009247600008250404000", "0152503500000000"); // 20-40%
  } else if (trainConfig == 603){
    cuts.AddCutPCM("16810a23", "00200009247600008250404000", "0152503500000000"); // 60-80%

    // eta cut variation. Standard: 0.9 (0)  00200009247600008250404000
    // 0.8 (d)
  } else if (trainConfig == 604){
    cuts.AddCutPCM("10110a13", "0d200009247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "0d200009247600008250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 605){
    cuts.AddCutPCM("16810a13", "0d200009247600008250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 606){ // with added particles
    cuts.AddCutPCM("10110a23", "0d200009247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "0d200009247600008250404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 607){
    cuts.AddCutPCM("16810a23", "0d200009247600008250404000", "0152501500000000"); // 60-80%

    // Chi2-PsiPair simultaneous variation. Standard: PsiPair 0.1 (5) Chi2 30 (2)
    // PsiPair 0.05 2D (6) Chi2 20 (8)
  } else if (trainConfig == 608){ //-------------__---------------------------------------
    cuts.AddCutPCM("10110a13", "00200009247600008860404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00200009247600008860404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 609){
    cuts.AddCutPCM("16810a13", "00200009247600008860404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 610){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600008860404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00200009247600008860404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 611){
    cuts.AddCutPCM("16810a23", "00200009247600008860404000", "0152501500000000"); // 60-80%
    // PsiPair 0.2 2D (8) Chi2 40 (c)
  } else if (trainConfig == 612){ //-------------__---------------------------------------
    cuts.AddCutPCM("10110a13", "00200009247600008c80404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a13", "00200009247600008c80404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 613){
    cuts.AddCutPCM("16810a13", "00200009247600008c80404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 614){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600008c80404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("12410a23", "00200009247600008c80404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 615){
    cuts.AddCutPCM("16810a23", "00200009247600008c80404000", "0152501500000000"); // 60-80%

    // TPC Chi2 cut
    // standard: no cut (9)
    // < 4 (a)
  } else if (trainConfig == 616){ //___-__________________
    cuts.AddCutPCM("10110a13", "0020000a247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "0020000a247600008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "0020000a247600008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "0020000a247600008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 617){
    cuts.AddCutPCM("14610a13", "0020000a247600008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "0020000a247600008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "0020000a247600008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "0020000a247600008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a13", "0020000a247600008250404000", "0152501500000000"); // 0-20%
  } else if (trainConfig == 618){ // with added particles
    cuts.AddCutPCM("10110a23", "0020000a247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "0020000a247600008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "0020000a247600008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "0020000a247600008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 619){
    cuts.AddCutPCM("14610a23", "0020000a247600008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "0020000a247600008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "0020000a247600008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "0020000a247600008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a23", "0020000a247600008250404000", "0152501500000000"); // 0-20%
    // < 3 (b)
  } else if (trainConfig == 620){ //___-__________________
    cuts.AddCutPCM("10110a13", "0020000b247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "0020000b247600008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "0020000b247600008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "0020000b247600008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 621){
    cuts.AddCutPCM("14610a13", "0020000b247600008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "0020000b247600008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "0020000b247600008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "0020000b247600008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a13", "0020000b247600008250404000", "0152501500000000"); // 0-20%
  } else if (trainConfig == 622){ // with added particles
    cuts.AddCutPCM("10110a23", "0020000b247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "0020000b247600008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "0020000b247600008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "0020000b247600008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 623){
    cuts.AddCutPCM("14610a23", "0020000b247600008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "0020000b247600008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "0020000b247600008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "0020000b247600008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a23", "0020000b247600008250404000", "0152501500000000"); // 0-20%
    // < 2.5 (c)
  } else if (trainConfig == 624){ //___-__________________
    cuts.AddCutPCM("10110a13", "0020000c247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "0020000c247600008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "0020000c247600008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "0020000c247600008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 625){
    cuts.AddCutPCM("14610a13", "0020000c247600008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "0020000c247600008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "0020000c247600008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "0020000c247600008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a13", "0020000c247600008250404000", "0152501500000000"); // 0-20%
  } else if (trainConfig == 626){ // with added particles
    cuts.AddCutPCM("10110a23", "0020000c247600008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "0020000c247600008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "0020000c247600008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "0020000c247600008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 627){
    cuts.AddCutPCM("14610a23", "0020000c247600008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "0020000c247600008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "0020000c247600008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "0020000c247600008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a23", "0020000c247600008250404000", "0152501500000000"); // 0-20%

    // dEdx nSigma around pion line rejection up to 100GeV instead of standard 2GeV
  } else if (trainConfig == 628){ //-------_--------------
    cuts.AddCutPCM("10110a13", "00200009247000008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200009247000008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200009247000008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200009247000008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 629){
    cuts.AddCutPCM("14610a13", "00200009247000008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200009247000008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200009247000008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200009247000008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a13", "00200009247000008250404000", "0152501500000000"); // 0-20%
  } else if (trainConfig == 630){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247000008250404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200009247000008250404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200009247000008250404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200009247000008250404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 631){
    cuts.AddCutPCM("14610a23", "00200009247000008250404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200009247000008250404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200009247000008250404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200009247000008250404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a23", "00200009247000008250404000", "0152501500000000"); // 0-20%

    // PsiPair < 0.1 and Chi2 < 20
  } else if (trainConfig == 632){ //-------------_----------------------------------------
    cuts.AddCutPCM("10110a13", "00200009247600008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200009247600008850404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200009247600008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200009247600008850404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 633){
    cuts.AddCutPCM("14610a13", "00200009247600008850404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200009247600008850404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200009247600008850404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200009247600008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a13", "00200009247600008850404000", "0152501500000000"); // 0-20%
  } else if (trainConfig == 634){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247600008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200009247600008850404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200009247600008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200009247600008850404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 635){
    cuts.AddCutPCM("14610a23", "00200009247600008850404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200009247600008850404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200009247600008850404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200009247600008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a23", "00200009247600008850404000", "0152501500000000"); // 0-20%

    // dEdx nSigma around pion line rejection up to 100GeV instead of standard 2GeV
    // AND Chi2 < 20
  } else if (trainConfig == 636){ //-------_-----_----------------------------------------
    cuts.AddCutPCM("10110a13", "00200009247000008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13", "00200009247000008850404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a13", "00200009247000008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a13", "00200009247000008850404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 637){
    cuts.AddCutPCM("14610a13", "00200009247000008850404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a13", "00200009247000008850404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a13", "00200009247000008850404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a13", "00200009247000008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a13", "00200009247000008850404000", "0152501500000000"); // 0-20%
  } else if (trainConfig == 638){ // with added particles
    cuts.AddCutPCM("10110a23", "00200009247000008850404000", "0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a23", "00200009247000008850404000", "0152501500000000"); // 10-20%
    cuts.AddCutPCM("12410a23", "00200009247000008850404000", "0152501500000000"); // 20-40%
    cuts.AddCutPCM("30110a23", "00200009247000008850404000", "0152501500000000"); // 0-5%
  } else if (trainConfig == 639){
    cuts.AddCutPCM("14610a23", "00200009247000008850404000", "0152501500000000"); // 40-60%
    cuts.AddCutPCM("12510a23", "00200009247000008850404000", "0152501500000000"); // 20-50%
    cuts.AddCutPCM("16810a23", "00200009247000008850404000", "0152501500000000"); // 60-80%
    cuts.AddCutPCM("31210a23", "00200009247000008850404000", "0152501500000000"); // 5-10%
    cuts.AddCutPCM("10210a23", "00200009247000008850404000", "0152501500000000"); // 0-20%

    // TPC dEdx pion rejection 3 sigma between 0.4 and 8 GeV, afterwards 2 sigma (b77)
  } else if (trainConfig == 640){ //-------___--------------
      cuts.AddCutPCM("10110a13", "002000092b7700008250404000", "0152501500000000"); // 0-10%
  } else if (trainConfig == 641){ // added particles
      cuts.AddCutPCM("10110a23", "002000092b7700008250404000", "0152501500000000"); // 0-10%

      // TOF nsigma -5,5 (2)     standard: no cut (0)
  } else if (trainConfig == 642){ //-----------_------------
      cuts.AddCutPCM("10110a13", "00200009247602008250404000", "0152501500000000"); // 0-10%
  } else if (trainConfig == 643){ // added particles
      cuts.AddCutPCM("10110a23", "00200009247602008250404000", "0152501500000000"); // 0-10%

      // Chi2 vs PsiPair 0.15 exp(-0.065 chi2) (fd)
      // standard:  triangular cut with chi2 < 30 and psipair < 0.1 (25)
  } else if (trainConfig == 644){ //---------------__-------
      cuts.AddCutPCM("10110a13", "00200009247600008fd0404000", "0152501500000000"); // 0-10%
  } else if (trainConfig == 645){ // added particles
      cuts.AddCutPCM("10110a23", "00200009247600008fd0404000", "0152501500000000"); // 0-10%

      // qT pT dependent, q max = 0.11 p , q max = 0.04, 2D qT-alpha with alpha < 0.95 (a)
      // standard: 2D qT-alpha with qT < 0.05 and alpha < 0.95, independent of pT (8)
  } else if (trainConfig == 646){ //--------------_---------
      cuts.AddCutPCM("10110a13", "0020000924760000a250404000", "0152501500000000"); // 0-10%
  } else if (trainConfig == 647){ // added particles
      cuts.AddCutPCM("10110a23", "0020000924760000a250404000", "0152501500000000"); // 0-10%

      // cospoint < 0.99 (9)     standard: < 0.85 (4)
  } else if (trainConfig == 648){ //------------------_-----
      cuts.AddCutPCM("10110a13", "00200009247600008250904000", "0152501500000000"); // 0-10%
  } else if (trainConfig == 649){ // added particles
      cuts.AddCutPCM("10110a23", "00200009247600008250904000", "0152501500000000"); // 0-10%


  //****************************************************************************************************
  // 5.02TeV Pb-Pb LHC18qr
  //****************************************************************************************************
  } else if (trainConfig == 650){
    cuts.AddCutPCM("10930013","00200009327000008250400000","0132501500000000"); // 0-90%
  } else if (trainConfig == 651){
    cuts.AddCutPCM("10130013","00200009327000008250400000","0132501500000000"); //
    cuts.AddCutPCM("11530013","00200009327000008250400000","0132501500000000"); //
    cuts.AddCutPCM("15930013","00200009327000008250400000","0132501500000000"); //
  } else if (trainConfig == 652){
    cuts.AddCutPCM("10130a13","00200009327000008250400000","0132501500000000"); //
    cuts.AddCutPCM("11230a13","00200009327000008250400000","0132501500000000"); //
    cuts.AddCutPCM("12430a13","00200009327000008250400000","0132501500000000"); //
    cuts.AddCutPCM("14630a13","00200009327000008250400000","0132501500000000"); //
    cuts.AddCutPCM("16830a13","00200009327000008250400000","0132501500000000"); //
  } else if (trainConfig == 653){
    cuts.AddCutPCM("10130a13","00200009327000008250400000","0132501500000000"); //
    cuts.AddCutPCM("11310a13","00200009327000008250400000","0132501500000000"); //
    cuts.AddCutPCM("13530a13","00200009327000008250400000","0132501500000000"); //
    cuts.AddCutPCM("15910a13","00200009327000008250400000","0132501500000000"); //
  } else if (trainConfig == 654){
    cuts.AddCutPCM("10130a13","00200009f9730000dge0400000","0133103100000000"); //
    cuts.AddCutPCM("11310a13","00200009f9730000dge0400000","0133103100000000"); //
    cuts.AddCutPCM("13530a13","00200009f9730000dge0400000","0133103100000000"); //
    cuts.AddCutPCM("15910a13","00200009f9730000dge0400000","0133103100000000"); //
  } else if (trainConfig == 655){
    cuts.AddCutPCM("30130a13","00200009f9730000dge0400000","0143103100000000"); //
    cuts.AddCutPCM("31230a13","00200009f9730000dge0400000","0143103100000000"); //
    cuts.AddCutPCM("11210a13","00200009f9730000dge0400000","0143103100000000"); //
    cuts.AddCutPCM("12310a13","00200009f9730000dge0400000","0143103100000000"); //
  } else if (trainConfig == 656){
    cuts.AddCutPCM("13430a13","00200009f9730000dge0400000","0143103100000000"); //
    cuts.AddCutPCM("14530a13","00200009f9730000dge0400000","0143103100000000"); //
    cuts.AddCutPCM("15610a13","00200009f9730000dge0400000","0143103100000000"); //
    cuts.AddCutPCM("16710a13","00200009f9730000dge0400000","0143103100000000"); //
    cuts.AddCutPCM("17810a13","00200009f9730000dge0400000","0143103100000000"); //
    cuts.AddCutPCM("18910a13","00200009f9730000dge0400000","0143103100000000"); //
  } else if (trainConfig == 657){
    cuts.AddCutPCM("30130a13","0dm00009f9730000dge0404000","0143105100000000"); //
    cuts.AddCutPCM("31230a13","0dm00009f9730000dge0404000","0143105100000000"); //
    cuts.AddCutPCM("11210a13","0dm00009f9730000dge0404000","0143105100000000"); //
    cuts.AddCutPCM("12310a13","0dm00009f9730000dge0404000","0143105100000000"); //
  } else if (trainConfig == 658){
    cuts.AddCutPCM("13430a13","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("14530a13","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("15610a13","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("16710a13","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("17810a13","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("18910a13","0dm00009f9730000dge0404000","0143103100000000"); //
  } else if (trainConfig == 659){
    cuts.AddCutPCM("30130a23","0dm00009f9730000dge0404000","0143105100000000"); //
    cuts.AddCutPCM("31230a23","0dm00009f9730000dge0404000","0143105100000000"); //
    cuts.AddCutPCM("11210a23","0dm00009f9730000dge0404000","0143105100000000"); //
    cuts.AddCutPCM("12310a23","0dm00009f9730000dge0404000","0143105100000000"); //
  } else if (trainConfig == 660){
    cuts.AddCutPCM("13430a23","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("14530a23","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("15610a23","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("16710a23","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("17810a23","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("18910a23","0dm00009f9730000dge0404000","0143103100000000"); //
  } else if (trainConfig == 661){ // 657 with SDDSSD pileup cut
    cuts.AddCutPCM("30130d13","0dm00009f9730000dge0404000","0143105100000000"); //
    cuts.AddCutPCM("31230d13","0dm00009f9730000dge0404000","0143105100000000"); //
    cuts.AddCutPCM("11210d13","0dm00009f9730000dge0404000","0143105100000000"); //
    cuts.AddCutPCM("12310d13","0dm00009f9730000dge0404000","0143105100000000"); //
  } else if (trainConfig == 662){ // 658 with SDDSSD pileup cut
    cuts.AddCutPCM("13430d13","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("14530d13","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("15610d13","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("16710d13","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("17810d13","0dm00009f9730000dge0404000","0143103100000000"); //
    cuts.AddCutPCM("18910d13","0dm00009f9730000dge0404000","0143103100000000"); //
  } else if (trainConfig == 663){ // to compare against 2015 5TeV, phot and mes cutno from 628
    cuts.AddCutPCM("30130a13","00200009247000008250404000","0152501500000000"); // 0-5%
    cuts.AddCutPCM("31230a13","00200009247000008250404000","0152501500000000"); // 5-10%
    cuts.AddCutPCM("10130a13","00200009247000008250404000","0152501500000000"); // 0-10%
    cuts.AddCutPCM("11210a13","00200009247000008250404000","0152501500000000"); // 10-20%
    cuts.AddCutPCM("12310a13","00200009247000008250404000","0152501500000000"); // 20-30%
  } else if (trainConfig == 664){
    cuts.AddCutPCM("13430a13","00200009247000008250404000","0152501500000000"); //
    cuts.AddCutPCM("14530a13","00200009247000008250404000","0152501500000000"); //
    cuts.AddCutPCM("15610a13","00200009247000008250404000","0152501500000000"); //
    cuts.AddCutPCM("16710a13","00200009247000008250404000","0152501500000000"); //
    cuts.AddCutPCM("17810a13","00200009247000008250404000","0152501500000000"); //
    cuts.AddCutPCM("18910a13","00200009247000008250404000","0152501500000000"); //

  } else if (trainConfig == 670){//additional highpthadron studies
    cuts.AddCutPCM("10130a13","0dm00009f9730000dge0404000","5143103100000000"); //
    cuts.AddCutPCM("11310a13","0dm00009f9730000dge0404000","5143103100000000"); //
    cuts.AddCutPCM("13530a13","0dm00009f9730000dge0404000","5143103100000000"); //
    cuts.AddCutPCM("15910a13","0dm00009f9730000dge0404000","5143103100000000"); //

   // **************************************
   //  RBins studies for 5.02TeV Pb-Pb 18qr
   // **************************************
  } else if (trainConfig == 701){ // central , a,b,c bins , V0-TPC pileup rejection
    cuts.AddCutPCM("10130a13","0d200009f9730200dge0404000", "0652501500000000"); //  0-10%
    cuts.AddCutPCM("10130a13","0da00009f9730200dge0404000", "0652501500000000"); //  0-10%  a
    cuts.AddCutPCM("10130a13","0db00009f9730200dge0404000", "0652501500000000"); //  0-10%  b
    cuts.AddCutPCM("10130a13","0dc00009f9730200dge0404000", "0652501500000000"); //  0-10%  c
  } else if (trainConfig == 702){ // semicentral, a,b,c bins, V0-TPC pileup rejection
    cuts.AddCutPCM("13530a13","0d200009f9730200dge0404000", "0652501500000000"); //  20-50%
    cuts.AddCutPCM("13530a13","0da00009f9730200dge0404000", "0652501500000000"); //  20-50%  a
    cuts.AddCutPCM("13530a13","0db00009f9730200dge0404000", "0652501500000000"); //  20-50%  b
    cuts.AddCutPCM("13530a13","0dc00009f9730200dge0404000", "0652501500000000"); //  20-50%  c

  } else if (trainConfig == 703){ // peripheral, a,b,c bins, V0-TPC pileup rejection
    cuts.AddCutPCM("15910a13","0d200009f9730200dge0404000", "0652501500000000"); //  50-90%
    cuts.AddCutPCM("15910a13","0da00009f9730200dge0404000", "0652501500000000"); //  50-90%  a
    cuts.AddCutPCM("15910a13","0db00009f9730200dge0404000", "0652501500000000"); //  50-90%  b
    cuts.AddCutPCM("15910a13","0dc00009f9730200dge0404000", "0652501500000000"); //  50-90%  c

    // To be used with MBW from 5TeV Nch

  } else if (trainConfig == 851){ // central , a,b,c bins , V0-TPC pileup rejection
    cuts.AddCutPCM("10130a13","0d200009f9730200dge0404000", "0652501500000000"); //  0-10%
    cuts.AddCutPCM("10130a13","0da00009f9730200dge0404000", "0652501500000000"); //  0-10%  a
    cuts.AddCutPCM("10130a13","0db00009f9730200dge0404000", "0652501500000000"); //  0-10%  b
    cuts.AddCutPCM("10130a13","0dc00009f9730200dge0404000", "0652501500000000"); //  0-10%  c
  } else if (trainConfig == 852){ // semicentral, a,b,c bins, V0-TPC pileup rejection
    cuts.AddCutPCM("13530a13","0d200009f9730200dge0404000", "0652501500000000"); //  20-50%
    cuts.AddCutPCM("13530a13","0da00009f9730200dge0404000", "0652501500000000"); //  20-50%  a
    cuts.AddCutPCM("13530a13","0db00009f9730200dge0404000", "0652501500000000"); //  20-50%  b
    cuts.AddCutPCM("13530a13","0dc00009f9730200dge0404000", "0652501500000000"); //  20-50%  c

  } else if (trainConfig == 853){ // peripheral, a,b,c bins, V0-TPC pileup rejection
    cuts.AddCutPCM("15910a13","0d200009f9730200dge0404000", "0652501500000000"); //  50-90%
    cuts.AddCutPCM("15910a13","0da00009f9730200dge0404000", "0652501500000000"); //  50-90%  a
    cuts.AddCutPCM("15910a13","0db00009f9730200dge0404000", "0652501500000000"); //  50-90%  b
    cuts.AddCutPCM("15910a13","0dc00009f9730200dge0404000", "0652501500000000"); //  50-90%  c

  //****************************************************************************************************
  // 5.02TeV Pb-Pb photon cuts
  //****************************************************************************************************
    // photon cuts adapted from LHC11h analysis
  } else if (trainConfig == 900){
      cuts.AddCutPCM("10110a13", "00200009247002008850404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a13", "00200009247002008850404000", "0152501500000000"); // 20-40%
  } else if (trainConfig == 901){ // added particles
      cuts.AddCutPCM("10110a23", "00200009247002008850404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a23", "00200009247002008850404000", "0152501500000000"); // 20-40%

    // LHC15o photon cuts
  } else if (trainConfig == 902){
      cuts.AddCutPCM("10110a13", "002000092b770200afd0904000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a13", "002000092b770200afd0904000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a13", "002000092b770200afd0904000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 903){ // added particles
      cuts.AddCutPCM("10110a23", "002000092b770200afd0904000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a23", "002000092b770200afd0904000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a23", "002000092b770200afd0904000", "0152501500000000"); // 60-80%

      // cut studies, variation wrt configs 568,569

      // eta cut (0.8 (d) instead of 0.9(0))
  } else if (trainConfig == 904){//_------------------------
      cuts.AddCutPCM("10110a13", "0d200009247600008250404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a13", "0d200009247600008250404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a13", "0d200009247600008250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 905){  // added particles
      cuts.AddCutPCM("10110a23", "0d200009247600008250404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a23", "0d200009247600008250404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a23", "0d200009247600008250404000", "0152501500000000"); // 60-80%

      // R bin cut (m)
  } else if (trainConfig == 906){//-_-----------------------
      cuts.AddCutPCM("10110a13", "00m00009247600008250404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a13", "00m00009247600008250404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a13", "00m00009247600008250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 907){  // added particles
      cuts.AddCutPCM("10110a23", "00m00009247600008250404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a23", "00m00009247600008250404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a23", "00m00009247600008250404000", "0152501500000000"); // 60-80%

      // TOF PID (-7,7 sigma (1))
  } else if (trainConfig == 908){//------------_-------------
      cuts.AddCutPCM("10110a13", "00200009247601008250404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a13", "00200009247601008250404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a13", "00200009247601008250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 909){  // added particles
      cuts.AddCutPCM("10110a23", "00200009247601008250404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a23", "00200009247601008250404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a23", "00200009247601008250404000", "0152501500000000"); // 60-80%
      // TOF PID (-10,6 sigma (a))
  } else if (trainConfig == 910){//------------_-------------
      cuts.AddCutPCM("10110a13", "0020000924760a008250404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a13", "0020000924760a008250404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a13", "0020000924760a008250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 911){  // added particles
      cuts.AddCutPCM("10110a23", "0020000924760a008250404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a23", "0020000924760a008250404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a23", "0020000924760a008250404000", "0152501500000000"); // 60-80%

      // TPC PID electron selection: -3,3 (a) instead of -3,5 (2)
  } else if (trainConfig == 912){//-------_-----------------
      cuts.AddCutPCM("10110a13", "00200009a47600008250404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a13", "00200009a47600008250404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a13", "00200009a47600008250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 913){  // added particles
      cuts.AddCutPCM("10110a23", "00200009a47600008250404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a23", "00200009a47600008250404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a23", "00200009a47600008250404000", "0152501500000000"); // 60-80%

      // stricter pT dependent qT cut with max 0.03 (l) instead of 0.05 (8)
  } else if (trainConfig == 914){//---------------_---------
      cuts.AddCutPCM("10110a13", "0020000924760000l250404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a13", "0020000924760000l250404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a13", "0020000924760000l250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 915){  // added particles
      cuts.AddCutPCM("10110a23", "0020000924760000l250404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a23", "0020000924760000l250404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a23", "0020000924760000l250404000", "0152501500000000"); // 60-80%

      // stricter pT dependent Chi2 PsiPair (ld instead of 25)
  } else if (trainConfig == 916){//----------------__-------
      cuts.AddCutPCM("10110a13", "00200009247600008ld0404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a13", "00200009247600008ld0404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a13", "00200009247600008ld0404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 917){  // added particles
      cuts.AddCutPCM("10110a23", "00200009247600008ld0404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a23", "00200009247600008ld0404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a23", "00200009247600008ld0404000", "0152501500000000"); // 60-80%

      // TOF PID (-4,4 sigma (b)) and minimum track momentum of 0.4GeV required
  } else if (trainConfig == 918){//------------_-------------
      cuts.AddCutPCM("10110a13", "0020000924760b008250404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a13", "0020000924760b008250404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a13", "0020000924760b008250404000", "0152501500000000"); // 60-80%
  } else if (trainConfig == 919){  // added particles
      cuts.AddCutPCM("10110a23", "0020000924760b008250404000", "0152501500000000"); // 0-10%
      cuts.AddCutPCM("12410a23", "0020000924760b008250404000", "0152501500000000"); // 20-40%
      cuts.AddCutPCM("16810a23", "0020000924760b008250404000", "0152501500000000"); // 60-80%

      // gamma eta < 0.8 ('d') and meson rapidity < 0.8 ('1')
      // for central and semi-central events in addition TOF 'b', Chi2-PsiPair 'ld', QT-alpha-pT 'a', TPC-e-PID 'a', TPC-pi-PID 'b77'
  } else if (trainConfig == 920){//_------____-_--___-------
      cuts.AddCutPCM("10110a13", "0d200009ab770b00ald0404000", "0152101500000000"); // 0-10%
      cuts.AddCutPCM("12410a13", "0d200009ab770b00ald0404000", "0152101500000000"); // 20-40%
      cuts.AddCutPCM("16810a13", "0d200009247600008250404000", "0152101500000000"); // 60-80%
  } else if (trainConfig == 921){  // added particles
      cuts.AddCutPCM("10110a23", "0d200009ab770b00ald0404000", "0152101500000000"); // 0-10%
      cuts.AddCutPCM("12410a23", "0d200009ab770b00ald0404000", "0152101500000000"); // 20-40%
      cuts.AddCutPCM("16810a23", "0d200009247600008250404000", "0152101500000000"); // 60-80%

  //****************************************************************************************************

  } else if (trainConfig == 1001){
    cuts.AddCutPCM("60100013", "04200009297002003220000000", "0152204500900000");
  } else if (trainConfig == 1002) {
    cuts.AddCutPCM("61200013", "04200009297002003220000000", "0152204500900000");
  } else if (trainConfig == 1003) {
    cuts.AddCutPCM("50100013", "04200009297002003220000000", "0152204500900000");
  } else if (trainConfig == 1004) {
    cuts.AddCutPCM("50200013", "04200009297002003220000000", "0152204500900000");
  } else if (trainConfig == 1005) {
    cuts.AddCutPCM("51200013", "04200009297002003220000000", "0152204500900000");
  } else if (trainConfig == 1006) {
    cuts.AddCutPCM("52400013", "04200009297002003220000000", "0152204500900000");
  } else if (trainConfig == 1007) {
    cuts.AddCutPCM("54600013", "04200009297002003220000000", "0152206500900000");
  } else if (trainConfig == 1008) {
    cuts.AddCutPCM("54800013", "04200009297002003220000000", "0152206500900000");
  } else if (trainConfig == 1009) {
    cuts.AddCutPCM("54500013", "04200009297002003220000000", "0152206500900000");
  } else if (trainConfig == 1010) {
    cuts.AddCutPCM("55600013", "04200009297002003220000000", "0152206500900000");
  } else if (trainConfig == 1011) {
    cuts.AddCutPCM("56800013", "04200009297002003220000000", "0152206500900000");
  } else if (trainConfig == 1012) {
    cuts.AddCutPCM("56700013", "04200009297002003220000000", "0152206500900000");
  } else if (trainConfig == 1013) {
    cuts.AddCutPCM("57800013", "04200009297002003220000000", "0152206500900000");
  } else if (trainConfig == 1014) {
    cuts.AddCutPCM("46900013", "04200009297002003220000000", "0152206500900000");
  } else if (trainConfig == 1015) {
    cuts.AddCutPCM("58900013", "04200009297002003220000000", "0152206500900000");
  } else  if (trainConfig == 1016){
    cuts.AddCutPCM("60100013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1017) {
    cuts.AddCutPCM("61200013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1018) {
    cuts.AddCutPCM("50100013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1019) {
    cuts.AddCutPCM("50200013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1020) {
    cuts.AddCutPCM("51200013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1021) {
    cuts.AddCutPCM("52400013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1022) {
    cuts.AddCutPCM("54600013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1023) {
    cuts.AddCutPCM("54800013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1024) {
    cuts.AddCutPCM("54500013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1025) {
    cuts.AddCutPCM("55600013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1026) {
    cuts.AddCutPCM("56800013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1027) {
    cuts.AddCutPCM("56700013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1028) {
    cuts.AddCutPCM("57800013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1029) {
    cuts.AddCutPCM("46900013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1030) {
    cuts.AddCutPCM("58900013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1031) {
    cuts.AddCutPCM("50800013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1032) {
    cuts.AddCutPCM("52500013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1033) {
    cuts.AddCutPCM("53500013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1034) {
    cuts.AddCutPCM("54500013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1035) {
    cuts.AddCutPCM("53400013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1036) {
    cuts.AddCutPCM("52300013", "00200009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1037) {
    cuts.AddCutPCM("52500013", "00700009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1038) {
    cuts.AddCutPCM("52400013", "00700009247602008250400000", "0152501500000000");
  } else if (trainConfig == 1039) { //latest std 11h cut
    cuts.AddCutPCM("52400013", "00200009247602008250404000", "0152501500000000");
  } else if (trainConfig == 1040) {
    cuts.AddCutPCM("52500013", "00200009247602008250404000", "0152501500000000");
  } else if (trainConfig == 1041) {
    cuts.AddCutPCM("50100013", "00200009247602008250404000", "0152501500000000");
    // Xe-Xe configs
  } else if (trainConfig == 1042) {
    cuts.AddCutPCM("10910113","00200009327000008250400000","0163103100000000"); // 0-90
  } else if (trainConfig == 1043) {
    cuts.AddCutPCM("10210113","00200009327000008250400000","0163103100000000"); // 0-20
  } else if (trainConfig == 1043) {
    cuts.AddCutPCM("12410113","00200009327000008250400000","0163103100000000"); // 20-40
  } else if (trainConfig == 1043) {
    cuts.AddCutPCM("10410113","00200009327000008250400000","0163103100000000"); // 0-40
  } else if (trainConfig == 1043) {
    cuts.AddCutPCM("14910113","00200009327000008250400000","0163103100000000"); // 40-90
    // Pb-Pb 5.02 TeV
  } else if (trainConfig == 1044) {
    cuts.AddCutPCM("52310013","00200009247602008250404000","0652501500000000"); // 20-30%
  } else if (trainConfig == 1045) {
    cuts.AddCutPCM("53410013","00200009247602008250404000","0652501500000000"); // 30-40%
  } else if (trainConfig == 1046) {
    cuts.AddCutPCM("54510013","00200009247602008250404000","0652501500000000"); // 40-50%
  } else if (trainConfig == 1047) {
    cuts.AddCutPCM("55610013","00200009247602008250404000","0652501500000000"); // 50-60%
  } else if (trainConfig == 1048) {
    cuts.AddCutPCM("56710013","00200009247602008250404000","0652501500000000"); // 60-70%
  } else if (trainConfig == 1049) {
    cuts.AddCutPCM("57810013","00200009247602008250404000","0652501500000000"); // 70-80%
  } else if (trainConfig == 1050) {
    cuts.AddCutPCM("58910013","00200009247602008250404000","0652501500000000"); // 80-90%
  } else if (trainConfig == 1051) {
    cuts.AddCutPCM("52310613","00200009247602008250404000","0652501500000000"); // 20-30% with PU cut
  } else if (trainConfig == 1052) {
    cuts.AddCutPCM("53410613","00200009247602008250404000","0652501500000000"); // 30-40% with PU cut
  } else if (trainConfig == 1053) {
    cuts.AddCutPCM("54610613","00200009247602008250404000","0652501500000000"); // 40-60% with PU cut
  } else if (trainConfig == 1054) {
    cuts.AddCutPCM("56810613","00200009247602008250404000","0652501500000000"); // 60-80% with PU cut
  } else if (trainConfig == 1055) {
    cuts.AddCutPCM("58910613","00200009247602008250404000","0652501500000000"); // 80-90% with PU cut

  } else if (trainConfig == 2000){
    cuts.AddCutPCM("10910a13","40200009327000008250400000","0163103100000000"); // BDT test


  } else {
    Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  if(!cuts.AreValid()){
    cout << "\n\n****************************************************" << endl;
    cout << "ERROR: No valid cuts stored in CutHandlerConv! Returning..." << endl;
    cout << "****************************************************\n\n" << endl;
    return;
  }

  Int_t numberOfCuts = cuts.GetNCuts();

  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();
  TList *MesonCutList = new TList();

  TList *HeaderList = new TList();
  if (generatorName.CompareTo("LHC13d2")==0){
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
//    TObjString *Header3 = new TObjString("eta_2");
//    HeaderList->Add(Header3);

  } else if (generatorName.CompareTo("LHC12a17x_fix")==0){
    TObjString *Header1 = new TObjString("PARAM");
    HeaderList->Add(Header1);
  } else if (generatorName.CompareTo("LHC14a1a")==0){
    if (doWeightingPart == 1){
      TObjString *Header1 = new TObjString("pi0_1");
      HeaderList->Add(Header1);
    } else if (doWeightingPart == 2){
      TObjString *Header1 = new TObjString("eta_2");
      HeaderList->Add(Header1);
    } else if (doWeightingPart == 3){
      TString nameHeaders[2]    = { "pi0_1", "eta_2" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (doWeightingPart == 4){
      TObjString *Header1 = new TObjString("pi0EMC_3");
      HeaderList->Add(Header1);
    } else if (doWeightingPart == 5){
      TObjString *Header1 = new TObjString("etaEMC_5");
      HeaderList->Add(Header1);
    } else if (doWeightingPart == 6){
      TString nameHeaders[2]    = { "pi0EMC_3", "etaEMC_5" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (doWeightingPart == 7){
      TString nameHeaders[4]    = { "pi0_1", "eta_2", "pi0EMC_3", "etaEMC_5" };
      for (Int_t iHead = 0; iHead < 4; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (doWeightingPart == 8){
      TObjString *Header1 = new TObjString("gEMCPhoton_7");
      HeaderList->Add(Header1);
    } else if (doWeightingPart == 9){
      TString nameHeaders[10]   = { "Pythia_Jets_PtHard_1_10", "Pythia_Jets_PtHard_2_10", "Pythia_Jets_PtHard_3_10", "Pythia_Jets_PtHard_4_10", "Pythia_Jets_PtHard_5_10",
        "Pythia_Jets_PtHard_6_10", "Pythia_Jets_PtHard_7_10", "Pythia_Jets_PtHard_8_10", "Pythia_Jets_PtHard_9_10", "Pythia_Jets_PtHard_10_10"
      };
      for (Int_t iHead = 0; iHead < 10; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (doWeightingPart == 10){
      TObjString *Header1 = new TObjString("pythia_bele_10_10");
      HeaderList->Add(Header1);
    } else if (doWeightingPart == 11){
      TString nameHeaders[34]   = { "Pythia_Jets_PtHard_1_10", "Pythia_Jets_PtHard_2_10", "Pythia_Jets_PtHard_3_10", "Pythia_Jets_PtHard_4_10", "Pythia_Jets_PtHard_5_10",
        "Pythia_Jets_PtHard_6_10", "Pythia_Jets_PtHard_7_10", "Pythia_Jets_PtHard_8_10", "Pythia_Jets_PtHard_9_10", "Pythia_Jets_PtHard_10_10",
        "gEMCPhoton_7", "flat pt kstar_8", "flat pt kstarbar_9", "pythia_cele_10_10", "pythia_cele_18_10",
        "pythia_cele_30_10", "pythia_cele_50_10", "pythia_ccbar_10_10", "pythia_ccbar_18_10", "pythia_ccbar_30_10",
        "pythia_ccbar_50_10", "pythia_bele_10_10", "pythia_bele_18_10", "pythia_bele_30_10", "pythia_bele_50_10",
        "pythia_bbbar_10_10", "pythia_bbbar_18_10", "pythia_bbbar_30_10", "pythia_bbbar_50_10", "pi0_1",
        "eta_2", "pi0EMC_3", "etaEMC_5", "hijing_0"
      };
      for (Int_t iHead = 0; iHead < 34; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else if (doWeightingPart == 12){
      TObjString *Header1 = new TObjString("pi0PHS_4");
      HeaderList->Add(Header1);
    } else if (doWeightingPart == 13){
      TObjString *Header1 = new TObjString("etaPHS_6");
      HeaderList->Add(Header1);
    } else if (doWeightingPart == 14){
      TString nameHeaders[2]    = { "pi0PHS_4", "etaPHS_6" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    } else {
      TString nameHeaders[2]    = { "pi0_1", "eta_2" };
      for (Int_t iHead = 0; iHead < 2; iHead++ ){
        TObjString *Header = new TObjString(nameHeaders[iHead]);
        HeaderList->Add(Header);
      }
    }
  } else if (generatorName.CompareTo("LHC14a1b")==0 || generatorName.CompareTo("LHC14a1c")==0){
    if (doWeightingPart == 1 || doWeightingPart == 2 || doWeightingPart == 3 ){
      TObjString *Header1 = new TObjString("BOX");
      HeaderList->Add(Header1);
    } if (doWeightingPart == 4 || doWeightingPart == 5 || doWeightingPart == 6 ){
      TObjString *Header1 = new TObjString("PARAM_EMC");
      HeaderList->Add(Header1);
    } if (doWeightingPart == 12 || doWeightingPart == 13 || doWeightingPart == 14 ){
      TObjString *Header1 = new TObjString("PARAM_PHOS");
      HeaderList->Add(Header1);
    }
  } else if (generatorName.CompareTo("LHC16h4")==0 || generatorName.CompareTo("LHC19h3")==0){
    if (doWeightingPart == 1){
      TObjString *Header1 = new TObjString("Injector (pi0)_1");
      HeaderList->Add(Header1);
    } else if (doWeightingPart == 2){
      TObjString *Header1 = new TObjString("Injector (eta)_2");
      HeaderList->Add(Header1);
    } else {
      TObjString *Header1 = new TObjString("Injector (pi0)_1");
      HeaderList->Add(Header1);
      TObjString *Header2 = new TObjString("Injector (eta)_2");
      HeaderList->Add(Header2);
    }
  }

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];
  Bool_t initializedMatBudWeigths_existing    = kFALSE;

  for(Int_t i = 0; i<numberOfCuts; i++){

    analysisEventCuts[i] = new AliConvEventCuts();

    if(generatorName.Contains("LHC11h") && enableFlattening){
      cout << "entering the cent. flattening loop -> searching for file: " << fileNameCentFlattening.Data() << endl;

      if( fileNameCentFlattening.Contains("FlatFile") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableFlattening, fileNameCentFlattening, "Cent");
      } else if( fileNameCentFlattening.Contains("Good") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableFlattening, fileNameCentFlattening, "CentGoodRuns");
      }else if( fileNameCentFlattening.Contains("SemiGood") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableFlattening, fileNameCentFlattening, "CentSemiGoodRuns");
      }else {
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(enableFlattening, fileNameCentFlattening, "CentTotalRuns");
      }
    }

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    if (enableMultiplicityWeighting){
      cout << "INFO enableling mult weighting" << endl;
      if(periodNameAnchor.CompareTo("LHC15o")==0 || periodNameAnchor.CompareTo("LHC18q")==0){
        TString cutNumber = cuts.GetEventCut(i);
        TString centCut = cutNumber(0,3);  // first three digits of event cut
        dataInputMultHisto = Form("%s_%s", periodNameAnchor.Data(), centCut.Data());
        mcInputMultHisto   = Form("%s_%s", periodNameV0Reader.Data(), centCut.Data());
        cout << "INFO read " << dataInputMultHisto.Data() << " and " <<  mcInputMultHisto.Data() << " from " << fileNameMultWeights.Data() << endl;
      } else {
        if(i == 0){
          dataInputMultHisto      = Form("%s_0005", periodNameAnchor.Data());
          mcInputMultHisto        = "LHC14a1a_0005";
        } else if(i == 1){
          dataInputMultHisto      = Form("%s_0510", periodNameAnchor.Data());
          mcInputMultHisto        = "LHC14a1a_0510";
        } else if(i == 2){
          dataInputMultHisto      = Form("%s_0010", periodNameAnchor.Data());
          mcInputMultHisto        = "LHC14a1a_0010";
        } else if(i == 3){
          dataInputMultHisto      = Form("%s_2040", periodNameAnchor.Data());
          mcInputMultHisto        = "LHC14a1b_2040";
        } else if(i == 4){
          dataInputMultHisto      = Form("%s_2050", periodNameAnchor.Data());
          mcInputMultHisto        = "LHC14a1b_2050";
        }
      }
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile(kTRUE, fileNameMultWeights, dataInputMultHisto, mcInputMultHisto );
    }

    TString histoNameMCPi0PT = "";  TString histoNameMCEtaPT = "";  TString histoNameMCK0sPT = "";    // pT spectra to be weighted
    TString fitNamePi0PT = "";      TString fitNameEtaPT = "";      TString fitNameK0sPT = "";        // fit to correct shape of pT spectra
    Bool_t weightPi0 = kFALSE;      Bool_t weightEta = kFALSE;      Bool_t weightK0s = kFALSE;
    if(enablePtWeighting){
      cout << "INFO enabeling pT weighting" << endl;
      if(periodNameAnchor.CompareTo("LHC15o")==0 || periodNameAnchor.CompareTo("LHC18q")==0){
        TString eventCutString  = cuts.GetEventCut(i);
        TString eventCutShort   = eventCutString(0,6);   // first six digits
        weightPi0         = kTRUE;
        histoNameMCPi0PT  = Form("Pi0_%s_5TeV_%s",   periodNameV0Reader.Data(), eventCutString.Data());  // MC
        fitNamePi0PT      = Form("Pi0_Data_5TeV_%s", eventCutShort.Data());                              // fit to data
        weightEta         = kTRUE;
        histoNameMCEtaPT  = Form("Eta_%s_5TeV_%s",   periodNameV0Reader.Data(), eventCutString.Data());
        fitNameEtaPT      = Form("Eta_Data_5TeV_%s", eventCutShort.Data());
      }
      analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(weightPi0, weightEta, weightK0s, fileNamePtWeights, histoNameMCPi0PT, histoNameMCEtaPT, histoNameMCK0sPT, fitNamePi0PT, fitNameEtaPT, fitNameK0sPT);
    }

    if (  trainConfig == 1   || trainConfig == 5   || trainConfig == 9   || trainConfig == 13   || trainConfig == 17   ||
        trainConfig == 21   || trainConfig == 25   || trainConfig == 29   || trainConfig == 33   || trainConfig == 37  ||
        trainConfig == 300 || trainConfig == 302 || trainConfig == 304){
      if (i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
      if (i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
      if (i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
      if (i == 3 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNamePtWeights, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      if (i == 4 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNamePtWeights, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
    } else if (  trainConfig == 2   || trainConfig == 6   || trainConfig == 10   || trainConfig == 14   || trainConfig == 18   ||
          trainConfig == 22   || trainConfig == 26   || trainConfig == 30   || trainConfig == 34   || trainConfig == 38  ||
          trainConfig == 301 || trainConfig == 303 || trainConfig == 305){
      if (i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      if (i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
      if (i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
      if (i == 3 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
    } else if ( trainConfig == 3   || trainConfig == 7    || trainConfig == 11   || trainConfig == 15  || trainConfig == 19   ||
          trainConfig == 23   || trainConfig == 27   || trainConfig == 31   || trainConfig == 35   || trainConfig == 39){
      if (i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
      if (i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
      if (i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
      if (i == 3 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      if (i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
    } else if (  trainConfig == 4   ||trainConfig == 8     || trainConfig == 12   || trainConfig == 16   || trainConfig == 20   ||
          trainConfig == 24   || trainConfig == 28   || trainConfig == 32   || trainConfig == 36   || trainConfig == 40){
      if (i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      if (i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
      if (i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
      if (i == 3 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNamePtWeights, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
    }

    if (trainConfig == 56 ){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
        if ( i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 3 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_1020TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_1020TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
        if ( i == 4 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0020TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0020TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M","Eta_Fit_Data_PbPb_2760GeV_0020V0M");
      }
    }
    if (trainConfig == 57 ){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
        if ( i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 3 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_1020TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_1020TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
        if ( i == 4 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0020TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0020TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M","Eta_Fit_Data_PbPb_2760GeV_0020V0M");
      }
    }
    if (trainConfig == 58 ){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_4060TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4060TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M","Eta_Fit_Data_PbPb_2760GeV_4060V0M");
        if ( i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_6080TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_6080TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M","Eta_Fit_Data_PbPb_2760GeV_6080V0M");
        if ( i == 3 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_4080TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4080TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
        if ( i == 4 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_3050TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_3050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3050V0M","Eta_Fit_Data_PbPb_2760GeV_3050V0M");
      }
    }
    if (trainConfig == 59 ){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4060TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4060TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M","Eta_Fit_Data_PbPb_2760GeV_4060V0M");
        if ( i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_6080TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_6080TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M","Eta_Fit_Data_PbPb_2760GeV_6080V0M");
        if ( i == 3 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4080TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4080TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
        if ( i == 4 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_3050TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_3050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3050V0M","Eta_Fit_Data_PbPb_2760GeV_3050V0M");
      }
    }
    if (trainConfig == 60 ){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2030TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2030TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2030V0M","Eta_Fit_Data_PbPb_2760GeV_2030V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_3040TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_3040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3040V0M","Eta_Fit_Data_PbPb_2760GeV_3040V0M");
        if ( i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_4050TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4050V0M","Eta_Fit_Data_PbPb_2760GeV_4050V0M");
        if ( i == 3 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_5060TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_5060TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_5060V0M","Eta_Fit_Data_PbPb_2760GeV_5060V0M");
        if ( i == 4 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }
    if (trainConfig == 61 ){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2030TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2030TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2030V0M","Eta_Fit_Data_PbPb_2760GeV_2030V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_3040TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_3040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3040V0M","Eta_Fit_Data_PbPb_2760GeV_3040V0M");
        if ( i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4050TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4050V0M","Eta_Fit_Data_PbPb_2760GeV_4050V0M");
        if ( i == 3 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_5060TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_5060TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_5060V0M","Eta_Fit_Data_PbPb_2760GeV_5060V0M");
        if ( i == 4 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (trainConfig == 70    || trainConfig == 72    || trainConfig == 74    || trainConfig == 76   || trainConfig == 78   ||
        trainConfig == 80    || trainConfig == 82    || trainConfig == 84    || trainConfig == 86   || trainConfig == 88   ||
        trainConfig == 90    || trainConfig == 92    || trainConfig == 94    || trainConfig == 96   || trainConfig == 98   ||
        trainConfig == 100   || trainConfig == 102   || trainConfig == 104   || trainConfig == 106  || trainConfig == 108  ||
        trainConfig == 110   || trainConfig == 112   || trainConfig == 114   || trainConfig == 116  || trainConfig == 118  ||
        trainConfig == 120   || trainConfig == 122   || trainConfig == 124   || trainConfig == 126  || trainConfig == 128  ||
        trainConfig == 130   || trainConfig == 132   || trainConfig == 134   || trainConfig == 136  || trainConfig == 138  ||
        trainConfig == 140   || trainConfig == 142   || trainConfig == 144   || trainConfig == 146  || trainConfig == 148  ||
        trainConfig == 150   || trainConfig == 152   || trainConfig == 154   || trainConfig == 156  || trainConfig == 158  ||
        trainConfig == 160   || trainConfig == 162   || trainConfig == 164   || trainConfig == 166  || trainConfig == 168  ||
        trainConfig == 170   || trainConfig == 172   || trainConfig == 174   || trainConfig == 182  || trainConfig == 184  ||
        trainConfig == 234   || trainConfig == 236   || trainConfig == 238   || trainConfig == 240  || trainConfig == 242  ||
        trainConfig == 244   || trainConfig == 357   || trainConfig == 358   ){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
        if ( i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 3 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 4 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (trainConfig == 71    || trainConfig == 73    || trainConfig == 75    || trainConfig == 77   || trainConfig == 79   ||
        trainConfig == 81    || trainConfig == 83    || trainConfig == 85    || trainConfig == 87   || trainConfig == 89   ||
        trainConfig == 91    || trainConfig == 93    || trainConfig == 95    || trainConfig == 97   || trainConfig == 99   ||
        trainConfig == 101   || trainConfig == 103   || trainConfig == 105   || trainConfig == 107  || trainConfig == 109  ||
        trainConfig == 111   || trainConfig == 113   || trainConfig == 115   || trainConfig == 117  || trainConfig == 119  ||
        trainConfig == 121   || trainConfig == 123   || trainConfig == 125   || trainConfig == 127  || trainConfig == 129  ||
        trainConfig == 131   || trainConfig == 133   || trainConfig == 135   || trainConfig == 137  || trainConfig == 139  ||
        trainConfig == 141   || trainConfig == 143   || trainConfig == 145   || trainConfig == 147  || trainConfig == 149  ||
        trainConfig == 151   || trainConfig == 153   || trainConfig == 155   || trainConfig == 157  || trainConfig == 159  ||
        trainConfig == 161   || trainConfig == 163   || trainConfig == 165   || trainConfig == 167  || trainConfig == 169  ||
        trainConfig == 171   || trainConfig == 173   || trainConfig == 175   || trainConfig == 183  || trainConfig == 185  ||
        trainConfig == 235   || trainConfig == 237   || trainConfig == 239   || trainConfig == 241  || trainConfig == 243  ||
        trainConfig == 245){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
        if ( i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 3 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 4 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (trainConfig == 186   || trainConfig == 187 || trainConfig == 190   || trainConfig == 191 || trainConfig == 359   || trainConfig == 360){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }
    if (trainConfig == 188   || trainConfig == 189 || trainConfig == 192   || trainConfig == 193 || trainConfig == 361   || trainConfig == 362){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (trainConfig == 194 || trainConfig == 195 || trainConfig == 198 || trainConfig == 199 || trainConfig == 202 || trainConfig == 203 ||
        trainConfig == 206 || trainConfig == 207 || trainConfig == 210 || trainConfig == 211 || trainConfig == 214 || trainConfig == 215 ||
        trainConfig == 218 || trainConfig == 219 || trainConfig == 222 || trainConfig == 223 || trainConfig == 226 || trainConfig == 227 ||
        trainConfig == 230 || trainConfig == 231){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }
    if (trainConfig == 196 || trainConfig == 197 || trainConfig == 200 || trainConfig == 201 || trainConfig == 204 || trainConfig == 205 ||
        trainConfig == 208 || trainConfig == 209 || trainConfig == 212 || trainConfig == 213 || trainConfig == 216 || trainConfig == 217 ||
        trainConfig == 220 || trainConfig == 221 || trainConfig == 224 || trainConfig == 225 || trainConfig == 228 || trainConfig == 229 ||
        trainConfig == 232 || trainConfig == 233){
    if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
        if ( i == 3 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (trainConfig == 313   || trainConfig == 314 || trainConfig == 317   || trainConfig == 318 || trainConfig == 363   || trainConfig == 364){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }
    if (trainConfig == 315   || trainConfig == 316 || trainConfig == 319   || trainConfig == 320 || trainConfig == 365   || trainConfig == 366){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (trainConfig == 321 || trainConfig == 322 || trainConfig == 325 || trainConfig == 326 || trainConfig == 329 || trainConfig == 330 ||
        trainConfig == 333 || trainConfig == 334 || trainConfig == 337 || trainConfig == 338 || trainConfig == 341 || trainConfig == 342 ||
        trainConfig == 345 || trainConfig == 346 || trainConfig == 349 || trainConfig == 350 || trainConfig == 353 || trainConfig == 354){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
    }
    }
    if (trainConfig == 323 || trainConfig == 324 || trainConfig == 327 || trainConfig == 328 || trainConfig == 331 || trainConfig == 332 ||
        trainConfig == 335 || trainConfig == 336 || trainConfig == 339 || trainConfig == 340 || trainConfig == 343 || trainConfig == 344 ||
        trainConfig == 347 || trainConfig == 348 || trainConfig == 351 || trainConfig == 352 || trainConfig == 355 || trainConfig == 356){
    if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
        if ( i == 3 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if(trainConfig == 367 || trainConfig == 368){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_2030TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2030TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2030V0M","Eta_Fit_Data_PbPb_2760GeV_2030V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_3040TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_3040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3040V0M","Eta_Fit_Data_PbPb_2760GeV_3040V0M");
        if ( i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_3050TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_3050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3050V0M","Eta_Fit_Data_PbPb_2760GeV_3050V0M");
      }
    }
    if (trainConfig == 369   || trainConfig == 370){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2030TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2030TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2030V0M","Eta_Fit_Data_PbPb_2760GeV_2030V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_3040TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_3040TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3040V0M","Eta_Fit_Data_PbPb_2760GeV_3040V0M");
        if ( i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_3050TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_3050TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3050V0M","Eta_Fit_Data_PbPb_2760GeV_3050V0M");
      }
    }
    if(trainConfig == 371 || trainConfig == 372){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
        if ( i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_PbPb_2760GeV_1020TPC",generatorName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_1020TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
      }
    }
    if (trainConfig == 373   || trainConfig == 374){
      if (generatorName.CompareTo("LHC14a1a") ==0 || generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
        if ( i == 1 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
        if ( i == 2 && enablePtWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNamePtWeights, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_1020TPC",generatorName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_1020TPC",generatorName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
      }
    }

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    if(fMinPtHardSet)
      analysisEventCuts[i]->SetMinFacPtHard(minFacPtHard);
    if(fMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    if(fSingleMaxPtHardSet)
      analysisEventCuts[i]->SetMaxFacPtHardSingleParticle(maxFacPtHardSingle);
    if(fJetFinderUsage)
      analysisEventCuts[i]->SetUseJetFinderForOutliers(kTRUE);
    if(fUsePtHardFromFile)
      analysisEventCuts[i]->SetUsePtHardBinFromFile(kTRUE);
    if(fUseAddOutlierRej)
      analysisEventCuts[i]->SetUseAdditionalOutlierRejection(kTRUE);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->SetLightOutput(enableLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (generatorName.CompareTo("LHC14a1b") ==0 || generatorName.CompareTo("LHC14a1c") ==0 ){
      if (doWeightingPart == 1 || doWeightingPart == 4 || doWeightingPart == 12 ) analysisEventCuts[i]->SetAddedSignalPDGCode(111);
      if (doWeightingPart == 2 || doWeightingPart == 5 || doWeightingPart == 13 ) analysisEventCuts[i]->SetAddedSignalPDGCode(221);
    }
    EventCutList->Add(analysisEventCuts[i]);
    if (trainConfig == 37 || trainConfig == 38){
      analysisEventCuts[i]->SelectSpecialTrigger(AliVEvent::kMB, "AliVEvent::kMB" );
    }
    if (trainConfig == 39 || trainConfig == 40){
      analysisEventCuts[i]->SelectSpecialTrigger(AliVEvent::kCentral,"AliVEvent::kCentral" );
    }
    if (trainConfig == 41 || trainConfig == 42){
      analysisEventCuts[i]->SelectSpecialTrigger(AliVEvent::kSemiCentral,"AliVEvent::kSemiCentral" );
    }
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);


    analysisCuts[i] = new AliConversionPhotonCuts();
    if (trainConfig == 300 || trainConfig == 301 || trainConfig == 302  || trainConfig == 303  || trainConfig == 304  || trainConfig == 305 ||
        trainConfig == 306 || trainConfig == 307 || trainConfig == 308  || trainConfig == 309  || trainConfig == 310  || trainConfig == 311 || trainConfig == 312)
      analysisCuts[i]->SetSwitchToKappaInsteadOfNSigdEdxTPC(kTRUE);
    if(trainConfig == 182 || trainConfig == 183 || trainConfig == 184 || trainConfig == 185)
      analysisCuts[i]->SetDodEdxSigmaCut(kFALSE);

    if (enableMatBudWeightsPi0 > 0){
        if (isMC > 0){
            if (analysisCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,fileNameMatBudWeights)){
                initializedMatBudWeigths_existing = kTRUE;}
            else {cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;}
        }
        else {cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;}
    }

    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetLightOutput(enableLightOutput);
    if (enableElecDeDxPostCalibration){
      if (isMC == 0){
        if(fileNamedEdxPostCalib.CompareTo("") != 0){
          analysisCuts[i]->SetElecDeDxPostCalibrationCustomFile(fileNamedEdxPostCalib);
          cout << "Setting custom dEdx recalibration file: " << fileNamedEdxPostCalib.Data() << endl;
        }
        analysisCuts[i]->SetDoElecDeDxPostCalibration(enableElecDeDxPostCalibration);
        cout << "Enabled TPC dEdx recalibration." << endl;
      } else{
        cout << "ERROR enableElecDeDxPostCalibration set to True even if MC file. Automatically reset to 0"<< endl;
        enableElecDeDxPostCalibration=kFALSE;
        analysisCuts[i]->SetDoElecDeDxPostCalibration(kFALSE);
      }
    }
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());

    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(enableLightOutput);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");

    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);//Attention new switch small for Photon QA
  task->SetDoChargedPrimary(enableChargedPrimary);
  task->SetDoPlotVsCentrality(kTRUE);
  task->SetDoTHnSparse(enableTHnSparse);
  task->SetDoCentFlattening(enableFlattening);
  if (initializedMatBudWeigths_existing) {
      task->SetDoMaterialBudgetWeightingOfGammasForTrueMesons(kTRUE);
  }

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);
  Int_t nContainer = 2;
  for(Int_t i = 0; i<numberOfCuts; i++){
    if(enableQAPhotonTask>1){
      if (initializedMatBudWeigths_existing) {
	mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s MBW Photon DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      }else{
	mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s Photon DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      }
      nContainer++;
    }
    if(enableQAMesonTask>1){
      if (initializedMatBudWeigths_existing) {
	mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s MBW Meson DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      }else{
	mgr->ConnectOutput(task,nContainer,mgr->CreateContainer(Form("%s_%s_%s Meson DCA tree",(cuts.GetEventCut(i)).Data(),(cuts.GetPhotonCut(i)).Data(),(cuts.GetMesonCut(i)).Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("GammaConvV1_%i.root",trainConfig)) );
      }
      nContainer++;
    }
  }

  return;

}
