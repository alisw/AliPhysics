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
//CutHandler contains all cuts for a certain analysis and trainconfig,
//it automatically checks length of cutStrings and takes care of the number of added cuts,
//no specification of the variable 'numberOfCuts' needed anymore.
//***************************************************************************************
class CutHandlerConv{
  public:
    CutHandlerConv(Int_t nMax=10){
      nCuts=0; nMaxCuts=nMax; validCuts = true;
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; mesonCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConv: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerConv: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut;
      nCuts++;
      return;
    }
    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerConv: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetMesonCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* mesonCutArray;
};

//***************************************************************************************
//main function
//***************************************************************************************
void AddTask_GammaConvV1_PbPb(  Int_t     trainConfig                     = 1,                              //change different set of cuts
                                Int_t     isMC                            = 0,                              //run MC
                                Int_t     enableQAMesonTask               = 0,                              //enable QA in AliAnalysisTaskGammaConvV1
                                Int_t     enableQAPhotonTask              = 0,                              // enable additional QA task
                                TString   fileNameInputForWeighting       = "MCSpectraInput.root",          // path to file for weigting input
                                Int_t     headerSelectionInt              = 0,                              // 1 pi0 header, 2 eta header, 3 both (only for "named" boxes)
                                TString   cutnumberAODBranch              = "100000006008400000001500000",  // cutnumber with which AODs have been filtered
                                TString   periodName                      = "LHC13d2",                      // name of the period for added signals and weighting
                                Bool_t    doWeighting                     = kFALSE,                         // enable Weighting
                                Bool_t    enableUseTHnSparse              = kTRUE,                          // enable THnSparse for mixed event BG
                                Int_t     enableV0EffiStudies             = 0,                              // enables V0finding efficiency histograms
                                TString   fileNameInputForCentFlattening  = "InterpValuesAndFlattening.root",
                                Int_t     doFlattening                    = 0,
                                Bool_t    enableChargedPrimary            = kFALSE,
                                Bool_t    enableTriggerMimicking          = kFALSE,                         // enable trigger mimicking
                                Bool_t    enableTriggerOverlapRej         = kFALSE,                         // enable trigger overlap rejection
                                Float_t   maxFacPtHard                    = 3.,                             // maximum factor between hardest jet and ptHard generated
                                TString   periodNameV0Reader              = "",
                                Bool_t    doMultiplicityWeighting         = kFALSE,                         //
                                TString   fileNameInputForMultWeighing    = "Multiplicity.root",            //
                                TString   periodNameAnchor                = "",
                                Bool_t    runLightOutput                  = kFALSE,                         // switch to run light output (only essential histograms for afterburner)
                                Int_t     enableMatBudWeightsPi0          = 0,                              // 1 = three radial bins, 2 = 10 radial bins
                                TString   filenameMatBudWeights           = "MCInputFileMaterialBudgetWeights.root",
                                Bool_t    processAODcheckForV0s           = kFALSE,                         // flag for AOD check if V0s contained in AliAODs.root and AliAODGammaConversion.root
                                TString   additionalTrainConfig           = "0"                             // additional counter for trainconfig, this has to be always the last parameter
                          )  {

  Int_t isHeavyIon = 1;
  if (additionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + additionalTrainConfig.Atoi();
  }

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaConvV1_%i",trainConfig), "No analysis manager found.");
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  Bool_t isMCForOtherSettings = 0;
  if (isMC > 0) isMCForOtherSettings = 1;
  //========= Add PID Reponse to ANALYSIS manager ====
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMCForOtherSettings);
  }

  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = "00000070000000000500004000";
  TString cutnumberEvent = "10000003";

  Bool_t enableV0findingEffi = kFALSE;
  if(enableV0EffiStudies > 0){
    enableV0findingEffi = kTRUE;
    cutnumberPhoton = "00000070000000000500004000";
    if(enableV0EffiStudies == 1){
      cutnumberEvent = "50100013";
    }else if(enableV0EffiStudies == 2){
      cutnumberEvent = "52500013";
    }
  }

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
      if (periodNameV0Reader.CompareTo("") != 0) fEventCuts->SetPeriodEnum(periodNameV0Reader);
      fEventCuts->SetLightOutput(runLightOutput);
      if(fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())){
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
      fCuts->SetProcessAODCheck(processAODcheckForV0s); // if processAODcheckForV0s is kTRUE, also check for V0s to be contained in AliAODs and AliAODGammaConversion.root
      if(trainConfig == 182 || trainConfig == 183 || trainConfig == 184 || trainConfig == 185){
        fCuts->SetDodEdxSigmaCut(kFALSE);
      }
      if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
        fV0ReaderV1->SetConversionCuts(fCuts);
        fCuts->SetFillCutHistograms("",kTRUE);
      }
    }

    if(inputHandler->IsA()==AliAODInputHandler::Class()){
    // AOD mode
      fV0ReaderV1->SetDeltaAODBranchName(Form("GammaConv_%s_gamma",cutnumberAODBranch.Data()));
    }
    fV0ReaderV1->Init();

    AliLog::SetGlobalLogLevel(AliLog::kInfo);

    //connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);

  }

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  AliAnalysisTaskGammaConvV1 *task=NULL;
  task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(runLightOutput);
  // Cut Numbers to use in Analysis

  CutHandlerConv cuts;

  if (trainConfig == 1){ // Standard cuts
    cuts.AddCut("60100013", "04200009297002003220000000", "0152204500900000"); // 0-5%
    cuts.AddCut("61200013", "04200009297002003220000000", "0152204500900000"); // 5-10%
    cuts.AddCut("50100013", "04200009297002003220000000", "0152204500900000"); // 0-10%
    cuts.AddCut("51200013", "04200009297002003220000000", "0152204500900000"); // 10-20%
    cuts.AddCut("50200013", "04200009297002003220000000", "0152204500900000"); // 0-20%
  } else if (trainConfig == 2) { // Standard cuts
    cuts.AddCut("52400013", "04200009297002003220000000", "0152204500900000"); // 20-40%
    cuts.AddCut("54600013", "04200009297002003220000000", "0152206500900000"); // 40-60%
    cuts.AddCut("56800013", "04200009297002003220000000", "0152206500900000"); // 60-80%
    cuts.AddCut("54800013", "04200009297002003220000000", "0152206500900000"); // 40-80%
    cuts.AddCut("54900013", "04200009297002003220000000", "0152206500900000"); // 40-90%
  } else if (trainConfig == 3) { // Standard cuts only added signals
    cuts.AddCut("60100023", "04200009297002003220000000", "0152204500900000"); // 0-5%
    cuts.AddCut("61200023", "04200009297002003220000000", "0152204500900000"); // 5-10%
    cuts.AddCut("50100023", "04200009297002003220000000", "0152204500900000"); // 0-10%
    cuts.AddCut("51200023", "04200009297002003220000000", "0152204500900000"); // 10-20%
    cuts.AddCut("50200023", "04200009297002003220000000", "0152204500900000"); // 0-20%
  } else if (trainConfig == 4) { // Standard cuts only added signals
    cuts.AddCut("52400023", "04200009297002003220000000", "0152204500900000"); // 20-40%
    cuts.AddCut("54600023", "04200009297002003220000000", "0152206500900000"); // 40-60%
    cuts.AddCut("56800023", "04200009297002003220000000", "0152206500900000"); // 60-80%
    cuts.AddCut("54800023", "04200009297002003220000000", "0152206500900000"); // 20-40%
    cuts.AddCut("54900023", "04200009297002003220000000", "0152206500900000"); // 40-90%
  } else if (trainConfig == 5){ // R-minCut 7.5 cm
    cuts.AddCut("60100013", "04900009297002003220000000", "0152204500900000"); // 0-5%
    cuts.AddCut("61200013", "04900009297002003220000000", "0152204500900000"); // 5-10%
    cuts.AddCut("50100013", "04900009297002003220000000", "0152204500900000"); // 0-10%
    cuts.AddCut("51200013", "04900009297002003220000000", "0152204500900000"); // 10-20%
    cuts.AddCut("50200013", "04900009297002003220000000", "0152204500900000"); // 0-20%
  } else if (trainConfig == 6) { // R-minCut 7.5 cm
    cuts.AddCut("52400013", "04900009297002003220000000", "0152204500900000"); // 20-40%
    cuts.AddCut("54600013", "04900009297002003220000000", "0152206500900000"); // 40-60%
    cuts.AddCut("56800013", "04900009297002003220000000", "0152206500900000"); // 60-80%
    cuts.AddCut("54800013", "04900009297002003220000000", "0152206500900000"); // 40-80%
    cuts.AddCut("54900013", "04900009297002003220000000", "0152206500900000"); // 40-90%
  } else if (trainConfig == 7) {// R-minCut 7.5 cm
    cuts.AddCut("60100023", "04900009297002003220000000", "0152204500900000"); // 0-5%
    cuts.AddCut("61200023", "04900009297002003220000000", "0152204500900000"); // 5-10%
    cuts.AddCut("50100023", "04900009297002003220000000", "0152204500900000"); // 0-10%
    cuts.AddCut("51200023", "04900009297002003220000000", "0152204500900000"); // 10-20%
    cuts.AddCut("50200023", "04900009297002003220000000", "0152204500900000"); // 0-20%
  } else if (trainConfig == 8) { // R-minCut 7.5 cm
    cuts.AddCut("52400023", "04900009297002003220000000", "0152204500900000"); // 20-40%
    cuts.AddCut("54600023", "04900009297002003220000000", "0152206500900000"); // 40-60%
    cuts.AddCut("56800023", "04900009297002003220000000", "0152206500900000"); // 60-80%
    cuts.AddCut("54800023", "04900009297002003220000000", "0152206500900000"); // 20-40%
    cuts.AddCut("54900023", "04900009297002003220000000", "0152206500900000"); // 40-90%
  } else if (trainConfig == 9){ // R-minCut 12.5 cm
    cuts.AddCut("60100013", "04800009297002003220000000", "0152204500900000"); // 0-5%
    cuts.AddCut("61200013", "04800009297002003220000000", "0152204500900000"); // 5-10%
    cuts.AddCut("50100013", "04800009297002003220000000", "0152204500900000"); // 0-10%
    cuts.AddCut("51200013", "04800009297002003220000000", "0152204500900000"); // 10-20%
    cuts.AddCut("50200013", "04800009297002003220000000", "0152204500900000"); // 0-20%
  } else if (trainConfig == 10) { // R-minCut 12.5 cm
    cuts.AddCut("52400013", "04800009297002003220000000", "0152204500900000"); // 20-40%
    cuts.AddCut("54600013", "04800009297002003220000000", "0152206500900000"); // 40-60%
    cuts.AddCut("56800013", "04800009297002003220000000", "0152206500900000"); // 60-80%
    cuts.AddCut("54800013", "04800009297002003220000000", "0152206500900000"); // 40-80%
    cuts.AddCut("54900013", "04800009297002003220000000", "0152206500900000"); // 40-90%
  } else if (trainConfig == 11) {// R-minCut 12.5 cm
    cuts.AddCut("60100023", "04800009297002003220000000", "0152204500900000"); // 0-5%
    cuts.AddCut("61200023", "04800009297002003220000000", "0152204500900000"); // 5-10%
    cuts.AddCut("50100023", "04800009297002003220000000", "0152204500900000"); // 0-10%
    cuts.AddCut("51200023", "04800009297002003220000000", "0152204500900000"); // 10-20%
    cuts.AddCut("50200023", "04800009297002003220000000", "0152204500900000"); // 0-20%
  } else if (trainConfig == 12) { // R-minCut 12.5 cm
    cuts.AddCut("52400023", "04800009297002003220000000", "0152204500900000"); // 20-40%
    cuts.AddCut("54600023", "04800009297002003220000000", "0152206500900000"); // 40-60%
    cuts.AddCut("56800023", "04800009297002003220000000", "0152206500900000"); // 60-80%
    cuts.AddCut("54800023", "04800009297002003220000000", "0152206500900000"); // 20-40%
    cuts.AddCut("54900023", "04800009297002003220000000", "0152206500900000"); // 40-90%
  }else  if (trainConfig == 13){ // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCut("60100013", "03200009297002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200013", "03200009297002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100013", "03200009297002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200013", "03200009297002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200013", "03200009297002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 14) {  // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCut("52400013", "03200009297002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600013", "03200009297002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800013", "03200009297002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800013", "03200009297002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCut("54900013", "03200009297002003220000000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 15) { // LHC10h standard, eta 0.65, y = 0.6  added signals
    cuts.AddCut("60100023", "03200009297002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200023", "03200009297002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100023", "03200009297002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200023", "03200009297002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200023", "03200009297002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 16) { // LHC10h standard, eta 0.65, y = 0.6  added signals
    cuts.AddCut("52400023", "03200009297002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600023", "03200009297002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800023", "03200009297002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800023", "03200009297002003220000000", "0152306500900000"); // 20-40%
    cuts.AddCut("54900023", "03200009297002003220000000", "0152306500900000"); // 40-90%
  }else  if (trainConfig == 17){ // LHC10h standard, eta 0.65, y = 0.6, photon quality 1
    cuts.AddCut("60100013", "03200009297002003220020000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200013", "03200009297002003220020000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100013", "03200009297002003220020000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200013", "03200009297002003220020000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200013", "03200009297002003220020000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 18) {  // LHC10h standard, eta 0.65, y = 0.6, photon quality 1
    cuts.AddCut("52400013", "03200009297002003220020000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600013", "03200009297002003220020000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800013", "03200009297002003220020000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800013", "03200009297002003220020000", "0152306500900000"); // 40-80%
    cuts.AddCut("54900013", "03200009297002003220020000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 19) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 1  added signal
    cuts.AddCut("60100023", "03200009297002003220020000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200023", "03200009297002003220020000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100023", "03200009297002003220020000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200023", "03200009297002003220020000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200023", "03200009297002003220020000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 20) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 1 added signal
    cuts.AddCut("52400023", "03200009297002003220020000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600023", "03200009297002003220020000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800023", "03200009297002003220020000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800023", "03200009297002003220020000", "0152306500900000"); // 20-40%
    cuts.AddCut("54900023", "03200009297002003220020000", "0152306500900000"); // 40-90%
  }else  if (trainConfig == 21){ // LHC10h standard, eta 0.65, y = 0.6, photon quality 2
    cuts.AddCut("60100013", "03200009297002003220030000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200013", "03200009297002003220030000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100013", "03200009297002003220030000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200013", "03200009297002003220030000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200013", "03200009297002003220030000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 22) {  // LHC10h standard, eta 0.65, y = 0.6, photon quality 2
    cuts.AddCut("52400013", "03200009297002003220030000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600013", "03200009297002003220030000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800013", "03200009297002003220030000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800013", "03200009297002003220030000", "0152306500900000"); // 40-80%
    cuts.AddCut("54900013", "03200009297002003220030000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 23) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 2  added signal
    cuts.AddCut("60100023", "03200009297002003220030000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200023", "03200009297002003220030000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100023", "03200009297002003220030000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200023", "03200009297002003220030000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200023", "03200009297002003220030000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 24) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 2 added signal
    cuts.AddCut("52400023", "03200009297002003220030000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600023", "03200009297002003220030000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800023", "03200009297002003220030000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800023", "03200009297002003220030000", "0152306500900000"); // 20-40%
    cuts.AddCut("54900023", "03200009297002003220030000", "0152306500900000"); // 40-90%
  }else  if (trainConfig == 25){ // LHC10h standard, eta 0.65, y = 0.6, photon quality 3
    cuts.AddCut("60100013", "03200009297002003220040000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200013", "03200009297002003220040000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100013", "03200009297002003220040000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200013", "03200009297002003220040000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200013", "03200009297002003220040000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 26) {  // LHC10h standard, eta 0.65, y = 0.6, photon quality 3
    cuts.AddCut("52400013", "03200009297002003220040000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600013", "03200009297002003220040000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800013", "03200009297002003220040000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800013", "03200009297002003220040000", "0152306500900000"); // 40-80%
    cuts.AddCut("54900013", "03200009297002003220040000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 27) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 3  added signal
    cuts.AddCut("60100023", "03200009297002003220040000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200023", "03200009297002003220040000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100023", "03200009297002003220040000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200023", "03200009297002003220040000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200023", "03200009297002003220040000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 28) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 3 added signal
    cuts.AddCut("52400023", "03200009297002003220040000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600023", "03200009297002003220040000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800023", "03200009297002003220040000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800023", "03200009297002003220040000", "0152306500900000"); // 20-40%
    cuts.AddCut("54900023", "03200009297002003220040000", "0152306500900000"); // 40-90%
  }else  if (trainConfig == 29){ // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm
    cuts.AddCut("60100013", "03700009297002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200013", "03700009297002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100013", "03700009297002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200013", "03700009297002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200013", "03700009297002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 30) {  // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm
    cuts.AddCut("52400013", "03700009297002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600013", "03700009297002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800013", "03700009297002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800013", "03700009297002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCut("54900013", "03700009297002003220000000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 31) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm  added signals
    cuts.AddCut("60100023", "03700009297002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200023", "03700009297002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100023", "03700009297002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200023", "03700009297002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200023", "03700009297002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 32) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm  added signals
    cuts.AddCut("52400023", "03700009297002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600023", "03700009297002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800023", "03700009297002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800023", "03700009297002003220000000", "0152306500900000"); // 20-40%
    cuts.AddCut("54900023", "03700009297002003220000000", "0152306500900000"); // 40-90%
  }else  if (trainConfig == 33){ // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1
    cuts.AddCut("60100013", "03700009297002003220020000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200013", "03700009297002003220020000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100013", "03700009297002003220020000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200013", "03700009297002003220020000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200013", "03700009297002003220020000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 34) {  // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1
    cuts.AddCut("52400013", "03700009297002003220020000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600013", "03700009297002003220020000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800013", "03700009297002003220020000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800013", "03700009297002003220020000", "0152306500900000"); // 40-80%
    cuts.AddCut("54900013", "03700009297002003220020000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 35) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1  added signals
    cuts.AddCut("60100023", "03700009297002003220020000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200023", "03700009297002003220020000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100023", "03700009297002003220020000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200023", "03700009297002003220020000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200023", "03700009297002003220020000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 36) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1  added signals
    cuts.AddCut("52400023", "03700009297002003220020000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600023", "03700009297002003220020000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800023", "03700009297002003220020000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800023", "03700009297002003220020000", "0152306500900000"); // 20-40%
    cuts.AddCut("54900023", "03700009297002003220020000", "0152306500900000"); // 40-90%
  }else  if (trainConfig == 37){ // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3
    cuts.AddCut("60100013", "03700009297002003220040000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200013", "03700009297002003220040000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100013", "03700009297002003220040000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200013", "03700009297002003220040000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200013", "03700009297002003220040000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 38) {  // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3
    cuts.AddCut("52400013", "03700009297002003220040000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600013", "03700009297002003220040000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800013", "03700009297002003220040000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800013", "03700009297002003220040000", "0152306500900000"); // 40-80%
    cuts.AddCut("54900013", "03700009297002003220040000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 39) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3  added signals
    cuts.AddCut("60100023", "03700009297002003220040000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200023", "03700009297002003220040000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100023", "03700009297002003220040000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200023", "03700009297002003220040000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200023", "03700009297002003220040000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 40) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3  added signals
    cuts.AddCut("52400023", "03700009297002003220040000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600023", "03700009297002003220040000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800023", "03700009297002003220040000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800023", "03700009297002003220040000", "0152306500900000"); // 20-40%
    cuts.AddCut("54900023", "03700009297002003220040000", "0152306500900000"); // 40-90%
  } else if (trainConfig == 41){ // Standard cuts, eta 0.9, only to be run on data :kSemiCentral
    cuts.AddCut("60160013", "00200009297002003220000000", "0152504500900000"); // 0-5%
    cuts.AddCut("61260013", "00200009297002003220000000", "0152504500900000"); // 5-10%
    cuts.AddCut("50160013", "00200009297002003220000000", "0152504500900000"); // 0-10%
    cuts.AddCut("61260013", "00200009297002003220000000", "0152504500900000"); // 10-20%
    cuts.AddCut("50260013", "00200009297002003220000000", "0152504500900000"); // 0-20%
  } else if (trainConfig == 42) { // Standard cuts, eta 0.9, only to be run on data :kSemiCentral
    cuts.AddCut("52360013", "00200009297002003220000000", "0152506500900000");
    cuts.AddCut("53460013", "00200009297002003220000000", "0152506500900000");
    cuts.AddCut("54560013", "00200009297002003220000000", "0152506500900000");
    cuts.AddCut("55660013", "00200009297002003220000000", "0152506500900000");
    cuts.AddCut("56760013", "00200009297002003220000000", "0152506500900000");
  } else if (trainConfig == 43){ // Standard cuts, eta 0.9, only to be run on data
    cuts.AddCut("52300013", "00200009297002003220000000", "0152504500900000");
    cuts.AddCut("53400013", "00200009297002003220000000", "0152506500900000");
    cuts.AddCut("54500013", "00200009297002003220000000", "0152506500900000");
    cuts.AddCut("55600013", "00200009297002003220000000", "0152506500900000");
    cuts.AddCut("52500013", "00200009297002003220000000", "0152506500900000");
  } else if ( trainConfig == 44){ // qt elipse cut 0.05
    cuts.AddCut("60100013", "00200009297002008220000000", "0152506500900000"); // 0-5%
    cuts.AddCut("61200013", "00200009297002008220000000", "0152506500900000"); // 5-10%
    cuts.AddCut("50100013", "00200009297002008220000000", "0152506500900000"); // 0-10%
    cuts.AddCut("51200013", "00200009297002008220000000", "0152506500900000"); // 10-20%
    cuts.AddCut("50200013", "00200009297002008220000000", "0152506500900000"); // 0-20%
  } else if ( trainConfig == 45) { // qt elipse cut 0.05
    cuts.AddCut("52400013", "00200009297002008220000000", "0152506500900000"); // 20-40%
    cuts.AddCut("54600013", "00200009297002008220000000", "0152506500900000"); // 40-60%
    cuts.AddCut("56800013", "00200009297002008220000000", "0152506500900000"); // 60-80%
    cuts.AddCut("54800013", "00200009297002008220000000", "0152506500900000"); // 40-80%
    cuts.AddCut("53500013", "00200009297002008220000000", "0152506500900000"); // 30-50%
  } else if ( trainConfig == 46){ // qt elipse cut 0.05
    cuts.AddCut("52300013", "00200009297002008220000000", "0152506500900000"); // 20-30%
    cuts.AddCut("53400013", "00200009297002008220000000", "0152506500900000"); // 30-40%
    cuts.AddCut("54500013", "00200009297002008220000000", "0152506500900000"); // 40-50%
    cuts.AddCut("55600013", "00200009297002008220000000", "0152506500900000"); // 50-60%
    cuts.AddCut("52500013", "00200009297002008220000000", "0152506500900000"); // 60-70%
  } else if ( trainConfig == 47){ // cos(theta_point) cut 0.85
    cuts.AddCut("60100013", "00200009297002003220400000", "0152506500900000"); // 0-5%
    cuts.AddCut("61200013", "00200009297002003220400000", "0152506500900000"); // 5-10%
    cuts.AddCut("50100013", "00200009297002003220400000", "0152506500900000"); // 0-10%
    cuts.AddCut("51200013", "00200009297002003220400000", "0152506500900000"); // 10-20%
    cuts.AddCut("50200013", "00200009297002003220400000", "0152506500900000"); // 0-20%
  } else if ( trainConfig == 48) { // cos(theta_point) cut 0.85
    cuts.AddCut("52400013", "00200009297002003220400000", "0152506500900000"); // 20-40%
    cuts.AddCut("54600013", "00200009297002003220400000", "0152506500900000"); // 40-60%
    cuts.AddCut("56800013", "00200009297002003220400000", "0152506500900000"); // 60-80%
    cuts.AddCut("54800013", "00200009297002003220400000", "0152506500900000"); // 40-80%
    cuts.AddCut("53500013", "00200009297002003220400000", "0152506500900000"); // 30-50%
  } else if ( trainConfig == 49){ // cos(theta_point) cut 0.85
    cuts.AddCut("52300013", "00200009297002003220400000", "0152506500900000"); // 20-30%
    cuts.AddCut("53400013", "00200009297002003220400000", "0152506500900000"); // 30-40%
    cuts.AddCut("54500013", "00200009297002003220400000", "0152506500900000"); // 40-50%
    cuts.AddCut("55600013", "00200009297002003220400000", "0152506500900000"); // 50-60%
    cuts.AddCut("52500013", "00200009297002003220400000", "0152506500900000"); // 60-70%
  } else if ( trainConfig == 50){ // psi pair 2D 0.05
    cuts.AddCut("60100013", "00200009297002003260000000", "0152506500900000"); // 0-5%
    cuts.AddCut("61200013", "00200009297002003260000000", "0152506500900000"); // 5-10%
    cuts.AddCut("50100013", "00200009297002003260000000", "0152506500900000"); // 0-10%
    cuts.AddCut("51200013", "00200009297002003260000000", "0152506500900000"); // 10-20%
    cuts.AddCut("50200013", "00200009297002003260000000", "0152506500900000"); // 0-20%
  } else if ( trainConfig == 51) { // psi pair 2D 0.05
    cuts.AddCut("52400013", "00200009297002003260000000", "0152506500900000"); // 20-40%
    cuts.AddCut("54600013", "00200009297002003260000000", "0152506500900000"); // 40-60%
    cuts.AddCut("56800013", "00200009297002003260000000", "0152506500900000"); // 60-80%
    cuts.AddCut("54800013", "00200009297002003260000000", "0152506500900000"); // 40-80%
    cuts.AddCut("53500013", "00200009297002003260000000", "0152506500900000"); // 30-50%
  } else if ( trainConfig == 52){ // psi pair 2D 0.05
    cuts.AddCut("52300013", "00200009297002003260000000", "0152506500900000"); // 20-30%
    cuts.AddCut("53400013", "00200009297002003260000000", "0152506500900000"); // 30-40%
    cuts.AddCut("54500013", "00200009297002003260000000", "0152506500900000"); // 40-50%
    cuts.AddCut("55600013", "00200009297002003260000000", "0152506500900000"); // 50-60%
    cuts.AddCut("52500013", "00200009297002003260000000", "0152506500900000"); // 60-70%
  } else if ( trainConfig == 53){ // psi pair 2D 0.1
    cuts.AddCut("60100013", "00200009297002003250000000", "0152506500900000"); // 0-5%
    cuts.AddCut("61200013", "00200009297002003250000000", "0152506500900000"); // 5-10%
    cuts.AddCut("50100013", "00200009297002003250000000", "0152506500900000"); // 0-10%
    cuts.AddCut("51200013", "00200009297002003250000000", "0152506500900000"); // 10-20%
    cuts.AddCut("50200013", "00200009297002003250000000", "0152506500900000"); // 0-20%
  } else if ( trainConfig == 54) { // psi pair 2D 0.1
    cuts.AddCut("52400013", "00200009297002003250000000", "0152506500900000"); // 20-40%
    cuts.AddCut("54600013", "00200009297002003250000000", "0152506500900000"); // 40-60%
    cuts.AddCut("56800013", "00200009297002003250000000", "0152506500900000"); // 60-80%
    cuts.AddCut("54800013", "00200009297002003250000000", "0152506500900000"); // 40-80%
    cuts.AddCut("53500013", "00200009297002003250000000", "0152506500900000"); // 30-50%
  } else if ( trainConfig == 55){ // psi pair 2D 0.1
    cuts.AddCut("52300013", "00200009297002003250000000", "0152506500900000"); // 20-30%
    cuts.AddCut("53400013", "00200009297002003250000000", "0152506500900000"); // 30-40%
    cuts.AddCut("54500013", "00200009297002003250000000", "0152506500900000"); // 40-50%
    cuts.AddCut("55600013", "00200009297002003250000000", "0152506500900000"); // 50-60%
    cuts.AddCut("52500013", "00200009297002003250000000", "0152506500900000"); // 60-70%
  } else if ( trainConfig == 56){ // cleaner cuts central classes
    cuts.AddCut("60100013", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCut("51200013", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCut("50200013", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 57){ // cleaner cuts central classes - added signal
    cuts.AddCut("60100023", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCut("51200023", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCut("50200023", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 58) { // cleaner cuts semicentral classes
    cuts.AddCut("52400013", "00200009247602008250400000", "0152506500000000"); // 20-40%
    cuts.AddCut("54600013", "00200009247602008250400000", "0152506500000000"); // 40-60%
    cuts.AddCut("56800013", "00200009247602008250400000", "0152506500000000"); // 60-80%
    cuts.AddCut("54800013", "00200009247602008250400000", "0152506500000000"); // 40-80%
    cuts.AddCut("53500013", "00200009247602008250400000", "0152506500000000"); // 30-50%
  } else if ( trainConfig == 59) { // cleaner cuts semicentral classes - added signal
    cuts.AddCut("52400023", "00200009247602008250400000", "0152506500000000"); // 20-40%
    cuts.AddCut("54600023", "00200009247602008250400000", "0152506500000000"); // 40-60%
    cuts.AddCut("56800023", "00200009247602008250400000", "0152506500000000"); // 60-80%
    cuts.AddCut("54800023", "00200009247602008250400000", "0152506500000000"); // 40-80%
    cuts.AddCut("53500023", "00200009247602008250400000", "0152506500000000"); // 30-50%
  } else if ( trainConfig == 60){ // cleaner cuts finer centralities
    cuts.AddCut("52300013", "00200009247602008250400000", "0152506500000000"); // 20-30%
    cuts.AddCut("53400013", "00200009247602008250400000", "0152506500000000"); // 30-40%
    cuts.AddCut("54500013", "00200009247602008250400000", "0152506500000000"); // 40-50%
    cuts.AddCut("55600013", "00200009247602008250400000", "0152506500000000"); // 50-60%
    cuts.AddCut("52500013", "00200009247602008250400000", "0152506500000000"); // 60-70%
  } else if ( trainConfig == 61){ // cleaner cuts finer centralities - added signal
    cuts.AddCut("52300023", "00200009247602008250400000", "0152506500000000"); // 20-30%
    cuts.AddCut("53400023", "00200009247602008250400000", "0152506500000000"); // 30-40%
    cuts.AddCut("54500023", "00200009247602008250400000", "0152506500000000"); // 40-50%
    cuts.AddCut("55600023", "00200009247602008250400000", "0152506500000000"); // 50-60%
    cuts.AddCut("52500023", "00200009247602008250400000", "0152506500000000"); // 60-70%
  } else if ( trainConfig == 62){ // cleaner cuts finer centralities
    cuts.AddCut("62300013", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCut("63400013", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCut("64500013", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCut("65600013", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCut("66700013", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 63){ // cleaner cuts finer centralities - added signal
    cuts.AddCut("62300023", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCut("63400023", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCut("64500023", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCut("65600023", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCut("66700023", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 64){ // cleaner cuts
    cuts.AddCut("67800013", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCut("68900013", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCut("56700013", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCut("57800013", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCut("58900013", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 65){ // cleaner cuts added signal
    cuts.AddCut("67800023", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCut("68900023", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCut("56700023", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCut("57800023", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCut("58900023", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 66){ // cleaner cuts
    cuts.AddCut("70100013", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCut("71200013", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCut("72300013", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCut("73400013", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCut("74500013", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 67){ // cleaner cuts added signal
    cuts.AddCut("70100023", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCut("71200023", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCut("72300023", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCut("73400023", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCut("74500023", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 68){ // cleaner cuts
    cuts.AddCut("75600013", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCut("76700013", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCut("77800013", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCut("78900013", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCut("70900013", "00200009247602008250400000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 69){ // cleaner cuts added signal
    cuts.AddCut("75600023", "00200009247602008250400000", "0152506500000000"); // 0-5%
    cuts.AddCut("76700023", "00200009247602008250400000", "0152506500000000"); // 5-10%
    cuts.AddCut("77800023", "00200009247602008250400000", "0152506500000000"); // 0-10%
    cuts.AddCut("78900023", "00200009247602008250400000", "0152506500000000"); // 10-20%
    cuts.AddCut("70900023", "00200009247602008250400000", "0152506500000000"); // 0-20%

    //----------- here start the syst var for LHC11h (and not only) -------------------
  } else if ( trainConfig == 70){ // min R = 35 cm
    cuts.AddCut("60100013", "00700009247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00700009247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00700009247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00700009247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCut("52500013", "00700009247602008250404000", "0152501500000000"); // 0-20%
  } else if ( trainConfig == 71){ // min R = 35 cm added signal
    cuts.AddCut("60100023", "00700009247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00700009247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00700009247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00700009247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCut("52500023", "00700009247602008250404000", "0152501500000000"); // 0-20%
  } else if ( trainConfig == 72){ // min R = 35 cm with phi cut
    cuts.AddCut("60100013", "00716609247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00716609247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00716609247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00716609247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCut("52500013", "00716609247602008250404000", "0152501500000000"); // 0-20%
  } else if ( trainConfig == 73){ // min R = 35 cm with phi cut - added signal
    cuts.AddCut("60100023", "00716609247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00716609247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00716609247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00716609247602008250404000", "0152501500000000"); // 10-20%
    cuts.AddCut("52500023", "00716609247602008250404000", "0152501500000000"); // 0-20%
  } else if ( trainConfig == 74){ // single pt 0.075 --------------------------------------
    cuts.AddCut("60100013", "00200049247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200049247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200049247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200049247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200049247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 75){ // single pt 0.075 - added signal
    cuts.AddCut("60100023", "00200049247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200049247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200049247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200049247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200049247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 76){ // single pt 0.075 with phi cut
    cuts.AddCut("60100013", "00216649247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216649247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216649247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216649247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216649247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 77){ // single pt 0.075 with phi cut - added signal
    cuts.AddCut("60100023", "00216649247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216649247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216649247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216649247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216649247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 78){ // single pt 0.1 ----------------------------------------
    cuts.AddCut("60100013", "00200019247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200019247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200019247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200019247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200019247602008250404000", "0152501500000000"); //20-50%
  } else if ( trainConfig == 79){ // single pt 0.1 - added signal
    cuts.AddCut("60100023", "00200019247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200019247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200019247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200019247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200019247602008250404000", "0152501500000000"); //20-50%
  } else if ( trainConfig == 80){ // single pt 0.1 with phi cut
    cuts.AddCut("60100013", "00216619247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216619247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216619247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216619247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216619247602008250404000", "0152501500000000"); //20-50%
  } else if ( trainConfig == 81){ // single pt 0.1  with phi cut - added signal
    cuts.AddCut("60100023", "00216619247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216619247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216619247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216619247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216619247602008250404000", "0152501500000000"); //20-50%
  } else if ( trainConfig == 82){ // variation TPC cls 0.7 --------------------------------
    cuts.AddCut("60100013", "00200006247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200006247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200006247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200006247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200006247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 83){ // variation TPC cls 0.7 added signal
    cuts.AddCut("60100023", "00200006247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200006247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200006247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200006247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200006247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 84){ // variation TPC cls 0.7 with phi cut
    cuts.AddCut("60100013", "00216606247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216606247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216606247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216606247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216606247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 85){ // variation TPC cls 0.7 with phi cut - added signal
    cuts.AddCut("60100023", "00216606247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216606247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216606247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216606247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216606247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 86){ // variation TPC cls 0.35 -------------------------------
    cuts.AddCut("60100013", "00200008247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200008247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200008247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200008247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200008247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 87){ // variation TPC cls 0.35 - added signal
    cuts.AddCut("60100023", "00200008247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200008247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200008247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200008247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200008247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 88){ // variation TPC cls 0.35 with phi cut
    cuts.AddCut("60100013", "00216608247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216608247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216608247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216608247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216608247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 89){ // variation TPC cls 0.35  with phi cut - added signal
    cuts.AddCut("60100023", "00216608247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216608247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216608247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216608247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216608247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 90){ // variation edEdx -4,5 ---------------------------------
    cuts.AddCut("60100013", "00200009347602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009347602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009347602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009347602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009347602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 91){ // variation edEdx  -4,5 - added signal
    cuts.AddCut("60100023", "00200009347602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009347602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009347602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009347602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009347602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 92){ // variation edEdx  -4,5 with phi cut
    cuts.AddCut("60100013", "00216609347602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609347602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609347602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609347602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609347602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 93){ // variation edEdx  -4,5 with phi cut - added signal
    cuts.AddCut("60100023", "00216609347602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609347602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609347602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609347602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609347602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 94){ // variation edEdx  -2.5,4 ------------------------------
    cuts.AddCut("60100013", "00200009647602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009647602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009647602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009647602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009647602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 95){ // variation edEdx  -2.5,4 - added signal
    cuts.AddCut("60100023", "00200009647602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009647602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009647602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009647602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009647602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 96){ // variation edEdx  -2.5,4 with phi cut
    cuts.AddCut("60100013", "00216609647602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609647602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609647602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609647602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609647602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 97){ // variation edEdx  -2.5,4 with phi cut - added signal
    cuts.AddCut("60100023", "00216609647602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609647602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609647602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609647602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609647602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 98){ // variation pdEdx 2.0 sigma, 1 sigma high pt -----------
    cuts.AddCut("60100013", "00200009287602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009287602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009287602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009287602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009287602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 99){ // variation pdEdx 2.0 sigma, 1 sigma high pt - added signal
    cuts.AddCut("60100023", "00200009287602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009287602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009287602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009287602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009287602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 100){ // variation pdEdx 2.0 sigma, 1 sigma high pt with phi cut
    cuts.AddCut("60100013", "00216609287602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609287602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609287602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609287602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609287602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 101){ // variation pdEdx 2.0 sigma, 1 sigma high pt with phi cut - added signal
    cuts.AddCut("60100023", "00216609287602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609287602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609287602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609287602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609287602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 102){ // variation pdEdx 3.0 sigma, no high pt ---------------
    cuts.AddCut("60100013", "00200009247002008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009247002008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009247002008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009247002008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009247002008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 103){ // variation pdEdx 3.0 sigma, no high pt - added signal
    cuts.AddCut("60100023", "00200009247002008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009247002008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009247002008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009247002008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009247002008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 104){ // variation pdEdx 3.0 sigma, no high pt with phi cut
    cuts.AddCut("60100013", "00216609247002008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609247002008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609247002008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609247002008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609247002008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 105){ // variation pdEdx 3.0 sigma, no high pt  with phi cut - added signal
    cuts.AddCut("60100023", "00216609247002008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609247002008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609247002008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609247002008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609247002008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 106){ // variation pdEdx std (3.0) sigma, but 0.3-3. ----------
    cuts.AddCut("60100013", "00200009245402008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009245402008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009245402008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009245402008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009245402008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 107){ //variation pdEdx 3.0 sigma,  0.3-3. - added signal
    cuts.AddCut("60100023", "00200009245402008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009245402008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009245402008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009245402008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009245402008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 108){ //variation pdEdx 3.0 sigma, 0.3-3. with phi cut
    cuts.AddCut("60100013", "00216609245402008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609245402008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609245402008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609245402008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609245402008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 109){ //variation pdEdx 3.0 sigma, 0.3-3. with phi cut - added signal
    cuts.AddCut("60100023", "00216609245402008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609245402008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609245402008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609245402008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609245402008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 110){ // TOF el. PID -3,5 --------------------------------------
    cuts.AddCut("60100013", "00200009247603008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009247603008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009247603008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009247603008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009247603008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 111){ // TOF el. PID -3,5 - added signal
    cuts.AddCut("60100023", "00200009247603008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009247603008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009247603008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009247603008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009247603008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 112){ // TOF el. PID -3,5 with phi cut
    cuts.AddCut("60100013", "00216609247603008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609247603008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609247603008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609247603008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609247603008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 113){ // TOF el. PID -3,5 with phi cut - added signal
    cuts.AddCut("60100023", "00216609247603008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609247603008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609247603008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609247603008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609247603008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 114){ // TOF el. PID -2,3 ---------------------------------------
    cuts.AddCut("60100013", "00200009247604008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009247604008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009247604008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009247604008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009247604008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 115){ // TOF el. PID -2,3 - added signal
    cuts.AddCut("60100023", "00200009247604008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009247604008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009247604008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009247604008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009247604008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 116){ // TOF el. PID -2,3 with phi cut
    cuts.AddCut("60100013", "00216609247604008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609247604008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609247604008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609247604008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609247604008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 117){ // TOF el. PID -2,3  with phi cut - added signal
    cuts.AddCut("60100023", "00216609247604008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609247604008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609247604008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609247604008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609247604008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 118){ // qt 0.03 2D ----------------------------------------
    cuts.AddCut("60100013", "00200009247602009250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009247602009250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009247602009250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009247602009250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009247602009250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 119){ // qt 0.03 2D - added signal
    cuts.AddCut("60100023", "00200009247602009250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009247602009250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009247602009250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009247602009250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009247602009250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 120){ // qt 0.03 2D with phi cut
    cuts.AddCut("60100013", "00216609247602009250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609247602009250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609247602009250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609247602009250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609247602009250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 121){ // qt 0.03 2D with phi cut - added signal
    cuts.AddCut("60100023", "00216609247602009250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609247602009250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609247602009250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609247602009250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609247602009250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 122){ // qt 0.06 2D ---------------------------------------
    cuts.AddCut("60100013", "00200009247602002250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009247602002250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009247602002250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009247602002250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009247602002250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 123){ // qt 0.06 2D - added signal
    cuts.AddCut("60100023", "00200009247602002250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009247602002250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009247602002250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009247602002250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009247602002250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 124){ // qt 0.06 2D with phi cut
    cuts.AddCut("60100013", "00216609247602002250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609247602002250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609247602002250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609247602002250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609247602002250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 125){ // qt 0.06 2D with phi cut - added signal
    cuts.AddCut("60100023", "00216609247602002250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609247602002250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609247602002250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609247602002250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609247602002250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 126){ // chi2  50. ----------------------------------------
    cuts.AddCut("60100013", "00200009247602008150404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009247602008150404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009247602008150404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009247602008150404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009247602008150404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 127){ // chi2  50.  added signal
    cuts.AddCut("60100023", "00200009247602008150404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009247602008150404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009247602008150404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009247602008150404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009247602008150404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 128){ // chi2  50. with phi cut
    cuts.AddCut("60100013", "00216609247602008150404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609247602008150404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609247602008150404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609247602008150404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609247602008150404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 129){ // chi2  50. with phi cut - added signal
    cuts.AddCut("60100023", "00216609247602008150404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609247602008150404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609247602008150404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609247602008150404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609247602008150404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 130){ // chi2  20. -----------------------------------------
    cuts.AddCut("60100013", "00200009247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 131){ // chi2  20.  added signal
    cuts.AddCut("60100023", "00200009247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 132){ // chi2  20. with phi cut
    cuts.AddCut("60100013", "00216609247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 133){ // chi2  20. with phi cut - added signal
    cuts.AddCut("60100023", "00216609247602008850404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609247602008850404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609247602008850404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609247602008850404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609247602008850404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 134){ // psi pair 0.05 2D ---------------------------------
    cuts.AddCut("60100013", "00200009247602008260404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009247602008260404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009247602008260404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009247602008260404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009247602008260404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 135){ // psi pair 0.05 2D - added signal
    cuts.AddCut("60100023", "00200009247602008260404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009247602008260404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009247602008260404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009247602008260404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009247602008260404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 136){ // psi pair 0.05 2D with phi cut
    cuts.AddCut("60100013", "00216609247602008260404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609247602008260404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609247602008260404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609247602008260404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609247602008260404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 137){ // psi pair 0.05 2D with phi cut - added signal
    cuts.AddCut("60100023", "00216609247602008260404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609247602008260404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609247602008260404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609247602008260404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609247602008260404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 138){ // psi pair 0.2 2D ----------------------------------
    cuts.AddCut("60100013", "00200009247602008280404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009247602008280404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009247602008280404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009247602008280404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009247602008280404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 139){ // psi pair 0.2 2D - added signal
    cuts.AddCut("60100023", "00200009247602008280404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009247602008280404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009247602008280404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009247602008280404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009247602008280404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 140){ // psi pair 0.2 2D with phi cut
    cuts.AddCut("60100013", "00216609247602008280404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609247602008280404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609247602008280404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609247602008280404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609247602008280404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 141){ // psi pair 0.2 2D with phi cut - added signal
    cuts.AddCut("60100023", "00216609247602008280404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609247602008280404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609247602008280404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609247602008280404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609247602008280404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 142){ // cosPA -1 ----------------------------------------
    cuts.AddCut("60100013", "00200009247602008250004000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009247602008250004000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009247602008250004000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009247602008250004000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009247602008250004000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 143){ // cosPA -1 - added signal
    cuts.AddCut("60100023", "00200009247602008250004000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009247602008250004000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009247602008250004000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009247602008250004000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009247602008250004000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 144){ // cosPA -1 with phi cut
    cuts.AddCut("60100013", "00216609247602008250004000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609247602008250004000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609247602008250004000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609247602008250004000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609247602008250004000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 145){ // cosPA -1 with phi cut - added signal
    cuts.AddCut("60100023", "00216609247602008250004000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609247602008250004000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609247602008250004000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609247602008250004000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609247602008250004000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 146){ // variation alpha 0.75 ------------------------------
    cuts.AddCut("60100013", "00200009247602008250404000", "0152505500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009247602008250404000", "0152505500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009247602008250404000", "0152505500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009247602008250404000", "0152505500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009247602008250404000", "0152505500000000"); // 20-50%
  } else if ( trainConfig == 147){ // variation alpha 0.75 - added signal
    cuts.AddCut("60100023", "00200009247602008250404000", "0152505500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009247602008250404000", "0152505500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009247602008250404000", "0152505500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009247602008250404000", "0152505500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009247602008250404000", "0152505500000000"); // 20-50%
  } else if ( trainConfig == 148){ // variation alpha 0.75 with phi cut
    cuts.AddCut("60100013", "00216609247602008250404000", "0152505500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609247602008250404000", "0152505500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609247602008250404000", "0152505500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609247602008250404000", "0152505500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609247602008250404000", "0152505500000000"); // 20-50%
  } else if ( trainConfig == 149){ // variation alpha 0.75 with phi cut - added signal
    cuts.AddCut("60100023", "00216609247602008250404000", "0152505500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609247602008250404000", "0152505500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609247602008250404000", "0152505500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609247602008250404000", "0152505500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609247602008250404000", "0152505500000000"); // 20-50%
  } else if ( trainConfig == 150){ // variation alpha 1. --------------------------------
    cuts.AddCut("60100013", "00200009247602008250404000", "0152503500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009247602008250404000", "0152503500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009247602008250404000", "0152503500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009247602008250404000", "0152503500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009247602008250404000", "0152503500000000"); // 20-50%
  } else if ( trainConfig == 151){ // variation alpha 1. - added signal
    cuts.AddCut("60100023", "00200009247602008250404000", "0152503500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009247602008250404000", "0152503500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009247602008250404000", "0152503500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009247602008250404000", "0152503500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009247602008250404000", "0152503500000000"); // 20-50%
  } else if ( trainConfig == 152){ // variation alpha 1. with phi cut
    cuts.AddCut("60100013", "00216609247602008250404000", "0152503500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609247602008250404000", "0152503500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609247602008250404000", "0152503500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609247602008250404000", "0152503500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609247602008250404000", "0152503500000000"); // 20-50%
  } else if ( trainConfig == 153){ // variation alpha 1. with phi cut - added signal
    cuts.AddCut("60100023", "00216609247602008250404000", "0152503500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609247602008250404000", "0152503500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609247602008250404000", "0152503500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609247602008250404000", "0152503500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609247602008250404000", "0152503500000000"); // 20-50%
  } else if ( trainConfig == 154){ // variation with phi cut at 2.0 - 4.0 --------------
    cuts.AddCut("60100013", "00215509247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00215509247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00215509247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00215509247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00215509247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 155){ // variation with phi cut at 2.0 - 4.0 - added signal
    cuts.AddCut("60100023", "00215509247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00215509247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00215509247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00215509247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00215509247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 156){ // variation with phi cut at 2.4 - 3.6 --------------
    cuts.AddCut("60100013", "00217709247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00217709247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00217709247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00217709247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00217709247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 157){ // variation with phi cut at 2.4 - 3.6 - added signal
    cuts.AddCut("60100023", "00217709247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00217709247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00217709247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00217709247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00217709247602008250404000", "0152501500000000"); // 20-50%
    //-------------------- end LHC11h syst cut variations ---------------------------

  } else if ( trainConfig == 158){ // standard LHC11h cut selection
    cuts.AddCut("60100013", "00200009247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 159){ // standard LHC11h cut selection - double rejec - added signal
    cuts.AddCut("60100023", "00200009247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 160){ // standard LHC11h cut selection - double rejec with phi cut
    cuts.AddCut("60100013", "00216609247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 161){ // standard LHC11h cut selection - double rejec with phi cut - added signal
    cuts.AddCut("60100023", "00216609247602008250404000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609247602008250404000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609247602008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609247602008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609247602008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 162){ // RP EM V0 mult background
    cuts.AddCut("60100013", "00200009247602008250404000", "0652501500000000"); //
    cuts.AddCut("61200013", "00200009247602008250404000", "0652501500000000"); //
    cuts.AddCut("50100013", "00200009247602008250404000", "0652501500000000"); //
    cuts.AddCut("52400013", "00200009247602008250404000", "0652501500000000"); //
    cuts.AddCut("52500013", "00200009247602008250404000", "0652501500000000"); //
  } else if ( trainConfig == 163){ // added signals
    cuts.AddCut("60100023", "00200009247602008250404000", "0652501500000000"); //
    cuts.AddCut("61200023", "00200009247602008250404000", "0652501500000000"); //
    cuts.AddCut("50100023", "00200009247602008250404000", "0652501500000000"); //
    cuts.AddCut("52400023", "00200009247602008250404000", "0652501500000000"); //
    cuts.AddCut("52500023", "00200009247602008250404000", "0652501500000000"); //
  } else if ( trainConfig == 164){ //  with phi cut
    cuts.AddCut("60100013", "00216609247602008250404000", "0652501500000000"); //
    cuts.AddCut("61200013", "00216609247602008250404000", "0652501500000000"); //
    cuts.AddCut("50100013", "00216609247602008250404000", "0652501500000000"); //
    cuts.AddCut("52400013", "00216609247602008250404000", "0652501500000000"); //
    cuts.AddCut("52500013", "00216609247602008250404000", "0652501500000000"); //
  } else if ( trainConfig == 165){ //  with phi cut - added signals
    cuts.AddCut("60100023", "00216609247602008250404000", "0652501500000000"); //
    cuts.AddCut("61200023", "00216609247602008250404000", "0652501500000000"); //
    cuts.AddCut("50100023", "00216609247602008250404000", "0652501500000000"); //
    cuts.AddCut("52400023", "00216609247602008250404000", "0652501500000000"); //
    cuts.AddCut("52500023", "00216609247602008250404000", "0652501500000000"); //

  } else if ( trainConfig == 166){ // cleaner cuts photon Quality 1
    cuts.AddCut("60100013", "00200009297002008250420000", "0152506500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009297002008250420000", "0152506500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009297002008250420000", "0152506500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009297002008250420000", "0152506500000000"); // 10-20%
    cuts.AddCut("52500013", "00200009297002008250420000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 167){ // cleaner cuts added signal photon Quality 1
    cuts.AddCut("60100023", "00200009297002008250420000", "0152506500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009297002008250420000", "0152506500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009297002008250420000", "0152506500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009297002008250420000", "0152506500000000"); // 10-20%
    cuts.AddCut("52500023", "00200009297002008250420000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 168){ // cleaner cuts photon Quality 2
    cuts.AddCut("60100013", "00200009297002008250430000", "0152506500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009297002008250430000", "0152506500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009297002008250430000", "0152506500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009297002008250430000", "0152506500000000"); // 10-20%
    cuts.AddCut("52500013", "00200009297002008250430000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 169){ // cleaner cuts added signal photon Quality 2
    cuts.AddCut("60100023", "00200009297002008250430000", "0152506500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009297002008250430000", "0152506500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009297002008250430000", "0152506500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009297002008250430000", "0152506500000000"); // 10-20%
    cuts.AddCut("52500023", "00200009297002008250430000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 170){ // cleaner cuts photon Quality 3
    cuts.AddCut("60100013", "00200009297002008250440000", "0152506500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009297002008250440000", "0152506500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009297002008250440000", "0152506500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009297002008250440000", "0152506500000000"); // 10-20%
    cuts.AddCut("52500013", "00200009297002008250440000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 171){ // cleaner cuts added signal photon Quality 3
    cuts.AddCut("60100023", "00200009297002008250440000", "0152506500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009297002008250440000", "0152506500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009297002008250440000", "0152506500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009297002008250440000", "0152506500000000"); // 10-20%
    cuts.AddCut("52500023", "00200009297002008250440000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 172){ // cleaner cuts, photon Quality 1, min R = 35 cm
    cuts.AddCut("60100013", "00700009297002008250420000", "0152506500000000"); // 0-5%
    cuts.AddCut("61200013", "00700009297002008250420000", "0152506500000000"); // 5-10%
    cuts.AddCut("50100013", "00700009297002008250420000", "0152506500000000"); // 0-10%
    cuts.AddCut("52400013", "00700009297002008250420000", "0152506500000000"); // 10-20%
    cuts.AddCut("52500013", "00700009297002008250420000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 173){ // cleaner cuts added signal, photon Quality 1, min R = 35 cm
    cuts.AddCut("60100023", "00700009297002008250420000", "0152506500000000"); // 0-5%
    cuts.AddCut("61200023", "00700009297002008250420000", "0152506500000000"); // 5-10%
    cuts.AddCut("50100023", "00700009297002008250420000", "0152506500000000"); // 0-10%
    cuts.AddCut("52400023", "00700009297002008250420000", "0152506500000000"); // 10-20%
    cuts.AddCut("52500023", "00700009297002008250420000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 174){ // cleaner cuts, photon Quality 3, min R = 35 cm
    cuts.AddCut("60100013", "00700009297002008250440000", "0152506500000000"); // 0-5%
    cuts.AddCut("61200013", "00700009297002008250440000", "0152506500000000"); // 5-10%
    cuts.AddCut("50100013", "00700009297002008250440000", "0152506500000000"); // 0-10%
    cuts.AddCut("52400013", "00700009297002008250440000", "0152506500000000"); // 10-20%
    cuts.AddCut("52500013", "00700009297002008250440000", "0152506500000000"); // 0-20%
  } else if ( trainConfig == 175){ // cleaner cuts added signal, photon Quality 3, min R = 35 cm
    cuts.AddCut("60100023", "00700009297002008250440000", "0152506500000000"); // 0-5%
    cuts.AddCut("61200023", "00700009297002008250440000", "0152506500000000"); // 5-10%
    cuts.AddCut("50100023", "00700009297002008250440000", "0152506500000000"); // 0-10%
    cuts.AddCut("52400023", "00700009297002008250440000", "0152506500000000"); // 10-20%
    cuts.AddCut("52500023", "00700009297002008250440000", "0152506500000000"); // 0-20%

  } else if ( trainConfig == 176){ // flow cuts with eta = 0.9, y = 0.85
    cuts.AddCut("60100013", "00200009297002008250400000", "0152506500000000");
    cuts.AddCut("61200013", "00200009297002008250400000", "0152506500000000");
    cuts.AddCut("51200013", "00200009297002008250400000", "0152506500000000");
    cuts.AddCut("52300013", "00200009297002008250400000", "0152506500000000");
    cuts.AddCut("53400013", "00200009297002008250400000", "0152506500000000");
    cuts.AddCut("54600013", "00200009297002008250400000", "0152506500000000");
    cuts.AddCut("56800013", "00200009297002008250400000", "0152506500000000");
  } else if ( trainConfig == 177){ // flow cuts with eta = 0.65, y = 0.6
    cuts.AddCut("60100013", "03200009297002008250400000", "0152306500000000");
    cuts.AddCut("61200013", "03200009297002008250400000", "0152306500000000");
    cuts.AddCut("51200013", "03200009297002008250400000", "0152306500000000");
    cuts.AddCut("52300013", "03200009297002008250400000", "0152306500000000");
    cuts.AddCut("53400013", "03200009297002008250400000", "0152306500000000");
    cuts.AddCut("54600013", "03200009297002008250400000", "0152306500000000");
    cuts.AddCut("56800013", "03200009297002008250400000", "0152306500000000");
  } else if ( trainConfig == 178){ // flow cuts with eta = 0.6, y = 0.5
    cuts.AddCut("60100013", "01200009297002008250400000", "0152406500000000");
    cuts.AddCut("61200013", "01200009297002008250400000", "0152406500000000");
    cuts.AddCut("51200013", "01200009297002008250400000", "0152406500000000");
    cuts.AddCut("52300013", "01200009297002008250400000", "0152406500000000");
    cuts.AddCut("53400013", "01200009297002008250400000", "0152406500000000");
    cuts.AddCut("54600013", "01200009297002008250400000", "0152406500000000");
    cuts.AddCut("56800013", "01200009297002008250400000", "0152406500000000");
  } else if ( trainConfig == 179){ // flow cuts with eta = 0.9, y = 0.85
    cuts.AddCut("60100013", "00200009297002208250400000", "0152506500000000");
    cuts.AddCut("61200013", "00200009297002208250400000", "0152506500000000");
    cuts.AddCut("51200013", "00200009297002208250400000", "0152506500000000");
    cuts.AddCut("52300013", "00200009297002208250400000", "0152506500000000");
    cuts.AddCut("53400013", "00200009297002208250400000", "0152506500000000");
    cuts.AddCut("54600013", "00200009297002208250400000", "0152506500000000");
    cuts.AddCut("56800013", "00200009297002208250400000", "0152506500000000");
  } else if ( trainConfig == 180){ // flow cuts with eta = 0.65, y = 0.6
    cuts.AddCut("60100013", "03200009297002208250400000", "0152306500000000");
    cuts.AddCut("61200013", "03200009297002208250400000", "0152306500000000");
    cuts.AddCut("51200013", "03200009297002208250400000", "0152306500000000");
    cuts.AddCut("52300013", "03200009297002208250400000", "0152306500000000");
    cuts.AddCut("53400013", "03200009297002208250400000", "0152306500000000");
    cuts.AddCut("54600013", "03200009297002208250400000", "0152306500000000");
    cuts.AddCut("56800013", "03200009297002208250400000", "0152306500000000");
  } else if ( trainConfig == 181){ // flow cuts with eta = 0.6, y = 0.5
    cuts.AddCut("60100013", "01200009297002208250400000", "0152406500000000");
    cuts.AddCut("61200013", "01200009297002208250400000", "0152406500000000");
    cuts.AddCut("51200013", "01200009297002208250400000", "0152406500000000");
    cuts.AddCut("52300013", "01200009297002208250400000", "0152406500000000");
    cuts.AddCut("53400013", "01200009297002208250400000", "0152406500000000");
    cuts.AddCut("54600013", "01200009297002208250400000", "0152406500000000");
    cuts.AddCut("56800013", "01200009297002208250400000", "0152406500000000");

  } else if ( trainConfig == 182){ // standard LHC11h cut selection - no dEdx
    cuts.AddCut("60100013", "00200009000002008250400000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00200009000002008250400000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00200009000002008250400000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009000002008250400000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009000002008250400000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 183){ // standard LHC11h cut selection - no dEdx add signals
    cuts.AddCut("60100023", "00200009000002008250400000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00200009000002008250400000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00200009000002008250400000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00200009000002008250400000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00200009000002008250400000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 184){ // standard LHC11h cut selection - no dEdx with phi cut
    cuts.AddCut("60100013", "00216609000002008250400000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200013", "00216609000002008250400000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100013", "00216609000002008250400000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609000002008250400000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609000002008250400000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 185){ // standard LHC11h cut selection - no dEdx with phi cut add signals
    cuts.AddCut("60100023", "00216609000002008250400000", "0152501500000000"); // 0-5%
    cuts.AddCut("61200023", "00216609000002008250400000", "0152501500000000"); // 5-10%
    cuts.AddCut("50100023", "00216609000002008250400000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400023", "00216609000002008250400000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500023", "00216609000002008250400000", "0152501500000000"); // 20-50%

  } else if ( trainConfig == 186){ // dir photons (tighter dEdx+TOF)
    cuts.AddCut("50100013", "00200009847005008750404000", "0152501500000000"); //
  } else if ( trainConfig == 187){ //
    cuts.AddCut("50100013", "00216609847005008750404000", "0152501500000000"); //
  } else if ( trainConfig == 188){ // dir photons (tighter dEdx+TOF)
    cuts.AddCut("52400013", "00200009847005008750404000", "0152501500000000"); //
  } else if ( trainConfig == 189){ //
    cuts.AddCut("52400013", "00216609847005008750404000", "0152501500000000"); //
  } else if ( trainConfig == 190){ // dir photons (tighter dEdx+TOF)
    cuts.AddCut("52500013", "00200009847005008750404000", "0152501500000000"); //
  } else if ( trainConfig == 191){ //
    cuts.AddCut("52500013", "00216609847005008750404000", "0152501500000000"); //

  } else if ( trainConfig == 192){ // phi variation
    cuts.AddCut("50100013", "00215509847005008750404000", "0152501500000000"); // 2 - 4 rad
    cuts.AddCut("50100013", "00217709847005008750404000", "0152501500000000"); // 2.6 - 3.4 rad
  } else if ( trainConfig == 193){ // phi variation
    cuts.AddCut("52400013", "00215509247002008750404000", "0152501500000000"); // 2 - 4 rad
    cuts.AddCut("52400013", "00217709847005008750404000", "0152501500000000"); // 2.6 - 3.4 rad
  } else if ( trainConfig == 194){ // single pt variation
    cuts.AddCut("50100013", "00200049847005008750404000", "0152501500000000"); // 0.075
    cuts.AddCut("50100013", "00200019847005008750404000", "0152501500000000"); // 0.1
  } else if ( trainConfig == 195){ //  single pt variation
    cuts.AddCut("52400013", "00200049847005008750404000", "0152501500000000"); // 0.075
    cuts.AddCut("52400013", "00200019847005008750404000", "0152501500000000"); // 0.1
  } else if ( trainConfig == 196){ // single pt variation - phi
    cuts.AddCut("50100013", "00216649847005008750404000", "0152501500000000"); // 0.075
    cuts.AddCut("50100013", "00216619847005008750404000", "0152501500000000"); // 0.1
  } else if ( trainConfig == 197){ //  single pt variation - phi
    cuts.AddCut("52400013", "00216649847005008750404000", "0152501500000000"); // 0.075
    cuts.AddCut("52400013", "00216619847005008750404000", "0152501500000000"); // 0.1
  } else if ( trainConfig == 198){ // TPC clusters
    cuts.AddCut("50100013", "00200008847005008750404000", "0152501500000000"); // 0.35
    cuts.AddCut("50100013", "00200006847005008750404000", "0152501500000000"); // 0.7
  } else if ( trainConfig == 199){ // TPC clusters
    cuts.AddCut("52400013", "00200008847005008750404000", "0152501500000000"); // 0.35
    cuts.AddCut("52400013", "00200006847005008750404000", "0152501500000000"); // 0.7
  } else if ( trainConfig == 200){ // TPC clusters - phi
    cuts.AddCut("50100013", "00216608847005008750404000", "0152501500000000"); // 0.35
    cuts.AddCut("50100013", "00216606847005008750404000", "0152501500000000"); // 0.7
  } else if ( trainConfig == 201){ // TPC clusters - phi
    cuts.AddCut("52400013", "00216608847005008750404000", "0152501500000000"); // 0.35
    cuts.AddCut("52400013", "00216606847005008750404000", "0152501500000000"); // 0.7
  } else if ( trainConfig == 202){ // edEdx
    cuts.AddCut("50100013", "00200009247005008750404000", "0152501500000000"); // -3, 5
    cuts.AddCut("50100013", "00200009647005008750404000", "0152501500000000"); // -2.5, 4
  } else if ( trainConfig == 203){ // edEdx
    cuts.AddCut("52400013", "00200009247005008750404000", "0152501500000000"); // -3, 5
    cuts.AddCut("52400013", "00200009647005008750404000", "0152501500000000"); // -2.5, 4
  } else if ( trainConfig == 204){ // edEdx - phi
    cuts.AddCut("50100013", "00216609247005008750404000", "0152501500000000"); // -3, 5
    cuts.AddCut("50100013", "00216609647005008750404000", "0152501500000000"); // -2.5, 4
  } else if ( trainConfig == 205){ // edEdx - phi
    cuts.AddCut("52400013", "00216609247005008750404000", "0152501500000000"); // -3, 5
    cuts.AddCut("52400013", "00216609647005008750404000", "0152501500000000"); // -2.5, 4
  } else if ( trainConfig == 206){ // pdEdx
    cuts.AddCut("50100013", "00200009840005008750404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCut("50100013", "00200009847605008750404000", "0152501500000000"); // 3&1sigma
  } else if ( trainConfig == 207){ // pdEdx
    cuts.AddCut("52400013", "00200009840005008750404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCut("52400013", "00200009847605008750404000", "0152501500000000"); // 3&1sigma
  } else if ( trainConfig == 208){ // pdEdx - phi
    cuts.AddCut("50100013", "00216609840005008750404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCut("50100013", "00216609847605008750404000", "0152501500000000"); // 3&1sigma
  } else if ( trainConfig == 209){ // pdEdx - phi
    cuts.AddCut("52400013", "00216609840005008750404000", "0152501500000000"); // low pt at 0.5
    cuts.AddCut("52400013", "00216609847605008750404000", "0152501500000000"); // 3&1sigma
  } else if ( trainConfig == 210){ // pdEdx & TOF
    cuts.AddCut("50100013", "00200009837005008750404000", "0152501500000000"); // 2.5sigma
    cuts.AddCut("50100013", "00200009847002008750404000", "0152501500000000"); // TOF -5, 5
  } else if ( trainConfig == 211){ // pdEdx & TOF
    cuts.AddCut("52400013", "00200009837005008750404000", "0152501500000000"); // 2.5sigma
    cuts.AddCut("52400013", "00200009847002008750404000", "0152501500000000"); // TOF -5, 5
  } else if ( trainConfig == 212){ // pdEdx & TOF - phi
    cuts.AddCut("50100013", "00216609837005008750404000", "0152501500000000"); // 2.5sigma
    cuts.AddCut("50100013", "00216609847002008750404000", "0152501500000000"); // TOF -5, 5
  } else if ( trainConfig == 213){ // pdEdx & TOF - phi
    cuts.AddCut("52400013", "00216609837005008750404000", "0152501500000000"); // 2.5sigma
    cuts.AddCut("52400013", "00216609847002008750404000", "0152501500000000"); // TOF -5, 5
  } else if ( trainConfig == 214){ // qT
    cuts.AddCut("50100013", "00200009847005002750404000", "0152501500000000"); // 0.06 2D
    cuts.AddCut("50100013", "00200009847005009750404000", "0152501500000000"); // 0.03 2D
  } else if ( trainConfig == 215){ // qT
    cuts.AddCut("52400013", "00200009847005002750404000", "0152501500000000"); // 0.06 2D
    cuts.AddCut("52400013", "00200009847005009750404000", "0152501500000000"); // 0.03 2D
  } else if ( trainConfig == 216){ // qT - phi
    cuts.AddCut("50100013", "00216609847005002750404000", "0152501500000000"); // 0.06 2D
    cuts.AddCut("50100013", "00216609847005009750404000", "0152501500000000"); // 0.03 2D
  } else if ( trainConfig == 217){ // qT - phi
    cuts.AddCut("52400013", "00216609847005002750404000", "0152501500000000"); // 0.06 2D
    cuts.AddCut("52400013", "00216609847005009750404000", "0152501500000000"); // 0.03 2D
  } else if ( trainConfig == 218){ // chi2 + psi pair
    cuts.AddCut("50100013", "00200009847005008710404000", "0152501500000000"); // 10+0.1 1D
    cuts.AddCut("50100013", "00200009847005008450404000", "0152501500000000"); // 7+0.1 2D
  } else if ( trainConfig == 219){ // chi2 + psi pair
    cuts.AddCut("52400013", "00200009847005008710404000", "0152501500000000"); // 10+0.1 1D
    cuts.AddCut("52400013", "00200009847005008450404000", "0152501500000000"); // 7+0.1 2D
  } else if ( trainConfig == 220){ // chi2 + psi pair - phi
    cuts.AddCut("50100013", "00216609847005008710404000", "0152501500000000"); // 10+0.1 1D
    cuts.AddCut("50100013", "00216609847005008450404000", "0152501500000000"); // 7+0.1 2D
  } else if ( trainConfig == 221){ // chi2 + psi pair - phi
    cuts.AddCut("52400013", "00216609847005008710404000", "0152501500000000"); // 10+0.1 1D
    cuts.AddCut("52400013", "00216609847005008450404000", "0152501500000000"); // 7+0.1 2D
  } else if ( trainConfig == 222){ // chi2 + psi pair
    cuts.AddCut("50100013", "00200009847005008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCut("50100013", "00200009847005008770404000", "0152501500000000"); // 10+0.07 2D
  } else if ( trainConfig == 223){ // chi2 + psi pair
    cuts.AddCut("52400013", "00200009847005008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCut("52400013", "00200009847005008770404000", "0152501500000000"); // 10+0.07 2D
  } else if ( trainConfig == 224){ // chi2 + psi pair - phi
    cuts.AddCut("50100013", "00216609847005008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCut("50100013", "00216609847005008770404000", "0152501500000000"); // 10+0.07 2D
  } else if ( trainConfig == 225){ // chi2 + psi pair - phi
    cuts.AddCut("52400013", "00216609847005008950404000", "0152501500000000"); // 15+0.1 2D
    cuts.AddCut("52400013", "00216609847005008770404000", "0152501500000000"); // 10+0.07 2D
  } else if ( trainConfig == 226){ // asym & cosPA
    cuts.AddCut("50100013", "00200009847005008757404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCut("50100013", "00200009847005008750004000", "0152501500000000"); // cosPA off
  } else if ( trainConfig == 227){ // asym & cosPA
    cuts.AddCut("52400013", "00200009847005008757404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCut("52400013", "00200009847005008750004000", "0152501500000000"); // cosPA off
  } else if ( trainConfig == 228){ // asym & cosPA - phi
    cuts.AddCut("50100013", "00216609847005008757404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCut("50100013", "00216609847005008750004000", "0152501500000000"); // cosPA off
  } else if ( trainConfig == 229){ // asym & cosPA - phi
    cuts.AddCut("52400013", "00216609847005008757404000", "0152501500000000"); // asym > 8GeV
    cuts.AddCut("52400013", "00216609847005008750004000", "0152501500000000"); // cosPA off
  } else if ( trainConfig == 230){ // asym 1 & 2
    cuts.AddCut("50100013", "00200009847005008756404000", "0152501500000000"); // asym > 6GeV
    cuts.AddCut("50100013", "00200009847005008754404000", "0152501500000000"); // asym 2D
  } else if ( trainConfig == 231){ // asym 1 & 2
    cuts.AddCut("52400013", "00200009847005008756404000", "0152501500000000"); // asym > 6GeV
    cuts.AddCut("52400013", "00200009847005008754404000", "0152501500000000"); // asym 2D
  } else if ( trainConfig == 232){ // asym 1 & 2 - phi
    cuts.AddCut("50100013", "00216609847005008756404000", "0152501500000000"); // asym > 6GeV
    cuts.AddCut("50100013", "00216609847005008754404000", "0152501500000000"); // asym 2D
  } else if ( trainConfig == 233){ // asym 1 & 2 - phi
    cuts.AddCut("52400013", "00216609847005008756404000", "0152501500000000"); // asym > 6GeV
    cuts.AddCut("52400013", "00216609847005008754404000", "0152501500000000"); // asym 2D

  } else if ( trainConfig == 237){ // K on the signal -3, 5
    cuts.AddCut("50100013", "00200009300000008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009300000008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009300000008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 238){ //
    cuts.AddCut("50100013", "00216609300000008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609300000008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609300000008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 239){ // -5, 10
    cuts.AddCut("50100013", "00200009500000008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009500000008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009500000008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 240){ //
    cuts.AddCut("50100013", "00216609500000008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609500000008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609500000008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 241){ // -3, 10
    cuts.AddCut("50100013", "00200009600000008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00200009600000008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00200009600000008250404000", "0152501500000000"); // 20-50%
  } else if ( trainConfig == 242){ //
    cuts.AddCut("50100013", "00216609600000008250404000", "0152501500000000"); // 0-10%
    cuts.AddCut("52400013", "00216609600000008250404000", "0152501500000000"); // 20-40%
    cuts.AddCut("52500013", "00216609600000008250404000", "0152501500000000"); // 20-50%

  } else if (trainConfig == 243){ // LHC15o: 11h photon and meson cuts with kINT7 trigger and PU cuts, V0M
    cuts.AddCut("30110113", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCut("31210113", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCut("12410113", "00200009247602008250404000", "0652501500000000"); // 20-40%
    cuts.AddCut("14610113", "00200009247602008250404000", "0652501500000000"); // 40-60%
    cuts.AddCut("16810113", "00200009247602008250404000", "0652501500000000"); // 60-80%
  } else if (trainConfig == 244){ // LHC15o: 11h photon and meson cuts with kINT7 trigger and PU cuts, V0M
    cuts.AddCut("10110113", "00200009247602008250404000", "0652501500000000"); //  0-10%
    cuts.AddCut("11210113", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCut("12510113", "00200009247602008250404000", "0652501500000000"); // 20-50%
    cuts.AddCut("15910113", "00200009247602008250404000", "0652501500000000"); // 50-90%
    cuts.AddCut("10010113", "00200009247602008250404000", "0652501500000000"); //  0-100%
  } else if (trainConfig == 245){ // LHC15o: 11h photon and meson cuts with kINT7 trigger, V0M
    cuts.AddCut("30110013", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCut("31210013", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCut("12410013", "00200009247602008250404000", "0652501500000000"); // 20-40%
    cuts.AddCut("14610013", "00200009247602008250404000", "0652501500000000"); // 40-60%
    cuts.AddCut("16810013", "00200009247602008250404000", "0652501500000000"); // 60-80%
  } else if (trainConfig == 246){ // LHC15o: 11h photon and meson cuts, with kINT7 trigger, V0M
    cuts.AddCut("10110013", "00200009247602008250404000", "0652501500000000"); //  0-10%
    cuts.AddCut("11210013", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCut("12510013", "00200009247602008250404000", "0652501500000000"); // 20-50%
    cuts.AddCut("15910013", "00200009247602008250404000", "0652501500000000"); // 50-90%
    cuts.AddCut("10010013", "00200009247602008250404000", "0652501500000000"); //  0-100%
  } else if (trainConfig == 247){ // LHC15o: 11h photon and meson cuts, with kINT7 trigger, V0M - user defined header!
    cuts.AddCut("10110023", "00200009247602008250404000", "0652501500000000"); //  0-10%
    cuts.AddCut("11210023", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCut("12510023", "00200009247602008250404000", "0652501500000000"); // 20-50%
    cuts.AddCut("15910023", "00200009247602008250404000", "0652501500000000"); // 50-90%
    cuts.AddCut("10010023", "00200009247602008250404000", "0652501500000000"); //  0-100%
  } else if (trainConfig == 248){ // LHC15o: 11h photon and meson cuts, with kINT7 trigger, V0M - user defined header!
    cuts.AddCut("30110023", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCut("31210023", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCut("12410023", "00200009247602008250404000", "0652501500000000"); // 20-40%
    cuts.AddCut("14610023", "00200009247602008250404000", "0652501500000000"); // 40-60%
    cuts.AddCut("16810023", "00200009247602008250404000", "0652501500000000"); // 60-80%
  } else if (trainConfig == 249){ // LHC15o: 11h photon and meson cuts, with kINT7 trigger and PU cuts, V0M - user defined header!
    cuts.AddCut("10110123", "00200009247602008250404000", "0652501500000000"); //  0-10%
    cuts.AddCut("11210123", "00200009247602008250404000", "0652501500000000"); // 10-20%
    cuts.AddCut("12510123", "00200009247602008250404000", "0652501500000000"); // 20-50%
    cuts.AddCut("15910123", "00200009247602008250404000", "0652501500000000"); // 50-90%
    cuts.AddCut("10010123", "00200009247602008250404000", "0652501500000000"); //  0-100%
  } else if (trainConfig == 250){ // LHC15o: 11h photon and meson cuts, with kINT7 trigger, V0M and PU cuts- user defined header!
    cuts.AddCut("30110123", "00200009247602008250404000", "0652501500000000"); //  0-5%
    cuts.AddCut("31210123", "00200009247602008250404000", "0652501500000000"); //  5-10%
    cuts.AddCut("12410123", "00200009247602008250404000", "0652501500000000"); // 20-40%
    cuts.AddCut("14610123", "00200009247602008250404000", "0652501500000000"); // 40-60%
    cuts.AddCut("16810123", "00200009247602008250404000", "0652501500000000"); // 60-80%

  } else  if (trainConfig == 300){ // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCut("60100013", "03200009300002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200013", "03200009300002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100013", "03200009300002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200013", "03200009300002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200013", "03200009300002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 301) {  // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCut("52400013", "03200009300002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600013", "03200009300002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800013", "03200009300002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800013", "03200009300002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCut("54900013", "03200009300002003220000000", "0152306500900000"); // 40-90%
  } else  if (trainConfig == 302){ // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCut("60100013", "03200009400002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200013", "03200009400002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100013", "03200009400002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200013", "03200009400002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200013", "03200009400002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 303) {  // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCut("52400013", "03200009400002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600013", "03200009400002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800013", "03200009400002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800013", "03200009400002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCut("54900013", "03200009400002003220000000", "0152306500900000"); // 40-90%
  } else  if (trainConfig == 304){ // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCut("60100013", "03200009500002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200013", "03200009500002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100013", "03200009500002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200013", "03200009500002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200013", "03200009500002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 305) {  // LHC10h standard, eta 0.65, y = 0.6
    cuts.AddCut("52400013", "03200009500002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600013", "03200009500002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800013", "03200009500002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800013", "03200009500002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCut("54900013", "03200009500002003220000000", "0152306500900000"); // 40-90%
  } else  if (trainConfig == 306){ // LHC10h standard, eta 0.9, y = 0.6
    cuts.AddCut("60100013", "00200009300002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200013", "00200009300002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100013", "00200009300002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200013", "00200009300002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200013", "00200009300002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 307) {  // LHC10h standard, eta 0.9, y = 0.6
    cuts.AddCut("52400013", "00200009300002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600013", "00200009300002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800013", "00200009300002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800013", "00200009300002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCut("54900013", "00200009300002003220000000", "0152306500900000"); // 40-90%
  } else  if (trainConfig == 308){ // LHC10h standard, eta 0.9, y = 0.6
    cuts.AddCut("60100013", "00200009400002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200013", "00200009400002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100013", "00200009400002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200013", "00200009400002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200013", "00200009400002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 309) { // LHC10h standard, eta 0.9, y = 0.6
    cuts.AddCut("52400013", "00200009400002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600013", "00200009400002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800013", "00200009400002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800013", "00200009400002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCut("54900013", "00200009400002003220000000", "0152306500900000"); // 40-90%
  } else  if (trainConfig == 310){ // LHC10h standard, eta 0.9, y = 0.6
    cuts.AddCut("60100013", "00200009500002003220000000", "0152304500900000"); // 0-5%
    cuts.AddCut("61200013", "00200009500002003220000000", "0152304500900000"); // 5-10%
    cuts.AddCut("50100013", "00200009500002003220000000", "0152304500900000"); // 0-10%
    cuts.AddCut("51200013", "00200009500002003220000000", "0152304500900000"); // 10-20%
    cuts.AddCut("50200013", "00200009500002003220000000", "0152304500900000"); // 0-20%
  } else if (trainConfig == 311) {  // LHC10h standard, eta 0.9, y = 0.6
    cuts.AddCut("52400013", "00200009500002003220000000", "0152304500900000"); // 20-40%
    cuts.AddCut("54600013", "00200009500002003220000000", "0152306500900000"); // 40-60%
    cuts.AddCut("56800013", "00200009500002003220000000", "0152306500900000"); // 60-80%
    cuts.AddCut("54800013", "00200009500002003220000000", "0152306500900000"); // 40-80%
    cuts.AddCut("54900013", "00200009500002003220000000", "0152306500900000"); // 40-90%
  } else  if (trainConfig == 312){ // LHC10h standard direct photon flow cuts
    cuts.AddCut("50200013", "00200009307000008250400000", "0152304500900000");
    cuts.AddCut("52400013", "00200009307000008250400000", "0152304500900000");
    cuts.AddCut("54800013", "00200009307000008250400000", "0152304500900000");

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
  if (periodName.CompareTo("LHC13d2")==0){
    TObjString *Header1 = new TObjString("pi0_1");
    HeaderList->Add(Header1);
//    TObjString *Header3 = new TObjString("eta_2");
//    HeaderList->Add(Header3);

  } else if (periodName.CompareTo("LHC12a17x_fix")==0){
    TObjString *Header1 = new TObjString("PARAM");
    HeaderList->Add(Header1);
  } else if (periodName.CompareTo("LHC14a1a")==0){
    if (headerSelectionInt == 1){
      TObjString *Header1 = new TObjString("pi0_1");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 2){
      TObjString *Header1 = new TObjString("eta_2");
      HeaderList->Add(Header1);
    }else {
      TObjString *Header1 = new TObjString("pi0_1");
      HeaderList->Add(Header1);
      TObjString *Header2 = new TObjString("eta_2");
      HeaderList->Add(Header2);
    }
  } else if (periodName.CompareTo("LHC14a1b")==0 || periodName.CompareTo("LHC14a1c")==0){
    TObjString *Header1 = new TObjString("BOX");
    HeaderList->Add(Header1);
  } else if (periodName.CompareTo("LHC16h4")==0){
    if (headerSelectionInt == 1){
      TObjString *Header1 = new TObjString("Injector (pi0)_1");
      HeaderList->Add(Header1);
    } else if (headerSelectionInt == 2){
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

    if(periodName.Contains("LHC11h") && doFlattening){
      cout << "entering the cent. flattening loop -> searching for file: " << fileNameInputForCentFlattening.Data() << endl;

      if( fileNameInputForCentFlattening.Contains("FlatFile") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "Cent");
      } else if( fileNameInputForCentFlattening.Contains("Good") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentGoodRuns");
      }else if( fileNameInputForCentFlattening.Contains("SemiGood") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentSemiGoodRuns");
      }else {
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentTotalRuns");
      }
    }

    TString dataInputMultHisto  = "";
    TString mcInputMultHisto    = "";
    if (doMultiplicityWeighting){
      cout << "enableling mult weighting" << endl;

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

      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile(kTRUE, fileNameInputForMultWeighing, dataInputMultHisto, mcInputMultHisto );

    }

    if (  trainConfig == 1   || trainConfig == 5   || trainConfig == 9   || trainConfig == 13   || trainConfig == 17   ||
        trainConfig == 21   || trainConfig == 25   || trainConfig == 29   || trainConfig == 33   || trainConfig == 37  ||
        trainConfig == 300 || trainConfig == 302 || trainConfig == 304){
      if (i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
      if (i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
      if (i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
      if (i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      if (i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
    } else if (  trainConfig == 2   || trainConfig == 6   || trainConfig == 10   || trainConfig == 14   || trainConfig == 18   ||
          trainConfig == 22   || trainConfig == 26   || trainConfig == 30   || trainConfig == 34   || trainConfig == 38  ||
          trainConfig == 301 || trainConfig == 303 || trainConfig == 305){
      if (i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      if (i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
      if (i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
      if (i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
    } else if ( trainConfig == 3   || trainConfig == 7    || trainConfig == 11   || trainConfig == 15  || trainConfig == 19   ||
          trainConfig == 23   || trainConfig == 27   || trainConfig == 31   || trainConfig == 35   || trainConfig == 39){
      if (i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
      if (i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
      if (i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
      if (i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      if (i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
    } else if (  trainConfig == 4   ||trainConfig == 8     || trainConfig == 12   || trainConfig == 16   || trainConfig == 20   ||
          trainConfig == 24   || trainConfig == 28   || trainConfig == 32   || trainConfig == 36   || trainConfig == 40){
      if (i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      if (i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
      if (i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
      if (i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
    }

    if (trainConfig == 56 ){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_1020TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_1020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
        if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0020TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M","Eta_Fit_Data_PbPb_2760GeV_0020V0M");
      }
    }
    if (trainConfig == 57 ){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_1020TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_1020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
        if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0020TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0020TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M","Eta_Fit_Data_PbPb_2760GeV_0020V0M");
      }
    }
    if (trainConfig == 58 ){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_4060TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M","Eta_Fit_Data_PbPb_2760GeV_4060V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_6080TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_6080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M","Eta_Fit_Data_PbPb_2760GeV_6080V0M");
        if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_4080TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
        if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_3050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_3050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3050V0M","Eta_Fit_Data_PbPb_2760GeV_3050V0M");
      }
    }
    if (trainConfig == 59 ){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4060TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M","Eta_Fit_Data_PbPb_2760GeV_4060V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_6080TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_6080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M","Eta_Fit_Data_PbPb_2760GeV_6080V0M");
        if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4080TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4080TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M","Eta_Fit_Data_PbPb_2760GeV_1020V0M");
        if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_3050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_3050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3050V0M","Eta_Fit_Data_PbPb_2760GeV_3050V0M");
      }
    }
    if (trainConfig == 60 ){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2030TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2030TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2030V0M","Eta_Fit_Data_PbPb_2760GeV_2030V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_3040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_3040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3040V0M","Eta_Fit_Data_PbPb_2760GeV_3040V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_4050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_4050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4050V0M","Eta_Fit_Data_PbPb_2760GeV_4050V0M");
        if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_5060TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_5060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_5060V0M","Eta_Fit_Data_PbPb_2760GeV_5060V0M");
        if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }
    if (trainConfig == 61 ){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2030TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2030TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2030V0M","Eta_Fit_Data_PbPb_2760GeV_2030V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_3040TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_3040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_3040V0M","Eta_Fit_Data_PbPb_2760GeV_3040V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_4050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_4050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_4050V0M","Eta_Fit_Data_PbPb_2760GeV_4050V0M");
        if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_5060TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_5060TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_5060V0M","Eta_Fit_Data_PbPb_2760GeV_5060V0M");
        if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
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
        trainConfig == 170   || trainConfig == 172   || trainConfig == 174   || trainConfig == 182  || trainConfig == 184  ){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
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
        trainConfig == 171   || trainConfig == 173   || trainConfig == 175   || trainConfig == 183  || trainConfig == 185  ){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (trainConfig == 186   || trainConfig == 187){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }
    if (trainConfig == 188   || trainConfig == 189){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
      }
    }
    if (trainConfig == 190   || trainConfig == 191){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    if (trainConfig == 192 || trainConfig == 194 || trainConfig == 196 || trainConfig == 198 || trainConfig == 200 || trainConfig == 202 ||
        trainConfig == 204 || trainConfig == 206 || trainConfig == 208 || trainConfig == 210 || trainConfig == 212 || trainConfig == 214 ||
        trainConfig == 216 || trainConfig == 218 || trainConfig == 220 || trainConfig == 222 || trainConfig == 224 || trainConfig == 226 ||
        trainConfig == 228 || trainConfig == 230 || trainConfig == 232){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }
    if (trainConfig == 193 || trainConfig == 195 || trainConfig == 197 || trainConfig == 199 || trainConfig == 201 || trainConfig == 203 ||
        trainConfig == 205 || trainConfig == 207 || trainConfig == 209 || trainConfig == 211 || trainConfig == 213 || trainConfig == 215 ||
        trainConfig == 217 || trainConfig == 219 || trainConfig == 221 || trainConfig == 223 || trainConfig == 225 || trainConfig == 227 ||
        trainConfig == 229 || trainConfig == 231 || trainConfig == 233){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
      }
    }

    if (trainConfig == 237   || trainConfig == 238   || trainConfig == 239   || trainConfig == 240  || trainConfig == 241  || trainConfig == 242){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());
    if (periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
      if (headerSelectionInt == 1) analysisEventCuts[i]->SetAddedSignalPDGCode(111);
      if (headerSelectionInt == 2) analysisEventCuts[i]->SetAddedSignalPDGCode(221);
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
        trainConfig == 306 || trainConfig == 307 || trainConfig == 308  || trainConfig == 309  || trainConfig == 310  || trainConfig == 311 || trainConfig == 312 ||
        trainConfig == 237 || trainConfig == 238 || trainConfig == 239  || trainConfig == 240  || trainConfig == 241  || trainConfig == 242)
      analysisCuts[i]->SetSwitchToKappaInsteadOfNSigdEdxTPC(kTRUE);
    if(trainConfig == 182 || trainConfig == 183 || trainConfig == 184 || trainConfig == 185)
      analysisCuts[i]->SetDodEdxSigmaCut(kFALSE);

    if (enableMatBudWeightsPi0 > 0){
        if (isMC > 0){
            if (analysisCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,filenameMatBudWeights)){
                initializedMatBudWeigths_existing = kTRUE;}
            else {cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;}
        }
        else {cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;}
    }

    analysisCuts[i]->SetIsHeavyIon(isHeavyIon);
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetLightOutput(runLightOutput);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());

    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(runLightOutput);
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
  task->SetDoTHnSparse(enableUseTHnSparse);
  task->SetDoCentFlattening(doFlattening);
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

  return;

}
