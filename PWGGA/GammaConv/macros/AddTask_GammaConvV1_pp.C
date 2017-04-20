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
//pp together with all supporting classes
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
      eventCutArray = new TString[nMaxCuts]; photonCutArray = new TString[nMaxCuts]; mesonCutArray = new TString[nMaxCuts]; clusterCutArray = new TString[nMaxCuts];
      for(Int_t i=0; i<nMaxCuts; i++) {eventCutArray[i] = ""; photonCutArray[i] = ""; mesonCutArray[i] = ""; clusterCutArray[i] = "";}
    }

    void AddCut(TString eventCut, TString photonCut, TString mesonCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConv: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 ) {cout << "ERROR in CutHandlerConv: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut; clusterCutArray[nCuts]="";
      nCuts++;
      return;
    }
    void AddCut(TString eventCut, TString photonCut, TString mesonCut, TString clusterCut){
      if(nCuts>=nMaxCuts) {cout << "ERROR in CutHandlerConv: Exceeded maximum number of cuts!" << endl; validCuts = false; return;}
      if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 || clusterCut.Length()!=19 ) {cout << "ERROR in CutHandlerConv: Incorrect length of cut string!" << endl; validCuts = false; return;}
      eventCutArray[nCuts]=eventCut; photonCutArray[nCuts]=photonCut; mesonCutArray[nCuts]=mesonCut; clusterCutArray[nCuts]=clusterCut;
      nCuts++;
      return;
    }

    Bool_t AreValid(){return validCuts;}
    Int_t GetNCuts(){if(validCuts) return nCuts; else return 0;}
    TString GetEventCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return eventCutArray[i]; else{cout << "ERROR in CutHandlerConv: GetEventCut wrong index i" << endl;return "";}}
    TString GetPhotonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return photonCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetPhotonCut wrong index i" << endl;return "";}}
    TString GetMesonCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return mesonCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetMesonCut wrong index i" << endl;return "";}}
    TString GetClusterCut(Int_t i){if(validCuts&&i<nMaxCuts&&i>=0) return clusterCutArray[i]; else {cout << "ERROR in CutHandlerConv: GetClusterCut wrong index i" << endl;return "";}}
  private:
    Bool_t validCuts;
    Int_t nCuts; Int_t nMaxCuts;
    TString* eventCutArray;
    TString* photonCutArray;
    TString* mesonCutArray;
    TString* clusterCutArray;
};

void AddTask_GammaConvV1_pp(  Int_t   trainConfig                     = 1,                    // change different set of cuts
                              Int_t   isMC                            = 0,                    // run MC
                              Int_t   enableQAMesonTask               = 0,                    // enable meson QA in AliAnalysisTaskGammaConvV1
                              Int_t   enableQAPhotonTask              = 0,                    // enable photon QA in AliAnalysisTaskGammaConvV1
                              TString fileNameInputForPartWeighting   = "MCSpectraInput.root",// path to file for weigting input
                              TString cutnumberAODBranch              = "000000006008400001001500000",  // cutnumber for AOD branch
                              TString periodname                      = "LHC12f1x",           // period name
                              Bool_t  doParticleWeighting             = kFALSE,               // enables weighting
                              Bool_t  enableV0findingEffi             = kFALSE,               // enables V0finding efficiency histograms
                              Bool_t  enableTriggerMimicking          = kFALSE,               // enable trigger mimicking
                              Bool_t  enableTriggerOverlapRej         = kFALSE,               // enable trigger overlap rejection
                              Float_t maxFacPtHard                    = 3.,                   // maximum factor between hardest jet and ptHard generated
                              TString periodNameV0Reader              = "",
                              Bool_t  doMultiplicityWeighting         = kFALSE,                  //
                              TString fileNameInputForMultWeighing    = "Multiplicity.root",    //
                              TString periodNameAnchor                = "",
                              Bool_t  runLightOutput                  = kFALSE,                // switch to run light output (only essential histograms for afterburner)
                              Bool_t  enableChargedPrimary            = kFALSE,
                              Bool_t    doSmear                       = kFALSE,                 // switches to run user defined smearing
                              Double_t  bremSmear                     = 1.,
                              Double_t  smearPar                      = 0.,                     // conv photon smearing params
                              Double_t  smearParConst                 = 0.,                      // conv photon smearing params
                              Int_t   enableMatBudWeightsPi0          = 0,                      // 1 = three radial bins, 2 = 10 radial bins
                              TString filenameMatBudWeights           = "MCInputFileMaterialBudgetWeights.root",
                              TString   additionalTrainConfig         = "0"                     // additional counter for trainconfig, this has to be always the last parameter
                            ) {

  Int_t isHeavyIon = 0; 
  //parse additionalTrainConfig flag
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){cout << "ERROR: AddTask_GammaConvV1_pp during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl; return;}
  TObjString* rAdditionalTrainConfig;
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(i==0){ rAdditionalTrainConfig = (TObjString*)rAddConfigArr->At(i);
    } else {
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      cout<< tempStr.Data()<<endl;
      
      if(tempStr.Contains("MaterialBudgetWeights") && enableMatBudWeightsPi0 > 0){
         if(tempStr.Contains("MaterialBudgetWeightsNONE")){
            enableMatBudWeightsPi0 = 0;
            cout << "INFO:  AddTask_GammaConvV1_pp materialBudgetWeights switched off signaled by additionalTrainConfigFlag" << endl;
         } else {
            TObjArray *fileNameMatBudWeightsArr = filenameMatBudWeights.Tokenize("/");
            if(fileNameMatBudWeightsArr->GetEntries()<1 ){
                cout<<"ERROR: AddTask_GammaConvV1_pp when reading material budget weights file name" << filenameMatBudWeights.Data()<< "'" << endl; 
                return;
            }  
            TObjString * oldMatObjStr = (TObjString*)fileNameMatBudWeightsArr->At( fileNameMatBudWeightsArr->GetEntries()-1);
            TString  oldfileName  = oldMatObjStr->GetString();
            TString  newFileName  = Form("MCInputFile%s.root",tempStr.Data());
            cout<<newFileName.Data()<<endl;
            if( oldfileName.EqualTo(newFileName.Data()) == 0 ){
              filenameMatBudWeights.ReplaceAll(oldfileName.Data(),newFileName.Data());
              cout << "INFO: AddTask_GammaConvV1_pp the material budget weights file has been change to " <<filenameMatBudWeights.Data()<<"'"<< endl;
          }
        }
      }
    }
  }
  
  TString sAdditionalTrainConfig = rAdditionalTrainConfig->GetString();
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
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
  TString cutnumberPhoton = "00200008400000002200000000";
  TString cutnumberEvent = "00000003"; 
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
    if(trainConfig == 101){
      fV0ReaderV1->SetProduceImpactParamHistograms(kTRUE);
      cutnumberPhoton="00200009227302008250400400";
      cutnumberEvent ="00010113";
    }
    if(trainConfig>=60 && trainConfig<80) fV0ReaderV1->SetImprovedPsiPair(0); //switch off for 8TeV as AODs are used for which improved psipair is not available

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
      if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
        fV0ReaderV1->SetConversionCuts(fCuts);
        fCuts->SetFillCutHistograms("",kTRUE);
      }
      if( trainConfig  == 73 || trainConfig  == 74 || (trainConfig  >= 80 && trainConfig  <= 87)   ){
        fCuts->SetDodEdxSigmaCut(kFALSE);
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
  //            find input container
  AliAnalysisTaskGammaConvV1 *task=NULL;
  task= new AliAnalysisTaskGammaConvV1(Form("GammaConvV1_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  task->SetV0ReaderName(V0ReaderName);
  task->SetLightOutput(runLightOutput);
  // Cut Numbers to use in Analysis
  
  CutHandlerConv cuts;
  
  //---------  standard configurations for 2.76TeV V00R without SDD --------------------------
  if(trainConfig == 1){    // various standard cuts
    cuts.AddCut("00000113", "00200009397300008250400000", "0163103100900000"); // new pi0/eta cut 2.76TeV 
    cuts.AddCut("00000113", "00200009397300008250400000", "0163103100000000"); // new pi0/eta cut 2.76TeV without MC smearing
    cuts.AddCut("00000113", "00200009366300003800000000", "0163103100900000"); // standard cut Pi0 pp 2.76TeV PbPb paper 2012
  } else if (trainConfig == 2) { // various standard cuts added signals
    cuts.AddCut("00000123", "00200009397300008250400000", "0163103100900000"); // new pi0/eta cut 2.76TeV 
    cuts.AddCut("00000123", "00200009397300008250400000", "0163103100000000"); // new pi0/eta cut 2.76TeV without MC smearing
    cuts.AddCut("00000123", "00200009366300003800000000", "0163103100900000"); // standard cut Pi0 pp 2.76TeV PbPb paper 2012
  } else if (trainConfig == 3) { // additional standards
    cuts.AddCut("00000113", "00200009297002008250400000", "0163103100900000"); // standard cut LHC11h pp 2.76TeV 
    cuts.AddCut("00000113", "00200009227302008250404000", "0163101500000000"); // Ana eta analysis prefered
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); // standard cut Pi0 pp 7TeV, all photon qualities
  } else if (trainConfig == 4) { // additional standards added signals
    cuts.AddCut("00000123", "00200009297002008250400000", "0163103100900000"); // standard cut LHC11h pp 2.76TeV 
    cuts.AddCut("00000123", "00200009227302008250404000", "0163101500000000"); // Ana eta analysis prefered
    cuts.AddCut("00000123", "00200008366300000200000000", "0163103100900000"); // standard cut Pi0 pp 7TeV, all photon qualities
  } else if (trainConfig == 5) { 
    // variations to different standards  
    cuts.AddCut("00000113", "00200009327000008250400000", "0163103100900000"); // go with 1sigma pi rejec to infty
    cuts.AddCut("00000113", "00200009317000008250400000", "0163103100900000"); // go with 0sigma pi rejec to infty
    cuts.AddCut("00000113", "00200009357000008250400000", "0163103100900000"); // go with 2sigma pi reject to infy
    // eta cut variation    
    cuts.AddCut("00000113", "03200009397300008250400000", "0163103100900000"); // eta 0.65
    cuts.AddCut("00000113", "04200009397300008250400000", "0163103100900000"); // eta 0.75
  } else if (trainConfig == 6) { 
    // variations to different standards added signals
    cuts.AddCut("00000123", "00200009327000008250400000", "0163103100900000"); // go with 1sigma pi rejec to infty
    cuts.AddCut("00000123", "00200009317000008250400000", "0163103100900000"); // go with 0sigma pi rejec to infty
    cuts.AddCut("00000123", "00200009357000008250400000", "0163103100900000"); // go with 2sigma pi reject to infy
    // eta cut variation added signals
    cuts.AddCut("00000123", "03200009397300008250400000", "0163103100900000"); // eta 0.65
    cuts.AddCut("00000123", "04200009397300008250400000", "0163103100900000"); // eta 0.75
  } else if (trainConfig == 7) {  //pion rejection variations
    cuts.AddCut("00000113", "00200009395300008250400000", "0163103100900000"); // min for pi dEdx 0.3GeV
    cuts.AddCut("00000113", "00200009390300008250400000", "0163103100900000"); // min for pi dEdx  0.5GeV
    cuts.AddCut("00000113", "00200009397400008250400000", "0163103100900000"); // max for pi dEdx  3 GeV
    cuts.AddCut("00000113", "00200009397200008250400000", "0163103100900000"); // new pi0/eta cut 4 GeV
  } else if (trainConfig == 8) { //pion rejection variations added signals
    cuts.AddCut("00000123", "00200009395300008250400000", "0163103100900000"); // min for pi dEdx 0.3GeV
    cuts.AddCut("00000123", "00200009390300008250400000", "0163103100900000"); // min for pi dEdx  0.5GeV
    cuts.AddCut("00000123", "00200009397400008250400000", "0163103100900000"); // max for pi dEdx  3 GeV
    cuts.AddCut("00000123", "00200009397200008250400000", "0163103100900000"); // new pi0/eta cut 4 GeV
  } else if (trainConfig == 9) {  //electron rejection variations
    cuts.AddCut("00000113", "00200009197300008250400000", "0163103100900000"); // dEdx e +-5 sigma 
    cuts.AddCut("00000113", "00200009497300008250400000", "0163103100900000"); // dEdx e -6,+7 sigma
    cuts.AddCut("00000113", "00200009297300008250400000", "0163103100900000"); // dEdx e -3,+5 sigma
    cuts.AddCut("00000113", "00200009597300008250400000", "0163103100900000"); // dEdx e +-4 sigma
  } else if (trainConfig == 10) { //electron rejection variations added signals
    cuts.AddCut("00000123", "00200009197300008250400000", "0163103100900000"); // dEdx e +-5 sigma 
    cuts.AddCut("00000123", "00200009497300008250400000", "0163103100900000"); // dEdx e -6,+7 sigma
    cuts.AddCut("00000123", "00200009297300008250400000", "0163103100900000"); // dEdx e -3,+5 sigma
    cuts.AddCut("00000123", "00200009597300008250400000", "0163103100900000"); // dEdx e +-4 sigma
  } else if (trainConfig == 11) { // single leg cuts
    cuts.AddCut("00000113", "00200049397300008250400000", "0163103100900000"); // variation pt 0.075 
    cuts.AddCut("00000113", "00200019397300008250400000", "0163103100900000"); // variation pt 0.1
    cuts.AddCut("00000113", "00200006397300008250400000", "0163103100900000"); // variation TPC cls 0.7
    cuts.AddCut("00000113", "00200008397300008250400000", "0163103100900000"); // variation TPC cls 0.35 
  } else if (trainConfig == 12) { // single leg cuts added signals
    cuts.AddCut("00000123", "00200049397300008250400000", "0163103100900000"); // variation pt 0.075 
    cuts.AddCut("00000123", "00200019397300008250400000", "0163103100900000"); // variation pt 0.1
    cuts.AddCut("00000123", "00200006397300008250400000", "0163103100900000"); // variation TPC cls 0.7
    cuts.AddCut("00000123", "00200008397300008250400000", "0163103100900000"); // variation TPC cls 0.35 
  } else if (trainConfig == 13) { // Qt variations 
    cuts.AddCut("00000113", "00200009397300009250400000", "0163103100900000"); // variation qt 0.03
    cuts.AddCut("00000113", "00200009397300002250400000", "0163103100900000"); // variation qt 0.06 
    cuts.AddCut("00000113", "00200009397300003250400000", "0163103100900000"); // variation qt 0.05 no 2D    
    cuts.AddCut("00000113", "00200009397300008250400000", "0163105100900000"); // tighter alpha meson 0.75
  } else if (trainConfig == 14) { // Qt variations added signals   
    cuts.AddCut("00000123", "00200009397300009250400000", "0163103100900000"); // variation qt 0.03
    cuts.AddCut("00000123", "00200009397300002250400000", "0163103100900000"); // variation qt 0.06 
    cuts.AddCut("00000123", "00200009397300003250400000", "0163103100900000"); // variation qt 0.05 no 2D
    cuts.AddCut("00000123", "00200009397300008250400000", "0163105100900000"); // tighter alpha meson 0.75
  } else if (trainConfig == 15) { // chi2 - Psi pair variations
    cuts.AddCut("00000113", "00200009397300008150400000", "0163103100900000"); // chi2 50 with psi pair 0.1
    cuts.AddCut("00000113", "00200009397300008850400000", "0163103100900000"); // chi2 20 with psi pair 0.1
    cuts.AddCut("00000113", "00200009397300008280400000", "0163103100900000"); // chi2 30 with psi pair 0.2
    cuts.AddCut("00000113", "00200009397300008260400000", "0163103100900000"); // chi2 30 with psi pair 0.05
  } else if (trainConfig == 16) { // chi2 - Psi pair variations
    cuts.AddCut("00000123", "00200009397300008150400000", "0163103100900000"); // chi2 50 with psi pair 0.1
    cuts.AddCut("00000123", "00200009397300008850400000", "0163103100900000"); // chi2 20 with psi pair 0.1
    cuts.AddCut("00000123", "00200009397300008280400000", "0163103100900000"); // chi2 30 with psi pair 0.2
    cuts.AddCut("00000123", "00200009397300008260400000", "0163103100900000"); // chi2 30 with psi pair 0.05
  } else if (trainConfig == 17) { // chi2 - Psi pair variations
    cuts.AddCut("00000113", "00200009397300008180400000", "0163103100900000"); // chi2 50 with psi pair 0.2
    cuts.AddCut("00000113", "00200009397300008860400000", "0163103100900000"); // chi2 20 with psi pair 0.05
    cuts.AddCut("00000113", "00200009397300008000400000", "0163103100900000"); // chi2 100, no psi pair
  } else if (trainConfig == 18) { // chi2 - Psi pair variations    
    cuts.AddCut("00000123", "00200009397300008180400000", "0163103100900000"); // chi2 50 with psi pair 0.2
    cuts.AddCut("00000123", "00200009397300008860400000", "0163103100900000"); // chi2 20 with psi pair 0.05
    cuts.AddCut("00000123", "00200009397300008000400000", "0163103100900000"); // chi2 100, no psi pair
  } else if (trainConfig == 19) { // photon quality
    cuts.AddCut("00000113", "00200009397300008250420000", "0163103100900000"); // photon quality 1
    cuts.AddCut("00000113", "00200009397300008250430000", "0163103100900000"); // photon quality 2
    cuts.AddCut("00000113", "00200009397300008250440000", "0163103100900000"); // photon quality 3
    cuts.AddCut("00000113", "00200009397300008250400000", "0163106100900000"); // tighter alpha meson 0.85
  } else if (trainConfig == 20) { // photon quality added signals
    cuts.AddCut("00000123", "00200009397300008250420000", "0163103100900000"); // photon quality 1
    cuts.AddCut("00000123", "00200009397300008250430000", "0163103100900000"); // photon quality 2
    cuts.AddCut("00000123", "00200009397300008250440000", "0163103100900000"); // photon quality 3
    cuts.AddCut("00000123", "00200009397300008250400000", "0163106100900000"); // tighter alpha meson 0.85
  } else if (trainConfig == 21) { // much stricter min R cut
    cuts.AddCut("00000113", "00700009397300008250400000", "0163103100900000"); // min R = 35 cm
    cuts.AddCut("00000113", "00700009397300008250420000", "0163103100900000"); // photon quality 1, min R = 35 cm
    cuts.AddCut("00000113", "00700009397300008250430000", "0163103100900000"); // photon quality 2, min R = 35 cm
    cuts.AddCut("00000113", "00700009397300008250440000", "0163103100900000"); // photon quality 3, min R = 35 cm
  } else if (trainConfig == 22) { // much stricter min R cut added signals
    cuts.AddCut("00000123", "00700009397300008250400000", "0163103100900000"); // min R = 35 cm
    cuts.AddCut("00000123", "00700009397300008250420000", "0163103100900000"); // photon quality 1, min R = 35 cm
    cuts.AddCut("00000123", "00700009397300008250430000", "0163103100900000"); // photon quality 2, min R = 35 cm
    cuts.AddCut("00000123", "00700009397300008250440000", "0163103100900000"); // photon quality 3, min R = 35 cm
  } else if (trainConfig == 23) { // meson cut variations
    cuts.AddCut("00000113", "00200009397300008250400000", "0163203100900000"); // y meson < 0.7
    cuts.AddCut("00000113", "00200009397300008250400000", "0163503100900000"); // y meson < 0.85
    cuts.AddCut("00000113", "00200009397300008250400000", "0263203100900000"); // mixed event with track mult
    cuts.AddCut("00000113", "00200009397300008250400000", "0063503100900000"); // BG with rotation
  } else if (trainConfig == 24) { // meson cut variations added signals
    cuts.AddCut("00000123", "00200009397300008250400000", "0163203100900000"); // y meson < 0.7
    cuts.AddCut("00000123", "00200009397300008250400000", "0163503100900000"); // y meson < 0.85
    cuts.AddCut("00000123", "00200009397300008250400000", "0263203100900000"); // mixed event with track mult
    cuts.AddCut("00000123", "00200009397300008250400000", "0063503100900000"); // BG with rotation

  //--------- testing triggers for EMC in 2.76TeV ----------------------------------------------
  } else if (trainConfig == 25) { //LHC11a
    cuts.AddCut("00003113", "00200009397300008250400000", "0163103100000000","1111121057032220000"); //WSDD test with triggers INT1
    cuts.AddCut("00003113", "00202209397300008250400000", "0163103100000000","1111121057032220000"); //WSDD test with triggers INT1 restricted wide EMC range
    cuts.AddCut("00003113", "00204409397300008250400000", "0163103100000000","1111121057032220000"); //WSDD test with triggers INT1 restricted tight EMC range
    cuts.AddCut("00051113", "00200009397300008250400000", "0163103100000000","1111121057032220000"); //WSDD test with triggers EMC1
    cuts.AddCut("00051113", "00202209397300008250400000", "0163103100000000","1111121057032220000"); //WSDD test with triggers EMC1 restricted wide EMC range
    cuts.AddCut("00051113", "00204409397300008250400000", "0163103100000000","1111121057032220000"); //WSDD test with triggers EMC1 restricted tight EMC range
  } else if (trainConfig == 26) { //LHC13g full acceptance
    cuts.AddCut("00010113", "00200009397300008250400000", "0163103100000000","1111121067032220000"); //INT7
    cuts.AddCut("00052113", "00200009397300008250400000", "0163103100000000","1111121067032220000"); //EMC7
    cuts.AddCut("00085113", "00200009397300008250400000", "0163103100000000","1111121067032220000"); //EG2
    cuts.AddCut("00083113", "00200009397300008250400000", "0163103100000000","1111121067032220000"); //EG1
  } else if (trainConfig == 27) { //LHC13g loose EMC acceptance for PCM photons
    cuts.AddCut("00010113", "00202209397300008250400000", "0163103100000000","1111121067032220000"); //INT7
    cuts.AddCut("00052113", "00202209397300008250400000", "0163103100000000","1111121067032220000"); //EMC7
    cuts.AddCut("00085113", "00202209397300008250400000", "0163103100000000","1111121067032220000"); //EG2
    cuts.AddCut("00083113", "00202209397300008250400000", "0163103100000000","1111121067032220000"); //EG1
  } else if (trainConfig == 28) { //LHC13g tight EMC acceptance for PCM photons
    cuts.AddCut("00010113", "00204409397300008250400000", "0163103100000000","1111121067032220000"); //INT7
    cuts.AddCut("00052113", "00204409397300008250400000", "0163103100000000","1111121067032220000"); //EMC7
    cuts.AddCut("00085113", "00204409397300008250400000", "0163103100000000","1111121067032220000"); //EG2
    cuts.AddCut("00083113", "00204409397300008250400000", "0163103100000000","1111121067032220000"); //EG1

  //---------  standard configurations for 2.76TeV V00R with SDD -----------------------------
  } else if(trainConfig == 30){    // various standard cuts
    cuts.AddCut("00003113", "00200009397300008250400000", "0163103100900000"); // new pi0/eta cut 2.76TeV 
    cuts.AddCut("00003113", "00200009397300008250400000", "0163103100000000"); // new pi0/eta cut 2.76TeV without MC smearing
    cuts.AddCut("00003113", "00200009366300003800000000", "0163103100900000"); // standard cut Pi0 pp 2.76TeV PbPb paper 2012
  } else if (trainConfig == 31) {
     cuts.AddCut("00003113", "00200009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD , only Minbias MC
  } else if (trainConfig == 32) {
     cuts.AddCut("00003123", "00200009366300003800000000", "0163103100900000"); //standard cut Pi0 pp 2.76TeV with SDD , only Boxes MC
    

  //--------- Ana marin: variations for eta reanlysis 2.76TeV 2015 ----------------------------
  } else if (trainConfig == 40) {
    cuts.AddCut("00000113", "00200009217302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCut("00000113", "00200009217302008250404000", "0152106500000000"); // alpha variation 0.8
    cuts.AddCut("00000113", "00200009217302008250404000", "0152103500000000"); // alpha variation  1
    cuts.AddCut("00000113", "00200049217302008250404000", "0152101500000000"); // single pT cut 0.075
  } else if (trainConfig == 41) {
    cuts.AddCut("00000123", "00200009217302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCut("00000123", "00200009217302008250404000", "0152106500000000"); // alpha variation  0.8
    cuts.AddCut("00000123", "00200009217302008250404000", "0152103500000000"); // alpha variation  1
    cuts.AddCut("00000123", "00200049217302008250404000", "0152101500000000"); // single pT cut 0.075
  } else if (trainConfig == 42) {
    cuts.AddCut("00000113", "00200019217302008250404000", "0152101500000000"); // single pT cut 0.1
    cuts.AddCut("00000113", "00200008217302008250404000", "0152101500000000"); // TPC cls 0.35
    cuts.AddCut("00000113", "00200006217302008250404000", "0152101500000000"); // TPC cls 0.7
    cuts.AddCut("00000113", "00200009217302009250404000", "0152101500000000"); // qT cut 0.03 2D
  } else if (trainConfig == 43) {
    cuts.AddCut("00000123", "00200019217302008250404000", "0152101500000000"); // single pT cut 0.1
    cuts.AddCut("00000123", "00200008217302008250404000", "0152101500000000"); // TPC cls 0.35
    cuts.AddCut("00000123", "00200006217302008250404000", "0152101500000000"); // TPC cls 0.7
    cuts.AddCut("00000123", "00200009217302009250404000", "0152101500000000"); // qT cut 0.03 2D
  } else if (trainConfig == 44) {
    cuts.AddCut("00000113", "00200009217302008210404000", "0152101500000000"); // variation chi2 30 psi pair 0.1 1D
    cuts.AddCut("00000113", "00200009217302008180404000", "0152101500000000"); // variation chi2 50 psi pair 0.2 2D
    cuts.AddCut("00000113", "00200009217302008860404000", "0152101500000000"); // variation chi2 20 psi pair 0.05 2D
    cuts.AddCut("00000113", "00200009217302009250404000", "0252101500000000"); // variation BG scheme track mult
  } else if (trainConfig == 45) {
    cuts.AddCut("00000123", "00200009217302008210404000", "0152101500000000"); // variation chi2 30 psi pair 0.1 1D
    cuts.AddCut("00000123", "00200009217302008180404000", "0152101500000000"); // variation chi2 50 psi pair 0.2 2D
    cuts.AddCut("00000123", "00200009217302008860404000", "0152101500000000"); // variation chi2 20 psi pair 0.05 2D
    cuts.AddCut("00000123", "00200009217302009250404000", "0252101500000000"); // variation BG scheme track mult
  } else if (trainConfig == 46) {
    cuts.AddCut("00000113", "00200009227302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCut("00000113", "00200009217302002250404000", "0152101500000000"); // qT
    cuts.AddCut("00000113", "00200009217302009250404000", "0152101500000000"); // qT
    cuts.AddCut("00000113", "00200009217302008250004000", "0152101500000000"); // cosPA
  } else if (trainConfig == 47) {
    cuts.AddCut("00000123", "00200009227302008250404000", "0152101500000000"); //New standard cut for eta: alpha pT dependent
    cuts.AddCut("00000123", "00200009217302002250404000", "0152101500000000"); // qT
    cuts.AddCut("00000123", "00200009217302009250404000", "0152101500000000"); // qT
    cuts.AddCut("00000123", "00200009217302008250004000", "0152101500000000"); // cosPA
  } else if (trainConfig == 48) {
    cuts.AddCut("00000113", "00200009317302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCut("00000113", "00200009617302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCut("00000113", "00200009215302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCut("00000113", "00200009210302008250404000", "0152101500000000"); //dEdx variation
  } else if (trainConfig == 49) {
    cuts.AddCut("00000123", "00200009317302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCut("00000123", "00200009617302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCut("00000123", "00200009215302008250404000", "0152101500000000"); //dEdx variation
    cuts.AddCut("00000123", "00200009210302008250404000", "0152101500000000"); //dEdx variation

    
    //---------configs for V0AND 8TeV --------------------------//
  } else if (trainConfig == 60) {
    cuts.AddCut("00010113", "00200009227300008250404000", "0152103500000000"); //New standard cut for 8TeV analysis V0AND with double counting cut, TOF removed
    cuts.AddCut("00010013", "00200009227300008250404000", "0152103500000000"); // no SPD pileup cut
    cuts.AddCut("00010113", "00100009227300008250404000", "0152103500000000"); // R cut 2.8 -180 cm
    cuts.AddCut("00010113", "00500009227300008250404000", "0152103500000000"); // R cut 10. -180 cm
  } else if (trainConfig == 61) {
    cuts.AddCut("00010113", "00200069227300008250404000", "0152103500000000"); // min pT 40 MeV
    cuts.AddCut("00010113", "00200049227300008250404000", "0152103500000000"); // min pT 75 MeV
    cuts.AddCut("00010113", "00200019227300008250404000", "0152103500000000"); // min pT 100MeV
  } else if (trainConfig == 62) {
    cuts.AddCut("00010113", "00200008227300008250404000", "0152103500000000"); // TPC cluster 35%
    cuts.AddCut("00010113", "00200006227300008250404000", "0152103500000000"); // TPC cluster 70%
    cuts.AddCut("00010113", "00200009227300008250604000", "0152103500000000"); // cosPA 0.9
    cuts.AddCut("00010113", "00200009227300008250304000", "0152103500000000"); // cosPA 0.75
  } else if (trainConfig == 63) {
    cuts.AddCut("00010113", "00200009327300008250404000", "0152103500000000"); // nsig electron   -4,5
    cuts.AddCut("00010113", "00200009627300008250404000", "0152103500000000"); // nsig electron -2.5,4
    cuts.AddCut("00010113", "00200009257300008250404000", "0152103500000000"); // nsig pion 2,-10
    cuts.AddCut("00010113", "00200009217300008250404000", "0152103500000000"); // nsig pion 0,-10
  } else if (trainConfig == 64) {
    cuts.AddCut("00010113", "00200009220300008250404000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCut("00010113", "00200009226300008250404000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCut("00010113", "00200009227600008250404000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCut("00010113", "00200009227100008250404000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 65) {
    cuts.AddCut("00010113", "00200009227300003250404000", "0152103500000000"); // qT max 0.05 1D
    cuts.AddCut("00010113", "00200009227300002250404000", "0152103500000000"); // qT max 0.06 2D
    cuts.AddCut("00010113", "00200009227300009250404000", "0152103500000000"); // qT max 0.03 2D
  } else if (trainConfig == 66) {
    cuts.AddCut("00010113", "00200009227300008150404000", "0152103500000000"); // chi2 50
    cuts.AddCut("00010113", "00200009227300008850404000", "0152103500000000"); // chi2 20
    cuts.AddCut("00010113", "00200009227300008250400000", "0152103500000000"); // no double counting
    cuts.AddCut("00010113", "00200009227300008250406000", "0152103500000000"); // double count with open angle 0.04
  } else if (trainConfig == 67) {
    cuts.AddCut("00010113", "00200009227300008210404000", "0152103500000000"); // Psi pair 0.1  1D
    cuts.AddCut("00010113", "00200009227300008260404000", "0152103500000000"); // Psi pair 0.05 2D
    cuts.AddCut("00010113", "00200009227300008280404000", "0152103500000000"); // Psi pair 0.2  2D
  } else if (trainConfig == 68) {
    cuts.AddCut("00010113", "00200009227300008250404000", "0252103500000000"); // variation BG scheme track mult
    cuts.AddCut("00010113", "00200009227300008250404000", "0152107500000000"); // alpha meson 0.85
    cuts.AddCut("00010113", "00200009227300008250404000", "0152105500000000"); // alpha meson 0.75
  } else if (trainConfig == 69) {
    cuts.AddCut("00010113", "00200009227300008250404000", "0152103500000000"); //New standard cut for 8TeV analysis V0AND with double counting cut, TOF removed
    
  //---------configs for 8TeV triggers --------------------------//  
  } else if (trainConfig == 70) { //pp 8TeV cuts with EMC triggers
    cuts.AddCut("00010113", "00200009227302008250400000", "0152103500000000","1111111067032220000"); // standard cut 8tev
    cuts.AddCut("00052113", "00200009227302008250400000", "0152103500000000","1111111067032220000"); // trigger kEMC7
    cuts.AddCut("00081113", "00200009227302008250400000", "0152103500000000","1111111067032220000"); // trigger kEGA7
  } else if (trainConfig == 71) { //pp 8TeV cuts with EMC triggers, restricted phi regio EMC tight
    cuts.AddCut("00010113", "00204409227302008250400000", "0152103500000000","1111111067032220000"); // standard cut 8tev
    cuts.AddCut("00052113", "00204409227302008250400000", "0152103500000000","1111111067032220000"); // trigger kEMC7 
    cuts.AddCut("00081113", "00204409227302008250400000", "0152103500000000","1111111067032220000"); // trigger kEGA7
  } else if (trainConfig == 72) { //pp 8TeV cuts with EMC triggers, restricted phi regio EMC wide
    cuts.AddCut("00010113", "00202209227302008250400000", "0152103500000000","1111111067032220000"); // standard cut 8tev
    cuts.AddCut("00052113", "00202209227302008250400000", "0152103500000000","1111111067032220000"); // trigger kEMC7 
    cuts.AddCut("00081113", "00202209227302008250400000", "0152103500000000","1111111067032220000"); // trigger kEGA7


    //---------systematic studies mesons and direct photon 2016 pp 7TeV------------------//
  } else if (trainConfig == 80) { 
    cuts.AddCut("00000113", "00200009227300008250404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed
    cuts.AddCut("00000013", "00200009227300008250404000", "0152103500000000"); // no SPD pileup cut
    cuts.AddCut("00000113", "00100009227300008250404000", "0152103500000000"); // R cut 2.8 -180 cm
    cuts.AddCut("00000113", "00500009227300008250404000", "0152103500000000"); // R cut 10. -180 cm
  } else if (trainConfig == 81) {   
    cuts.AddCut("00000113", "00200069227300008250404000", "0152103500000000"); // min pT 40 MeV
    cuts.AddCut("00000113", "00200049227300008250404000", "0152103500000000"); // min pT 75 MeV
    cuts.AddCut("00000113", "00200019227300008250404000", "0152103500000000"); // min pT 100MeV
  } else if (trainConfig == 82) {  
    cuts.AddCut("00000113", "00200008227300008250404000", "0152103500000000"); // TPC cluster 35%
    cuts.AddCut("00000113", "00200006227300008250404000", "0152103500000000"); // TPC cluster 70%
    cuts.AddCut("00000113", "00200009227300008250604000", "0152103500000000"); // cosPA 0.9
    cuts.AddCut("00000113", "00200009227300008250304000", "0152103500000000"); // cosPA 0.75
  } else if (trainConfig == 83) {   
    cuts.AddCut("00000113", "00200009327300008250404000", "0152103500000000"); // nsig electron   -4,5
    cuts.AddCut("00000113", "00200009627300008250404000", "0152103500000000"); // nsig electron -2.5,4
    cuts.AddCut("00000113", "00200009257300008250404000", "0152103500000000"); // nsig pion 2,-10
    cuts.AddCut("00000113", "00200009217300008250404000", "0152103500000000"); // nsig pion 0,-10
  } else if (trainConfig == 84) {  
    cuts.AddCut("00000113", "00200009220300008250404000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCut("00000113", "00200009226300008250404000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCut("00000113", "00200009227600008250404000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCut("00000113", "00200009227100008250404000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 85) {  
    cuts.AddCut("00000113", "00200009227300003250404000", "0152103500000000"); // qT max 0.05 1D
    cuts.AddCut("00000113", "00200009227300002250404000", "0152103500000000"); // qT max 0.06 2D
    cuts.AddCut("00000113", "00200009227300009250404000", "0152103500000000"); // qT max 0.03 2D
  } else if (trainConfig == 86) {   
    cuts.AddCut("00000113", "00200009227300008150404000", "0152103500000000"); // chi2 50
    cuts.AddCut("00000113", "00200009227300008850404000", "0152103500000000"); // chi2 20
    cuts.AddCut("00000113", "00200009227300008250400000", "0152103500000000"); // no double counting
    cuts.AddCut("00000113", "00200009227300008250406000", "0152103500000000"); // double count with open angle 0.04
  } else if (trainConfig == 87) { 
    cuts.AddCut("00000113", "00200009227300008210404000", "0152103500000000"); // Psi pair 0.1  1D
    cuts.AddCut("00000113", "00200009227300008260404000", "0152103500000000"); // Psi pair 0.05 2D
    cuts.AddCut("00000113", "00200009227300008280404000", "0152103500000000"); // Psi pair 0.2  2D
  } else if (trainConfig == 88) { 
    cuts.AddCut("00000113", "00200009227300008250404000", "0252103500000000"); // variation BG scheme track mult
    cuts.AddCut("00000113", "00200009227300008250404000", "0152107500000000"); // alpha meson 0.85
    cuts.AddCut("00000113", "00200009227300008250404000", "0152105500000000"); // alpha meson 0.75

   //--------- pp7TeV purity studies (kappa cut)    -------------------------//
  } else if (trainConfig == 89) {
    cuts.AddCut("00000113", "00200009300000008250404000", "0152103500000000"); // -3 < kappa <  5
    cuts.AddCut("00000113", "00200009500000008250404000", "0152103500000000"); // -5 < kappa < 10
    cuts.AddCut("00000113", "00200009600000008250404000", "0152103500000000"); // -3 < kappa < 10
    cuts.AddCut("00000113", "00200009700000008250404000", "0152103500000000"); //  0 < kappa < 10

  // --------- testing validity of old 7 TeV result -------------------------//
  } else if (trainConfig == 90) {
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 91) {
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200009227302008254400000", "0152103500000000"); //asym pt dep
    cuts.AddCut("00000113", "00200009227302008255400000", "0152103500000000"); //asym tight pt dep
  } else if (trainConfig == 92) {
    cuts.AddCut("00000113", "00200009227302008250404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut
    cuts.AddCut("00000113", "00200009227302008250404000", "0152503500000000"); //y < 0.85
    cuts.AddCut("00000113", "00200009227302008250404000", "0152303500000000"); //y < 0.60
  } else if (trainConfig == 93) {
    cuts.AddCut("00000113", "00200009227300008250404000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR with double counting cut, TOF removed
    

  // ------------------------- run 2 High mult triggers --------------------------------------
  } else if (trainConfig == 100) {  
    cuts.AddCut("00074113", "00200009227302008250404000", "0152103500000000"); // for V0 High-Mult trigger
    cuts.AddCut("00075113", "00200009227302008250404000", "0152103500000000"); // for SPD High-Mult trigger
    cuts.AddCut("00074013", "00200009227302008250400000", "0152103500000000"); // check # of entries w/ pileup rejection cut for V0HM
    cuts.AddCut("00075013", "00200009227302008250400000", "0152103500000000"); // check # of entries w/o pileup rejection cut for SPHM
    
    
  // ------------------------- run 2 configurations ------------------------------------------  
  } else if (trainConfig == 111) {  // standard config for run 2
    cuts.AddCut("00010113", "00200009227302008250400000", "0152103500000000"); //New standard cut for eta analysis
    cuts.AddCut("00010113", "00200009227302008250400000", "0152101500000000"); //variation alpha pT dependent
    cuts.AddCut("00010113", "00200009227302008250400000", "0152109500000000"); //variation alpha
    cuts.AddCut("00010113", "00200009227302008250400000", "0152101500000002"); //variation alpha opan max
  } else if (trainConfig == 112) {  // standard config for run 2 with 
    cuts.AddCut("00010113", "00200009227302008250400000", "0152103500900000"); //New standard cut for eta analysis
    cuts.AddCut("00010113", "00200009227302008250400000", "0152101500900000"); //variation alpha pT dependent
    cuts.AddCut("00010113", "00200009227302008250400000", "0152109500900000"); //variation alpha
    cuts.AddCut("00010113", "00200009227302008250400000", "0152101500900002"); //variation alpha opan max
  } else if (trainConfig == 113) { 
    cuts.AddCut("00010113", "00200009227302008250400000", "0152101500000000"); 
    cuts.AddCut("00010113", "00200009227302008250404000", "0152101500000000"); // double counting cut
    cuts.AddCut("00010113", "00200009327302008250400000", "0152101500000000"); // dEdx 4 sigma below e
    cuts.AddCut("00010113", "00200009327302008250404000", "0152101500000000");
  } else if (trainConfig == 114) { 
    cuts.AddCut("00010113", "00200009227300008254404000", "0152101500000000"); // 13TeV with asymmetry and pT dep alpha cut
    cuts.AddCut("00010113", "00200009227300008254404000", "0152103500000000"); // 13TeV with asymmetry
    cuts.AddCut("00010113", "00200009227300008250404000", "0152101500000000"); // 13TeV pT dep alpha cut
    cuts.AddCut("00010113", "00200009227300008250404000", "0152103500000000"); // 13TeV
  } else if (trainConfig == 115) { 
    cuts.AddCut("00010113", "00200009227302008250404000", "0152101500000000");
    cuts.AddCut("00010113", "00200009327302008250404000", "0152101500000000"); // dEdx 4 sigma below e
    cuts.AddCut("00010113", "00200079227302008250404000", "0152101500000000"); // pT cut at 0
    cuts.AddCut("00010113", "00200079327302008250404000", "0152101500000000");
  } else if (trainConfig == 116) {
    cuts.AddCut("00010113", "00200009227302008254404000", "0152101500000000"); // standard cut Gamma pp 13TeV
    cuts.AddCut("30110113", "00200009227302008254404000", "0152101500000000"); // mult.: 0-5%
    cuts.AddCut("31310113", "00200009227302008254404000", "0152101500000000"); // mult.: 5-15%
    cuts.AddCut("33610113", "00200009227302008254404000", "0152101500000000"); // mult.: 15-30%
  } else if (trainConfig == 117) {
    cuts.AddCut("13510113", "00200009227302008254404000", "0152101500000000"); // mult.: 30-50%
    cuts.AddCut("15010113", "00200009227302008254404000", "0152101500000000"); // mult.: 50-100%
    cuts.AddCut("10110113", "00200009227302008254404000", "0152101500000000"); // mult.: 0-10%
    cuts.AddCut("11010113", "00200009227302008254404000", "0152101500000000"); // mult.: 10-100%
  } else if (trainConfig == 118){ 
    cuts.AddCut("00010113", "00200009227300008250404000", "0152103500000000"); // New standard cut for pp 5 TeV analysis VAND
    cuts.AddCut("00010113", "00200079227300008250404000", "0152103500000000"); // min pT no cut
    cuts.AddCut("00010113", "00200019227300008250404000", "0152103500000000"); // min pT 100 MeV
    cuts.AddCut("00010113", "00200008227300008250404000", "0152103500000000"); // TPC cluster 35% 
    cuts.AddCut("00010113", "00200006227300008250404000", "0152103500000000"); // TPC cluster 70%
  } else if (trainConfig == 119){ 
    cuts.AddCut("00010113", "00200009327300008250404000", "0152103500000000"); // edEdx -4,5
    cuts.AddCut("00010113", "00200009627300008250404000", "0152103500000000"); // edEdx -2.5,4
    cuts.AddCut("00010113", "00200009257300008250404000", "0152103500000000"); // pidEdx 2,-10
    cuts.AddCut("00010113", "00200009217300008250404000", "0152103500000000"); // pidEdx 0,-10
  } else if (trainConfig == 120){ 
    cuts.AddCut("00010113", "00200009220300008250404000", "0152103500000000"); // pion nsig min mom 0.50 GeV/c
    cuts.AddCut("00010113", "00200009226300008250404000", "0152103500000000"); // pion nsig min mom 0.25 GeV/c
    cuts.AddCut("00010113", "00200009227600008250404000", "0152103500000000"); // pion nsig max mom 2.00 GeV/c
    cuts.AddCut("00010113", "00200009227100008250404000", "0152103500000000"); // pion nsig max mom 5.00 GeV/c
  } else if (trainConfig == 121){ 
    cuts.AddCut("00010113", "00200009227300003250404000", "0152103500000000"); // qT max 0.05 1D
    cuts.AddCut("00010113", "00200009227300002250404000", "0152103500000000"); // qT max 0.06 2D
    cuts.AddCut("00010113", "00200009227300009250404000", "0152103500000000"); // qT max 0.03 2D    
    cuts.AddCut("00010113", "00200009227300008210404000", "0152103500000000"); // Psi pair 0.1  1D 
    cuts.AddCut("00010113", "00200009227300008260404000", "0152103500000000"); // Psi pair 0.05 2D
    cuts.AddCut("00010113", "00200009227300008280404000", "0152103500000000"); // Psi pair 0.2  2D
  } else if (trainConfig == 122){ 
    cuts.AddCut("00010113", "00200009227300008210404000", "0152101500000000"); // variation chi2 30 psi pair 0.1 1D
    cuts.AddCut("00010113", "00200009227300008860404000", "0152101500000000"); // variation chi2 20 psi pair 0.05 2D
    cuts.AddCut("00010113", "00200009227300008180404000", "0152101500000000"); // variation chi2 50 psi pair 0.2 2D
    cuts.AddCut("00010113", "00200009227300008254404000", "0152103500000000"); // Photon Asymmetry Cut
  } else if (trainConfig == 123){ 
    cuts.AddCut("00010113", "00200009227300008250604000", "0152103500000000"); // CosPA 0.9
    cuts.AddCut("00010113", "00200009227300008250304000", "0152103500000000"); // CosPA 0.75
    cuts.AddCut("00010113", "00200009227300008250400000", "0152103500000000"); // no double counting
    cuts.AddCut("00010113", "00200009227300008250404000", "0152105500000000"); // meson alpha < 0.75
    cuts.AddCut("00010113", "00200009227300008250404000", "0152107500000000"); // meson alpha < 0.85
    // -------- Material weights -several configs with same sets for running with various mat weights ----------
  } else if (trainConfig == 160) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 161) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 162) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 163) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 164) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 165) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 166) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 167) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 168) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz
  } else if (trainConfig == 169) { // like last last two in 70 and dalitz standard 7TeV
    cuts.AddCut("00000113", "00200009227302008250400000", "0152103500000000"); //New standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200008366300000200000000", "0163103100900000"); //Old standard cut for 7TeV analysis V0OR
    cuts.AddCut("00000113", "00200009360300007800004000", "0263103100900000"); //dalitz: New Standard Only MB, standard pp7Tev cut dalitz

    
  // -------------------------- mult cut studies --------------------------------------------
  } else if (trainConfig == 200) { // kMB
    cuts.AddCut("00100113", "00200009227302008250400000", "0152103500000000"); // 0 -2
    cuts.AddCut("01200113", "00200009227302008250400000", "0152103500000000"); // 2 -5
    cuts.AddCut("02300113", "00200009227302008250400000", "0152103500000000"); // 5 -10
    cuts.AddCut("03400113", "00200009227302008250400000", "0152103500000000"); // 10-30
    cuts.AddCut("04500113", "00200009227302008250400000", "0152103500000000"); // 30-100
    
  // -------------------------A. Marin,   2016 pp, open cuts --------------------------------------

 // Min Bias
  } else if (trainConfig == 300) {  
    cuts.AddCut("00010113", "00200009297302001280004000", "0152103500000000"); // Min Bias
    cuts.AddCut("00010113", "00200009297302001280004000", "0152101500000000"); // alpha pT dependent

 // High Mult V0
  } else if (trainConfig == 301) {  
    cuts.AddCut("00074113", "00200009297302001280004000", "0152103500000000"); // for V0 High-Mult trigger
    cuts.AddCut("00074013", "00200009297302001280004000", "0152103500000000"); // check # of entries w/ pileup rejection cut for V0HM

 // Low B Field
  } else if (trainConfig == 302) {  
    cuts.AddCut("00010113", "00200089297302001280004000", "0152103500000000"); // Min Bias
    cuts.AddCut("00010113", "00200089397302001280004000", "0152103500000000"); // Open dEdx
   

    
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
  TList *ClusterCutList = new TList();

  TList *HeaderList = new TList();
  if (periodname.Contains("LHC12i3")){
    TObjString *Header2 = new TObjString("BOX");
    HeaderList->Add(Header2);
  } else if (periodname.CompareTo("LHC14e2b")==0){
    TObjString *Header2 = new TObjString("pi0_1");
    HeaderList->Add(Header2);
    TObjString *Header3 = new TObjString("eta_2");
    HeaderList->Add(Header3);
  }
  
  TString energy = "";
  TString mcName = "";
  TString mcNameAdd = "";
  if (periodname.Contains("WOSDD")){
    mcNameAdd = "_WOSDD";
  } else if (periodname.Contains("WSDD")){
    mcNameAdd = "_WSDD";
  } 
  if (periodname.Contains("LHC12i3")){
    energy = "2760GeV";
    mcName = "Pythia8_LHC12i3";
  } else if (periodname.Contains("LHC12f1a")){
    energy = "2760GeV";
    mcName = "Pythia8_LHC12f1a";
  } else if (periodname.Contains("LHC12f1b")){
    energy = "2760GeV";
    mcName = "Phojet_LHC12f1b";
  } else if (periodname.Contains("LHC14e2a")){
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2a";
  } else if (periodname.Contains("LHC14e2b")){
    energy = "8TeV";
    mcName = "Pythia8_LHC14e2b";
  } else if (periodname.Contains("LHC14e2c")){
    energy = "8TeV";
    mcName = "Phojet_LHC14e2c";
  }
  
  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts        = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts      = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts  = new AliConversionMesonCuts*[numberOfCuts];
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts     = new AliCaloPhotonCuts*[numberOfCuts];
  Bool_t enableClustersForTrigger             = kFALSE;
  Bool_t initializedMatBudWeigths_existing    = kFALSE;
  
  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i]          = new AliConvEventCuts();
    TString fitNamePi0            = Form("Pi0_Fit_Data_%s",energy.Data());
    TString fitNameEta            = Form("Eta_Fit_Data_%s",energy.Data());
    TString fAddedSignalString    = (cuts.GetEventCut(i)).Data();
    fAddedSignalString            = fAddedSignalString(6,1);
    Bool_t fAddedSignal           = kFALSE;
    if (fAddedSignalString.CompareTo("2") == 0) 
      fAddedSignal                = kTRUE;

    TString mcInputNamePi0        = "";
    TString mcInputNameEta        = "";
    if (fAddedSignal && (periodname.Contains("LHC12i3") || periodname.CompareTo("LHC14e2b")==0)){
      mcInputNamePi0              = Form("Pi0_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta              = Form("Eta_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    } else {
      mcInputNamePi0              = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta              = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    }
    //if (doParticleWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kFALSE, kFALSE, kFALSE, fileNameInputForPartWeighting, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);
    if (doParticleWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForPartWeighting, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);

    TString dataInputMultHisto    = "";
    TString mcInputMultHisto      = "";
    TString triggerString         = (cuts.GetEventCut(i)).Data();
    triggerString                 = triggerString(3,2);
    if (triggerString.CompareTo("03")==0) 
      triggerString               = "00";
    if (periodNameAnchor.CompareTo("LHC13g") == 0 && triggerString.CompareTo("10")== 0 )
      triggerString               = "00";

    dataInputMultHisto            = Form("%s_%s", periodNameAnchor.Data(), triggerString.Data());
    mcInputMultHisto              = Form("%s_%s", periodNameV0Reader.Data(), triggerString.Data());
   
    if (doMultiplicityWeighting){
      cout << "enabling mult weighting" << endl;
      analysisEventCuts[i]->SetUseWeightMultiplicityFromFile( kTRUE, fileNameInputForMultWeighing, dataInputMultHisto, mcInputMultHisto );
    }

    
    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->SetV0ReaderName(V0ReaderName);
    if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts[i]->SetPeriodEnum(periodNameV0Reader);
    analysisEventCuts[i]->SetLightOutput(runLightOutput);
    analysisEventCuts[i]->InitializeCutsFromCutString((cuts.GetEventCut(i)).Data());

    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);

    if ( (trainConfig > 24 && trainConfig < 29) || ( trainConfig > 69 && trainConfig < 73 ) ){
        TString caloCutPos = cuts.GetClusterCut(i);
        caloCutPos.Resize(1);
        TString TrackMatcherName = Form("CaloTrackMatcher_%s",caloCutPos.Data());
        if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherName.Data()) ){
          AliCaloTrackMatcher* fTrackMatcher = new AliCaloTrackMatcher(TrackMatcherName.Data(),caloCutPos.Atoi());
          fTrackMatcher->SetV0ReaderName(V0ReaderName);
          mgr->AddTask(fTrackMatcher);
          mgr->ConnectInput(fTrackMatcher,0,cinput);
        }

        enableClustersForTrigger  = kTRUE;
        analysisClusterCuts[i]    = new AliCaloPhotonCuts();
        analysisClusterCuts[i]->SetV0ReaderName(V0ReaderName);
        analysisClusterCuts[i]->SetLightOutput(runLightOutput);
        analysisClusterCuts[i]->InitializeCutsFromCutString((cuts.GetClusterCut(i)).Data());
        ClusterCutList->Add(analysisClusterCuts[i]);
        analysisClusterCuts[i]->SetFillCutHistograms("");   
    }  
    
    analysisCuts[i]               = new AliConversionPhotonCuts();
    if ( trainConfig == 89 ){
      analysisCuts[i]->SetSwitchToKappaInsteadOfNSigdEdxTPC(kTRUE);
    }
    if (enableMatBudWeightsPi0 > 0){
        if (isMC > 0){
            if (analysisCuts[i]->InitializeMaterialBudgetWeights(enableMatBudWeightsPi0,filenameMatBudWeights)){
                initializedMatBudWeigths_existing = kTRUE;}
            else {cout << "ERROR The initialization of the materialBudgetWeights did not work out." << endl;}
        }
        else {cout << "ERROR 'enableMatBudWeightsPi0'-flag was set > 0 even though this is not a MC task. It was automatically reset to 0." << endl;}
    }
    analysisCuts[i]->SetV0ReaderName(V0ReaderName);
    analysisCuts[i]->SetLightOutput(runLightOutput);
    analysisCuts[i]->InitializeCutsFromCutString((cuts.GetPhotonCut(i)).Data());
    
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);

    analysisMesonCuts[i]          = new AliConversionMesonCuts();
    analysisMesonCuts[i]->SetLightOutput(runLightOutput);
    analysisMesonCuts[i]->InitializeCutsFromCutString((cuts.GetMesonCut(i)).Data());
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
    if(doSmear) analysisMesonCuts[i]->SetDefaultSmearing(bremSmear,smearPar,smearParConst);
  
  }

  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetConversionCutList(numberOfCuts,ConvCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetMoveParticleAccordingToVertex(kTRUE);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA
  task->SetDoChargedPrimary(enableChargedPrimary);
  if (enableClustersForTrigger){
    task->SetDoClusterSelectionForTriggerNorm(enableClustersForTrigger);
    task->SetClusterCutList(numberOfCuts,ClusterCutList);
  }  
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
