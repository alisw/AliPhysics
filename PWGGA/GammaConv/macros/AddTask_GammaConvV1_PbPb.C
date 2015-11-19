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
                                Bool_t    enableChargedPrimary            = kTRUE,
                                Bool_t    enableTriggerMimicking          = kFALSE,                         // enable trigger mimicking
                                Bool_t    enableTriggerOverlapRej         = kFALSE,                         // enable trigger overlap rejection
                                Float_t   maxFacPtHard                    = 3.,                             // maximum factor between hardest jet and ptHard generated
                                TString   periodNameV0Reader              = ""
                          ) {

  // ================= Load Librariers =================================
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCDB");
  gSystem->Load("libSTEER");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libTender");
  gSystem->Load("libTenderSupplies");
  gSystem->Load("libPWGflowBase");
  gSystem->Load("libPWGflowTasks");
  gSystem->Load("libPWGGAGammaConv");

  Int_t isHeavyIon = 1;

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
  if( !(AliV0ReaderV1*)mgr->GetTask("V0ReaderV1") ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");
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
      if( trainConfig >= 198 && trainConfig < 220){
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
  // Cut Numbers to use in Analysis
  Int_t numberOfCuts;
  if (trainConfig == 176 || trainConfig == 177 || trainConfig == 178 || trainConfig == 179 || trainConfig == 180 || trainConfig == 181) numberOfCuts = 7;
  else if( trainConfig >= 198 && trainConfig < 226 || trainConfig == 234) numberOfCuts = 1;
  else numberOfCuts = 5; 

  TString *eventCutArray = new TString[numberOfCuts];
  TString *photonCutArray = new TString[numberOfCuts];
  TString *mesonCutArray = new TString[numberOfCuts];

  if (trainConfig == 1){ // Standard cuts
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "04200009297002003220000000"; mesonCutArray[ 1] = "0152204500900000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "04200009297002003220000000"; mesonCutArray[ 2] = "0152204500900000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "04200009297002003220000000"; mesonCutArray[ 3] = "0152204500900000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "04200009297002003220000000"; mesonCutArray[ 4] = "0152204500900000"; // 0-20%
  } else if (trainConfig == 2) { // Standard cuts
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "04200009297002003220000000"; mesonCutArray[ 1] = "0152206500900000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "04200009297002003220000000"; mesonCutArray[ 2] = "0152206500900000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "04200009297002003220000000"; mesonCutArray[ 3] = "0152206500900000"; // 40-80%
    eventCutArray[ 4] = "54900013"; photonCutArray[ 4] = "04200009297002003220000000"; mesonCutArray[ 4] = "0152206500900000"; // 40-90%
  } else if (trainConfig == 3) { // Standard cuts only added signals
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "04200009297002003220000000"; mesonCutArray[ 1] = "0152204500900000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "04200009297002003220000000"; mesonCutArray[ 2] = "0152204500900000"; // 0-10%
    eventCutArray[ 3] = "51200023"; photonCutArray[ 3] = "04200009297002003220000000"; mesonCutArray[ 3] = "0152204500900000"; // 10-20%
    eventCutArray[ 4] = "50200023"; photonCutArray[ 4] = "04200009297002003220000000"; mesonCutArray[ 4] = "0152204500900000"; // 0-20%
  } else if (trainConfig == 4) { // Standard cuts only added signals
    eventCutArray[ 0] = "52400023"; photonCutArray[ 0] = "04200009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; // 20-40%
    eventCutArray[ 1] = "54600023"; photonCutArray[ 1] = "04200009297002003220000000"; mesonCutArray[ 1] = "0152206500900000"; // 40-60%
    eventCutArray[ 2] = "56800023"; photonCutArray[ 2] = "04200009297002003220000000"; mesonCutArray[ 2] = "0152206500900000"; // 60-80%
    eventCutArray[ 3] = "54800023"; photonCutArray[ 3] = "04200009297002003220000000"; mesonCutArray[ 3] = "0152206500900000"; // 20-40% 
    eventCutArray[ 4] = "54900023"; photonCutArray[ 4] = "04200009297002003220000000"; mesonCutArray[ 4] = "0152206500900000"; // 40-90%
  } else if (trainConfig == 5){ // R-minCut 7.5 cm
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "04900009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "04900009297002003220000000"; mesonCutArray[ 1] = "0152204500900000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "04900009297002003220000000"; mesonCutArray[ 2] = "0152204500900000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "04900009297002003220000000"; mesonCutArray[ 3] = "0152204500900000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "04900009297002003220000000"; mesonCutArray[ 4] = "0152204500900000"; // 0-20%
  } else if (trainConfig == 6) { // R-minCut 7.5 cm
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "04900009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "04900009297002003220000000"; mesonCutArray[ 1] = "0152206500900000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "04900009297002003220000000"; mesonCutArray[ 2] = "0152206500900000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "04900009297002003220000000"; mesonCutArray[ 3] = "0152206500900000"; // 40-80%
    eventCutArray[ 4] = "54900013"; photonCutArray[ 4] = "04900009297002003220000000"; mesonCutArray[ 4] = "0152206500900000"; // 40-90%
  } else if (trainConfig == 7) {// R-minCut 7.5 cm
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "04900009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "04900009297002003220000000"; mesonCutArray[ 1] = "0152204500900000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "04900009297002003220000000"; mesonCutArray[ 2] = "0152204500900000"; // 0-10%
    eventCutArray[ 3] = "51200023"; photonCutArray[ 3] = "04900009297002003220000000"; mesonCutArray[ 3] = "0152204500900000"; // 10-20%
    eventCutArray[ 4] = "50200023"; photonCutArray[ 4] = "04900009297002003220000000"; mesonCutArray[ 4] = "0152204500900000"; // 0-20%
  } else if (trainConfig == 8) { // R-minCut 7.5 cm
    eventCutArray[ 0] = "52400023"; photonCutArray[ 0] = "04900009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; // 20-40%
    eventCutArray[ 1] = "54600023"; photonCutArray[ 1] = "04900009297002003220000000"; mesonCutArray[ 1] = "0152206500900000"; // 40-60%
    eventCutArray[ 2] = "56800023"; photonCutArray[ 2] = "04900009297002003220000000"; mesonCutArray[ 2] = "0152206500900000"; // 60-80%
    eventCutArray[ 3] = "54800023"; photonCutArray[ 3] = "04900009297002003220000000"; mesonCutArray[ 3] = "0152206500900000"; // 20-40% 
    eventCutArray[ 4] = "54900023"; photonCutArray[ 4] = "04900009297002003220000000"; mesonCutArray[ 4] = "0152206500900000"; // 40-90% 
  } else if (trainConfig == 9){ // R-minCut 12.5 cm
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "04800009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "04800009297002003220000000"; mesonCutArray[ 1] = "0152204500900000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "04800009297002003220000000"; mesonCutArray[ 2] = "0152204500900000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "04800009297002003220000000"; mesonCutArray[ 3] = "0152204500900000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "04800009297002003220000000"; mesonCutArray[ 4] = "0152204500900000"; // 0-20%
  } else if (trainConfig == 10) { // R-minCut 12.5 cm
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "04800009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "04800009297002003220000000"; mesonCutArray[ 1] = "0152206500900000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "04800009297002003220000000"; mesonCutArray[ 2] = "0152206500900000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "04800009297002003220000000"; mesonCutArray[ 3] = "0152206500900000"; // 40-80%
    eventCutArray[ 4] = "54900013"; photonCutArray[ 4] = "04800009297002003220000000"; mesonCutArray[ 4] = "0152206500900000"; // 40-90%
  } else if (trainConfig == 11) {// R-minCut 12.5 cm
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "04800009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "04800009297002003220000000"; mesonCutArray[ 1] = "0152204500900000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "04800009297002003220000000"; mesonCutArray[ 2] = "0152204500900000"; // 0-10%
    eventCutArray[ 3] = "51200023"; photonCutArray[ 3] = "04800009297002003220000000"; mesonCutArray[ 3] = "0152204500900000"; // 10-20%
    eventCutArray[ 4] = "50200023"; photonCutArray[ 4] = "04800009297002003220000000"; mesonCutArray[ 4] = "0152204500900000"; // 0-20%
  } else if (trainConfig == 12) { // R-minCut 12.5 cm
    eventCutArray[ 0] = "52400023"; photonCutArray[ 0] = "04800009297002003220000000"; mesonCutArray[ 0] = "0152204500900000"; // 20-40%
    eventCutArray[ 1] = "54600023"; photonCutArray[ 1] = "04800009297002003220000000"; mesonCutArray[ 1] = "0152206500900000"; // 40-60%
    eventCutArray[ 2] = "56800023"; photonCutArray[ 2] = "04800009297002003220000000"; mesonCutArray[ 2] = "0152206500900000"; // 60-80%
    eventCutArray[ 3] = "54800023"; photonCutArray[ 3] = "04800009297002003220000000"; mesonCutArray[ 3] = "0152206500900000"; // 20-40% 
    eventCutArray[ 4] = "54900023"; photonCutArray[ 4] = "04800009297002003220000000"; mesonCutArray[ 4] = "0152206500900000"; // 40-90%
  }else  if (trainConfig == 13){ // LHC10h standard, eta 0.65, y = 0.6 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "03200009297002003220000000"; mesonCutArray[ 0] = "0152304500900000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "03200009297002003220000000"; mesonCutArray[ 1] = "0152304500900000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "03200009297002003220000000"; mesonCutArray[ 2] = "0152304500900000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "03200009297002003220000000"; mesonCutArray[ 3] = "0152304500900000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "03200009297002003220000000"; mesonCutArray[ 4] = "0152304500900000"; // 0-20%
  } else if (trainConfig == 14) {  // LHC10h standard, eta 0.65, y = 0.6 
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "03200009297002003220000000"; mesonCutArray[ 0] = "0152304500900000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "03200009297002003220000000"; mesonCutArray[ 1] = "0152306500900000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "03200009297002003220000000"; mesonCutArray[ 2] = "0152306500900000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "03200009297002003220000000"; mesonCutArray[ 3] = "0152306500900000"; // 40-80%
    eventCutArray[ 4] = "54900013"; photonCutArray[ 4] = "03200009297002003220000000"; mesonCutArray[ 4] = "0152306500900000"; // 40-90%
  } else if (trainConfig == 15) { // LHC10h standard, eta 0.65, y = 0.6  added signals
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "03200009297002003220000000"; mesonCutArray[ 0] = "0152304500900000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "03200009297002003220000000"; mesonCutArray[ 1] = "0152304500900000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "03200009297002003220000000"; mesonCutArray[ 2] = "0152304500900000"; // 0-10%
    eventCutArray[ 3] = "51200023"; photonCutArray[ 3] = "03200009297002003220000000"; mesonCutArray[ 3] = "0152304500900000"; // 10-20%
    eventCutArray[ 4] = "50200023"; photonCutArray[ 4] = "03200009297002003220000000"; mesonCutArray[ 4] = "0152304500900000"; // 0-20%
  } else if (trainConfig == 16) { // LHC10h standard, eta 0.65, y = 0.6  added signals
    eventCutArray[ 0] = "52400023"; photonCutArray[ 0] = "03200009297002003220000000"; mesonCutArray[ 0] = "0152304500900000"; // 20-40%
    eventCutArray[ 1] = "54600023"; photonCutArray[ 1] = "03200009297002003220000000"; mesonCutArray[ 1] = "0152306500900000"; // 40-60%
    eventCutArray[ 2] = "56800023"; photonCutArray[ 2] = "03200009297002003220000000"; mesonCutArray[ 2] = "0152306500900000"; // 60-80%
    eventCutArray[ 3] = "54800023"; photonCutArray[ 3] = "03200009297002003220000000"; mesonCutArray[ 3] = "0152306500900000"; // 20-40% 
    eventCutArray[ 4] = "54900023"; photonCutArray[ 4] = "03200009297002003220000000"; mesonCutArray[ 4] = "0152306500900000"; // 40-90%
  }else  if (trainConfig == 17){ // LHC10h standard, eta 0.65, y = 0.6, photon quality 1 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "03200009297002003220020000"; mesonCutArray[ 0] = "0152304500900000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "03200009297002003220020000"; mesonCutArray[ 1] = "0152304500900000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "03200009297002003220020000"; mesonCutArray[ 2] = "0152304500900000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "03200009297002003220020000"; mesonCutArray[ 3] = "0152304500900000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "03200009297002003220020000"; mesonCutArray[ 4] = "0152304500900000"; // 0-20%
  } else if (trainConfig == 18) {  // LHC10h standard, eta 0.65, y = 0.6, photon quality 1 
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "03200009297002003220020000"; mesonCutArray[ 0] = "0152304500900000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "03200009297002003220020000"; mesonCutArray[ 1] = "0152306500900000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "03200009297002003220020000"; mesonCutArray[ 2] = "0152306500900000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "03200009297002003220020000"; mesonCutArray[ 3] = "0152306500900000"; // 40-80%
    eventCutArray[ 4] = "54900013"; photonCutArray[ 4] = "03200009297002003220020000"; mesonCutArray[ 4] = "0152306500900000"; // 40-90%
  } else if (trainConfig == 19) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 1  added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "03200009297002003220020000"; mesonCutArray[ 0] = "0152304500900000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "03200009297002003220020000"; mesonCutArray[ 1] = "0152304500900000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "03200009297002003220020000"; mesonCutArray[ 2] = "0152304500900000"; // 0-10%
    eventCutArray[ 3] = "51200023"; photonCutArray[ 3] = "03200009297002003220020000"; mesonCutArray[ 3] = "0152304500900000"; // 10-20%
    eventCutArray[ 4] = "50200023"; photonCutArray[ 4] = "03200009297002003220020000"; mesonCutArray[ 4] = "0152304500900000"; // 0-20%
  } else if (trainConfig == 20) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 1 added signal
    eventCutArray[ 0] = "52400023"; photonCutArray[ 0] = "03200009297002003220020000"; mesonCutArray[ 0] = "0152304500900000"; // 20-40%
    eventCutArray[ 1] = "54600023"; photonCutArray[ 1] = "03200009297002003220020000"; mesonCutArray[ 1] = "0152306500900000"; // 40-60%
    eventCutArray[ 2] = "56800023"; photonCutArray[ 2] = "03200009297002003220020000"; mesonCutArray[ 2] = "0152306500900000"; // 60-80%
    eventCutArray[ 3] = "54800023"; photonCutArray[ 3] = "03200009297002003220020000"; mesonCutArray[ 3] = "0152306500900000"; // 20-40% 
    eventCutArray[ 4] = "54900023"; photonCutArray[ 4] = "03200009297002003220020000"; mesonCutArray[ 4] = "0152306500900000"; // 40-90%
  }else  if (trainConfig == 21){ // LHC10h standard, eta 0.65, y = 0.6, photon quality 2 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "03200009297002003220030000"; mesonCutArray[ 0] = "0152304500900000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "03200009297002003220030000"; mesonCutArray[ 1] = "0152304500900000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "03200009297002003220030000"; mesonCutArray[ 2] = "0152304500900000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "03200009297002003220030000"; mesonCutArray[ 3] = "0152304500900000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "03200009297002003220030000"; mesonCutArray[ 4] = "0152304500900000"; // 0-20%
  } else if (trainConfig == 22) {  // LHC10h standard, eta 0.65, y = 0.6, photon quality 2 
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "03200009297002003220030000"; mesonCutArray[ 0] = "0152304500900000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "03200009297002003220030000"; mesonCutArray[ 1] = "0152306500900000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "03200009297002003220030000"; mesonCutArray[ 2] = "0152306500900000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "03200009297002003220030000"; mesonCutArray[ 3] = "0152306500900000"; // 40-80%
    eventCutArray[ 4] = "54900013"; photonCutArray[ 4] = "03200009297002003220030000"; mesonCutArray[ 4] = "0152306500900000"; // 40-90%
  } else if (trainConfig == 23) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 2  added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "03200009297002003220030000"; mesonCutArray[ 0] = "0152304500900000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "03200009297002003220030000"; mesonCutArray[ 1] = "0152304500900000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "03200009297002003220030000"; mesonCutArray[ 2] = "0152304500900000"; // 0-10%
    eventCutArray[ 3] = "51200023"; photonCutArray[ 3] = "03200009297002003220030000"; mesonCutArray[ 3] = "0152304500900000"; // 10-20%
    eventCutArray[ 4] = "50200023"; photonCutArray[ 4] = "03200009297002003220030000"; mesonCutArray[ 4] = "0152304500900000"; // 0-20%
  } else if (trainConfig == 24) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 2 added signal
    eventCutArray[ 0] = "52400023"; photonCutArray[ 0] = "03200009297002003220030000"; mesonCutArray[ 0] = "0152304500900000"; // 20-40%
    eventCutArray[ 1] = "54600023"; photonCutArray[ 1] = "03200009297002003220030000"; mesonCutArray[ 1] = "0152306500900000"; // 40-60%
    eventCutArray[ 2] = "56800023"; photonCutArray[ 2] = "03200009297002003220030000"; mesonCutArray[ 2] = "0152306500900000"; // 60-80%
    eventCutArray[ 3] = "54800023"; photonCutArray[ 3] = "03200009297002003220030000"; mesonCutArray[ 3] = "0152306500900000"; // 20-40% 
    eventCutArray[ 4] = "54900023"; photonCutArray[ 4] = "03200009297002003220030000"; mesonCutArray[ 4] = "0152306500900000"; // 40-90%
  }else  if (trainConfig == 25){ // LHC10h standard, eta 0.65, y = 0.6, photon quality 3 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "03200009297002003220040000"; mesonCutArray[ 0] = "0152304500900000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "03200009297002003220040000"; mesonCutArray[ 1] = "0152304500900000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "03200009297002003220040000"; mesonCutArray[ 2] = "0152304500900000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "03200009297002003220040000"; mesonCutArray[ 3] = "0152304500900000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "03200009297002003220040000"; mesonCutArray[ 4] = "0152304500900000"; // 0-20%
  } else if (trainConfig == 26) {  // LHC10h standard, eta 0.65, y = 0.6, photon quality 3 
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "03200009297002003220040000"; mesonCutArray[ 0] = "0152304500900000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "03200009297002003220040000"; mesonCutArray[ 1] = "0152306500900000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "03200009297002003220040000"; mesonCutArray[ 2] = "0152306500900000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "03200009297002003220040000"; mesonCutArray[ 3] = "0152306500900000"; // 40-80%
    eventCutArray[ 4] = "54900013"; photonCutArray[ 4] = "03200009297002003220040000"; mesonCutArray[ 4] = "0152306500900000"; // 40-90%
  } else if (trainConfig == 27) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 3  added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "03200009297002003220040000"; mesonCutArray[ 0] = "0152304500900000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "03200009297002003220040000"; mesonCutArray[ 1] = "0152304500900000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "03200009297002003220040000"; mesonCutArray[ 2] = "0152304500900000"; // 0-10%
    eventCutArray[ 3] = "51200023"; photonCutArray[ 3] = "03200009297002003220040000"; mesonCutArray[ 3] = "0152304500900000"; // 10-20%
    eventCutArray[ 4] = "50200023"; photonCutArray[ 4] = "03200009297002003220040000"; mesonCutArray[ 4] = "0152304500900000"; // 0-20%
  } else if (trainConfig == 28) { // LHC10h standard, eta 0.65, y = 0.6, photon quality 3 added signal
    eventCutArray[ 0] = "52400023"; photonCutArray[ 0] = "03200009297002003220040000"; mesonCutArray[ 0] = "0152304500900000"; // 20-40%
    eventCutArray[ 1] = "54600023"; photonCutArray[ 1] = "03200009297002003220040000"; mesonCutArray[ 1] = "0152306500900000"; // 40-60%
    eventCutArray[ 2] = "56800023"; photonCutArray[ 2] = "03200009297002003220040000"; mesonCutArray[ 2] = "0152306500900000"; // 60-80%
    eventCutArray[ 3] = "54800023"; photonCutArray[ 3] = "03200009297002003220040000"; mesonCutArray[ 3] = "0152306500900000"; // 20-40% 
    eventCutArray[ 4] = "54900023"; photonCutArray[ 4] = "03200009297002003220040000"; mesonCutArray[ 4] = "0152306500900000"; // 40-90%
  }else  if (trainConfig == 29){ // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "03700009297002003220000000"; mesonCutArray[ 0] = "0152304500900000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "03700009297002003220000000"; mesonCutArray[ 1] = "0152304500900000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "03700009297002003220000000"; mesonCutArray[ 2] = "0152304500900000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "03700009297002003220000000"; mesonCutArray[ 3] = "0152304500900000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "03700009297002003220000000"; mesonCutArray[ 4] = "0152304500900000"; // 0-20%
  } else if (trainConfig == 30) {  // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm 
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "03700009297002003220000000"; mesonCutArray[ 0] = "0152304500900000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "03700009297002003220000000"; mesonCutArray[ 1] = "0152306500900000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "03700009297002003220000000"; mesonCutArray[ 2] = "0152306500900000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "03700009297002003220000000"; mesonCutArray[ 3] = "0152306500900000"; // 40-80%
    eventCutArray[ 4] = "54900013"; photonCutArray[ 4] = "03700009297002003220000000"; mesonCutArray[ 4] = "0152306500900000"; // 40-90%
  } else if (trainConfig == 31) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm  added signals
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "03700009297002003220000000"; mesonCutArray[ 0] = "0152304500900000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "03700009297002003220000000"; mesonCutArray[ 1] = "0152304500900000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "03700009297002003220000000"; mesonCutArray[ 2] = "0152304500900000"; // 0-10%
    eventCutArray[ 3] = "51200023"; photonCutArray[ 3] = "03700009297002003220000000"; mesonCutArray[ 3] = "0152304500900000"; // 10-20%
    eventCutArray[ 4] = "50200023"; photonCutArray[ 4] = "03700009297002003220000000"; mesonCutArray[ 4] = "0152304500900000"; // 0-20%
  } else if (trainConfig == 32) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm  added signals
    eventCutArray[ 0] = "52400023"; photonCutArray[ 0] = "03700009297002003220000000"; mesonCutArray[ 0] = "0152304500900000"; // 20-40%
    eventCutArray[ 1] = "54600023"; photonCutArray[ 1] = "03700009297002003220000000"; mesonCutArray[ 1] = "0152306500900000"; // 40-60%
    eventCutArray[ 2] = "56800023"; photonCutArray[ 2] = "03700009297002003220000000"; mesonCutArray[ 2] = "0152306500900000"; // 60-80%
    eventCutArray[ 3] = "54800023"; photonCutArray[ 3] = "03700009297002003220000000"; mesonCutArray[ 3] = "0152306500900000"; // 20-40% 
    eventCutArray[ 4] = "54900023"; photonCutArray[ 4] = "03700009297002003220000000"; mesonCutArray[ 4] = "0152306500900000"; // 40-90%
  }else  if (trainConfig == 33){ // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "03700009297002003220020000"; mesonCutArray[ 0] = "0152304500900000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "03700009297002003220020000"; mesonCutArray[ 1] = "0152304500900000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "03700009297002003220020000"; mesonCutArray[ 2] = "0152304500900000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "03700009297002003220020000"; mesonCutArray[ 3] = "0152304500900000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "03700009297002003220020000"; mesonCutArray[ 4] = "0152304500900000"; // 0-20%
  } else if (trainConfig == 34) {  // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1  
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "03700009297002003220020000"; mesonCutArray[ 0] = "0152304500900000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "03700009297002003220020000"; mesonCutArray[ 1] = "0152306500900000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "03700009297002003220020000"; mesonCutArray[ 2] = "0152306500900000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "03700009297002003220020000"; mesonCutArray[ 3] = "0152306500900000"; // 40-80%
    eventCutArray[ 4] = "54900013"; photonCutArray[ 4] = "03700009297002003220020000"; mesonCutArray[ 4] = "0152306500900000"; // 40-90%
  } else if (trainConfig == 35) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1  added signals
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "03700009297002003220020000"; mesonCutArray[ 0] = "0152304500900000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "03700009297002003220020000"; mesonCutArray[ 1] = "0152304500900000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "03700009297002003220020000"; mesonCutArray[ 2] = "0152304500900000"; // 0-10%
    eventCutArray[ 3] = "51200023"; photonCutArray[ 3] = "03700009297002003220020000"; mesonCutArray[ 3] = "0152304500900000"; // 10-20%
    eventCutArray[ 4] = "50200023"; photonCutArray[ 4] = "03700009297002003220020000"; mesonCutArray[ 4] = "0152304500900000"; // 0-20%
  } else if (trainConfig == 36) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 1  added signals
    eventCutArray[ 0] = "52400023"; photonCutArray[ 0] = "03700009297002003220020000"; mesonCutArray[ 0] = "0152304500900000"; // 20-40%
    eventCutArray[ 1] = "54600023"; photonCutArray[ 1] = "03700009297002003220020000"; mesonCutArray[ 1] = "0152306500900000"; // 40-60%
    eventCutArray[ 2] = "56800023"; photonCutArray[ 2] = "03700009297002003220020000"; mesonCutArray[ 2] = "0152306500900000"; // 60-80%
    eventCutArray[ 3] = "54800023"; photonCutArray[ 3] = "03700009297002003220020000"; mesonCutArray[ 3] = "0152306500900000"; // 20-40% 
    eventCutArray[ 4] = "54900023"; photonCutArray[ 4] = "03700009297002003220020000"; mesonCutArray[ 4] = "0152306500900000"; // 40-90%
  }else  if (trainConfig == 37){ // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "03700009297002003220040000"; mesonCutArray[ 0] = "0152304500900000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "03700009297002003220040000"; mesonCutArray[ 1] = "0152304500900000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "03700009297002003220040000"; mesonCutArray[ 2] = "0152304500900000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "03700009297002003220040000"; mesonCutArray[ 3] = "0152304500900000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "03700009297002003220040000"; mesonCutArray[ 4] = "0152304500900000"; // 0-20%
  } else if (trainConfig == 38) {  // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3  
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "03700009297002003220040000"; mesonCutArray[ 0] = "0152304500900000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "03700009297002003220040000"; mesonCutArray[ 1] = "0152306500900000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "03700009297002003220040000"; mesonCutArray[ 2] = "0152306500900000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "03700009297002003220040000"; mesonCutArray[ 3] = "0152306500900000"; // 40-80%
    eventCutArray[ 4] = "54900013"; photonCutArray[ 4] = "03700009297002003220040000"; mesonCutArray[ 4] = "0152306500900000"; // 40-90%
  } else if (trainConfig == 39) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3  added signals
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "03700009297002003220040000"; mesonCutArray[ 0] = "0152304500900000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "03700009297002003220040000"; mesonCutArray[ 1] = "0152304500900000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "03700009297002003220040000"; mesonCutArray[ 2] = "0152304500900000"; // 0-10%
    eventCutArray[ 3] = "51200023"; photonCutArray[ 3] = "03700009297002003220040000"; mesonCutArray[ 3] = "0152304500900000"; // 10-20%
    eventCutArray[ 4] = "50200023"; photonCutArray[ 4] = "03700009297002003220040000"; mesonCutArray[ 4] = "0152304500900000"; // 0-20%
  } else if (trainConfig == 40) { // LHC10h standard, eta 0.65, y = 0.6, min R = 35 cm, photon quality 3  added signals
    eventCutArray[ 0] = "52400023"; photonCutArray[ 0] = "03700009297002003220040000"; mesonCutArray[ 0] = "0152304500900000"; // 20-40%
    eventCutArray[ 1] = "54600023"; photonCutArray[ 1] = "03700009297002003220040000"; mesonCutArray[ 1] = "0152306500900000"; // 40-60%
    eventCutArray[ 2] = "56800023"; photonCutArray[ 2] = "03700009297002003220040000"; mesonCutArray[ 2] = "0152306500900000"; // 60-80%
    eventCutArray[ 3] = "54800023"; photonCutArray[ 3] = "03700009297002003220040000"; mesonCutArray[ 3] = "0152306500900000"; // 20-40% 
    eventCutArray[ 4] = "54900023"; photonCutArray[ 4] = "03700009297002003220040000"; mesonCutArray[ 4] = "0152306500900000"; // 40-90%
  } else if (trainConfig == 41){ // Standard cuts, eta 0.9, only to be run on data :kSemiCentral
    eventCutArray[ 0] = "60160013"; photonCutArray[ 0] = "00200009297002003220000000"; mesonCutArray[ 0] = "0152504500900000"; // 0-5%
    eventCutArray[ 1] = "61260013"; photonCutArray[ 1] = "00200009297002003220000000"; mesonCutArray[ 1] = "0152504500900000"; // 5-10%
    eventCutArray[ 2] = "50160013"; photonCutArray[ 2] = "00200009297002003220000000"; mesonCutArray[ 2] = "0152504500900000"; // 0-10%
    eventCutArray[ 3] = "61260013"; photonCutArray[ 3] = "00200009297002003220000000"; mesonCutArray[ 3] = "0152504500900000"; // 10-20%
    eventCutArray[ 4] = "50260013"; photonCutArray[ 4] = "00200009297002003220000000"; mesonCutArray[ 4] = "0152504500900000"; // 0-20%
  } else if (trainConfig == 42) { // Standard cuts, eta 0.9, only to be run on data :kSemiCentral
    eventCutArray[ 0] = "52360013"; photonCutArray[ 0] = "00200009297002003220000000"; mesonCutArray[ 0] = "0152506500900000"; 
    eventCutArray[ 1] = "53460013"; photonCutArray[ 1] = "00200009297002003220000000"; mesonCutArray[ 1] = "0152506500900000"; 
    eventCutArray[ 2] = "54560013"; photonCutArray[ 2] = "00200009297002003220000000"; mesonCutArray[ 2] = "0152506500900000"; 
    eventCutArray[ 3] = "55660013"; photonCutArray[ 3] = "00200009297002003220000000"; mesonCutArray[ 3] = "0152506500900000"; 
    eventCutArray[ 4] = "56760013"; photonCutArray[ 4] = "00200009297002003220000000"; mesonCutArray[ 4] = "0152506500900000"; 
  } else if (trainConfig == 43){ // Standard cuts, eta 0.9, only to be run on data
    eventCutArray[ 0] = "52300013"; photonCutArray[ 0] = "00200009297002003220000000"; mesonCutArray[ 0] = "0152504500900000"; 
    eventCutArray[ 1] = "53400013"; photonCutArray[ 1] = "00200009297002003220000000"; mesonCutArray[ 1] = "0152506500900000"; 
    eventCutArray[ 2] = "54500013"; photonCutArray[ 2] = "00200009297002003220000000"; mesonCutArray[ 2] = "0152506500900000"; 
    eventCutArray[ 3] = "55600013"; photonCutArray[ 3] = "00200009297002003220000000"; mesonCutArray[ 3] = "0152506500900000"; 
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002003220000000"; mesonCutArray[ 4] = "0152506500900000"; 
  } else if ( trainConfig == 44){ // qt elipse cut 0.05
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002008220000000"; mesonCutArray[ 0] = "0152506500900000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002008220000000"; mesonCutArray[ 1] = "0152506500900000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002008220000000"; mesonCutArray[ 2] = "0152506500900000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "00200009297002008220000000"; mesonCutArray[ 3] = "0152506500900000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "00200009297002008220000000"; mesonCutArray[ 4] = "0152506500900000"; // 0-20%
  } else if ( trainConfig == 45) { // qt elipse cut 0.05
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "00200009297002008220000000"; mesonCutArray[ 0] = "0152506500900000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "00200009297002008220000000"; mesonCutArray[ 1] = "0152506500900000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "00200009297002008220000000"; mesonCutArray[ 2] = "0152506500900000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "00200009297002008220000000"; mesonCutArray[ 3] = "0152506500900000"; // 40-80%
    eventCutArray[ 4] = "53500013"; photonCutArray[ 4] = "00200009297002008220000000"; mesonCutArray[ 4] = "0152506500900000"; // 30-50%
  } else if ( trainConfig == 46){ // qt elipse cut 0.05
    eventCutArray[ 0] = "52300013"; photonCutArray[ 0] = "00200009297002008220000000"; mesonCutArray[ 0] = "0152506500900000"; // 20-30% 
    eventCutArray[ 1] = "53400013"; photonCutArray[ 1] = "00200009297002008220000000"; mesonCutArray[ 1] = "0152506500900000"; // 30-40% 
    eventCutArray[ 2] = "54500013"; photonCutArray[ 2] = "00200009297002008220000000"; mesonCutArray[ 2] = "0152506500900000"; // 40-50% 
    eventCutArray[ 3] = "55600013"; photonCutArray[ 3] = "00200009297002008220000000"; mesonCutArray[ 3] = "0152506500900000"; // 50-60% 
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002008220000000"; mesonCutArray[ 4] = "0152506500900000"; // 60-70% 
  } else if ( trainConfig == 47){ // cos(theta_point) cut 0.85
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002003220400000"; mesonCutArray[ 0] = "0152506500900000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002003220400000"; mesonCutArray[ 1] = "0152506500900000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002003220400000"; mesonCutArray[ 2] = "0152506500900000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "00200009297002003220400000"; mesonCutArray[ 3] = "0152506500900000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "00200009297002003220400000"; mesonCutArray[ 4] = "0152506500900000"; // 0-20%
  } else if ( trainConfig == 48) { // cos(theta_point) cut 0.85
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "00200009297002003220400000"; mesonCutArray[ 0] = "0152506500900000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "00200009297002003220400000"; mesonCutArray[ 1] = "0152506500900000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "00200009297002003220400000"; mesonCutArray[ 2] = "0152506500900000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "00200009297002003220400000"; mesonCutArray[ 3] = "0152506500900000"; // 40-80%
    eventCutArray[ 4] = "53500013"; photonCutArray[ 4] = "00200009297002003220400000"; mesonCutArray[ 4] = "0152506500900000"; // 30-50%
  } else if ( trainConfig == 49){ // cos(theta_point) cut 0.85
    eventCutArray[ 0] = "52300013"; photonCutArray[ 0] = "00200009297002003220400000"; mesonCutArray[ 0] = "0152506500900000"; // 20-30% 
    eventCutArray[ 1] = "53400013"; photonCutArray[ 1] = "00200009297002003220400000"; mesonCutArray[ 1] = "0152506500900000"; // 30-40% 
    eventCutArray[ 2] = "54500013"; photonCutArray[ 2] = "00200009297002003220400000"; mesonCutArray[ 2] = "0152506500900000"; // 40-50% 
    eventCutArray[ 3] = "55600013"; photonCutArray[ 3] = "00200009297002003220400000"; mesonCutArray[ 3] = "0152506500900000"; // 50-60% 
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002003220400000"; mesonCutArray[ 4] = "0152506500900000"; // 60-70% 
  } else if ( trainConfig == 50){ // psi pair 2D 0.05
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002003260000000"; mesonCutArray[ 0] = "0152506500900000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002003260000000"; mesonCutArray[ 1] = "0152506500900000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002003260000000"; mesonCutArray[ 2] = "0152506500900000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "00200009297002003260000000"; mesonCutArray[ 3] = "0152506500900000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "00200009297002003260000000"; mesonCutArray[ 4] = "0152506500900000"; // 0-20%
  } else if ( trainConfig == 51) { // psi pair 2D 0.05
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "00200009297002003260000000"; mesonCutArray[ 0] = "0152506500900000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "00200009297002003260000000"; mesonCutArray[ 1] = "0152506500900000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "00200009297002003260000000"; mesonCutArray[ 2] = "0152506500900000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "00200009297002003260000000"; mesonCutArray[ 3] = "0152506500900000"; // 40-80%
    eventCutArray[ 4] = "53500013"; photonCutArray[ 4] = "00200009297002003260000000"; mesonCutArray[ 4] = "0152506500900000"; // 30-50%
  } else if ( trainConfig == 52){ // psi pair 2D 0.05
    eventCutArray[ 0] = "52300013"; photonCutArray[ 0] = "00200009297002003260000000"; mesonCutArray[ 0] = "0152506500900000"; // 20-30% 
    eventCutArray[ 1] = "53400013"; photonCutArray[ 1] = "00200009297002003260000000"; mesonCutArray[ 1] = "0152506500900000"; // 30-40% 
    eventCutArray[ 2] = "54500013"; photonCutArray[ 2] = "00200009297002003260000000"; mesonCutArray[ 2] = "0152506500900000"; // 40-50% 
    eventCutArray[ 3] = "55600013"; photonCutArray[ 3] = "00200009297002003260000000"; mesonCutArray[ 3] = "0152506500900000"; // 50-60% 
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002003260000000"; mesonCutArray[ 4] = "0152506500900000"; // 60-70% 
  } else if ( trainConfig == 53){ // psi pair 2D 0.1
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002003250000000"; mesonCutArray[ 0] = "0152506500900000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002003250000000"; mesonCutArray[ 1] = "0152506500900000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002003250000000"; mesonCutArray[ 2] = "0152506500900000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "00200009297002003250000000"; mesonCutArray[ 3] = "0152506500900000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "00200009297002003250000000"; mesonCutArray[ 4] = "0152506500900000"; // 0-20%
  } else if ( trainConfig == 54) { // psi pair 2D 0.1
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "00200009297002003250000000"; mesonCutArray[ 0] = "0152506500900000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "00200009297002003250000000"; mesonCutArray[ 1] = "0152506500900000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "00200009297002003250000000"; mesonCutArray[ 2] = "0152506500900000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "00200009297002003250000000"; mesonCutArray[ 3] = "0152506500900000"; // 40-80%
    eventCutArray[ 4] = "53500013"; photonCutArray[ 4] = "00200009297002003250000000"; mesonCutArray[ 4] = "0152506500900000"; // 30-50%
  } else if ( trainConfig == 55){ // psi pair 2D 0.1
    eventCutArray[ 0] = "52300013"; photonCutArray[ 0] = "00200009297002003250000000"; mesonCutArray[ 0] = "0152506500900000"; // 20-30% 
    eventCutArray[ 1] = "53400013"; photonCutArray[ 1] = "00200009297002003250000000"; mesonCutArray[ 1] = "0152506500900000"; // 30-40% 
    eventCutArray[ 2] = "54500013"; photonCutArray[ 2] = "00200009297002003250000000"; mesonCutArray[ 2] = "0152506500900000"; // 40-50% 
    eventCutArray[ 3] = "55600013"; photonCutArray[ 3] = "00200009297002003250000000"; mesonCutArray[ 3] = "0152506500900000"; // 50-60% 
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002003250000000"; mesonCutArray[ 4] = "0152506500900000"; // 60-70% 
  } else if ( trainConfig == 56){ // cleaner cuts central classes 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "50200013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20%
  } else if ( trainConfig == 57){ // cleaner cuts central classes - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "51200023"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "50200023"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20% 
  } else if ( trainConfig == 58) { // cleaner cuts semicentral classes
    eventCutArray[ 0] = "52400013"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; // 20-40%
    eventCutArray[ 1] = "54600013"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152506500000000"; // 40-60%
    eventCutArray[ 2] = "56800013"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152506500000000"; // 60-80%
    eventCutArray[ 3] = "54800013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152506500000000"; // 40-80%
    eventCutArray[ 4] = "53500013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152506500000000"; // 30-50%
  } else if ( trainConfig == 59) { // cleaner cuts semicentral classes - added signal
    eventCutArray[ 0] = "52400023"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; // 20-40%
    eventCutArray[ 1] = "54600023"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152506500000000"; // 40-60%
    eventCutArray[ 2] = "56800023"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152506500000000"; // 60-80%
    eventCutArray[ 3] = "54800023"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152506500000000"; // 40-80%
    eventCutArray[ 4] = "53500023"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152506500000000"; // 30-50%
  } else if ( trainConfig == 60){ // cleaner cuts finer centralities
    eventCutArray[ 0] = "52300013"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; // 20-30% 
    eventCutArray[ 1] = "53400013"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152506500000000"; // 30-40% 
    eventCutArray[ 2] = "54500013"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152506500000000"; // 40-50% 
    eventCutArray[ 3] = "55600013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152506500000000"; // 50-60% 
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152506500000000"; // 60-70%
  } else if ( trainConfig == 61){ // cleaner cuts finer centralities - added signal
    eventCutArray[ 0] = "52300023"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; // 20-30% 
    eventCutArray[ 1] = "53400023"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152506500000000"; // 30-40% 
    eventCutArray[ 2] = "54500023"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152506500000000"; // 40-50% 
    eventCutArray[ 3] = "55600023"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152506500000000"; // 50-60% 
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152506500000000"; // 60-70% 
  } else if ( trainConfig == 62){ // cleaner cuts finer centralities
    eventCutArray[ 0] = "62300013"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "63400013"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "64500013"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "65600013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "66700013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20%
  } else if ( trainConfig == 63){ // cleaner cuts finer centralities - added signal
    eventCutArray[ 0] = "62300023"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "63400023"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "64500023"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "65600023"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "66700023"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20%
  } else if ( trainConfig == 64){ // cleaner cuts
    eventCutArray[ 0] = "67800013"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "68900013"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "56700013"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "57800013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "58900013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20%
  } else if ( trainConfig == 65){ // cleaner cuts added signal
    eventCutArray[ 0] = "67800023"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "68900023"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "56700023"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "57800023"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "58900023"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20%
  } else if ( trainConfig == 66){ // cleaner cuts
    eventCutArray[ 0] = "70100013"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "71200013"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "72300013"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "73400013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "74500013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20%
  } else if ( trainConfig == 67){ // cleaner cuts added signal
    eventCutArray[ 0] = "70100023"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "71200023"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "72300023"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "73400023"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "74500023"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20%
  } else if ( trainConfig == 68){ // cleaner cuts
    eventCutArray[ 0] = "75600013"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0]= "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "76700013"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1]= "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "77800013"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2]= "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "78900013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3]= "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "70900013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4]= "0152506500000000"; // 0-20%
  } else if ( trainConfig == 69){ // cleaner cuts added signal
    eventCutArray[ 0] = "75600023"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0]= "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "76700023"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1]= "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "77800023"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2]= "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "78900023"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3]= "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "70900023"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4]= "0152506500000000"; // 0-20%
  } else if ( trainConfig == 70){ // variation eta  0.65 ----------- here start the syst var for LHC11h (and not only) -------------------------
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "03200009247602008250400000"; mesonCutArray[ 0] = "0152301500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "03200009247602008250400000"; mesonCutArray[ 1] = "0152301500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "03200009247602008250400000"; mesonCutArray[ 2] = "0152301500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "03200009247602008250400000"; mesonCutArray[ 3] = "0152301500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "03200009247602008250400000"; mesonCutArray[ 4] = "0152301500000000"; // 20-50% 
  } else if ( trainConfig == 71){ // variation eta  0.65 - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "03200009247602008250400000"; mesonCutArray[ 0] = "0152301500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "03200009247602008250400000"; mesonCutArray[ 1] = "0152301500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "03200009247602008250400000"; mesonCutArray[ 2] = "0152301500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "03200009247602008250400000"; mesonCutArray[ 3] = "0152301500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "03200009247602008250400000"; mesonCutArray[ 4] = "0152301500000000"; // 20-50% 
  } else if ( trainConfig == 72){ // variation eta  0.65 with phi cut 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "03216609247602008250400000"; mesonCutArray[ 0] = "0152301500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "03216609247602008250400000"; mesonCutArray[ 1] = "0152301500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "03216609247602008250400000"; mesonCutArray[ 2] = "0152301500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "03216609247602008250400000"; mesonCutArray[ 3] = "0152301500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "03216609247602008250400000"; mesonCutArray[ 4] = "0152301500000000"; // 20-50% 
  } else if ( trainConfig == 73){ // variation eta  0.65 with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "03216609247602008250400000"; mesonCutArray[ 0] = "0152301500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "03216609247602008250400000"; mesonCutArray[ 1] = "0152301500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "03216609247602008250400000"; mesonCutArray[ 2] = "0152301500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "03216609247602008250400000"; mesonCutArray[ 3] = "0152301500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "03216609247602008250400000"; mesonCutArray[ 4] = "0152301500000000"; // 20-50% 
  } else if ( trainConfig == 74){ // variation eta  0.75
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "04200009247602008250400000"; mesonCutArray[ 0] = "0152201500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "04200009247602008250400000"; mesonCutArray[ 1] = "0152201500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "04200009247602008250400000"; mesonCutArray[ 2] = "0152201500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "04200009247602008250400000"; mesonCutArray[ 3] = "0152201500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "04200009247602008250400000"; mesonCutArray[ 4] = "0152201500000000"; // 20-50% 
  } else if ( trainConfig == 75){ // variation eta  0.75 added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "04200009247602008250400000"; mesonCutArray[ 0] = "0152201500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "04200009247602008250400000"; mesonCutArray[ 1] = "0152201500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "04200009247602008250400000"; mesonCutArray[ 2] = "0152201500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "04200009247602008250400000"; mesonCutArray[ 3] = "0152201500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "04200009247602008250400000"; mesonCutArray[ 4] = "0152201500000000"; // 20-50% 
  } else if ( trainConfig == 76){ // variation eta  0.75 with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "04216609247602008250400000"; mesonCutArray[ 0] = "0152201500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "04216609247602008250400000"; mesonCutArray[ 1] = "0152201500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "04216609247602008250400000"; mesonCutArray[ 2] = "0152201500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "04216609247602008250400000"; mesonCutArray[ 3] = "0152201500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "04216609247602008250400000"; mesonCutArray[ 4] = "0152201500000000"; // 20-50% 
  } else if ( trainConfig == 77){ // variation eta  0.75 with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "04216609247602008250400000"; mesonCutArray[ 0] = "0152201500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "04216609247602008250400000"; mesonCutArray[ 1] = "0152201500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "04216609247602008250400000"; mesonCutArray[ 2] = "0152201500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "04216609247602008250400000"; mesonCutArray[ 3] = "0152201500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "04216609247602008250400000"; mesonCutArray[ 4] = "0152201500000000"; // 20-50% 
  } else if ( trainConfig == 78){ // min R = 35 cm
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00700009247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00700009247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00700009247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00700009247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 10-20%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00700009247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 0-20%
  } else if ( trainConfig == 79){ // min R = 35 cm added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00700009247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00700009247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00700009247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00700009247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 10-20%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00700009247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 0-20%
  } else if ( trainConfig == 80){ // min R = 35 cm with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00716609247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00716609247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00716609247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00716609247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 10-20%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00716609247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 0-20%
  } else if ( trainConfig == 81){ // min R = 35 cm with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00716609247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00716609247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00716609247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00716609247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 10-20%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00716609247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 0-20%
  } else if ( trainConfig == 82){ // single pt 0.075
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200049247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200049247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200049247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200049247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200049247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 83){ // single pt 0.075 - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200049247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200049247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200049247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200049247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200049247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 84){ // single pt 0.075 with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216649247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216649247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216649247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216649247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216649247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 85){ // single pt 0.075 with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216649247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216649247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216649247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216649247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216649247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 86){ // single pt 0.1
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200019247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200019247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200019247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200019247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200019247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; //20-50%
  } else if ( trainConfig == 87){ // single pt 0.1 - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200019247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200019247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200019247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200019247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200019247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; //20-50%
  } else if ( trainConfig == 88){ // single pt 0.1 with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216619247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216619247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216619247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216619247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216619247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; //20-50%
  } else if ( trainConfig == 89){ // single pt 0.1  with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216619247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216619247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216619247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216619247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216619247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; //20-50%
  } else if ( trainConfig == 90){ // variation TPC cls 0.7
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200006247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200006247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200006247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200006247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200006247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 91){ // variation TPC cls 0.7 added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200006247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200006247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200006247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200006247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200006247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 92){ // variation TPC cls 0.7 with phi cut 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216606247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216606247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216606247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216606247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216606247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 93){ // variation TPC cls 0.7 with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216606247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216606247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216606247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216606247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216606247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 94){ // variation TPC cls 0.35
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200008247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200008247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200008247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200008247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200008247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 95){ // variation TPC cls 0.35 - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200008247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200008247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200008247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200008247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200008247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 96){ // variation TPC cls 0.35 with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216608247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216608247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216608247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216608247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216608247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 97){ // variation TPC cls 0.35  with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216608247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216608247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216608247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216608247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216608247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 98){ // variation edEdx  -4,5
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009347602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009347602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009347602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009347602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009347602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 99){ // variation edEdx  -4,5 - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009347602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009347602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009347602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009347602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009347602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 100){ // variation edEdx  -4,5 with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609347602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609347602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609347602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609347602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609347602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 101){ // variation edEdx  -4,5 with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609347602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609347602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609347602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609347602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609347602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 102){ // variation edEdx  -2.5,4
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009647602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009647602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009647602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009647602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009647602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 103){ // variation edEdx  -2.5,4 - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009647602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009647602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009647602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009647602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009647602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 104){ // variation edEdx  -2.5,4 with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609647602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609647602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609647602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609647602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609647602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 105){ // variation edEdx  -2.5,4 with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609647602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609647602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609647602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609647602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609647602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 106){ //variation pion p dEdx 2.0 sigma, 1 sigma high pt cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009287602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009287602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009287602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009287602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009287602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 107){ //variation pion p dEdx 2.0 sigma, 1 sigma high pt cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009287602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009287602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009287602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009287602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009287602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 108){ //variation pion p dEdx 2.0 sigma, 1 sigma high pt cut with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609287602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609287602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609287602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609287602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609287602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 109){ //variation pion p dEdx 2.0 sigma, 1 sigma high pt cut with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609287602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609287602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609287602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609287602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609287602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 110){ //variation pion p dEdx 2.5 sigma, no high pt cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009237002008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009237002008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009237002008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009237002008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009237002008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 111){ //variation pion p dEdx  2.5 sigma, no high pt cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009237002008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009237002008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009237002008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009237002008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009237002008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 112){ //variation pion p dEdx 2.5 sigma, no high pt cut with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609237002008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609237002008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609237002008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609237002008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609237002008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 113){ //variation pion p dEdx  2.5 sigma, no high pt cut  with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609237002008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609237002008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609237002008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609237002008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609237002008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 114){ //variation pion p dEdx 3.0 sigma, 0.3-3.
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009245402008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009245402008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009245402008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009245402008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009245402008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 115){ //variation pion p dEdx 3.0 sigma,  0.3-3. - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009245402008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009245402008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009245402008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009245402008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009245402008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 116){ //variation pion p dEdx 3.0 sigma, 0.3-3. with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609245402008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609245402008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609245402008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609245402008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609245402008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 117){ //variation pion p dEdx 3.0 sigma,  0.3-3. with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609245402008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609245402008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609245402008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609245402008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609245402008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 118){ // TOF el. PID -3,5
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247603008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247603008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247603008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247603008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247603008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 119){ // TOF el. PID -3,5 - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247603008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247603008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247603008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247603008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247603008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 120){ // TOF el. PID -3,5 with phi cut 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247603008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247603008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247603008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247603008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247603008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 121){ // TOF el. PID -3,5 with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247603008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247603008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247603008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247603008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247603008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 122){ // TOF el. PID -2,3
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247604008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247604008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247604008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247604008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247604008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 123){ // TOF el. PID -2,3 - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247604008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247604008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247604008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247604008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247604008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 124){ // TOF el. PID -2,3 with phi cut 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247604008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247604008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247604008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247604008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247604008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 125){ // TOF el. PID -2,3  with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247604008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247604008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247604008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247604008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247604008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 126){ // qt 0.03 2D
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602009250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602009250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602009250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602009250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602009250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 127){ // qt 0.03 2D - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602009250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602009250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602009250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247602009250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247602009250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 128){ // qt 0.03 2D with phi cut 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247602009250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247602009250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602009250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247602009250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247602009250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 129){ // qt 0.03 2D with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247602009250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247602009250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602009250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247602009250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247602009250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 130){ // qt 0.07 no2D
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602002250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602002250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602002250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602002250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602002250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 131){ // qt 0.07 no2D - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602002250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602002250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602002250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247602002250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247602002250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 132){ // qt 0.07 no2D with phi cut 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247602002250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247602002250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602002250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247602002250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247602002250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 133){ // qt 0.07 no2D with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247602002250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247602002250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602002250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247602002250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247602002250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 134){ // chi2  50.
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602008150400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602008150400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602008150400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008150400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008150400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 135){ // chi2  50.  added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602008150400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602008150400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602008150400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247602008150400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247602008150400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 136){ // chi2  50. with phi cut 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247602008150400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247602008150400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602008150400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247602008150400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247602008150400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 137){ // chi2  50. with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247602008150400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247602008150400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602008150400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247602008150400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247602008150400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 138){ // chi2  20.
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602008850400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602008850400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602008850400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008850400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008850400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 139){ // chi2  20.  added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602008850400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602008850400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602008850400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247602008850400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247602008850400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 140){ // chi2  20. with phi cut 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247602008850400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247602008850400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602008850400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247602008850400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247602008850400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 141){ // chi2  20. with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247602008850400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247602008850400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602008850400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247602008850400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247602008850400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 142){ // psi pair 0.05 2D
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602008260400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602008260400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602008260400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008260400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008260400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 143){ // psi pair 0.05 2D - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602008260400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602008260400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602008260400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247602008260400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247602008260400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 144){ // psi pair 0.05 2D with phi cut 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247602008260400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247602008260400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602008260400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247602008260400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247602008260400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 145){ // psi pair 0.05 2D with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247602008260400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247602008260400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602008260400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247602008260400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247602008260400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 146){ // psi pair 0.2 2D
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602008280400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602008280400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602008280400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008280400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008280400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 147){ // psi pair 0.2 2D - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602008280400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602008280400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602008280400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247602008280400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247602008280400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 148){ // psi pair 0.2 2D with phi cut 
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247602008280400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247602008280400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602008280400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247602008280400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247602008280400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 149){ // psi pair 0.2 2D with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247602008280400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247602008280400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602008280400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247602008280400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247602008280400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 150){ // cosPA -1
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602008250000000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602008250000000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602008250000000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008250000000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008250000000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 151){ // cosPA -1 - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602008250000000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602008250000000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602008250000000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247602008250000000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247602008250000000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 152){ // cosPA -1 with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247602008250000000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247602008250000000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602008250000000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247602008250000000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247602008250000000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 153){ // cosPA -1 with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247602008250000000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247602008250000000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602008250000000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247602008250000000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247602008250000000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 154){ // variation alpha 0.75
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0]= "0152505500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1]= "0152505500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2]= "0152505500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3]= "0152505500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4]= "0152505500000000"; // 20-50%
  } else if ( trainConfig == 155){ // variation alpha 0.75 - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0]= "0152505500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1]= "0152505500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2]= "0152505500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3]= "0152505500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4]= "0152505500000000"; // 20-50%
  } else if ( trainConfig == 156){ // variation alpha 0.75 with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247602008250400000"; mesonCutArray[ 0]= "0152505500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247602008250400000"; mesonCutArray[ 1]= "0152505500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602008250400000"; mesonCutArray[ 2]= "0152505500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247602008250400000"; mesonCutArray[ 3]= "0152505500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247602008250400000"; mesonCutArray[ 4]= "0152505500000000"; // 20-50%
  } else if ( trainConfig == 157){ // variation alpha 0.75 with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247602008250400000"; mesonCutArray[ 0]= "0152505500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247602008250400000"; mesonCutArray[ 1]= "0152505500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602008250400000"; mesonCutArray[ 2]= "0152505500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247602008250400000"; mesonCutArray[ 3]= "0152505500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247602008250400000"; mesonCutArray[ 4]= "0152505500000000"; // 20-50%
  } else if ( trainConfig == 158){ // variation alpha 1.
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0]= "0152503500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1]= "0152503500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2]= "0152503500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3]= "0152503500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4]= "0152503500000000"; // 20-50%
  } else if ( trainConfig == 159){ // variation alpha 1. - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0]= "0152503500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1]= "0152503500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2]= "0152503500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3]= "0152503500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4]= "0152503500000000"; // 20-50%
  } else if ( trainConfig == 160){ // variation alpha 1. with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247602008250400000"; mesonCutArray[ 0]= "0152503500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247602008250400000"; mesonCutArray[ 1]= "0152503500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602008250400000"; mesonCutArray[ 2]= "0152503500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247602008250400000"; mesonCutArray[ 3]= "0152503500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247602008250400000"; mesonCutArray[ 4]= "0152503500000000"; // 20-50%
  } else if ( trainConfig == 161){ // variation alpha 1. with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247602008250400000"; mesonCutArray[ 0]= "0152503500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247602008250400000"; mesonCutArray[ 1]= "0152503500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602008250400000"; mesonCutArray[ 2]= "0152503500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247602008250400000"; mesonCutArray[ 3]= "0152503500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247602008250400000"; mesonCutArray[ 4]= "0152503500000000"; // 20-50%
  } else if ( trainConfig == 162){ // standard LHC11h cut selection
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 163){ // standard LHC11h cut selection - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 164){ // standard LHC11h cut selection with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247602008250400000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247602008250400000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602008250400000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247602008250400000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247602008250400000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 165){ // standard LHC11h cut selection with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247602008250400000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247602008250400000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602008250400000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247602008250400000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247602008250400000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 166){ // cleaner cuts photon Quality 1
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002008250420000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002008250420000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002008250420000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009297002008250420000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002008250420000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20%
  } else if ( trainConfig == 167){ // cleaner cuts added signal photon Quality 1
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009297002008250420000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009297002008250420000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009297002008250420000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009297002008250420000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009297002008250420000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20% 
  } else if ( trainConfig == 168){ // cleaner cuts photon Quality 2
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002008250430000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002008250430000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002008250430000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009297002008250430000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002008250430000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20%
  } else if ( trainConfig == 169){ // cleaner cuts added signal photon Quality 2
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009297002008250430000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009297002008250430000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009297002008250430000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009297002008250430000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009297002008250430000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20% 
  } else if ( trainConfig == 170){ // cleaner cuts photon Quality 3
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002008250440000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002008250440000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009297002008250440000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009297002008250440000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009297002008250440000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20%
  } else if ( trainConfig == 171){ // cleaner cuts added signal photon Quality 3
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009297002008250440000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009297002008250440000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009297002008250440000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009297002008250440000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009297002008250440000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20% 
  } else if ( trainConfig == 172){ // cleaner cuts, photon Quality 1, min R = 35 cm
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00700009297002008250420000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00700009297002008250420000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00700009297002008250420000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00700009297002008250420000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00700009297002008250420000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20%
  } else if ( trainConfig == 173){ // cleaner cuts added signal, photon Quality 1, min R = 35 cm
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00700009297002008250420000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00700009297002008250420000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00700009297002008250420000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00700009297002008250420000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00700009297002008250420000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20% 
  } else if ( trainConfig == 174){ // cleaner cuts, photon Quality 3, min R = 35 cm
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00700009297002008250440000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00700009297002008250440000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00700009297002008250440000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00700009297002008250440000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00700009297002008250440000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20%
  } else if ( trainConfig == 175){ // cleaner cuts added signal, photon Quality 3, min R = 35 cm
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00700009297002008250440000"; mesonCutArray[ 0] = "0152506500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00700009297002008250440000"; mesonCutArray[ 1] = "0152506500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00700009297002008250440000"; mesonCutArray[ 2] = "0152506500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00700009297002008250440000"; mesonCutArray[ 3] = "0152506500000000"; // 10-20%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00700009297002008250440000"; mesonCutArray[ 4] = "0152506500000000"; // 0-20% 
  } else if ( trainConfig == 176){ // flow cuts with eta = 0.9, y = 0.85
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[ 0] = "0152506500000000";
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002008250400000"; mesonCutArray[ 1] = "0152506500000000";
    eventCutArray[ 2] = "51200013"; photonCutArray[ 2] = "00200009297002008250400000"; mesonCutArray[ 2] = "0152506500000000";
    eventCutArray[ 3] = "52300013"; photonCutArray[ 3] = "00200009297002008250400000"; mesonCutArray[ 3] = "0152506500000000";
    eventCutArray[ 4] = "53400013"; photonCutArray[ 4] = "00200009297002008250400000"; mesonCutArray[ 4] = "0152506500000000";
    eventCutArray[ 5] = "54600013"; photonCutArray[ 5] = "00200009297002008250400000"; mesonCutArray[ 5] = "0152506500000000";
    eventCutArray[ 6] = "56800013"; photonCutArray[ 6] = "00200009297002008250400000"; mesonCutArray[ 6] = "0152506500000000";
  } else if ( trainConfig == 177){ // flow cuts with eta = 0.65, y = 0.6
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "03200009297002008250400000"; mesonCutArray[ 0] = "0152306500000000";
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "03200009297002008250400000"; mesonCutArray[ 1] = "0152306500000000";
    eventCutArray[ 2] = "51200013"; photonCutArray[ 2] = "03200009297002008250400000"; mesonCutArray[ 2] = "0152306500000000";
    eventCutArray[ 3] = "52300013"; photonCutArray[ 3] = "03200009297002008250400000"; mesonCutArray[ 3] = "0152306500000000";
    eventCutArray[ 4] = "53400013"; photonCutArray[ 4] = "03200009297002008250400000"; mesonCutArray[ 4] = "0152306500000000";
    eventCutArray[ 5] = "54600013"; photonCutArray[ 5] = "03200009297002008250400000"; mesonCutArray[ 5] = "0152306500000000";
    eventCutArray[ 6] = "56800013"; photonCutArray[ 6] = "03200009297002008250400000"; mesonCutArray[ 6] = "0152306500000000";
  } else if ( trainConfig == 178){ // flow cuts with eta = 0.6, y = 0.5
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "01200009297002008250400000"; mesonCutArray[ 0] = "0152406500000000";
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "01200009297002008250400000"; mesonCutArray[ 1] = "0152406500000000";
    eventCutArray[ 2] = "51200013"; photonCutArray[ 2] = "01200009297002008250400000"; mesonCutArray[ 2] = "0152406500000000";
    eventCutArray[ 3] = "52300013"; photonCutArray[ 3] = "01200009297002008250400000"; mesonCutArray[ 3] = "0152406500000000";
    eventCutArray[ 4] = "53400013"; photonCutArray[ 4] = "01200009297002008250400000"; mesonCutArray[ 4] = "0152406500000000";
    eventCutArray[ 5] = "54600013"; photonCutArray[ 5] = "01200009297002008250400000"; mesonCutArray[ 5] = "0152406500000000";
    eventCutArray[ 6] = "56800013"; photonCutArray[ 6] = "01200009297002008250400000"; mesonCutArray[ 6] = "0152406500000000";
  } else if ( trainConfig == 179){ // flow cuts with eta = 0.9, y = 0.85
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009297002208250400000"; mesonCutArray[ 0] = "0152506500000000";
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009297002208250400000"; mesonCutArray[ 1] = "0152506500000000";
    eventCutArray[ 2] = "51200013"; photonCutArray[ 2] = "00200009297002208250400000"; mesonCutArray[ 2] = "0152506500000000";
    eventCutArray[ 3] = "52300013"; photonCutArray[ 3] = "00200009297002208250400000"; mesonCutArray[ 3] = "0152506500000000";
    eventCutArray[ 4] = "53400013"; photonCutArray[ 4] = "00200009297002208250400000"; mesonCutArray[ 4] = "0152506500000000";
    eventCutArray[ 5] = "54600013"; photonCutArray[ 5] = "00200009297002208250400000"; mesonCutArray[ 5] = "0152506500000000";
    eventCutArray[ 6] = "56800013"; photonCutArray[ 6] = "00200009297002208250400000"; mesonCutArray[ 6] = "0152506500000000";
  } else if ( trainConfig == 180){ // flow cuts with eta = 0.65, y = 0.6
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "03200009297002208250400000"; mesonCutArray[ 0] = "0152306500000000";
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "03200009297002208250400000"; mesonCutArray[ 1] = "0152306500000000";
    eventCutArray[ 2] = "51200013"; photonCutArray[ 2] = "03200009297002208250400000"; mesonCutArray[ 2] = "0152306500000000";
    eventCutArray[ 3] = "52300013"; photonCutArray[ 3] = "03200009297002208250400000"; mesonCutArray[ 3] = "0152306500000000";
    eventCutArray[ 4] = "53400013"; photonCutArray[ 4] = "03200009297002208250400000"; mesonCutArray[ 4] = "0152306500000000";
    eventCutArray[ 5] = "54600013"; photonCutArray[ 5] = "03200009297002208250400000"; mesonCutArray[ 5] = "0152306500000000";
    eventCutArray[ 6] = "56800013"; photonCutArray[ 6] = "03200009297002208250400000"; mesonCutArray[ 6] = "0152306500000000";
  } else if ( trainConfig == 181){ // flow cuts with eta = 0.6, y = 0.5
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "01200009297002208250400000"; mesonCutArray[ 0] = "0152406500000000";
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "01200009297002208250400000"; mesonCutArray[ 1] = "0152406500000000";
    eventCutArray[ 2] = "51200013"; photonCutArray[ 2] = "01200009297002208250400000"; mesonCutArray[ 2] = "0152406500000000";
    eventCutArray[ 3] = "52300013"; photonCutArray[ 3] = "01200009297002208250400000"; mesonCutArray[ 3] = "0152406500000000";
    eventCutArray[ 4] = "53400013"; photonCutArray[ 4] = "01200009297002208250400000"; mesonCutArray[ 4] = "0152406500000000";
    eventCutArray[ 5] = "54600013"; photonCutArray[ 5] = "01200009297002208250400000"; mesonCutArray[ 5] = "0152406500000000";
    eventCutArray[ 6] = "56800013"; photonCutArray[ 6] = "01200009297002208250400000"; mesonCutArray[ 6] = "0152406500000000";
  } else if ( trainConfig == 182){ // standard LHC11h cut selection -> for centr. flattening
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
//     eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; //with fDoCentralityFlat = 2
//     eventCutArray[ 4] = "50800013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; //with fDoCentralityFlat = 8
  } else if ( trainConfig == 183){ // standard LHC11h cut selection - added signal -> for centr. flattening
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
//     eventCutArray[ 3] = "51200023"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; //with fDoCentralityFlat = 2
//     eventCutArray[ 4] = "50800023"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; //with fDoCentralityFlat = 8
  } else if ( trainConfig == 184){ // standard LHC11h cut selection with phi cut -> for centr. flattening
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
//     eventCutArray[ 3] = "51200013"; photonCutArray[ 3] = "00216609247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; //with fDoCentralityFlat = 2
//     eventCutArray[ 4] = "50800013"; photonCutArray[ 4] = "00216609247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; //with fDoCentralityFlat = 8
  } else if ( trainConfig == 185){ // standard LHC11h cut selection with phi cut - added signal -> for centr. flattening
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
//     eventCutArray[ 3] = "51200023"; photonCutArray[ 3] = "00216609247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; //with fDoCentralityFlat = 2
//     eventCutArray[ 4] = "50800023"; photonCutArray[ 4] = "00216609247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; //with fDoCentralityFlat = 8
  } else if ( trainConfig == 186){ // variation with phi cut at 2.0 - 4.0
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00215509247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00215509247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00215509247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00215509247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00215509247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 187){ // variation with phi cut at 2.0 - 4.0 - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00215509247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00215509247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00215509247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00215509247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00215509247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 188){ // variation with phi cut at 2.4 - 3.6
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00217709247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00217709247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00217709247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00217709247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00217709247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 189){ // variation with phi cut at 2.4 - 3.6 - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00217709247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00217709247602008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00217709247602008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00217709247602008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00217709247602008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 190){ // open an. cut studies 
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // std: 0 - 3.14 (pi)
    eventCutArray[ 1] = "50100013"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152501500000002"; // 0 - pt dep
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152501500000022"; // pt dep - pt dep
    eventCutArray[ 3] = "50100013"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152501500000001"; // 0 - pt dep
    eventCutArray[ 4] = "50100013"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152501500000021"; // pt dep - pt dep
  } else if ( trainConfig == 191){ // open an. cut studies - added signal
    eventCutArray[ 0] = "50100023"; photonCutArray[ 0] = "00200009247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // std: 0 - 3.14 (pi)
    eventCutArray[ 1] = "50100023"; photonCutArray[ 1] = "00200009247602008250400000"; mesonCutArray[ 1] = "0152501500000002"; // 0 - pt dep
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602008250400000"; mesonCutArray[ 2] = "0152501500000022"; // pt dep - pt dep
    eventCutArray[ 3] = "50100023"; photonCutArray[ 3] = "00200009247602008250400000"; mesonCutArray[ 3] = "0152501500000001"; // 0 - pt dep
    eventCutArray[ 4] = "50100023"; photonCutArray[ 4] = "00200009247602008250400000"; mesonCutArray[ 4] = "0152501500000021"; // pt dep - pt dep
  } else if ( trainConfig == 192){ // open an. cut studies - with phi cut 
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00216609247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // std: 0 - 3.14 (pi)
    eventCutArray[ 1] = "50100013"; photonCutArray[ 1] = "00216609247602008250400000"; mesonCutArray[ 1] = "0152501500000002"; // 0 - pt dep
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602008250400000"; mesonCutArray[ 2] = "0152501500000022"; // pt dep - pt dep
    eventCutArray[ 3] = "50100013"; photonCutArray[ 3] = "00216609247602008250400000"; mesonCutArray[ 3] = "0152501500000001"; // 0 - pt dep
    eventCutArray[ 4] = "50100013"; photonCutArray[ 4] = "00216609247602008250400000"; mesonCutArray[ 4] = "0152501500000021"; // pt dep - pt dep
  } else if ( trainConfig == 193){ // open an. cut studies - with phi cut added signals
    eventCutArray[ 0] = "50100023"; photonCutArray[ 0] = "00216609247602008250400000"; mesonCutArray[ 0] = "0152501500000000"; // std: 0 - 3.14 (pi)
    eventCutArray[ 1] = "50100023"; photonCutArray[ 1] = "00216609247602008250400000"; mesonCutArray[ 1] = "0152501500000002"; // 0 - pt dep
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602008250400000"; mesonCutArray[ 2] = "0152501500000022"; // pt dep - pt dep
    eventCutArray[ 3] = "50100023"; photonCutArray[ 3] = "00216609247602008250400000"; mesonCutArray[ 3] = "0152501500000001"; // 0 - pt dep
    eventCutArray[ 4] = "50100023"; photonCutArray[ 4] = "00216609247602008250400000"; mesonCutArray[ 4] = "0152501500000021"; // pt dep - pt dep
  } else if ( trainConfig == 194){ // standard LHC11h cut selection - no dEdx
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009000002008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009000002008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009000002008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009000002008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009000002008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 195){ // standard LHC11h cut selection - no dEdx add signals
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009000002008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009000002008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009000002008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009000002008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009000002008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 196){ // standard LHC11h cut selection - no dEdx with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609000002008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609000002008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609000002008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609000002008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609000002008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 197){ // standard LHC11h cut selection - no dEdx with phi cut add signals
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609000002008250400000"; mesonCutArray[ 0] = "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609000002008250400000"; mesonCutArray[ 1] = "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609000002008250400000"; mesonCutArray[ 2] = "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609000002008250400000"; mesonCutArray[ 3] = "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609000002008250400000"; mesonCutArray[ 4] = "0152501500000000"; // 20-50%
  } else if ( trainConfig == 198){ // V0 f. eff., open cut//00200009247602008250400000
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00000070000000000500004000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 199){ // V0 f. eff., open cut
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00000070000000000500004000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 200){ // V0 f. eff., el. dEdx
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00000070200000000500004000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 201){ // V0 f. eff., el. dEdx
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00000070200000000500004000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 202){ // V0 f. eff., chi2/psi pair
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00000070000000000250004000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 203){ // V0 f. eff., chi2/psi pair
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00000070000000000250004000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 204){ // V0 f. eff., el. dEdx and chi2/psi pair
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00000070200000000250004000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 205){ // V0 f. eff., el. dEdx and chi2/psi pair
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00000070200000000250004000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 206){ // V0 f. eff., single pt, el. dEdx and chi2/psi pair
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00000000200000000250004000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 207){ // V0 f. eff., single pt, el. dEdx and chi2/psi pair
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00000000200000000250004000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 208){ // V0 f. eff., single pt, el. dEdx and chi2/psi pair, cosPA
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00000000200000000250404000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 209){ // V0 f. eff., single pt, el. dEdx and chi2/psi pair, cosPA
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00000000200000000250404000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 210){ // V0 f. eff., ... plus qt cut  
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00000000200000008250404000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 211){ // V0 f. eff., ... plus qt cut 
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00000000200000008250404000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 212){ // V0 f. eff., ... plus pion dedx
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00000000247600008250404000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 213){ // V0 f. eff., ... plus pion dedx 
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00000000247600008250404000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 214){ // V0 f. eff., ... plus R min
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00200000247600008250404000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 215){ // V0 f. eff., ... plus R min 
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00200000247600008250404000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 216){ // V0 f. eff., ... plus phi cut
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00216600247600008250404000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 217){ // V0 f. eff., ... plus phi cut 
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00216600247600008250404000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 218){ // V0 f. eff., only phi cut
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00016670000000000500004000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 219){ // V0 f. eff., only phi cut 
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00016670000000000500004000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 220){ // V0 f. eff., all but psipairchi2 cut
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00200009247602008500400000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 221){ // V0 f. eff., all but psipairchi2 cut 
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00200009247602008500400000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 222){ // V0 f. eff., all but qt cut
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00200009247602000250400000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 223){ // V0 f. eff., all but qt cut 
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00200009247602000250400000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 224){ // V0 f. eff., all but psipairchi2 and qt cut
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00200009247602000500400000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 225){ // V0 f. eff., all but psipairchi2 and qt cut 
    eventCutArray[ 0] = "52500013"; photonCutArray[ 0] = "00200009247602000500400000"; mesonCutArray[ 0] = "0142501500000000";
  } else if ( trainConfig == 226){ // standard LHC11h cut selection - double rejec
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602008250404000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602008250404000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602008250404000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008250404000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008250404000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 227){ // standard LHC11h cut selection - double rejec - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602008250404000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602008250404000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602008250404000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247602008250404000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247602008250404000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 228){ // standard LHC11h cut selection - double rejec with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247602008250404000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247602008250404000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602008250404000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247602008250404000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247602008250404000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 229){ // standard LHC11h cut selection - double rejec with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247602008250404000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247602008250404000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602008250404000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247602008250404000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247602008250404000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 230){ // standard LHC11h cut selection - double rejec, no Psi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602008200404000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602008200404000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602008200404000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008200404000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008200404000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 231){ // standard LHC11h cut selection - double rejec, no Psi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602008200404000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602008200404000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602008200404000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247602008200404000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247602008200404000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 232){ // standard LHC11h cut selection - double rejec, no Psi cut with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247602008200404000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247602008200404000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602008200404000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247602008200404000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247602008200404000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 233){ // standard LHC11h cut selection - double rejec, no Psi cut with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247602008200404000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247602008200404000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602008200404000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247602008200404000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247602008200404000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 234){ // standard LHC11h cut selection - double rejec, no Psi cut with phi cut - added signal
    eventCutArray[ 0] = "56800013"; photonCutArray[ 0] = "00200009247602008250404000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
  } else if ( trainConfig == 235){ // standard LHC11h cut selection - double rejec, no chi2
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00200009247602008510404000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00200009247602008510404000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009247602008510404000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00200009247602008510404000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00200009247602008510404000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 236){ // standard LHC11h cut selection - double rejec,  no chi2 - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00200009247602008510404000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00200009247602008510404000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00200009247602008510404000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00200009247602008510404000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00200009247602008510404000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 237){ // standard LHC11h cut selection - double rejec,  no chi2 cut with phi cut
    eventCutArray[ 0] = "60100013"; photonCutArray[ 0] = "00216609247602008510404000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200013"; photonCutArray[ 1] = "00216609247602008510404000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609247602008510404000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400013"; photonCutArray[ 3] = "00216609247602008510404000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500013"; photonCutArray[ 4] = "00216609247602008510404000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 238){ // standard LHC11h cut selection - double rejec,  no chi2 cut with phi cut - added signal
    eventCutArray[ 0] = "60100023"; photonCutArray[ 0] = "00216609247602008510404000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "61200023"; photonCutArray[ 1] = "00216609247602008510404000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100023"; photonCutArray[ 2] = "00216609247602008510404000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "52400023"; photonCutArray[ 3] = "00216609247602008510404000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "52500023"; photonCutArray[ 4] = "00216609247602008510404000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 239){ // standard LHC11h cut selection - double rejec
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00200009247602008250404000"; mesonCutArray[ 0]= "0152501500000000"; // 0-5%
    eventCutArray[ 1] = "50100013"; photonCutArray[ 1] = "00200009245602008250404000"; mesonCutArray[ 1]= "0152501500000000"; // 5-10%
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00200009246602008250404000"; mesonCutArray[ 2]= "0152501500000000"; // 0-10%
    eventCutArray[ 3] = "50100013"; photonCutArray[ 3] = "00200009240602008250404000"; mesonCutArray[ 3]= "0152501500000000"; // 20-40%
    eventCutArray[ 4] = "50100013"; photonCutArray[ 4] = "00200009248602008250404000"; mesonCutArray[ 4]= "0152501500000000"; // 20-50%
  } else if ( trainConfig == 240){ // standard LHC11h cut selection - double rejec with phi cut
    eventCutArray[ 0] = "50100013"; photonCutArray[ 0] = "00216609247602008250404000"; mesonCutArray[ 0]= "0152501500000000"; // min mom 0.4
    eventCutArray[ 1] = "50100013"; photonCutArray[ 1] = "00216609245602008250404000"; mesonCutArray[ 1]= "0152501500000000"; // 0.3
    eventCutArray[ 2] = "50100013"; photonCutArray[ 2] = "00216609246602008250404000"; mesonCutArray[ 2]= "0152501500000000"; // 0.25
    eventCutArray[ 3] = "50100013"; photonCutArray[ 3] = "002166092407602008250404000"; mesonCutArray[ 3]= "0152501500000000"; // 0.5
    eventCutArray[ 4] = "50100013"; photonCutArray[ 4] = "00216609248602008250404000"; mesonCutArray[ 4]= "0152501500000000"; // 0.2
  }else {
    Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

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
  }

  EventCutList->SetOwner(kTRUE);
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];


  for(Int_t i = 0; i<numberOfCuts; i++){

    analysisEventCuts[i] = new AliConvEventCuts();
    
    if(periodName.CompareTo("LHC11h") && (doFlattening > 0)){
      cout << "entering the flattening loop -> searching for file: " << fileNameInputForCentFlattening.Data() << endl;
      
      if( fileNameInputForCentFlattening.Contains("Low") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentLowRange");
      }else if( fileNameInputForCentFlattening.Contains("Middle") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentMiddleRange");
      }else if( fileNameInputForCentFlattening.Contains("High") ){
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "CentHighRange");
      }else {
        analysisEventCuts[i]->SetUseWeightFlatCentralityFromFile(doFlattening, fileNameInputForCentFlattening, "Cent");
      }
    }
    
    if (  trainConfig == 1   || trainConfig == 5   || trainConfig == 9   || trainConfig == 13   || trainConfig == 17   || 
        trainConfig == 21   || trainConfig == 25   || trainConfig == 29   || trainConfig == 33   || trainConfig == 37   ||
        trainConfig == 202  ){ // || trainConfig == 41 
      if (i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
      if (i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
      if (i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
      if (i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      if (i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
    } else if (  trainConfig == 2   || trainConfig == 6   || trainConfig == 10   || trainConfig == 14   || trainConfig == 18   ||
          trainConfig == 22   || trainConfig == 26   || trainConfig == 30   || trainConfig == 34   || trainConfig == 38   ||
          trainConfig == 203  ){ // || trainConfig == 42
      if (i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_2040TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M");
      if (i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4060TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4060V0M");
      if (i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_6080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_6080V0M");
      if (i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_PbPb_2760GeV_4080TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_4080V0M");
    } else if ( trainConfig == 3   || trainConfig == 7    || trainConfig == 11   || trainConfig == 15  || trainConfig == 19   || 
          trainConfig == 23   || trainConfig == 27   || trainConfig == 31   || trainConfig == 35   || trainConfig == 39   || 
          trainConfig == 204  ){ //|| trainConfig == 43 
      if (i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0005TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M");
      if (i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0510TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M");
      if (i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0010TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M");
      if (i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE,fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_1020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_1020V0M");
      if (i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kFALSE, kFALSE, fileNameInputForWeighting, "Pi0_Hijing_LHC13d2_addSig_PbPb_2760GeV_0020TPC", "", "","Pi0_Fit_Data_PbPb_2760GeV_0020V0M");
    } else if (  trainConfig == 4   ||trainConfig == 8     || trainConfig == 12   || trainConfig == 16   || trainConfig == 20   || 
          trainConfig == 24   || trainConfig == 28   || trainConfig == 32   || trainConfig == 36   || trainConfig == 40   || 
          trainConfig == 205){ // || trainConfig == 44 
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
    
    if (   trainConfig == 70   || trainConfig == 72    || trainConfig == 74    || trainConfig == 76   || trainConfig == 78   || 
        trainConfig == 80    || trainConfig == 82   || trainConfig == 84   || trainConfig == 86   || trainConfig == 88   ||
        trainConfig == 90   || trainConfig == 92  || trainConfig == 94   || trainConfig == 96   || trainConfig == 98   || 
        trainConfig == 100   || trainConfig == 102   || trainConfig == 104   || trainConfig == 106   || trainConfig == 108   || 
        trainConfig == 110   || trainConfig == 112   || trainConfig == 114   || trainConfig == 116   || trainConfig == 118   ||
        trainConfig == 120   || trainConfig == 122   || trainConfig == 124   || trainConfig == 126   || trainConfig == 128   || 
        trainConfig == 130   || trainConfig == 132   || trainConfig == 134  || trainConfig == 136   || trainConfig == 138   || 
        trainConfig == 140   || trainConfig == 142   || trainConfig == 144   || trainConfig == 146   || trainConfig == 148   ||
        trainConfig == 150   || trainConfig == 152   || trainConfig == 154   || trainConfig == 156   || trainConfig == 158   || 
        trainConfig == 160   || trainConfig == 162   || trainConfig == 164   || trainConfig == 166   || trainConfig == 168   ||
        trainConfig == 170   || trainConfig == 172   || trainConfig == 174   || trainConfig == 182   || trainConfig == 184   || 
        trainConfig == 186   || trainConfig == 188   || trainConfig == 194   || trainConfig == 196   || trainConfig == 226   || 
        trainConfig == 228  || trainConfig == 230 || trainConfig == 232		 || trainConfig == 235 || trainConfig == 237){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }
        
    if (   trainConfig == 71   || trainConfig == 73    || trainConfig == 75    || trainConfig == 77    || trainConfig == 79    ||
        trainConfig == 81    || trainConfig == 83   || trainConfig == 85   || trainConfig == 87    || trainConfig == 89    || 
        trainConfig == 91   || trainConfig == 93   || trainConfig == 95   || trainConfig == 97    || trainConfig == 99   ||
        trainConfig == 101   || trainConfig == 103    || trainConfig == 105   || trainConfig == 107   || trainConfig == 109   || 
        trainConfig == 111   || trainConfig == 113   || trainConfig == 115   || trainConfig == 117   || trainConfig == 119   ||
        trainConfig == 121   || trainConfig == 123   || trainConfig == 125   || trainConfig == 127   || trainConfig == 129   ||
        trainConfig == 131   || trainConfig == 133   || trainConfig == 135   || trainConfig == 137   || trainConfig == 139   ||
        trainConfig == 141   || trainConfig == 143   || trainConfig == 145   || trainConfig == 147   || trainConfig == 149   || 
        trainConfig == 151   || trainConfig == 153   || trainConfig == 155   || trainConfig == 157   || trainConfig == 159   || 
        trainConfig == 161   || trainConfig == 163   || trainConfig == 165   || trainConfig == 167   || trainConfig == 169   || 
        trainConfig == 171   || trainConfig == 173   || trainConfig == 175   || trainConfig == 183   || trainConfig == 185   ||
        trainConfig == 187   || trainConfig == 189   || trainConfig == 195   || trainConfig == 197  || trainConfig == 227   || 
        trainConfig == 229  || trainConfig == 231 || trainConfig == 233		 || trainConfig == 236 || trainConfig == 238){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2040TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2040V0M","Eta_Fit_Data_PbPb_2760GeV_2040V0M");
        if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_2050TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_2050V0M","Eta_Fit_Data_PbPb_2760GeV_2050V0M");
      }
    }
      
//     if ( trainConfig == 182 || trainConfig == 184){
//       if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
//         if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
//         if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
//         if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
//       }
//     }
// 
//     if ( trainConfig == 183 || trainConfig == 185){
//       if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
//         if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0005TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0005V0M","Eta_Fit_Data_PbPb_2760GeV_0005V0M");
//         if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0510TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0510V0M","Eta_Fit_Data_PbPb_2760GeV_0510V0M");
//         if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
//       }
//     }
    
    if(trainConfig == 190 || trainConfig == 192 || trainConfig == 239 || trainConfig == 240 ){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }

    if(trainConfig == 191 || trainConfig == 193){
      if (periodName.CompareTo("LHC14a1a") ==0 || periodName.CompareTo("LHC14a1b") ==0 || periodName.CompareTo("LHC14a1c") ==0 ){
        if ( i == 0 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 1 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 2 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 3 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
        if ( i == 4 && doWeighting)  analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE,fileNameInputForWeighting, Form("Pi0_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), Form("Eta_Hijing_%s_addSig_PbPb_2760GeV_0010TPC",periodName.Data()), "","Pi0_Fit_Data_PbPb_2760GeV_0010V0M","Eta_Fit_Data_PbPb_2760GeV_0010V0M");
      }
    }

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
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
    analysisCuts[i]->InitializeCutsFromCutString(photonCutArray[i].Data());
    if( trainConfig == 198 || trainConfig == 199 || trainConfig == 202 || trainConfig == 203){
      analysisCuts[i]->SetDodEdxSigmaCut(kFALSE);
    }
    ConvCutList->Add(analysisCuts[i]);
    analysisCuts[i]->SetFillCutHistograms("",kFALSE);
    


    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
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
    
  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
