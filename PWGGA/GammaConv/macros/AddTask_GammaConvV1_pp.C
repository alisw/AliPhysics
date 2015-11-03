void AddTask_GammaConvV1_pp(  Int_t   trainConfig                 = 1,                    // change different set of cuts
                              Int_t   isMC                        = 0,                    // run MC
                              Int_t   enableQAMesonTask           = 0,                    // enable meson QA in AliAnalysisTaskGammaConvV1
                              Int_t   enableQAPhotonTask          = 0,                    // enable photon QA in AliAnalysisTaskGammaConvV1
                              TString fileNameInputForWeighting 	= "MCSpectraInput.root",// path to file for weigting input
                              TString cutnumberAODBranch          = "000000006008400001001500000",  // cutnumber for AOD branch
                              TString periodname                  = "LHC12f1x",           // period name
                              Bool_t  doWeighting                 = kFALSE,               // enables weighting
                              Bool_t  enableV0findingEffi         = kFALSE,               // enables V0finding efficiency histograms
                              Bool_t  enableTriggerMimicking      = kFALSE,               // enable trigger mimicking
                              Bool_t  enableTriggerOverlapRej     = kFALSE,               // enable trigger overlap rejection
                              Float_t maxFacPtHard                = 3.,                   // maximum factor between hardest jet and ptHard generated
                              TString periodNameV0Reader          = ""
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

  Int_t isHeavyIon = 0;
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
  // Cut Numbers to use in Analysis
  Int_t numberOfCuts = 4;
  if ( trainConfig == 70) numberOfCuts = 2;
  if ( trainConfig == 100) numberOfCuts = 5;
  
  TString *eventCutArray = new TString[numberOfCuts];
  TString *photonCutArray = new TString[numberOfCuts];
  TString *mesonCutArray = new TString[numberOfCuts];
  
  if(trainConfig == 1){
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD , only Minbias MC
    eventCutArray[ 1] = "00012113"; photonCutArray[ 1] = "00200009366300003800000000"; mesonCutArray[1] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD, V0AND
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009326000003800000000"; mesonCutArray[2] = "0163103100900000"; //standard cut Gamma pp 2-76TeV
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009366000003800000000"; mesonCutArray[3] = "0163103100900000"; //standard cut Gamma pp 2-76TeV
  } else if (trainConfig == 2) {
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD , only boxes
    eventCutArray[ 1] = "00012123"; photonCutArray[ 1] = "00200009366300003800000000"; mesonCutArray[1] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD, V0AND , only boxes
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009326000003800000000"; mesonCutArray[2] = "0163103100900000"; //standard cut Gamma pp 2-76TeV , only boxes
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009366000003800000000"; mesonCutArray[3] = "0163103100900000"; //standard cut Gamma pp 2-76TeV 
  } else if (trainConfig == 3) {
    eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD , only Minbias MC
    eventCutArray[ 1] = "00013113"; photonCutArray[ 1] = "00200009366300003800000000"; mesonCutArray[1] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD, V0AND , only Minbias MC
    eventCutArray[ 2] = "00003123"; photonCutArray[ 2] = "00200009366300003800000000"; mesonCutArray[2] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD , only Boxes MC
    eventCutArray[ 3] = "00013123"; photonCutArray[ 3] = "00200009366300003800000000"; mesonCutArray[3] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD, V0AND, only Boxes MC
  } else if (trainConfig == 4) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD , all photon qualities
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009366300003800020000"; mesonCutArray[1] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 1
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009366300003800030000"; mesonCutArray[2] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 2
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009366300003800040000"; mesonCutArray[3] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 3
  } else if (trainConfig == 5) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00700009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD , all photon qualities, min R = 35 cm
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00700009366300003800020000"; mesonCutArray[1] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 1, min R = 35 cm
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00700009366300003800030000"; mesonCutArray[2] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 2, min R = 35 cm
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00700009366300003800040000"; mesonCutArray[3] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV without SDD , photon quality 3, min R = 35 cm
  } else if (trainConfig == 6) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200008366300003200000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, all photon qualities
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200008366300003200020000"; mesonCutArray[1] = "0163103100900000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 1
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200008366300003200030000"; mesonCutArray[2] = "0163103100900000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 2
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200008366300003200040000"; mesonCutArray[3] = "0163103100900000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 3   
  } else if (trainConfig == 7) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00700008366300003200000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, all photon qualities, min R = 35 cm
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00700008366300003200020000"; mesonCutArray[1] = "0163103100900000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 1, min R = 35 cm
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00700008366300003200030000"; mesonCutArray[2] = "0163103100900000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 2, min R = 35 cm
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00700008366300003200040000"; mesonCutArray[3] = "0163103100900000"; //standard cut Pi0 pp 7TeV, with qt cut 0.05, photon quality 3, min R = 35 cm
  } else if (trainConfig == 8) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200008366300000200000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 7TeV, all photon qualities
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200008366300000200020000"; mesonCutArray[1] = "0163103100900000"; //standard cut Pi0 pp 7TeV, photon quality 1
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200008366300000200030000"; mesonCutArray[2] = "0163103100900000"; //standard cut Pi0 pp 7TeV, photon quality 2
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200008366300000200040000"; mesonCutArray[3] = "0163103100900000"; //standard cut Pi0 pp 7TeV, photon quality 3   
  } else if (trainConfig == 9) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00700008366300000200000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 7TeV, all photon qualities, min R = 35 cm
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00700008366300000200020000"; mesonCutArray[1] = "0163103100900000"; //standard cut Pi0 pp 7TeV, photon quality 1, min R = 35 cm
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00700008366300000200030000"; mesonCutArray[2] = "0163103100900000"; //standard cut Pi0 pp 7TeV, photon quality 2, min R = 35 cm
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00700008366300000200040000"; mesonCutArray[3] = "0163103100900000"; //standard cut Pi0 pp 7TeV, photon quality 3, min R = 35 cm	   
  } else if (trainConfig == 10) {
    eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD , all photon qualities
    eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00200009366300003800020000"; mesonCutArray[1] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 1
    eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00200009366300003800030000"; mesonCutArray[2] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 2
    eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00200009366300003800040000"; mesonCutArray[3] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 3
  } else if (trainConfig == 11) {
    eventCutArray[ 0] = "00003113"; photonCutArray[ 0] = "00700009366300003800000000"; mesonCutArray[0] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD , all photon qualities, min R = 35 cm
    eventCutArray[ 1] = "00003113"; photonCutArray[ 1] = "00700009366300003800020000"; mesonCutArray[1] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 1, min R = 35 cm
    eventCutArray[ 2] = "00003113"; photonCutArray[ 2] = "00700009366300003800030000"; mesonCutArray[2] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 2, min R = 35 cm
    eventCutArray[ 3] = "00003113"; photonCutArray[ 3] = "00700009366300003800040000"; mesonCutArray[3] = "0163103100900000"; //standard cut Pi0 pp 2.76TeV with SDD , photon quality 3, min R = 35 cm
  } else if (trainConfig == 12) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[0] = "0152506500000000"; //standard cut LHC11h pp 2.76TeV 
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "03200009297002008250400000"; mesonCutArray[1] = "0152506500000000"; //variation eta 0.65
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "04200009297002008250400000"; mesonCutArray[2] = "0152506500000000"; //variation eta 0.75
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009295002008250400000"; mesonCutArray[3] = "0152506500000000"; //variation pion p dEdx 0.3-5.
  } else if (trainConfig == 13) { //added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[0] = "0152506500000000"; //standard cut LHC11h pp 2.76TeV 
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "03200009297002008250400000"; mesonCutArray[1] = "0152506500000000"; //variation eta 0.65
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "04200009297002008250400000"; mesonCutArray[2] = "0152506500000000"; //variation eta 0.75
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009295002008250400000"; mesonCutArray[3] = "0152506500000000"; //variation pion p dEdx 0.3-5.
  } else if (trainConfig == 14) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200049297002008250400000"; mesonCutArray[0] = "0152506500000000"; //variation pt 0.075 
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200019297002008250400000"; mesonCutArray[1] = "0152506500000000"; //variation pt 0.1
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200006297002008250400000"; mesonCutArray[2] = "0152506500000000"; //variation TPC cls 0.7
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200008297002008250400000"; mesonCutArray[3] = "0152506500000000"; //variation TPC cls 0.35 
  } else if (trainConfig == 15) { //added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200049297002008250400000"; mesonCutArray[0] = "0152506500000000"; //variation pt 0.075 
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200019297002008250400000"; mesonCutArray[1] = "0152506500000000"; //variation pt 0.1
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200006297002008250400000"; mesonCutArray[2] = "0152506500000000"; //variation TPC cls 0.7
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200008297002008250400000"; mesonCutArray[3] = "0152506500000000"; //variation TPC cls 0.35 
  } else if (trainConfig == 16) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009397002008250400000"; mesonCutArray[0] = "0152506500000000"; //variation edEdx -4,5
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009697002008250400000"; mesonCutArray[1] = "0152506500000000"; //variation edEdx -2.5,4
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009297003008250400000"; mesonCutArray[2] = "0152506500000000"; //variation TOF el. PID -3,5
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009297004008250400000"; mesonCutArray[3] = "0152506500000000"; //variation TOF el. PID -2,3
  } else if (trainConfig == 17) { //added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009397002008250400000"; mesonCutArray[0] = "0152506500000000"; //variation edEdx -4,5
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009697002008250400000"; mesonCutArray[1] = "0152506500000000"; //variation edEdx -2.5,4
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009297003008250400000"; mesonCutArray[2] = "0152506500000000"; //variation TOF el. PID -3,5
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009297004008250400000"; mesonCutArray[3] = "0152506500000000"; //variation TOF el. PID -2,3
  } else if (trainConfig == 18) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009297002009250400000"; mesonCutArray[0] = "0152506500000000"; //variation qt 0.03
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009297002002250400000"; mesonCutArray[1] = "0152506500000000"; //variation qt 0.07 no2D
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009297002008150400000"; mesonCutArray[2] = "0152506500000000"; //variation chi2 50.
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009297002008850400000"; mesonCutArray[3] = "0152506500000000"; //variation chi2 20.
  } else if (trainConfig == 19) { //added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009297002009250400000"; mesonCutArray[0] = "0152506500000000"; //variation qt 0.03
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009297002002250400000"; mesonCutArray[1] = "0152506500000000"; //variation qt 0.07 no2D
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009297002008150400000"; mesonCutArray[2] = "0152506500000000"; //variation chi2 50.
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009297002008850400000"; mesonCutArray[3] = "0152506500000000"; //variation chi2 20.
  } else if (trainConfig == 20) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009297002008260400000"; mesonCutArray[0] = "0152506500000000"; //variation psi pair 0.05
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009297002008280400000"; mesonCutArray[1] = "0152506500000000"; //variation psi pair 0.2
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009297002008250000000"; mesonCutArray[2] = "0152506500000000"; //variation cosPA -1
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009297002008250400000"; mesonCutArray[3] = "0152505500000000"; //variation alpha 0.75
  } else if (trainConfig == 21) { //added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009297002008260400000"; mesonCutArray[0] = "0152506500000000"; //variation psi pair 0.05
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009297002008280400000"; mesonCutArray[1] = "0152506500000000"; //variation psi pair 0.2
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009297002008250000000"; mesonCutArray[2] = "0152506500000000"; //variation cosPA -1
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009297002008250400000"; mesonCutArray[3] = "0152505500000000"; //variation alpha 0.75
  } else if (trainConfig == 22) {
    eventCutArray[ 0] = "00040113"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[0] = "0152506500000000"; // trigger kTRD
    eventCutArray[ 1] = "00050113"; photonCutArray[ 1] = "00200009297002008250400000"; mesonCutArray[1] = "0152506500000000"; // trigger kEMC
    eventCutArray[ 2] = "00060113"; photonCutArray[ 2] = "00200009297002008250400000"; mesonCutArray[2] = "0152506500000000"; // trigger kPHI
    eventCutArray[ 3] = "00070113"; photonCutArray[ 3] = "00200009297002008250400000"; mesonCutArray[3] = "0152506500000000"; // trigger kHighMult
  } else if (trainConfig == 23) {
    eventCutArray[ 0] = "00080113"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[0] = "0152506500000000"; // trigger kEMCEGA
    eventCutArray[ 1] = "00090113"; photonCutArray[ 1] = "00200009297002008250400000"; mesonCutArray[1] = "0152506500000000"; // trigger kEMCEJE
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009297002008250400000"; mesonCutArray[2] = "0152506500000000"; // minimum bias
    eventCutArray[ 3] = "00011113"; photonCutArray[ 3] = "00200009297002008250400000"; mesonCutArray[3] = "0152506500000000"; // trigger kINT8
  } else if (trainConfig == 24) {
    eventCutArray[ 0] = "00042113"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[0] = "0152506500000000"; // trigger kTRD CINT8 HEE
    eventCutArray[ 1] = "00044113"; photonCutArray[ 1] = "00200009297002008250400000"; mesonCutArray[1] = "0152506500000000"; // trigger kTRD CINT8 HSE
    eventCutArray[ 2] = "00046113"; photonCutArray[ 2] = "00200009297002008250400000"; mesonCutArray[2] = "0152506500000000"; // trigger kTRD CINT8 HJE
    eventCutArray[ 3] = "00048113"; photonCutArray[ 3] = "00200009297002008250400000"; mesonCutArray[3] = "0152506500000000"; // trigger kTRD CINT8 HQU
  } else if (trainConfig == 25) {
    eventCutArray[ 0] = "00041113"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[0] = "0152506500000000"; // trigger kTRD CINT7 HEE
    eventCutArray[ 1] = "00043113"; photonCutArray[ 1] = "00200009297002008250400000"; mesonCutArray[1] = "0152506500000000"; // trigger kTRD CINT7 HSE
    eventCutArray[ 2] = "00045113"; photonCutArray[ 2] = "00200009297002008250400000"; mesonCutArray[2] = "0152506500000000"; // trigger kTRD CINT7 HJE
    eventCutArray[ 3] = "00047113"; photonCutArray[ 3] = "00200009297002008250400000"; mesonCutArray[3] = "0152506500000000"; // trigger kTRD CINT7 HQU
  } else if (trainConfig == 26) {
    eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[0] = "0152506500000000"; // trigger kEMC7
    eventCutArray[ 1] = "00053113"; photonCutArray[ 1] = "00200009297002008250400000"; mesonCutArray[1] = "0152506500000000"; // trigger kEMC8
    eventCutArray[ 2] = "00062113"; photonCutArray[ 2] = "00200009297002008250400000"; mesonCutArray[2] = "0152506500000000"; // trigger kPHI7
    eventCutArray[ 3] = "00063113"; photonCutArray[ 3] = "00200009297002008250400000"; mesonCutArray[3] = "0152506500000000"; // trigger kPHI8
  } else if (trainConfig == 27) {
    eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[0] = "0152506500000000"; // trigger kEMC1
    eventCutArray[ 1] = "00071113"; photonCutArray[ 1] = "00200009297002008250400000"; mesonCutArray[1] = "0152506500000000"; // trigger kSHM1
    eventCutArray[ 2] = "00072113"; photonCutArray[ 2] = "00200009297002008250400000"; mesonCutArray[2] = "0152506500000000"; // trigger kSHM7
    eventCutArray[ 3] = "00073113"; photonCutArray[ 3] = "00200009297002008250400000"; mesonCutArray[3] = "0152506500000000"; // trigger kSHM8
  } else if (trainConfig == 28) {
    eventCutArray[ 0] = "00081113"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[0] = "0152506500000000"; // trigger kEMCEGA + CINT7
    eventCutArray[ 1] = "00082113"; photonCutArray[ 1] = "00200009297002008250400000"; mesonCutArray[1] = "0152506500000000"; // trigger kEMCEGA + CINT8
    eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009297002008250400000"; mesonCutArray[2] = "0152506500000000"; // trigger kEMCEG1 + CINT7
    eventCutArray[ 3] = "00084113"; photonCutArray[ 3] = "00200009297002008250400000"; mesonCutArray[3] = "0152506500000000"; // trigger kEMCEG1 + CINT8
  } else if (trainConfig == 29) {
    eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[0] = "0152506500000000"; // trigger kEMCEG2 + CINT7
    eventCutArray[ 1] = "00086113"; photonCutArray[ 1] = "00200009297002008250400000"; mesonCutArray[1] = "0152506500000000"; // trigger kEMCEG2 + CINT8
    eventCutArray[ 2] = "00091113"; photonCutArray[ 2] = "00200009297002008250400000"; mesonCutArray[2] = "0152506500000000"; // trigger kEMCEJE + CINT7
    eventCutArray[ 3] = "00092113"; photonCutArray[ 3] = "00200009297002008250400000"; mesonCutArray[3] = "0152506500000000"; // trigger kEMCEJE + CINT8
  } else if (trainConfig == 30) {
    eventCutArray[ 0] = "00093113"; photonCutArray[ 0] = "00200009297002008250400000"; mesonCutArray[0] = "0152506500000000"; // trigger kEMCEJ1 + CINT7
    eventCutArray[ 1] = "00094113"; photonCutArray[ 1] = "00200009297002008250400000"; mesonCutArray[1] = "0152506500000000"; // trigger kEMCEJ1 + CINT8
    eventCutArray[ 2] = "00095113"; photonCutArray[ 2] = "00200009297002008250400000"; mesonCutArray[2] = "0152506500000000"; // trigger kEMCEJ2 + CINT7
    eventCutArray[ 3] = "00096113"; photonCutArray[ 3] = "00200009297002008250400000"; mesonCutArray[3] = "0152506500000000"; // trigger kEMCEJ2 + CINT8		
  } else if (trainConfig == 31) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009257002008250400000"; mesonCutArray[0] = "0152106500000000"; //new standard cut for pp 8 TeV
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009357002008250400000"; mesonCutArray[1] = "0152106500000000"; //variation edEdx -4,5
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009657002008250400000"; mesonCutArray[2] = "0152106500000000"; //variation edEdx -2.5,4
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009255002008250400000"; mesonCutArray[3] = "0152106500000000"; //variation pion p dEdx 0.3-5.
  } else if (trainConfig == 32) { //added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009257002008250400000"; mesonCutArray[0] = "0152106500000000"; //new standard cut for pp 8 TeV
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009357002008250400000"; mesonCutArray[1] = "0152106500000000"; //variation edEdx -4,5
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009657002008250400000"; mesonCutArray[2] = "0152106500000000"; //variation edEdx -2.5,4
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009255002008250400000"; mesonCutArray[3] = "0152106500000000"; //variation pion p dEdx 0.3-5.
  } else if (trainConfig == 33) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200049257002008250400000"; mesonCutArray[0] = "0152106500000000"; //variation pt 0.075 
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200019257002008250400000"; mesonCutArray[1] = "0152106500000000"; //variation pt 0.1
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200006257002008250400000"; mesonCutArray[2] = "0152106500000000"; //variation TPC cls 0.7
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200008257002008250400000"; mesonCutArray[3] = "0152106500000000"; //variation TPC cls 0.35 
  } else if (trainConfig == 34) { //added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200049257002008250400000"; mesonCutArray[0] = "0152106500000000"; //variation pt 0.075 
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200019257002008250400000"; mesonCutArray[1] = "0152106500000000"; //variation pt 0.1
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200006257002008250400000"; mesonCutArray[2] = "0152106500000000"; //variation TPC cls 0.7
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200008257002008250400000"; mesonCutArray[3] = "0152106500000000"; //variation TPC cls 0.35 
  } else if (trainConfig == 35) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227002008250400000"; mesonCutArray[0] = "0152106500000000"; //variation pidEdx 1,-10
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009237002008250400000"; mesonCutArray[1] = "0152106500000000"; //variation pidEdx 2.5,-10
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009297002008250400000"; mesonCutArray[2] = "0152106500000000"; //variation pidEdx 3,-10
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009250002008250400000"; mesonCutArray[3] = "0152106500000000"; //variation pion p dEdx 0.5-5
  } else if (trainConfig == 36) { //added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009227002008250400000"; mesonCutArray[0] = "0152106500000000"; //variation pidEdx 1,-10
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009237002008250400000"; mesonCutArray[1] = "0152106500000000"; //variation pidEdx 2.5,-10
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009297002008250400000"; mesonCutArray[2] = "0152106500000000"; //variation pidEdx 3,-10
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009250002008250400000"; mesonCutArray[3] = "0152106500000000"; //variation pion p dEdx 0.5-5
  } else if (trainConfig == 37) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009257002009250400000"; mesonCutArray[0] = "0152106500000000"; //variation qt 0.03
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009257002002250400000"; mesonCutArray[1] = "0152106500000000"; //variation qt 0.07 no2D
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009257002008150400000"; mesonCutArray[2] = "0152106500000000"; //variation chi2 50.
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009257002008850400000"; mesonCutArray[3] = "0152106500000000"; //variation chi2 20.
  } else if (trainConfig == 38) { //added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009257002009250400000"; mesonCutArray[0] = "0152106500000000"; //variation qt 0.03
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009257002002250400000"; mesonCutArray[1] = "0152106500000000"; //variation qt 0.07 no2D
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009257002008150400000"; mesonCutArray[2] = "0152106500000000"; //variation chi2 50.
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009257002008850400000"; mesonCutArray[3] = "0152106500000000"; //variation chi2 20.
  } else if (trainConfig == 39) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009257002008260400000"; mesonCutArray[0] = "0152106500000000"; //variation psi pair 0.05
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009257002008280400000"; mesonCutArray[1] = "0152106500000000"; //variation psi pair 0.2
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009257002008250000000"; mesonCutArray[2] = "0152106500000000"; //variation cosPA -1
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009257002008250600000"; mesonCutArray[3] = "0152106500000000"; //variation cosPA 0.9
  } else if (trainConfig == 40) { //added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009257002008260400000"; mesonCutArray[0] = "0152106500000000"; //variation psi pair 0.05
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009257002008280400000"; mesonCutArray[1] = "0152106500000000"; //variation psi pair 0.2
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009257002008250000000"; mesonCutArray[2] = "0152106500000000"; //variation cosPA -1
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009257002008250600000"; mesonCutArray[3] = "0152106500000000"; //variation cosPA 0.9
  } else if (trainConfig == 41) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009257002008950400000"; mesonCutArray[0] = "0152106500000000"; //variation chi2 15
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009257002008230400000"; mesonCutArray[1] = "0152106500000000"; //variation psi pair 0.035
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009257002008250400000"; mesonCutArray[2] = "0152105500000000"; //variation alpha 0.75
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009257002008250400000"; mesonCutArray[3] = "0152107500000000"; //variation alpha 0.85
  } else if (trainConfig == 42) { //added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009257002008950400000"; mesonCutArray[0] = "0152106500000000"; //variation chi2 15
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009257002008230400000"; mesonCutArray[1] = "0152106500000000"; //variation psi pair 0.035
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009257002008250400000"; mesonCutArray[2] = "0152105500000000"; //variation alpha 0.75
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009257002008250400000"; mesonCutArray[3] = "0152107500000000"; //variation alpha 0.85
  } else if (trainConfig == 43) {
    eventCutArray[ 0] = "00040113"; photonCutArray[ 0] = "00200009257002008250400000"; mesonCutArray[0] = "0152106500000000"; // trigger kTRD with y 0.8
    eventCutArray[ 1] = "00050113"; photonCutArray[ 1] = "00200009257002008250400000"; mesonCutArray[1] = "0152106500000000"; // trigger kEMC with y 0.8
    eventCutArray[ 2] = "00060113"; photonCutArray[ 2] = "00200009257002008250400000"; mesonCutArray[2] = "0152106500000000"; // trigger kPHI with y 0.8
    eventCutArray[ 3] = "00070113"; photonCutArray[ 3] = "00200009257002008250400000"; mesonCutArray[3] = "0152106500000000"; // trigger kHighMult with y 0.8
  } else if (trainConfig == 44) {
    eventCutArray[ 0] = "00080113"; photonCutArray[ 0] = "00200009257002008250400000"; mesonCutArray[0] = "0152106500000000"; // trigger kEMCEGA with y 0.8
    eventCutArray[ 1] = "00090113"; photonCutArray[ 1] = "00200009257002008250400000"; mesonCutArray[1] = "0152106500000000"; // trigger kEMCEJE with y 0.8
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009257002008250400000"; mesonCutArray[2] = "0152106500000000"; // minimum bias with y 0.8
    eventCutArray[ 3] = "00011113"; photonCutArray[ 3] = "00200009257002008250400000"; mesonCutArray[3] = "0152106500000000"; // trigger kINT8 with y 0.8
  } else if (trainConfig == 45) {
    eventCutArray[ 0] = "00042113"; photonCutArray[ 0] = "00200009257002008250400000"; mesonCutArray[0] = "0152106500000000"; // trigger kTRD CINT8 HEE
    eventCutArray[ 1] = "00044113"; photonCutArray[ 1] = "00200009257002008250400000"; mesonCutArray[1] = "0152106500000000"; // trigger kTRD CINT8 HSE
    eventCutArray[ 2] = "00046113"; photonCutArray[ 2] = "00200009257002008250400000"; mesonCutArray[2] = "0152106500000000"; // trigger kTRD CINT8 HJE
    eventCutArray[ 3] = "00048113"; photonCutArray[ 3] = "00200009257002008250400000"; mesonCutArray[3] = "0152106500000000"; // trigger kTRD CINT8 HQU
  } else if (trainConfig == 46) {
    eventCutArray[ 0] = "00041113"; photonCutArray[ 0] = "00200009257002008250400000"; mesonCutArray[0] = "0152106500000000"; // trigger kTRD CINT7 HEE
    eventCutArray[ 1] = "00043113"; photonCutArray[ 1] = "00200009257002008250400000"; mesonCutArray[1] = "0152106500000000"; // trigger kTRD CINT7 HSE
    eventCutArray[ 2] = "00045113"; photonCutArray[ 2] = "00200009257002008250400000"; mesonCutArray[2] = "0152106500000000"; // trigger kTRD CINT7 HJE
    eventCutArray[ 3] = "00047113"; photonCutArray[ 3] = "00200009257002008250400000"; mesonCutArray[3] = "0152106500000000"; // trigger kTRD CINT7 HQU
  } else if (trainConfig == 47) {
    eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009257002008250400000"; mesonCutArray[0] = "0152106500000000"; // trigger kEMC7
    eventCutArray[ 1] = "00053113"; photonCutArray[ 1] = "00200009257002008250400000"; mesonCutArray[1] = "0152106500000000"; // trigger kEMC8
    eventCutArray[ 2] = "00062113"; photonCutArray[ 2] = "00200009257002008250400000"; mesonCutArray[2] = "0152106500000000"; // trigger kPHI7
    eventCutArray[ 3] = "00063113"; photonCutArray[ 3] = "00200009257002008250400000"; mesonCutArray[3] = "0152106500000000"; // trigger kPHI8
  } else if (trainConfig == 48) {
    eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200009257002008250400000"; mesonCutArray[0] = "0152106500000000"; // trigger kEMC1
    eventCutArray[ 1] = "00071113"; photonCutArray[ 1] = "00200009257002008250400000"; mesonCutArray[1] = "0152106500000000"; // trigger kSHM1
    eventCutArray[ 2] = "00072113"; photonCutArray[ 2] = "00200009257002008250400000"; mesonCutArray[2] = "0152106500000000"; // trigger kSHM7
    eventCutArray[ 3] = "00073113"; photonCutArray[ 3] = "00200009257002008250400000"; mesonCutArray[3] = "0152106500000000"; // trigger kSHM8
  } else if (trainConfig == 49) {
    eventCutArray[ 0] = "00081113"; photonCutArray[ 0] = "00200009257002008250400000"; mesonCutArray[0] = "0152106500000000"; // trigger kEMCEGA + CINT7
    eventCutArray[ 1] = "00082113"; photonCutArray[ 1] = "00200009257002008250400000"; mesonCutArray[1] = "0152106500000000"; // trigger kEMCEGA + CINT8
    eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009257002008250400000"; mesonCutArray[2] = "0152106500000000"; // trigger kEMCEG1 + CINT7
    eventCutArray[ 3] = "00084113"; photonCutArray[ 3] = "00200009257002008250400000"; mesonCutArray[3] = "0152106500000000"; // trigger kEMCEG1 + CINT8
  } else if (trainConfig == 50) {
    eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200009257002008250400000"; mesonCutArray[0] = "0152106500000000"; // trigger kEMCEG2 + CINT7
    eventCutArray[ 1] = "00086113"; photonCutArray[ 1] = "00200009257002008250400000"; mesonCutArray[1] = "0152106500000000"; // trigger kEMCEG2 + CINT8
    eventCutArray[ 2] = "00091113"; photonCutArray[ 2] = "00200009257002008250400000"; mesonCutArray[2] = "0152106500000000"; // trigger kEMCEJE + CINT7
    eventCutArray[ 3] = "00092113"; photonCutArray[ 3] = "00200009257002008250400000"; mesonCutArray[3] = "0152106500000000"; // trigger kEMCEJE + CINT8
  } else if (trainConfig == 51) {
    eventCutArray[ 0] = "00093113"; photonCutArray[ 0] = "00200009257002008250400000"; mesonCutArray[0] = "0152106500000000"; // trigger kEMCEJ1 + CINT7
    eventCutArray[ 1] = "00094113"; photonCutArray[ 1] = "00200009257002008250400000"; mesonCutArray[1] = "0152106500000000"; // trigger kEMCEJ1 + CINT8
    eventCutArray[ 2] = "00095113"; photonCutArray[ 2] = "00200009257002008250400000"; mesonCutArray[2] = "0152106500000000"; // trigger kEMCEJ2 + CINT7
    eventCutArray[ 3] = "00096113"; photonCutArray[ 3] = "00200009257002008250400000"; mesonCutArray[3] = "0152106500000000"; // trigger kEMCEJ2 + CINT8		
  } else if (trainConfig == 52) { //pp 2.76TeV cuts
    eventCutArray[ 0] = "00042113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; // trigger kTRD CINT8 HEE
    eventCutArray[ 1] = "00044113"; photonCutArray[ 1] = "00200009366300003800000000"; mesonCutArray[1] = "0163103100900000"; // trigger kTRD CINT8 HSE
    eventCutArray[ 2] = "00046113"; photonCutArray[ 2] = "00200009366300003800000000"; mesonCutArray[2] = "0163103100900000"; // trigger kTRD CINT8 HJE
    eventCutArray[ 3] = "00048113"; photonCutArray[ 3] = "00200009366300003800000000"; mesonCutArray[3] = "0163103100900000"; // trigger kTRD CINT8 HQU
  } else if (trainConfig == 53) { //pp 2.76TeV cuts
    eventCutArray[ 0] = "00041113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; // trigger kTRD CINT7 HEE
    eventCutArray[ 1] = "00043113"; photonCutArray[ 1] = "00200009366300003800000000"; mesonCutArray[1] = "0163103100900000"; // trigger kTRD CINT7 HSE
    eventCutArray[ 2] = "00045113"; photonCutArray[ 2] = "00200009366300003800000000"; mesonCutArray[2] = "0163103100900000"; // trigger kTRD CINT7 HJE
    eventCutArray[ 3] = "00047113"; photonCutArray[ 3] = "00200009366300003800000000"; mesonCutArray[3] = "0163103100900000"; // trigger kTRD CINT7 HQU
  } else if (trainConfig == 54) { //pp 2.76TeV cuts
    eventCutArray[ 0] = "00052113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; // trigger kEMC7
    eventCutArray[ 1] = "00053113"; photonCutArray[ 1] = "00200009366300003800000000"; mesonCutArray[1] = "0163103100900000"; // trigger kEMC8
    eventCutArray[ 2] = "00062113"; photonCutArray[ 2] = "00200009366300003800000000"; mesonCutArray[2] = "0163103100900000"; // trigger kPHI7
    eventCutArray[ 3] = "00063113"; photonCutArray[ 3] = "00200009366300003800000000"; mesonCutArray[3] = "0163103100900000"; // trigger kPHI8
  } else if (trainConfig == 55) { //pp 2.76TeV cuts
    eventCutArray[ 0] = "00051113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; // trigger kEMC1
    eventCutArray[ 1] = "00071113"; photonCutArray[ 1] = "00200009366300003800000000"; mesonCutArray[1] = "0163103100900000"; // trigger kSHM1
    eventCutArray[ 2] = "00072113"; photonCutArray[ 2] = "00200009366300003800000000"; mesonCutArray[2] = "0163103100900000"; // trigger kSHM7
    eventCutArray[ 3] = "00073113"; photonCutArray[ 3] = "00200009366300003800000000"; mesonCutArray[3] = "0163103100900000"; // trigger kSHM8
  } else if (trainConfig == 56) { //pp 2.76TeV cuts
    eventCutArray[ 0] = "00081113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; // trigger kEMCEGA + CINT7
    eventCutArray[ 1] = "00082113"; photonCutArray[ 1] = "00200009366300003800000000"; mesonCutArray[1] = "0163103100900000"; // trigger kEMCEGA + CINT8
    eventCutArray[ 2] = "00083113"; photonCutArray[ 2] = "00200009366300003800000000"; mesonCutArray[2] = "0163103100900000"; // trigger kEMCEG1 + CINT7
    eventCutArray[ 3] = "00084113"; photonCutArray[ 3] = "00200009366300003800000000"; mesonCutArray[3] = "0163103100900000"; // trigger kEMCEG1 + CINT8
  } else if (trainConfig == 57) { //pp 2.76TeV cuts
    eventCutArray[ 0] = "00085113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; // trigger kEMCEG2 + CINT7
    eventCutArray[ 1] = "00086113"; photonCutArray[ 1] = "00200009366300003800000000"; mesonCutArray[1] = "0163103100900000"; // trigger kEMCEG2 + CINT8
    eventCutArray[ 2] = "00091113"; photonCutArray[ 2] = "00200009366300003800000000"; mesonCutArray[2] = "0163103100900000"; // trigger kEMCEJE + CINT7
    eventCutArray[ 3] = "00092113"; photonCutArray[ 3] = "00200009366300003800000000"; mesonCutArray[3] = "0163103100900000"; // trigger kEMCEJE + CINT8
  } else if (trainConfig == 58) { //pp 2.76TeV cuts
    eventCutArray[ 0] = "00093113"; photonCutArray[ 0] = "00200009366300003800000000"; mesonCutArray[0] = "0163103100900000"; // trigger kEMCEJ1 + CINT7
    eventCutArray[ 1] = "00094113"; photonCutArray[ 1] = "00200009366300003800000000"; mesonCutArray[1] = "0163103100900000"; // trigger kEMCEJ1 + CINT8
    eventCutArray[ 2] = "00095113"; photonCutArray[ 2] = "00200009366300003800000000"; mesonCutArray[2] = "0163103100900000"; // trigger kEMCEJ2 + CINT7
    eventCutArray[ 3] = "00096113"; photonCutArray[ 3] = "00200009366300003800000000"; mesonCutArray[3] = "0163103100900000"; // trigger kEMCEJ2 + CINT8		
  } else if (trainConfig == 59) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009257002008250400000"; mesonCutArray[0] = "0152103500000000"; //alpha meson 1
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009217302008250400000"; mesonCutArray[1] = "0152106500000000"; //pion 0-sigma cut for 0.4GeV<p<3.5GeV above -10-sigma
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009227302008250400000"; mesonCutArray[2] = "0152106500000000"; //pion 1-sigma cut for 0.4GeV<p<3.5GeV above -10-sigma
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009287302008250400000"; mesonCutArray[3] = "0152106500000000"; //pion 2-sigma cut for 0.4GeV<p<3.5GeV above   1-sigma
  } else if (trainConfig == 60) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152103500000000"; //New standard cut for eta analysis
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009227302008250400000"; mesonCutArray[1] = "0152107500000000"; //variation alpha 0.85
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200049227302008250400000"; mesonCutArray[2] = "0152103500000000"; //variation single pt 0.075
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200019227302008250400000"; mesonCutArray[3] = "0152103500000000"; //variation single pt 0.1
  } else if (trainConfig == 61) {//added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152103500000000"; //New standard cut for eta analysis
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009227302008250400000"; mesonCutArray[1] = "0152107500000000"; //variation alpha 0.85
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200049227302008250400000"; mesonCutArray[2] = "0152103500000000"; //variation single pt 0.075
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200019227302008250400000"; mesonCutArray[3] = "0152103500000000"; //variation single pt 0.1
  } else if (trainConfig == 62) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200006227302008250400000"; mesonCutArray[0] = "0152103500000000"; //variation 70% TPC Cluster
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200008227302008250400000"; mesonCutArray[1] = "0152103500000000"; //variation 35% TPC Cluster
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009327302008250400000"; mesonCutArray[2] = "0152103500000000"; //variation dE/dx e -4<sigma<5
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009627302008250400000"; mesonCutArray[3] = "0152103500000000"; //variation dE/dx e -2.5<sigma<4
  } else if (trainConfig == 63) {//added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200006227302008250400000"; mesonCutArray[0] = "0152103500000000"; //variation 70% TPC Cluster
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200008227302008250400000"; mesonCutArray[1] = "0152103500000000"; //variation 35% TPC Cluster
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009327302008250400000"; mesonCutArray[2] = "0152103500000000"; //variation dE/dx e -4<sigma<5
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009627302008250400000"; mesonCutArray[3] = "0152103500000000"; //variation dE/dx e -2.5<sigma<4
  } else if (trainConfig == 64) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009225302008250400000"; mesonCutArray[0] = "0152103500000000"; //variation min p 0.3GeV pion line
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009220302008250400000"; mesonCutArray[1] = "0152103500000000"; //variation min p 0.5GeV pion line
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009227102008250400000"; mesonCutArray[2] = "0152103500000000"; //variation max p 5GeV pion line
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009217302008250400000"; mesonCutArray[3] = "0152103500000000"; //variation nsigma>0 pion line
  } else if (trainConfig == 65) {//added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009225302008250400000"; mesonCutArray[0] = "0152103500000000"; //variation min p 0.3GeV pion line
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009220302008250400000"; mesonCutArray[1] = "0152103500000000"; //variation min p 0.5GeV pion line
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009227102008250400000"; mesonCutArray[2] = "0152103500000000"; //variation max p 5GeV pion line
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009217302008250400000"; mesonCutArray[3] = "0152103500000000"; //variation nsigma>0 pion line
  } else if (trainConfig == 66) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009257302008250400000"; mesonCutArray[0] = "0152103500000000"; //variation nsigma>2 pion line
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009227302002250400000"; mesonCutArray[1] = "0152103500000000"; //variation max qt 0.07 1D
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009227302009250400000"; mesonCutArray[2] = "0152103500000000"; //variation max qt 0.03 2D
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009227302003250400000"; mesonCutArray[3] = "0152103500000000"; //variation max qt 0.05 1D
  } else if (trainConfig == 67) {//added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009257302008250400000"; mesonCutArray[0] = "0152103500000000"; //variation nsigma>2 pion line
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009227302002250400000"; mesonCutArray[1] = "0152103500000000"; //variation max qt 0.07 1D
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009227302009250400000"; mesonCutArray[2] = "0152103500000000"; //variation max qt 0.03 2D
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009227302003250400000"; mesonCutArray[3] = "0152103500000000"; //variation max qt 0.05 1D
  } else if (trainConfig == 68) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008180400000"; mesonCutArray[0] = "0152103500000000"; //variation chi2 50 psi pair 0.2 2D
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009227302008210400000"; mesonCutArray[1] = "0152103500000000"; //variation chi2 30 psi pair 0.1 1D
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009227302008860400000"; mesonCutArray[2] = "0152103500000000"; //variation chi2 20 psi pair 0.05 2D
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009227302008250400000"; mesonCutArray[3] = "0252103500000000"; //variation BG scheme track mult
  } else if (trainConfig == 69) {//added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009227302008180400000"; mesonCutArray[0] = "0152103500000000"; //variation chi2 50 psi pair 0.2 2D
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009227302008210400000"; mesonCutArray[1] = "0152103500000000"; //variation chi2 30 psi pair 0.1 1D
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009227302008860400000"; mesonCutArray[2] = "0152103500000000"; //variation chi2 20 psi pair 0.05 2D
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009227302008250400000"; mesonCutArray[3] = "0252103500000000"; //variation BG scheme track mult
  } else if (trainConfig == 70) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152103500000000"; //New standard cut for eta analysis
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009227302008250400000"; mesonCutArray[1] = "0152103500000000"; //New standard cut for eta analysis
  } else if (trainConfig == 71) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152103500000000"; //New standard cut for eta analysis
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009227302008250400000"; mesonCutArray[1] = "0152101500000000"; //variation alpha pT dependent
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009227302008250400000"; mesonCutArray[2] = "0152109500000000"; //variation alpha
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009227302008250400000"; mesonCutArray[3] = "0152101500000002"; //variation alpha/opan Max
  } else if (trainConfig == 72) {//added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152103500000000"; //New standard cut for eta analysis
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009227302008250400000"; mesonCutArray[1] = "0152101500000000"; //variation alpha 0.85
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009227302008250400000"; mesonCutArray[2] = "0152109500000000"; //variation alpha
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009227302008250400000"; mesonCutArray[3] = "0152101500000002"; //variation alpha opan max				
  } else if (trainConfig == 73) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152103500000000"; //New standard cut for eta analysis
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009227302008250404000"; mesonCutArray[1] = "0152101500000000"; //variation alpha pT dependent
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009227302008250404000"; mesonCutArray[2] = "0152109500000000"; //variation alpha
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009227302008250400000"; mesonCutArray[3] = "0152101500000000"; //double counting
  } else if (trainConfig == 74) {//added signals
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152103500000000"; //New standard cut for eta analysis
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009227302008250404000"; mesonCutArray[1] = "0152101500000000"; //variation alpha 0.85
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009227302008250404000"; mesonCutArray[2] = "0152109500000000"; //variation alpha
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009227302008250400000"; mesonCutArray[3] = "0152101500000000"; //double counting	
    } else if (trainConfig == 75) { //pp 8TeV cuts
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152103500000000"; //New standard cut for eta/pi0 analysis
    eventCutArray[ 1] = "00005211"; photonCutArray[ 1] = "00200009227302008250400000"; mesonCutArray[1] = "0152103500000000"; // trigger kEMC7
    eventCutArray[ 2] = "00006211"; photonCutArray[ 2] = "00200009227302008250400000"; mesonCutArray[2] = "0152103500000000"; // trigger kPHI7
    eventCutArray[ 3] = "00008111"; photonCutArray[ 3] = "00200009227302008250400000"; mesonCutArray[3] = "0152103500000000"; // trigger kEMCEGA + CINT7
    //---------systematic studies July 2015--------------------------//
  } else if (trainConfig == 80) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008250404000"; mesonCutArray[0] = "0152101500000000"; //New standard cut for eta: alpha pT dependent
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009227302008250404000"; mesonCutArray[1] = "0152106500000000"; // alpha variation 0.8		
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009227302008250404000"; mesonCutArray[2] = "0152103500000000"; // alpha variation  1
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200049227302008250404000"; mesonCutArray[3] = "0152101500000000"; // single pT cut 0.075
  } else if (trainConfig == 81) {
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009227302008250404000"; mesonCutArray[0] = "0152101500000000"; //New standard cut for eta: alpha pT dependent
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009227302008250404000"; mesonCutArray[1] = "0152106500000000"; // alpha variation  0.8
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009227302008250404000"; mesonCutArray[2] = "0152103500000000"; // alpha variation  1
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200049227302008250404000"; mesonCutArray[3] = "0152101500000000"; // single pT cut 0.075
  } else if (trainConfig == 82) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200019227302008250404000"; mesonCutArray[0] = "0152101500000000"; // single pT cut 0.1
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200008227302008250404000"; mesonCutArray[1] = "0152101500000000"; // TPC cls 0.35	
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200006227302008250404000"; mesonCutArray[2] = "0152101500000000"; // TPC cls 0.7
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009227302009250404000"; mesonCutArray[3] = "0152101500000000"; // qT cut 0.03 2D
  } else if (trainConfig == 83) {
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200019227302008250404000"; mesonCutArray[0] = "0152101500000000"; // single pT cut 0.1
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200008227302008250404000"; mesonCutArray[1] = "0152101500000000"; // TPC cls 0.35
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200006227302008250404000"; mesonCutArray[2] = "0152101500000000"; // TPC cls 0.7
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009227302009250404000"; mesonCutArray[3] = "0152101500000000"; // qT cut 0.03 2D
  } else if (trainConfig == 84) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008210404000"; mesonCutArray[0] = "0152101500000000"; // variation chi2 30 psi pair 0.1 1D
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009227302008180404000"; mesonCutArray[1] = "0152101500000000"; // variation chi2 50 psi pair 0.2 2D	
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009227302008860404000"; mesonCutArray[2] = "0152101500000000"; // variation chi2 20 psi pair 0.05 2D
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009227302009250404000"; mesonCutArray[3] = "0252101500000000"; // variation BG scheme track mult
  } else if (trainConfig == 85) {
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009227302008210404000"; mesonCutArray[0] = "0152101500000000"; // variation chi2 30 psi pair 0.1 1D
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200008227302008180404000"; mesonCutArray[1] = "0152101500000000"; // variation chi2 50 psi pair 0.2 2D
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200006227302008860404000"; mesonCutArray[2] = "0152101500000000"; // variation chi2 20 psi pair 0.05 2D
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009227302009250404000"; mesonCutArray[3] = "0252101500000000"; // variation BG scheme track mult
  } else if (trainConfig == 86) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008250404000"; mesonCutArray[0] = "0152101500000000"; //New standard cut for eta: alpha pT dependent
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009227302002250404000"; mesonCutArray[1] = "0152101500000000"; // qT		
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009227302009250404000"; mesonCutArray[2] = "0152101500000000"; // qT
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009227302008250004000"; mesonCutArray[3] = "0152101500000000"; // cosPA
  } else if (trainConfig == 87) {
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009227302008250404000"; mesonCutArray[0] = "0152101500000000"; //New standard cut for eta: alpha pT dependent
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009227302002250404000"; mesonCutArray[1] = "0152101500000000"; // qT
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009227302009250404000"; mesonCutArray[2] = "0152101500000000"; // qT
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009227302008250004000"; mesonCutArray[3] = "0152101500000000"; // cosPA
    //---------systematic studies Agost 2015 with dEdx--------------------------//
  } else if (trainConfig == 90) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009217302008250404000"; mesonCutArray[0] = "0152101500000000"; //New standard cut for eta: alpha pT dependent
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009217302008250404000"; mesonCutArray[1] = "0152106500000000"; // alpha variation 0.8		
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009217302008250404000"; mesonCutArray[2] = "0152103500000000"; // alpha variation  1
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200049217302008250404000"; mesonCutArray[3] = "0152101500000000"; // single pT cut 0.075
  } else if (trainConfig == 91) {
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009217302008250404000"; mesonCutArray[0] = "0152101500000000"; //New standard cut for eta: alpha pT dependent
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009217302008250404000"; mesonCutArray[1] = "0152106500000000"; // alpha variation  0.8
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009217302008250404000"; mesonCutArray[2] = "0152103500000000"; // alpha variation  1
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200049217302008250404000"; mesonCutArray[3] = "0152101500000000"; // single pT cut 0.075
  } else if (trainConfig == 92) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200019217302008250404000"; mesonCutArray[0] = "0152101500000000"; // single pT cut 0.1
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200008217302008250404000"; mesonCutArray[1] = "0152101500000000"; // TPC cls 0.35	
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200006217302008250404000"; mesonCutArray[2] = "0152101500000000"; // TPC cls 0.7
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009217302009250404000"; mesonCutArray[3] = "0152101500000000"; // qT cut 0.03 2D
  } else if (trainConfig == 93) {
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200019217302008250404000"; mesonCutArray[0] = "0152101500000000"; // single pT cut 0.1
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200008217302008250404000"; mesonCutArray[1] = "0152101500000000"; // TPC cls 0.35
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200006217302008250404000"; mesonCutArray[2] = "0152101500000000"; // TPC cls 0.7
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009217302009250404000"; mesonCutArray[3] = "0152101500000000"; // qT cut 0.03 2D
  } else if (trainConfig == 94) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009217302008210404000"; mesonCutArray[0] = "0152101500000000"; // variation chi2 30 psi pair 0.1 1D
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009217302008180404000"; mesonCutArray[1] = "0152101500000000"; // variation chi2 50 psi pair 0.2 2D	
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009217302008860404000"; mesonCutArray[2] = "0152101500000000"; // variation chi2 20 psi pair 0.05 2D
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009217302009250404000"; mesonCutArray[3] = "0252101500000000"; // variation BG scheme track mult
  } else if (trainConfig == 95) {
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009217302008210404000"; mesonCutArray[0] = "0152101500000000"; // variation chi2 30 psi pair 0.1 1D
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009217302008180404000"; mesonCutArray[1] = "0152101500000000"; // variation chi2 50 psi pair 0.2 2D
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009217302008860404000"; mesonCutArray[2] = "0152101500000000"; // variation chi2 20 psi pair 0.05 2D
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009217302009250404000"; mesonCutArray[3] = "0252101500000000"; // variation BG scheme track mult
  } else if (trainConfig == 96) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009227302008250404000"; mesonCutArray[0] = "0152101500000000"; //New standard cut for eta: alpha pT dependent
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009217302002250404000"; mesonCutArray[1] = "0152101500000000"; // qT		
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009217302009250404000"; mesonCutArray[2] = "0152101500000000"; // qT
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009217302008250004000"; mesonCutArray[3] = "0152101500000000"; // cosPA
  } else if (trainConfig == 97) {
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009227302008250404000"; mesonCutArray[0] = "0152101500000000"; //New standard cut for eta: alpha pT dependent
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009217302002250404000"; mesonCutArray[1] = "0152101500000000"; // qT
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009217302009250404000"; mesonCutArray[2] = "0152101500000000"; // qT
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009217302008250004000"; mesonCutArray[3] = "0152101500000000"; // cosPA
  } else if (trainConfig == 98) {
    eventCutArray[ 0] = "00000113"; photonCutArray[ 0] = "00200009317302008250404000"; mesonCutArray[0] = "0152101500000000"; //dEdx variation
    eventCutArray[ 1] = "00000113"; photonCutArray[ 1] = "00200009617302008250404000"; mesonCutArray[1] = "0152101500000000"; //dEdx variation
    eventCutArray[ 2] = "00000113"; photonCutArray[ 2] = "00200009215302008250404000"; mesonCutArray[2] = "0152101500000000"; //dEdx variation
    eventCutArray[ 3] = "00000113"; photonCutArray[ 3] = "00200009210302008250404000"; mesonCutArray[3] = "0152101500000000"; //dEdx variation
  } else if (trainConfig == 99) {
    eventCutArray[ 0] = "00000123"; photonCutArray[ 0] = "00200009317302008250404000"; mesonCutArray[0] = "0152101500000000"; //dEdx variation
    eventCutArray[ 1] = "00000123"; photonCutArray[ 1] = "00200009617302008250404000"; mesonCutArray[1] = "0152101500000000"; //dEdx variation
    eventCutArray[ 2] = "00000123"; photonCutArray[ 2] = "00200009215302008250404000"; mesonCutArray[2] = "0152101500000000"; //dEdx variation
    eventCutArray[ 3] = "00000123"; photonCutArray[ 3] = "00200009210302008250404000"; mesonCutArray[3] = "0152101500000000"; //dEdx variation
  } else if (trainConfig == 100) { // MB
    eventCutArray[ 0] = "00100113"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152103500000000"; // 0 -2
    eventCutArray[ 1] = "01200113"; photonCutArray[ 1] = "00200009227302008250400000"; mesonCutArray[1] = "0152103500000000"; // 2 -5
    eventCutArray[ 2] = "02300113"; photonCutArray[ 2] = "00200009227302008250400000"; mesonCutArray[2] = "0152103500000000"; // 5 -10
    eventCutArray[ 3] = "03400113"; photonCutArray[ 3] = "00200009227302008250400000"; mesonCutArray[3] = "0152103500000000"; // 10-30
    eventCutArray[ 4] = "04500113"; photonCutArray[ 4] = "00200009227302008250400000"; mesonCutArray[4] = "0152103500000000"; // 30-100
  } else if (trainConfig == 101) {  // like trainConfig 71, just with special trigger kINT7
    eventCutArray[ 0] = "00010113"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152103500000000"; //New standard cut for eta analysis
    eventCutArray[ 1] = "00010113"; photonCutArray[ 1] = "00200009227302008250400000"; mesonCutArray[1] = "0152101500000000"; //variation alpha pT dependent
    eventCutArray[ 2] = "00010113"; photonCutArray[ 2] = "00200009227302008250400000"; mesonCutArray[2] = "0152109500000000"; //variation alpha
    eventCutArray[ 3] = "00010113"; photonCutArray[ 3] = "00200009227302008250400000"; mesonCutArray[3] = "0152101500000002"; ///variation alpha opan max 
  } else if (trainConfig == 102) {  // like trainConfig 71, except with smearing added to meson cut
    eventCutArray[ 0] = "00010113"; photonCutArray[ 0] = "00200009227302008250400000"; mesonCutArray[0] = "0152103500900000"; //New standard cut for eta analysis
    eventCutArray[ 1] = "00010113"; photonCutArray[ 1] = "00200009227302008250400000"; mesonCutArray[1] = "0152101500900000"; //variation alpha pT dependent
    eventCutArray[ 2] = "00010113"; photonCutArray[ 2] = "00200009227302008250400000"; mesonCutArray[2] = "0152109500900000"; //variation alpha
    eventCutArray[ 3] = "00010113"; photonCutArray[ 3] = "00200009227302008250400000"; mesonCutArray[3] = "0152101500900002"; ///variation alpha opan max    
  }	else {
    Error(Form("GammaConvV1_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  TList *EventCutList = new TList();
  TList *ConvCutList = new TList();
  TList *MesonCutList = new TList();

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
  AliConvEventCuts **analysisEventCuts = new AliConvEventCuts*[numberOfCuts];
  ConvCutList->SetOwner(kTRUE);
  AliConversionPhotonCuts **analysisCuts = new AliConversionPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();
    TString fitNamePi0 = Form("Pi0_Fit_Data_%s",energy.Data());
    TString fitNameEta = Form("Eta_Fit_Data_%s",energy.Data());
    TString fAddedSignalString = eventCutArray[i];
    fAddedSignalString = fAddedSignalString(6,1);
    Bool_t fAddedSignal = kFALSE;
    if (fAddedSignalString.CompareTo("2") == 0) fAddedSignal = kTRUE;

    TString mcInputNamePi0 = "";
    TString mcInputNameEta = "";
    if (fAddedSignal && (periodname.Contains("LHC12i3") || periodname.CompareTo("LHC14e2b")==0)){
      mcInputNamePi0 = Form("Pi0_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta = Form("Eta_%s%s_addSig_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    } else {
      mcInputNamePi0 = Form("Pi0_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
      mcInputNameEta = Form("Eta_%s%s_%s", mcName.Data(), mcNameAdd.Data(), energy.Data() );
    }	
    //		if (doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kFALSE, kFALSE, kFALSE, fileNameInputForWeighting, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);
    if (doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);

    analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisCuts[i] = new AliConversionPhotonCuts();
    analysisCuts[i]->InitializeCutsFromCutString(photonCutArray[i].Data());
    
    if( trainConfig  == 73 || trainConfig  == 74 || (trainConfig  >= 80 && trainConfig  <= 87) ){
      
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
  task->SetDoPhotonQA(enableQAPhotonTask);  //Attention new switch small for Photon QA

  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaConvV1_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaConvV1_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
