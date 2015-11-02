void AddTask_GammaCalo_pp(  Int_t     trainConfig               = 1,                  // change different set of cuts
                            Int_t     isMC                      = 0,                // run MC
                            Int_t     enableQAMesonTask         = 0,                 // enable QA in AliAnalysisTaskGammaCalo
                            Int_t     enableQAClusterTask       = 0,                 // enable additional QA task
                            TString   fileNameInputForWeighting = "MCSpectraInput.root",       // path to file for weigting input
                            TString   cutnumberAODBranch        = "000000006008400001001500000",
                            TString   periodname                = "LHC12f1x",             // period name
                            Bool_t    doWeighting               = kFALSE,              // enables weighting
                            Bool_t    isUsingTHnSparse          = kTRUE,              // enable or disable usage of THnSparses for background estimation
                            Int_t     enableExtMatchAndQA       = 0,                // enable QA(3), extMatch+QA(2), extMatch(1), disabled (0)
                            Bool_t    enableTriggerMimicking    = kFALSE,              // enable trigger mimicking
                            Bool_t    enableTriggerOverlapRej   = kFALSE,              // enable trigger overlap rejection
                            Float_t   maxFacPtHard              = 3.,                // maximum factor between hardest jet and ptHard generated
                            TString   periodNameV0Reader        = ""
) {
  
  Int_t isHeavyIon = 0;
  
  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error(Form("AddTask_GammaCalo_pp_%i",trainConfig), "No analysis manager found.");
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
  
  Printf("here \n");
  
  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton = "00000008400100001500000000";
  TString cutnumberEvent = "00000003";
  Bool_t doEtaShift = kFALSE;
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  if( !(AliV0ReaderV1*)mgr->GetTask("V0ReaderV1") ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");
    if (periodNameV0Reader.CompareTo("") != 0) fV0ReaderV1->SetPeriodName(periodNameV0Reader);
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
    }
    AliConvEventCuts *fEventCuts=NULL;
    if(cutnumberEvent!=""){
      fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
      fEventCuts->SetPreSelectionCutFlag(kTRUE);
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

    AliLog::SetGlobalLogLevel(AliLog::kFatal);

    //connect input V0Reader
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);
  }

  //================================================
  //========= Add task to the ANALYSIS manager =====
  //================================================
  AliAnalysisTaskGammaCalo *task=NULL;
  task= new AliAnalysisTaskGammaCalo(Form("GammaCalo_%i",trainConfig));
  task->SetIsHeavyIon(isHeavyIon);
  task->SetIsMC(isMC);
  // Cut Numbers to use in Analysis
  Int_t numberOfCuts = 2;
  if (trainConfig == 32   || trainConfig == 101 || trainConfig == 201 || trainConfig == 63)
      numberOfCuts = 1;
  if (trainConfig == 111  || trainConfig == 114 || trainConfig == 117 || trainConfig == 120   || trainConfig == 121   || 
      trainConfig == 122  || trainConfig == 123 || trainConfig == 124 || trainConfig == 109)
      numberOfCuts = 3;
  if (trainConfig == 31   || trainConfig == 4   || trainConfig == 5   || trainConfig == 6   || trainConfig == 7   ||
      trainConfig == 10   || trainConfig == 51  || trainConfig == 52  || trainConfig == 53  || trainConfig == 54  ||
      trainConfig == 55   || trainConfig == 56  || trainConfig == 57  || trainConfig == 58  || 
      trainConfig == 60   || trainConfig == 61  || trainConfig == 62  || trainConfig == 64  || trainConfig == 70  ||
      trainConfig == 71   || trainConfig == 72  || trainConfig == 73  || trainConfig == 74  || trainConfig == 75  ||
      trainConfig == 76   || trainConfig == 77  || trainConfig == 78  || trainConfig == 79  || trainConfig == 80  ||
      trainConfig == 81   || trainConfig == 102 || trainConfig == 110 || trainConfig == 8)
      numberOfCuts = 4;
  if (trainConfig == 2    || trainConfig == 3   || trainConfig == 84  || trainConfig == 85  || trainConfig == 86  ||
      trainConfig == 87   || trainConfig == 88  || trainConfig == 89  || trainConfig == 98  || trainConfig == 99  ||
      trainConfig == 103  || trainConfig == 104  )
      numberOfCuts = 5;
  if (trainConfig == 65   || trainConfig == 66  || trainConfig == 67  || trainConfig == 68  || trainConfig == 82  ||
      trainConfig == 83   || trainConfig == 105)
      numberOfCuts = 6;
  if (trainConfig == 59)
      numberOfCuts = 8;

  
  TString *eventCutArray = new TString[numberOfCuts];
  TString *clusterCutArray = new TString[numberOfCuts];
  TString *mesonCutArray = new TString[numberOfCuts];

  // cluster cuts
  // 0 "ClusterType",  1 "EtaMin", 2 "EtaMax", 3 "PhiMin", 4 "PhiMax", 5 "DistanceToBadChannel", 6 "Timing", 7 "TrackMatching", 8 "ExoticCell",
  // 9 "MinEnergy", 10 "MinNCells", 11 "MinM02", 12 "MaxM02", 13 "MinM20", 14 "MaxM20", 15 "MaximumDispersion", 16 "NLM"
  
  // ************************************* EMCAL cuts ****************************************************
  // LHC11a
  if (trainConfig == 1){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000050"; // 400 MeV cluster min energy
    eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111121050032220000"; mesonCutArray[1] = "0163103100000050"; // 400 MeV cluster min energy
  } else if (trainConfig == 2){ //EMCAL minEnergy variation
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121050012220000"; mesonCutArray[0] = "0163103100000050"; //0.2 GeV/c
    eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111121050022220000"; mesonCutArray[1] = "0163103100000050"; //0.3 GeV/c
    eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111121050032220000"; mesonCutArray[2] = "0163103100000050"; //0.4 GeV/c default
    eventCutArray[ 3] = "00003113"; clusterCutArray[3] = "1111121050042220000"; mesonCutArray[3] = "0163103100000050"; //0.5 GeV/c
    eventCutArray[ 4] = "00003113"; clusterCutArray[4] = "1111121050052220000"; mesonCutArray[4] = "0163103100000050"; //0.6 GeV/c
  } else if (trainConfig == 3){ //EMCAL minNCells variation
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121050031220000"; mesonCutArray[0] = "0163103100000050"; //n cells >= 1
    eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111121050033220000"; mesonCutArray[1] = "0163103100000050"; //n cells >= 3
    eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111121050032000000"; mesonCutArray[2] = "0163103100000050"; //no M02 cut
    eventCutArray[ 3] = "00003113"; clusterCutArray[3] = "1113121050032220000"; mesonCutArray[3] = "0163103100000050"; //only modules with TRD infront
    eventCutArray[ 4] = "00003113"; clusterCutArray[4] = "1111221050032220000"; mesonCutArray[4] = "0163103100000050"; //no modules with TRD infront
  } else if (trainConfig == 4 ){ // EMCAL clusters 2.76 TeV NonLinearity
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111101050032220000"; mesonCutArray[0] = "0163103100000050"; // NonLinearity kSDMv5
    eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111121050032220000"; mesonCutArray[1] = "0163103100000050"; // NonLinearity LHC11a ConvCalo
    eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111122050032220000"; mesonCutArray[2] = "0163103100000050"; // NonLinearity LHC11a Calo
    eventCutArray[ 3] = "00003113"; clusterCutArray[3] = "1111100050032220000"; mesonCutArray[3] = "0163103100000050"; // NonLinearity none
  } else if (trainConfig == 5 ){ // EMCAL clusters open angle variation
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000050"; // min open angle - 0.0202 
    eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111121050032220000"; mesonCutArray[1] = "0163103100000030"; // min open angle - 0.01 
    eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111121050032220000"; mesonCutArray[2] = "0163103100000040"; // min open angle - 0.0152 
    eventCutArray[ 3] = "00003113"; clusterCutArray[3] = "1111121050032220000"; mesonCutArray[3] = "0163103100000060"; // min open angle - 0.0404 
  } else if (trainConfig == 6){ // EMCAL clusters 2.76 TeV NonLinearity
    eventCutArray[ 0] = "00051113"; clusterCutArray[0] = "1111101050032220000"; mesonCutArray[0] = "0163103100000050"; // NonLinearity kSDMv5
    eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111121050032220000"; mesonCutArray[1] = "0163103100000050"; // NonLinearity LHC11a ConvCalo
    eventCutArray[ 2] = "00051113"; clusterCutArray[2] = "1111122050032220000"; mesonCutArray[2] = "0163103100000050"; // NonLinearity LHC11a Calo
    eventCutArray[ 3] = "00051113"; clusterCutArray[3] = "1111100050032220000"; mesonCutArray[3] = "0163103100000050"; // NonLinearity none
  } else if (trainConfig == 7 ){ // EMCAL clusters open angle variation
    eventCutArray[ 0] = "00051113"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000050"; // min open angle - 0.0202 
    eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111121050032220000"; mesonCutArray[1] = "0163103100000030"; // min open angle - 0.01 
    eventCutArray[ 2] = "00051113"; clusterCutArray[2] = "1111121050032220000"; mesonCutArray[2] = "0163103100000040"; // min open angle - 0.0152 
    eventCutArray[ 3] = "00051113"; clusterCutArray[3] = "1111121050032220000"; mesonCutArray[3] = "0163103100000060"; // min open angle - 0.0404 
  } else if (trainConfig == 8 ){ // EMCAL clusters 2.76 TeV additional NonLinearity variations
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111123050032220000"; mesonCutArray[0] = "0163103100000050"; // NonLinearity kTestBeamv2 + ConvCalo
    eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111124050032220000"; mesonCutArray[1] = "0163103100000050"; // NonLinearity kTestBeamv2 + Calo
    eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111125050032220000"; mesonCutArray[2] = "0163103100000050"; // NonLinearity LHC11a kSDM ConvCalo
    eventCutArray[ 3] = "00003113"; clusterCutArray[3] = "1111126050032220000"; mesonCutArray[3] = "0163103100000050"; // NonLinearity LHC11a kSDM Calo
    // LHC13g  
  } else if (trainConfig == 10){  // EMCAL clusters, EMCEGA triggers
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000050"; // INT7
    eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111121060032220000"; mesonCutArray[1] = "0163103100000050"; // EMC7
    eventCutArray[ 2] = "00083113"; clusterCutArray[2] = "1111121060032220000"; mesonCutArray[2] = "0163103100000050"; // EMCEG1,
    eventCutArray[ 3] = "00085113"; clusterCutArray[3] = "1111121060032220000"; mesonCutArray[3] = "0163103100000050"; // EMCEG2,
  
  } else if (trainConfig == 12){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0) without and with added signals
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111100050032220000"; mesonCutArray[0] = "0163103100000050"; // 400 MeV cluster min energy
    eventCutArray[ 1] = "00003123"; clusterCutArray[1] = "1111100050032220000"; mesonCutArray[1] = "0163103100000050"; // 400 MeV cluster min energy
  } else if (trainConfig == 13){ // EMCAL clusters 2.76 TeV LHC11a, with SDD (0), kEMC1 (1)
    eventCutArray[ 0] = "00003123"; clusterCutArray[0] = "1111100050032230000"; mesonCutArray[0] = "0163103100000050"; // 400 MeV cluster min energy
    eventCutArray[ 1] = "00051123"; clusterCutArray[1] = "1111100050032230000"; mesonCutArray[1] = "0163103100000050"; // 400 MeV cluster min energy
  // ************************************* PHOS cuts ****************************************************
  } else if (trainConfig == 31) { //PHOS clusters
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "2444400040033200000"; mesonCutArray[0] = "0163103100000050"; //pp LHC11a with SDD, PHOS
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "2444400040033200000"; mesonCutArray[1] = "0163103100000050"; //pp LHC13g default MB
    eventCutArray[ 2] = "00061113"; clusterCutArray[2] = "2444400040033200000"; mesonCutArray[2] = "0163103100000050"; //pp LHC11a PHI1
    eventCutArray[ 3] = "00062113"; clusterCutArray[3] = "2444400040033200000"; mesonCutArray[3] = "0163103100000050"; //pp LHC11a PHI7
  } else if (trainConfig == 32){ // Validation PHOS
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "2444400040033200000"; mesonCutArray[0] = "0163003100900050";
  } else if (trainConfig == 33){ // PHOS clusters, without and with added signals
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "2444400040033200000"; mesonCutArray[0] = "0163003100900050";
    eventCutArray[ 1] = "00003123"; clusterCutArray[1] = "2444400040033200000"; mesonCutArray[1] = "0163003100900050";

    
    // here is the order of the cluster cut string
    // usually for EMCal we start with 11111: default values for            "ClusterType", "EtaMin", "EtaMax", "PhiMin", "PhiMax"
    // then two numbers for nonlinearity, e.g. 21: this is                  "NonLinearity1", "NonLinearity2"
    // Then some cuts on the clusters, e.g. 06003222: this is               "DistanceToBadChannel", "Timing", "TrackMatching", "ExoticCell", "MinEnergy", "MinNCells", "MinM02", "MaxM02"
    // finally some for now unused cuts, usually 0000: this is              "MinM20", "MaxM20", "MaximumDispersion", "NLM"
    
    // and the meson cut string
    // it starts with general criteria, 01631031: this is       "MesonKind", "BackgroundScheme", "NumberOfBgEvents", "DegreesForRotationMethod", "RapidityMesonCut", "RCut", "AlphaMesonCut", "SelectionWindow"
    // followed usually by 000000 (cuts for PCM!): this is      "SharedElectronCuts", "RejectToCloseV0s", "UseMCPSmearing", "DcaGammaGamma", "DcaRPrimVtx", "DcaZPrimVtx"
    // concluded by e.g 50 for the angle: this is               "MinOpanMesonCut", "MaxOpanMesonCut"
    
    // LHC13g cut studies
  } else if (trainConfig == 51){  // EMCAL clusters, EMCEG1 trigger
    eventCutArray[ 0] = "00083113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000050"; // EMCEG1, 400 MeV min energy, NCells >=2, M02 default cut, 35ns timing, 1 cell diagonal
    eventCutArray[ 1] = "00083113"; clusterCutArray[1] = "1111121060052220000"; mesonCutArray[1] = "0163103100000050"; // EMCEG1, 600 MeV min energy
    eventCutArray[ 2] = "00083113"; clusterCutArray[2] = "1111121060031220000"; mesonCutArray[2] = "0163103100000050"; // EMCEG1,                     NCells >=1
    eventCutArray[ 3] = "00083113"; clusterCutArray[3] = "1111121060033220000"; mesonCutArray[3] = "0163103100000050"; // EMCEG1,                     NCells >=3
  } else if (trainConfig == 52){  // EMCAL clusters, EMCEG1 trigger
    eventCutArray[ 0] = "00083113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000030"; // EMCEG1,                                                               0.01 opening
    eventCutArray[ 1] = "00083113"; clusterCutArray[1] = "1111121060032220000"; mesonCutArray[1] = "0163103100000040"; // EMCEG1,                                                               0.75 cell diagonal
    eventCutArray[ 2] = "00083113"; clusterCutArray[2] = "1111121060032220000"; mesonCutArray[2] = "0163103100000060"; // EMCEG1,                                                               2 cell diagonals
    eventCutArray[ 3] = "00083113"; clusterCutArray[3] = "1111121060032000000"; mesonCutArray[3] = "0163103100000050"; // EMCEG1,                                 no M02 cut
  } else if (trainConfig == 53){  // EMCAL clusters, EMCEG1 trigger
    eventCutArray[ 0] = "00083113"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000050"; // EMCEG1,                                                   50ns timing
    eventCutArray[ 1] = "00083113"; clusterCutArray[1] = "1111101060032220000"; mesonCutArray[1] = "0163103100000050"; // EMCEG1,                standard kSDMv5
    eventCutArray[ 2] = "00083113"; clusterCutArray[2] = "1111122060032220000"; mesonCutArray[2] = "0163103100000050"; // EMCEG1,              NonLinearity LHC11a Calo
    eventCutArray[ 3] = "00083113"; clusterCutArray[3] = "1111100060032220000"; mesonCutArray[3] = "0163103100000050"; // EMCEG1,              NonLinearity none
  
  } else if (trainConfig == 54){  // EMCAL clusters, EMCEG2 trigger
    eventCutArray[ 0] = "00085113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000050"; // EMCEG2, 400 MeV min energy, NCells >=2, M02 default cut, 35ns timing, 1 cell diagonal
    eventCutArray[ 1] = "00085113"; clusterCutArray[1] = "1111121060052220000"; mesonCutArray[1] = "0163103100000050"; // EMCEG2, 600 MeV min energy
    eventCutArray[ 2] = "00085113"; clusterCutArray[2] = "1111121060031220000"; mesonCutArray[2] = "0163103100000050"; // EMCEG2,                     NCells >=1
    eventCutArray[ 3] = "00085113"; clusterCutArray[3] = "1111121060033220000"; mesonCutArray[3] = "0163103100000050"; // EMCEG2,                     NCells >=3
  } else if (trainConfig == 55){  // EMCAL clusters, EMCEG2 trigger
    eventCutArray[ 0] = "00085113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000030"; // EMCEG2,                                                               0.01 opening
    eventCutArray[ 1] = "00085113"; clusterCutArray[1] = "1111121060032220000"; mesonCutArray[1] = "0163103100000040"; // EMCEG2,                                                               0.75 cell diagonal
    eventCutArray[ 2] = "00085113"; clusterCutArray[2] = "1111121060032220000"; mesonCutArray[2] = "0163103100000060"; // EMCEG2,                                                               2 cell diagonals
    eventCutArray[ 3] = "00085113"; clusterCutArray[3] = "1111121060032000000"; mesonCutArray[3] = "0163103100000050"; // EMCEG2,                                 no M02 cut
  } else if (trainConfig == 56){  // EMCAL clusters, EMCEG2 trigger
    eventCutArray[ 0] = "00085113"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000050"; // EMCEG2,                                                   50ns timing
    eventCutArray[ 1] = "00085113"; clusterCutArray[1] = "1111101060032220000"; mesonCutArray[1] = "0163103100000050"; // EMCEG2,                standard kSDMv5
    eventCutArray[ 2] = "00085113"; clusterCutArray[2] = "1111122060032220000"; mesonCutArray[2] = "0163103100000050"; // EMCEG2,              NonLinearity LHC11a Calo
    eventCutArray[ 3] = "00085113"; clusterCutArray[3] = "1111100060032220000"; mesonCutArray[3] = "0163103100000050"; // EMCEG2,              NonLinearity none
    
  } else if (trainConfig == 57){  // EMCAL clusters, INT7 trigger
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000050"; // INT7, 400 MeV min energy, NCells >=2, M02 default cut, 35ns timing, 1 cell diagonal
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111121060052220000"; mesonCutArray[1] = "0163103100000050"; // INT7, 600 MeV min energy
    eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111121060031220000"; mesonCutArray[2] = "0163103100000050"; // INT7,                     NCells >=1
    eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1111121060033220000"; mesonCutArray[3] = "0163103100000050"; // INT7,                     NCells >=3
  } else if (trainConfig == 58){  // EMCAL clusters, INT7 trigger
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000030"; // INT7,                                                               0.01 opening
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111121060032220000"; mesonCutArray[1] = "0163103100000040"; // INT7,                                                               0.75 cell diagonal
    eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111121060032220000"; mesonCutArray[2] = "0163103100000060"; // INT7,                                                               2 cell diagonals
    eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1111121060032000000"; mesonCutArray[3] = "0163103100000050"; // INT7,                                 no M02 cut
  } else if (trainConfig == 59){  // EMCAL clusters, INT7 trigger
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000050"; // INT7, 50ns timing
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111101060032220000"; mesonCutArray[1] = "0163103100000050"; // INT7, standard kSDMv5
    eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111122060032220000"; mesonCutArray[2] = "0163103100000050"; // INT7, NonLinearity LHC11a Calo
    eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1111100060032220000"; mesonCutArray[3] = "0163103100000050"; // INT7, NonLinearity none
    eventCutArray[ 4] = "00000113"; clusterCutArray[4] = "1111123060032220000"; mesonCutArray[4] = "0163103100000050"; // INT7, NonLinearity TB+ConvCalo
    eventCutArray[ 5] = "00000113"; clusterCutArray[5] = "1111124060032220000"; mesonCutArray[5] = "0163103100000050"; // INT7,              TB+Calo
    eventCutArray[ 6] = "00000113"; clusterCutArray[6] = "1111125060032220000"; mesonCutArray[6] = "0163103100000050"; // INT7,              kPi0MC+ConvCalo (replay of Jasons with ConvCalo)
    eventCutArray[ 7] = "00000113"; clusterCutArray[7] = "1111126060032220000"; mesonCutArray[7] = "0163103100000050"; // INT7,              kPi0MC+Calo (replay of Jasons with Calo)
    
  } else if (trainConfig == 60){  // EMCAL clusters, EMC7 trigger
    eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000050"; // EMC7, 400 MeV min energy, NCells >=2, M02 default cut, 35ns timing, 1 cell diagonal
    eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111121060052220000"; mesonCutArray[1] = "0163103100000050"; // EMC7, 600 MeV min energy
    eventCutArray[ 2] = "00052113"; clusterCutArray[2] = "1111121060031220000"; mesonCutArray[2] = "0163103100000050"; // EMC7,                     NCells >=1
    eventCutArray[ 3] = "00052113"; clusterCutArray[3] = "1111121060033220000"; mesonCutArray[3] = "0163103100000050"; // EMC7,                     NCells >=3
  } else if (trainConfig == 61){  // EMCAL clusters, EMC7 trigger
    eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000030"; // EMC7,                                                               0.01 opening
    eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111121060032220000"; mesonCutArray[1] = "0163103100000040"; // EMC7,                                                               0.75 cell diagonal
    eventCutArray[ 2] = "00052113"; clusterCutArray[2] = "1111121060032220000"; mesonCutArray[2] = "0163103100000060"; // EMC7,                                                               2 cell diagonals
    eventCutArray[ 3] = "00052113"; clusterCutArray[3] = "1111121060032000000"; mesonCutArray[3] = "0163103100000050"; // EMC7,                                 no M02 cut
  } else if (trainConfig == 62){  // EMCAL clusters, EMC7 trigger
    eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000050"; // EMC7,                                                   50ns timing
    eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111101060032220000"; mesonCutArray[1] = "0163103100000050"; // EMC7,                standard kSDMv5
    eventCutArray[ 2] = "00052113"; clusterCutArray[2] = "1111122060032220000"; mesonCutArray[2] = "0163103100000050"; // EMC7,              NonLinearity LHC11a Calo
    eventCutArray[ 3] = "00052113"; clusterCutArray[3] = "1111100060032220000"; mesonCutArray[3] = "0163103100000050"; // EMC7,              NonLinearity none
    
  } else if (trainConfig == 63){  // EMCAL clusters, INT7 trigger
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000050"; // INT7, 400 MeV min energy, NCells >=2, M02 default cut
  } else if (trainConfig == 64){  // no non NonLinearity
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111100060032220000"; mesonCutArray[0] = "0163103100000050"; // INT7
    eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111100060032220000"; mesonCutArray[1] = "0163103100000050"; // EMC7
    eventCutArray[ 2] = "00085113"; clusterCutArray[2] = "1111100060032220000"; mesonCutArray[2] = "0163103100000050"; // EG2
    eventCutArray[ 3] = "00083113"; clusterCutArray[3] = "1111100060032220000"; mesonCutArray[3] = "0163103100000050"; // EG1

  } else if (trainConfig == 65){  // trackMatching variations
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111121061032220000"; mesonCutArray[0] = "0163103100000050"; // INT7
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111121062032220000"; mesonCutArray[1] = "0163103100000050"; //
    eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111121063032220000"; mesonCutArray[2] = "0163103100000050"; //
    eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1111121064032220000"; mesonCutArray[3] = "0163103100000050"; //
    eventCutArray[ 4] = "00000113"; clusterCutArray[4] = "1111121065032220000"; mesonCutArray[4] = "0163103100000050"; //
    eventCutArray[ 5] = "00000113"; clusterCutArray[5] = "1111121066032220000"; mesonCutArray[5] = "0163103100000050"; //
  } else if (trainConfig == 66){  // trackMatching variations
    eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111121061032220000"; mesonCutArray[0] = "0163103100000050"; // EMC7
    eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111121062032220000"; mesonCutArray[1] = "0163103100000050"; //
    eventCutArray[ 2] = "00052113"; clusterCutArray[2] = "1111121063032220000"; mesonCutArray[2] = "0163103100000050"; //
    eventCutArray[ 3] = "00052113"; clusterCutArray[3] = "1111121064032220000"; mesonCutArray[3] = "0163103100000050"; //
    eventCutArray[ 4] = "00052113"; clusterCutArray[4] = "1111121065032220000"; mesonCutArray[4] = "0163103100000050"; //
    eventCutArray[ 5] = "00052113"; clusterCutArray[5] = "1111121066032220000"; mesonCutArray[5] = "0163103100000050"; //
  } else if (trainConfig == 67){  // trackMatching variations
    eventCutArray[ 0] = "00085113"; clusterCutArray[0] = "1111121061032220000"; mesonCutArray[0] = "0163103100000050"; // EG2
    eventCutArray[ 1] = "00085113"; clusterCutArray[1] = "1111121062032220000"; mesonCutArray[1] = "0163103100000050"; //
    eventCutArray[ 2] = "00085113"; clusterCutArray[2] = "1111121063032220000"; mesonCutArray[2] = "0163103100000050"; //
    eventCutArray[ 3] = "00085113"; clusterCutArray[3] = "1111121064032220000"; mesonCutArray[3] = "0163103100000050"; //
    eventCutArray[ 4] = "00085113"; clusterCutArray[4] = "1111121065032220000"; mesonCutArray[4] = "0163103100000050"; //
    eventCutArray[ 5] = "00085113"; clusterCutArray[5] = "1111121066032220000"; mesonCutArray[5] = "0163103100000050"; //
  } else if (trainConfig == 68){  // trackMatching variations
    eventCutArray[ 0] = "00083113"; clusterCutArray[0] = "1111121061032220000"; mesonCutArray[0] = "0163103100000050"; // EG1
    eventCutArray[ 1] = "00083113"; clusterCutArray[1] = "1111121062032220000"; mesonCutArray[1] = "0163103100000050"; //
    eventCutArray[ 2] = "00083113"; clusterCutArray[2] = "1111121063032220000"; mesonCutArray[2] = "0163103100000050"; //
    eventCutArray[ 3] = "00083113"; clusterCutArray[3] = "1111121064032220000"; mesonCutArray[3] = "0163103100000050"; //
    eventCutArray[ 4] = "00083113"; clusterCutArray[4] = "1111121065032220000"; mesonCutArray[4] = "0163103100000050"; //
    eventCutArray[ 5] = "00083113"; clusterCutArray[5] = "1111121066032220000"; mesonCutArray[5] = "0163103100000050"; //
    
    // LHC11a cut studies
  } else if (trainConfig == 70){  // EMCAL clusters, MB (INT1) trigger
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000050"; // MB, 400 MeV min energy, NCells >=2, M02 default cut, 50ns timing, 1 cell diagonal
    eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111121050052220000"; mesonCutArray[1] = "0163103100000050"; // MB, 600 MeV min energy
    eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111121050031220000"; mesonCutArray[2] = "0163103100000050"; // MB,                     NCells >=1
    eventCutArray[ 3] = "00003113"; clusterCutArray[3] = "1111121050033220000"; mesonCutArray[3] = "0163103100000050"; // MB,                     NCells >=3
  } else if (trainConfig == 71){  // EMCAL clusters, MB (INT1) trigger
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000030"; // MB,                                                               0.01 opening
    eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111121050032220000"; mesonCutArray[1] = "0163103100000040"; // MB,                                                               0.75 cell diagonal
    eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111121050032220000"; mesonCutArray[2] = "0163103100000060"; // MB,                                                               2 cell diagonals
    eventCutArray[ 3] = "00003113"; clusterCutArray[3] = "1111121050032000000"; mesonCutArray[3] = "0163103100000050"; // MB,                                 no M02 cut
  } else if (trainConfig == 72){  // EMCAL clusters, MB (INT1) trigger
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000050"; // MB,                                                   35ns timing
    eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111101050032220000"; mesonCutArray[1] = "0163103100000050"; // MB,                standard kSDMv5
    eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111122050032220000"; mesonCutArray[2] = "0163103100000050"; // MB,              NonLinearity LHC11a Calo
    eventCutArray[ 3] = "00003113"; clusterCutArray[3] = "1111100050032220000"; mesonCutArray[3] = "0163103100000050"; // MB,              NonLinearity none

  } else if (trainConfig == 73){  // EMCAL clusters, EMC1 trigger
    eventCutArray[ 0] = "00051113"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000050"; // EMC1, 400 MeV min energy, NCells >=2, M02 default cut, 50ns timing, 1 cell diagonal
    eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111121050052220000"; mesonCutArray[1] = "0163103100000050"; // EMC1, 600 MeV min energy
    eventCutArray[ 2] = "00051113"; clusterCutArray[2] = "1111121050031220000"; mesonCutArray[2] = "0163103100000050"; // EMC1,                     NCells >=1
    eventCutArray[ 3] = "00051113"; clusterCutArray[3] = "1111121050033220000"; mesonCutArray[3] = "0163103100000050"; // EMC1,                     NCells >=3
  } else if (trainConfig == 74){  // EMCAL clusters, EMC1 trigger
    eventCutArray[ 0] = "00051113"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000030"; // EMC1,                                                               0.01 opening
    eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111121050032220000"; mesonCutArray[1] = "0163103100000040"; // EMC1,                                                               0.75 cell diagonal
    eventCutArray[ 2] = "00051113"; clusterCutArray[2] = "1111121050032220000"; mesonCutArray[2] = "0163103100000060"; // EMC1,                                                               2 cell diagonals
    eventCutArray[ 3] = "00051113"; clusterCutArray[3] = "1111121050032000000"; mesonCutArray[3] = "0163103100000050"; // EMC1,                                 no M02 cut
  } else if (trainConfig == 75){  // EMCAL clusters, EMC1 trigger
    eventCutArray[ 0] = "00051113"; clusterCutArray[0] = "1111121060032330000"; mesonCutArray[0] = "0163103100000050"; // EMC1,                                                   35ns timing
    eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111101050032220000"; mesonCutArray[1] = "0163103100000050"; // EMC1,                standard kSDMv5
    eventCutArray[ 2] = "00051113"; clusterCutArray[2] = "1111122050032220000"; mesonCutArray[2] = "0163103100000050"; // EMC1,              NonLinearity LHC11a Calo
    eventCutArray[ 3] = "00051113"; clusterCutArray[3] = "1111100050032220000"; mesonCutArray[3] = "0163103100000050"; // EMC1,              NonLinearity none
    
  } else if (trainConfig == 76){  // EMCAL clusters, MB (INT1) trigger, for added signals
    eventCutArray[ 0] = "00003123"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000050"; // MB, 400 MeV min energy, NCells >=2, M02 default cut, 50ns timing, 1 cell diagonal
    eventCutArray[ 1] = "00003123"; clusterCutArray[1] = "1111121050052220000"; mesonCutArray[1] = "0163103100000050"; // MB, 600 MeV min energy
    eventCutArray[ 2] = "00003123"; clusterCutArray[2] = "1111121050031220000"; mesonCutArray[2] = "0163103100000050"; // MB,                     NCells >=1
    eventCutArray[ 3] = "00003123"; clusterCutArray[3] = "1111121050033220000"; mesonCutArray[3] = "0163103100000050"; // MB,                     NCells >=3
  } else if (trainConfig == 77){  // EMCAL clusters, MB (INT1) trigger, for added signals
    eventCutArray[ 0] = "00003123"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000030"; // MB,                                                               0.01 opening
    eventCutArray[ 1] = "00003123"; clusterCutArray[1] = "1111121050032220000"; mesonCutArray[1] = "0163103100000040"; // MB,                                                               0.75 cell diagonal
    eventCutArray[ 2] = "00003123"; clusterCutArray[2] = "1111121050032220000"; mesonCutArray[2] = "0163103100000060"; // MB,                                                               2 cell diagonals
    eventCutArray[ 3] = "00003123"; clusterCutArray[3] = "1111121050032000000"; mesonCutArray[3] = "0163103100000050"; // MB,                                 no M02 cut
  } else if (trainConfig == 78){  // EMCAL clusters, MB (INT1) trigger, for added signals
    eventCutArray[ 0] = "00003123"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000050"; // MB,                                                   35ns timing
    eventCutArray[ 1] = "00003123"; clusterCutArray[1] = "1111101050032220000"; mesonCutArray[1] = "0163103100000050"; // MB,                standard kSDMv5
    eventCutArray[ 2] = "00003123"; clusterCutArray[2] = "1111122050032220000"; mesonCutArray[2] = "0163103100000050"; // MB,              NonLinearity LHC11a Calo
    eventCutArray[ 3] = "00003123"; clusterCutArray[3] = "1111100050032220000"; mesonCutArray[3] = "0163103100000050"; // MB,              NonLinearity none
    
  } else if (trainConfig == 79){  // EMCAL clusters, EMC1 trigger, for added signals
    eventCutArray[ 0] = "00051123"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000050"; // EMC1, 400 MeV min energy, NCells >=2, M02 default cut, 50ns timing, 1 cell diagonal
    eventCutArray[ 1] = "00051123"; clusterCutArray[1] = "1111121050052220000"; mesonCutArray[1] = "0163103100000050"; // EMC1, 600 MeV min energy
    eventCutArray[ 2] = "00051123"; clusterCutArray[2] = "1111121050031220000"; mesonCutArray[2] = "0163103100000050"; // EMC1,                     NCells >=1
    eventCutArray[ 3] = "00051123"; clusterCutArray[3] = "1111121050033220000"; mesonCutArray[3] = "0163103100000050"; // EMC1,                     NCells >=3
  } else if (trainConfig == 80){  // EMCAL clusters, EMC1 trigger, for added signals
    eventCutArray[ 0] = "00051123"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000030"; // EMC1,                                                               0.01 opening
    eventCutArray[ 1] = "00051123"; clusterCutArray[1] = "1111121050032220000"; mesonCutArray[1] = "0163103100000040"; // EMC1,                                                               0.75 cell diagonal
    eventCutArray[ 2] = "00051123"; clusterCutArray[2] = "1111121050032220000"; mesonCutArray[2] = "0163103100000060"; // EMC1,                                                               2 cell diagonals
    eventCutArray[ 3] = "00051123"; clusterCutArray[3] = "1111121050032000000"; mesonCutArray[3] = "0163103100000050"; // EMC1,                                 no M02 cut
  } else if (trainConfig == 81){  // EMCAL clusters, EMC1 trigger, for added signals
    eventCutArray[ 0] = "00051123"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000050"; // EMC1,                                                   35ns timing
    eventCutArray[ 1] = "00051123"; clusterCutArray[1] = "1111101050032220000"; mesonCutArray[1] = "0163103100000050"; // EMC1,                standard kSDMv5
    eventCutArray[ 2] = "00051123"; clusterCutArray[2] = "1111122050032220000"; mesonCutArray[2] = "0163103100000050"; // EMC1,              NonLinearity LHC11a Calo
    eventCutArray[ 3] = "00051123"; clusterCutArray[3] = "1111100050032220000"; mesonCutArray[3] = "0163103100000050"; // EMC1,              NonLinearity none

  } else if (trainConfig == 82){  // trackMatching variations
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121051032220000"; mesonCutArray[0] = "0163103100000050"; // MB
    eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111121052032220000"; mesonCutArray[1] = "0163103100000050"; //
    eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111121053032220000"; mesonCutArray[2] = "0163103100000050"; //
    eventCutArray[ 3] = "00003113"; clusterCutArray[3] = "1111121054032220000"; mesonCutArray[3] = "0163103100000050"; //
    eventCutArray[ 4] = "00003113"; clusterCutArray[4] = "1111121055032220000"; mesonCutArray[4] = "0163103100000050"; //
    eventCutArray[ 5] = "00003113"; clusterCutArray[5] = "1111121056032220000"; mesonCutArray[5] = "0163103100000050"; //
  } else if (trainConfig == 83){  // trackMatching variations
    eventCutArray[ 0] = "00051113"; clusterCutArray[0] = "1111121051032220000"; mesonCutArray[0] = "0163103100000050"; // EMC1
    eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111121052032220000"; mesonCutArray[1] = "0163103100000050"; //
    eventCutArray[ 2] = "00051113"; clusterCutArray[2] = "1111121053032220000"; mesonCutArray[2] = "0163103100000050"; //
    eventCutArray[ 3] = "00051113"; clusterCutArray[3] = "1111121054032220000"; mesonCutArray[3] = "0163103100000050"; //
    eventCutArray[ 4] = "00051113"; clusterCutArray[4] = "1111121055032220000"; mesonCutArray[4] = "0163103100000050"; //
    eventCutArray[ 5] = "00051113"; clusterCutArray[5] = "1111121056032220000"; mesonCutArray[5] = "0163103100000050"; //
  } else if (trainConfig == 84){  // min Energy MB
    eventCutArray[ 0] = "00003113"; clusterCutArray[0] = "1111121050012220000"; mesonCutArray[0] = "0163103100000050"; //0.2 GeV/c
    eventCutArray[ 1] = "00003113"; clusterCutArray[1] = "1111121050022220000"; mesonCutArray[1] = "0163103100000050"; //0.3 GeV/c
    eventCutArray[ 2] = "00003113"; clusterCutArray[2] = "1111121050032220000"; mesonCutArray[2] = "0163103100000050"; //0.4 GeV/c default
    eventCutArray[ 3] = "00003113"; clusterCutArray[3] = "1111121050042220000"; mesonCutArray[3] = "0163103100000050"; //0.5 GeV/c
    eventCutArray[ 4] = "00003113"; clusterCutArray[4] = "1111121050052220000"; mesonCutArray[4] = "0163103100000050"; //0.6 GeV/c
  } else if (trainConfig == 85){  // min Energy EMC1
    eventCutArray[ 0] = "00051113"; clusterCutArray[0] = "1111121050012220000"; mesonCutArray[0] = "0163103100000050"; //0.2 GeV/c
    eventCutArray[ 1] = "00051113"; clusterCutArray[1] = "1111121050022220000"; mesonCutArray[1] = "0163103100000050"; //0.3 GeV/c
    eventCutArray[ 2] = "00051113"; clusterCutArray[2] = "1111121050032220000"; mesonCutArray[2] = "0163103100000050"; //0.4 GeV/c default
    eventCutArray[ 3] = "00051113"; clusterCutArray[3] = "1111121050042220000"; mesonCutArray[3] = "0163103100000050"; //0.5 GeV/c
    eventCutArray[ 4] = "00051113"; clusterCutArray[4] = "1111121050052220000"; mesonCutArray[4] = "0163103100000050"; //0.6 GeV/c
  } else if (trainConfig == 86){  // min Energy INT7
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111121050012220000"; mesonCutArray[0] = "0163103100000050"; //0.2 GeV/c
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111121050022220000"; mesonCutArray[1] = "0163103100000050"; //0.3 GeV/c
    eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111121050032220000"; mesonCutArray[2] = "0163103100000050"; //0.4 GeV/c default
    eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1111121050042220000"; mesonCutArray[3] = "0163103100000050"; //0.5 GeV/c
    eventCutArray[ 4] = "00000113"; clusterCutArray[4] = "1111121050052220000"; mesonCutArray[4] = "0163103100000050"; //0.6 GeV/c
  } else if (trainConfig == 87){  // min Energy EMC7
    eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111121050012220000"; mesonCutArray[0] = "0163103100000050"; //0.2 GeV/c
    eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111121050022220000"; mesonCutArray[1] = "0163103100000050"; //0.3 GeV/c
    eventCutArray[ 2] = "00052113"; clusterCutArray[2] = "1111121050032220000"; mesonCutArray[2] = "0163103100000050"; //0.4 GeV/c default
    eventCutArray[ 3] = "00052113"; clusterCutArray[3] = "1111121050042220000"; mesonCutArray[3] = "0163103100000050"; //0.5 GeV/c
    eventCutArray[ 4] = "00052113"; clusterCutArray[4] = "1111121050052220000"; mesonCutArray[4] = "0163103100000050"; //0.6 GeV/c
  } else if (trainConfig == 88){  // min Energy EG2
    eventCutArray[ 0] = "00085113"; clusterCutArray[0] = "1111121050012220000"; mesonCutArray[0] = "0163103100000050"; //0.2 GeV/c
    eventCutArray[ 1] = "00085113"; clusterCutArray[1] = "1111121050022220000"; mesonCutArray[1] = "0163103100000050"; //0.3 GeV/c
    eventCutArray[ 2] = "00085113"; clusterCutArray[2] = "1111121050032220000"; mesonCutArray[2] = "0163103100000050"; //0.4 GeV/c default
    eventCutArray[ 3] = "00085113"; clusterCutArray[3] = "1111121050042220000"; mesonCutArray[3] = "0163103100000050"; //0.5 GeV/c
    eventCutArray[ 4] = "00085113"; clusterCutArray[4] = "1111121050052220000"; mesonCutArray[4] = "0163103100000050"; //0.6 GeV/c
  } else if (trainConfig == 89){  // min Energy EG1
    eventCutArray[ 0] = "00083113"; clusterCutArray[0] = "1111121050012220000"; mesonCutArray[0] = "0163103100000050"; //0.2 GeV/c
    eventCutArray[ 1] = "00083113"; clusterCutArray[1] = "1111121050022220000"; mesonCutArray[1] = "0163103100000050"; //0.3 GeV/c
    eventCutArray[ 2] = "00083113"; clusterCutArray[2] = "1111121050032220000"; mesonCutArray[2] = "0163103100000050"; //0.4 GeV/c default
    eventCutArray[ 3] = "00083113"; clusterCutArray[3] = "1111121050042220000"; mesonCutArray[3] = "0163103100000050"; //0.5 GeV/c
    eventCutArray[ 4] = "00083113"; clusterCutArray[4] = "1111121050052220000"; mesonCutArray[4] = "0163103100000050"; //0.6 GeV/c
    
  } else if (trainConfig == 98){ // MB - with multiplicity bins
    eventCutArray[ 0] = "00103113"; clusterCutArray[0] = "1111121050032220000"; mesonCutArray[0] = "0163103100000050"; // 0 -2
    eventCutArray[ 1] = "01203113"; clusterCutArray[1] = "1111121050032220000"; mesonCutArray[1] = "0163103100000050"; // 2 -5
    eventCutArray[ 2] = "02303113"; clusterCutArray[2] = "1111121050032220000"; mesonCutArray[2] = "0163103100000050"; // 5 -10
    eventCutArray[ 3] = "03503113"; clusterCutArray[3] = "1111121050032220000"; mesonCutArray[3] = "0163103100000050"; // 10 -30
    eventCutArray[ 4] = "05703113"; clusterCutArray[4] = "1111121050032220000"; mesonCutArray[4] = "0163103100000050"; // 30 -100
  } else if (trainConfig == 99){ // INT7 - with multiplicity bins
    eventCutArray[ 0] = "00100113"; clusterCutArray[0] = "1111121060032220000"; mesonCutArray[0] = "0163103100000050"; // 0 -2
    eventCutArray[ 1] = "01200113"; clusterCutArray[1] = "1111121060032220000"; mesonCutArray[1] = "0163103100000050"; // 2 -5
    eventCutArray[ 2] = "02300113"; clusterCutArray[2] = "1111121060032220000"; mesonCutArray[2] = "0163103100000050"; // 5 -10
    eventCutArray[ 3] = "03500113"; clusterCutArray[3] = "1111121060032220000"; mesonCutArray[3] = "0163103100000050"; // 10 -30
    eventCutArray[ 4] = "05700113"; clusterCutArray[4] = "1111121060032220000"; mesonCutArray[4] = "0163103100000050"; // 30 -100
    
    
// 8 TeV configs

    // here is the order of the cluster cut string
    // usually for EMCal we start with 11111: default values for            "ClusterType", "EtaMin", "EtaMax", "PhiMin", "PhiMax"
    // then two numbers for nonlinearity, e.g. 21: this is                  "NonLinearity1", "NonLinearity2"
    // Then some cuts on the clusters, e.g. 06003222: this is               "DistanceToBadChannel", "Timing", "TrackMatching", "ExoticCell", "MinEnergy", "MinNCells", "MinM02", "MaxM02"
    // finally some for now unused cuts, usually 0000: this is              "MinM20", "MaxM20", "MaximumDispersion", "NLM"

    //standard cut
  } else if (trainConfig == 101){ // EMCAL clusters pp 8 TeV 
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060032230000"; mesonCutArray[0] = "0163103100000050"; // 400 MeV cluster min energy
    // 8 TeV variations
  } else if (trainConfig == 102){ // EMCAL clusters pp 8 TeV, timing variation
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111050032230000"; mesonCutArray[0] = "0163103100000050"; // time -50ns_50ns
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111060032230000"; mesonCutArray[1] = "0163103100000050"; // time -30ns_35ns - standard
    eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111111070032230000"; mesonCutArray[2] = "0163103100000050"; // time -30ns_30ns
    eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1111111080032230000"; mesonCutArray[3] = "0163103100000050"; // time -20ns_30ns
  } else if (trainConfig == 103){ //EMCAL minEnergy variation
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060012230000"; mesonCutArray[0] = "0163103100000050"; //0.2 GeV/c
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111060022230000"; mesonCutArray[1] = "0163103100000050"; //0.3 GeV/c
    eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111111060032230000"; mesonCutArray[2] = "0163103100000050"; //0.4 GeV/c default
    eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1111111060042230000"; mesonCutArray[3] = "0163103100000050"; //0.5 GeV/c
    eventCutArray[ 4] = "00000113"; clusterCutArray[4] = "1111111060052230000"; mesonCutArray[4] = "0163103100000050"; //0.6 GeV/c
  } else if (trainConfig == 104){ //EMCAL minNCells, M02, with/without TRD variation
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060031230000"; mesonCutArray[0] = "0163103100000050"; //n cells >= 1
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111060033230000"; mesonCutArray[1] = "0163103100000050"; //n cells >= 3
    eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111111060032000000"; mesonCutArray[2] = "0163103100000050"; //no M02 cut
    eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1113111060032230000"; mesonCutArray[3] = "0163103100000050"; //only modules with TRD infront
    eventCutArray[ 4] = "00000113"; clusterCutArray[4] = "1111211060032230000"; mesonCutArray[4] = "0163103100000050"; //no modules with TRD infront
  } else if (trainConfig == 105){  // trackMatching variations
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111061032230000"; mesonCutArray[0] = "0163103100000050"; //
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111062032230000"; mesonCutArray[1] = "0163103100000050"; //
    eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111111063032230000"; mesonCutArray[2] = "0163103100000050"; //
    eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1111111064032230000"; mesonCutArray[3] = "0163103100000050"; //
    eventCutArray[ 4] = "00000113"; clusterCutArray[4] = "1111111065032230000"; mesonCutArray[4] = "0163103100000050"; //
    eventCutArray[ 5] = "00000113"; clusterCutArray[5] = "1111111066032230000"; mesonCutArray[5] = "0163103100000050"; //

  } else if (trainConfig == 109){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111100060032230000"; mesonCutArray[0] = "0163103100000050"; // NonLinearity none
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111113060032230000"; mesonCutArray[1] = "0163103100000050"; // NonLinearity kTestBeamv2 + LHC12 ConvCalo
    eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111114060032230000"; mesonCutArray[2] = "0163103100000050"; // NonLinearity kTestBeamv2 + LHC12 Calo
  } else if (trainConfig == 110){ // EMCAL clusters pp 8 TeV, Different NonLinearities
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111101060032230000"; mesonCutArray[0] = "0163103100000050"; // NonLinearity kSDMv5
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111060032230000"; mesonCutArray[1] = "0163103100000050"; // NonLinearity LHC12 ConvCalo
    eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111112060032230000"; mesonCutArray[2] = "0163103100000050"; // NonLinearity LHC12 Calo
    eventCutArray[ 3] = "00000113"; clusterCutArray[3] = "1111100060032230000"; mesonCutArray[3] = "0163103100000050"; // NonLinearity none

   // LHC12fa-i and MC
    // default with three cuts
  } else if (trainConfig == 111){  // EMCAL clusters, different triggers no NonLinearity
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111100060032230000"; mesonCutArray[0] = "0163103100000050"; // INT7
    eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111100060032230000"; mesonCutArray[1] = "0163103100000050"; // EMC7
    eventCutArray[ 2] = "00081113"; clusterCutArray[2] = "1111100060032230000"; mesonCutArray[2] = "0163103100000050"; // EMCEG1,
    
    // all the cut variations
  } else if (trainConfig == 112){  // EMCAL clusters, EMCEG1 trigger
    eventCutArray[ 0] = "00081113"; clusterCutArray[0] = "1111111060032230000"; mesonCutArray[0] = "0163103100000050"; // EMCEGA, 400 MeV min energy, NCells >=2, M02 default cut
    eventCutArray[ 1] = "00081113"; clusterCutArray[1] = "1111111060052230000"; mesonCutArray[1] = "0163103100000050"; // EMCEGA, 600 MeV min energy
  } else if (trainConfig == 113){  // EMCAL clusters, EMCEG1 trigger
    eventCutArray[ 0] = "00081113"; clusterCutArray[0] = "1111111060031230000"; mesonCutArray[0] = "0163103100000050"; // EMCEGA,                     NCells >=1
    eventCutArray[ 1] = "00081113"; clusterCutArray[1] = "1111111060033230000"; mesonCutArray[1] = "0163103100000050"; // EMCEGA,                     NCells >=3
  } else if (trainConfig == 114){  // EMCAL clusters, EMCEG1 trigger
    eventCutArray[ 0] = "00081113"; clusterCutArray[0] = "1111111060032000000"; mesonCutArray[0] = "0163103100000050"; // EMCEGA,                                 no M02 cut
    eventCutArray[ 1] = "00081113"; clusterCutArray[1] = "1111111020032230000"; mesonCutArray[1] = "0163103100000050"; // EMCEGA,                                                 500ns timing
    eventCutArray[ 2] = "00085113"; clusterCutArray[2] = "1111111040032230000"; mesonCutArray[2] = "0163103100000050"; // EMCEGA,                                                 100ns timing
  } else if (trainConfig == 115){  // EMCAL clusters, INT7 trigger
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060032230000"; mesonCutArray[0] = "0163103100000050"; // INT7, 400 MeV min energy, NCells >=2, M02 default cut
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111060052230000"; mesonCutArray[1] = "0163103100000050"; // INT7, 600 MeV min energy
  } else if (trainConfig == 116){  // EMCAL clusters, INT7 trigger
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060031230000"; mesonCutArray[0] = "0163103100000050"; // INT7,                       NCells >=1
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111060033230000"; mesonCutArray[1] = "0163103100000050"; // INT7,                       NCells >=3
  } else if (trainConfig == 117){  // EMCAL clusters, INT7 trigger
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060032000000"; mesonCutArray[0] = "0163103100000050"; // INT7,                                   no M02 cut
    eventCutArray[ 1] = "00000113"; clusterCutArray[1] = "1111111020032230000"; mesonCutArray[1] = "0163103100000050"; // INT7,                                                   500ns timing
    eventCutArray[ 2] = "00000113"; clusterCutArray[2] = "1111111040032230000"; mesonCutArray[2] = "0163103100000050"; // INT7,                                                   100ns timing
  } else if (trainConfig == 118){  // EMCAL clusters, EMC7 trigger
    eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111111060032230000"; mesonCutArray[0] = "0163103100000050"; // EMC7, 400 MeV min energy, NCells >=2, M02 default cut
    eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111111060052230000"; mesonCutArray[1] = "0163103100000050"; // EMC7, 600 MeV min energy
  } else if (trainConfig == 119){  // EMCAL clusters, EMC7 trigger
    eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111111060031230000"; mesonCutArray[0] = "0163103100000050"; // EMC7,                     NCells >=1
    eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111111060033230000"; mesonCutArray[1] = "0163103100000050"; // EMC7,                     NCells >=3
  } else if (trainConfig == 120){  // EMCAL clusters, EMC7 trigger
    eventCutArray[ 0] = "00052113"; clusterCutArray[0] = "1111111060032000000"; mesonCutArray[0] = "0163103100000050"; // EMC7,                                 no M02 cut
    eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111111020032230000"; mesonCutArray[1] = "0163103100000050"; // EMC7,                                                 500ns timing
    eventCutArray[ 2] = "00052113"; clusterCutArray[2] = "1111111040032230000"; mesonCutArray[2] = "0163103100000050"; // EMC7,                                                 100ns timing

  }else if (trainConfig == 121){ // EMCAL clusters, different special triggers, different conv calo non lin
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111111060032230000"; mesonCutArray[0] = "0163103100000050"; // INT7
    eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111111060032230000"; mesonCutArray[1] = "0163103100000050"; // EMC7
    eventCutArray[ 2] = "00081113"; clusterCutArray[2] = "1111111060032230000"; mesonCutArray[2] = "0163103100000050"; // EMCEG1,
  }else if (trainConfig == 122){ // EMCAL clusters, different special triggers, different kSDMv5
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111101060032230000"; mesonCutArray[0] = "0163103100000050"; // INT7
    eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111101060032230000"; mesonCutArray[1] = "0163103100000050"; // EMC7
    eventCutArray[ 2] = "00081113"; clusterCutArray[2] = "1111101060032230000"; mesonCutArray[2] = "0163103100000050"; // EMCEG1,
  }else if (trainConfig == 123){ // EMCAL clusters, different special triggers, different conv calo non lin
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111112060032230000"; mesonCutArray[0] = "0163103100000050"; // INT7
    eventCutArray[ 1] = "00052113"; clusterCutArray[1] = "1111112060032230000"; mesonCutArray[1] = "0163103100000050"; // EMC7
    eventCutArray[ 2] = "00081113"; clusterCutArray[2] = "1111112060032230000"; mesonCutArray[2] = "0163103100000050"; // EMCEG1,


    // 7 TeV
  } else if (trainConfig == 201){ // EMCAL clusters pp 7 TeV
    eventCutArray[ 0] = "00000113"; clusterCutArray[0] = "1111100010032230000"; mesonCutArray[0] = "0163103100000050"; // 1000ns timing cut

  } else {
    Error(Form("GammaCalo_%i",trainConfig), "wrong trainConfig variable no cuts have been specified for the configuration");
    return;
  }

  TList *EventCutList = new TList();
  TList *ClusterCutList = new TList();
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
  ClusterCutList->SetOwner(kTRUE);
  AliCaloPhotonCuts **analysisClusterCuts = new AliCaloPhotonCuts*[numberOfCuts];
  MesonCutList->SetOwner(kTRUE);
  AliConversionMesonCuts **analysisMesonCuts = new AliConversionMesonCuts*[numberOfCuts];

  for(Int_t i = 0; i<numberOfCuts; i++){
    analysisEventCuts[i] = new AliConvEventCuts();   
    
    // definition of weighting input
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
    
    if (doWeighting) analysisEventCuts[i]->SetUseReweightingWithHistogramFromFile(kTRUE, kTRUE, kFALSE, fileNameInputForWeighting, mcInputNamePi0, mcInputNameEta, "",fitNamePi0,fitNameEta);

    analysisEventCuts[i]->SetTriggerMimicking(enableTriggerMimicking);
    analysisEventCuts[i]->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
    analysisEventCuts[i]->SetMaxFacPtHard(maxFacPtHard);
    analysisEventCuts[i]->InitializeCutsFromCutString(eventCutArray[i].Data());
    EventCutList->Add(analysisEventCuts[i]);
    analysisEventCuts[i]->SetFillCutHistograms("",kFALSE);
    
    analysisClusterCuts[i] = new AliCaloPhotonCuts((isMC==2));
    analysisClusterCuts[i]->SetIsPureCaloCut(2);
    analysisClusterCuts[i]->InitializeCutsFromCutString(clusterCutArray[i].Data());
    ClusterCutList->Add(analysisClusterCuts[i]);
    analysisClusterCuts[i]->SetExtendedQA(enableExtMatchAndQA);
    analysisClusterCuts[i]->SetFillCutHistograms("");
    
    analysisMesonCuts[i] = new AliConversionMesonCuts();
    analysisMesonCuts[i]->InitializeCutsFromCutString(mesonCutArray[i].Data());
    analysisMesonCuts[i]->SetIsMergedClusterCut(2);
    MesonCutList->Add(analysisMesonCuts[i]);
    analysisMesonCuts[i]->SetFillCutHistograms("");
    analysisEventCuts[i]->SetAcceptedHeader(HeaderList);
  }
  task->SetEventCutList(numberOfCuts,EventCutList);
  task->SetCaloCutList(numberOfCuts,ClusterCutList);
  task->SetMesonCutList(numberOfCuts,MesonCutList);
  task->SetDoMesonAnalysis(kTRUE);
  task->SetDoMesonQA(enableQAMesonTask); //Attention new switch for Pi0 QA
  task->SetDoClusterQA(enableQAClusterTask);  //Attention new switch small for Cluster QA
  task->SetDoTHnSparse(isUsingTHnSparse);
  if(enableExtMatchAndQA == 2 || enableExtMatchAndQA == 3){ task->SetPlotHistsExtQA(kTRUE);}
  
  //connect containers
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("GammaCalo_%i",trainConfig), TList::Class(),
              AliAnalysisManager::kOutputContainer,Form("GammaCalo_%i.root",trainConfig));

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,1,coutput);

  return;

}
