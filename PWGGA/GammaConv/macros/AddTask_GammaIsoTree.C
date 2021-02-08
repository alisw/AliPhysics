void AddTask_GammaIsoTree(
  Int_t     trainConfig                   = 1,
  Int_t     isMC                          = 0,
  Int_t     IsHeavyIon                    = 0,
  TString   photonCutNumberV0Reader       = "00200009327000008250400000",
  TString   periodNameV0Reader            = "",
  Bool_t    useHistograms                 = kFALSE, // if activated, analysis will be performed hist based instead of cut based
  Int_t     enableExtMatchAndQA           = 0,
  Bool_t    enableTriggerOverlapRej       = kTRUE,
  Int_t     enableTriggerMimicking        = 0,        // enable trigger mimicking
  TString   settingMaxFacPtHard           = "3.",       // maximum factor between hardest jet and ptHard generated
  Bool_t    makeAdditionalHistos          = kFALSE,
  Bool_t    storeTracks                   = kTRUE,
  Bool_t    storeEMCalCluster             = kTRUE,
  Bool_t    storePHOSCluster              = kTRUE,
  Bool_t    storeConversions              = kTRUE,
  Bool_t    doIsolation                   = kTRUE,
  Bool_t    doOwnTrackMatching            = kFALSE,
  // subwagon config
  TString   fileNameExternalInputs        = "",
  Double_t  genPtCut                      = 0, // only save particles from stack with gen pt > genPtCut
  Double_t  recPtCut                      = 0, // only save clusters with rec pt > recPtCut
  TString   additionalTrainConfig         = "0"       // additional counter for trainconfig
  ){

  //
  // ─── SET CONFIG ─────────────────────────────────────────────────────────────────
  //
  AliCutHandlerPCM cuts(13); // only for tokenize
  TString addTaskName = "AddTask_GammaIsoTree";

  // Default
  TString   TaskEventCutnumber                = "00010113";
  TString   TaskClusterCutnumberEMC           = "1111100010022700000";
  TString   TaskTMCut                         = "111110001f022700000"; // should be like normal cut, but only handles the track
                                                                      // matching for book keeping
  TString   TaskClusterCutnumberIsolationEMC = "111110001f022700000";
  TString   TaskClusterCutnumberTaggingEMC = "111110001f022700000";
  TString   TaskClusterCutnumberPHOS          = "2444411044013300000";
  TString   TaskConvCutnumber                 = "0dm0000922700000dge0404000";


  vector<Float_t> trackIsoR = {0.2,0.3,0.4};
  vector<Double_t> trackIsoE = {0.5,1.5,2.5};
  vector<Float_t> neutralIsoR = {0.2,0.3,0.4};
  vector<Double_t> neutralIsoE = {0.5,1.5,2.5};
  Double_t minSignalM02 = 0.1;
  Double_t maxSignalM02 = 0.5;
  Int_t trackMatcherRunningMode = 0; // CaloTrackMatcher running mode
  Bool_t backgroundTrackMatching = kFALSE; // obsolete
  Bool_t doNeutralIso            = kTRUE;
  Bool_t doChargedIso            = kTRUE;
  Bool_t doCellIso               = kTRUE;
  Bool_t doTagging               = kTRUE;

  TString sAdditionalTrainConfig      = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "", "", addTaskName);
  if (sAdditionalTrainConfig.Atoi() > 0){
    trainConfig = trainConfig + sAdditionalTrainConfig.Atoi();
    cout << "INFO: running additionalTrainConfig '" << sAdditionalTrainConfig.Atoi() << "', train config: '" << trainConfig << "'" << endl;
  }

  TString fileNamePtWeights           = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FPTW:");
  TString fileNameMultWeights         = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FMUW:");
  TString fileNameCustomTriggerMimicOADB   = cuts.GetSpecialFileNameFromString (fileNameExternalInputs, "FTRM:");

  TString corrTaskSetting             = cuts.GetSpecialSettingFromAddConfig(additionalTrainConfig, "CF", "", addTaskName);
  if(corrTaskSetting.CompareTo(""))
    cout << "corrTaskSetting: " << corrTaskSetting.Data() << endl;
  
  
  // pp 8 TeV
  // ────────────────────────────────────────────────────────────────────────────────
  Int_t                       fMinClsTPC = 70;  
  Double_t                    fChi2PerClsTPC = 5;   
  Int_t                       fMinClsITS = 0;  
  Double_t                    fEtaCut = 0.9;  
  Double_t                    fPtCut= 0.1;  
  Double_t                    fYMCCut = 9999;  

  Double_t                    fAntiIsolation[2] = {5.,10};
  if(trainConfig == 1){ 
      TaskEventCutnumber                = "00010113";
      TaskClusterCutnumberEMC           = "1111132060032230000";
      TaskClusterCutnumberIsolationEMC = "1111100060022700000";
      TaskClusterCutnumberTaggingEMC = "1111100060022700000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
  } else if(trainConfig == 2){ 
      TaskEventCutnumber                = "00052113";
      TaskClusterCutnumberEMC           = "1111132060032230000";
      TaskClusterCutnumberIsolationEMC = "1111100060022700000";
      TaskClusterCutnumberTaggingEMC = "1111100060022700000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";
  } else if(trainConfig == 3){ 
      TaskEventCutnumber                = "00081113";
      TaskClusterCutnumberEMC           = "1111132010032230000";
      TaskClusterCutnumberIsolationEMC = "1111100060022700000";
      TaskClusterCutnumberTaggingEMC = "1111100060022700000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

  } else if(trainConfig == 4){  // min bias loose cluster cuts
      TaskEventCutnumber                = "00010113";
      TaskClusterCutnumberEMC           = "111113200f000000000";
      TaskClusterCutnumberIsolationEMC  = "111113206f022700000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kTRUE;
      doChargedIso = kTRUE;
      doTagging = kTRUE;
      doCellIso = kTRUE;
  } else if(trainConfig == 5){  // min bias loose cluster cuts
      TaskEventCutnumber                = "00052103";
      TaskClusterCutnumberEMC           = "111113200f000000000";
      TaskClusterCutnumberIsolationEMC = "111113206f022700000";
      TaskClusterCutnumberTaggingEMC = "111113206f022700000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kTRUE;
      doChargedIso = kTRUE;
      doTagging = kTRUE;
      doCellIso = kTRUE;

  // cut based study
  } else if(trainConfig == 6){  // min bias (cuts from PCMEMC 84 + loose iso)
      TaskEventCutnumber                = "00010103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;
  } else if(trainConfig == 7){  // trigger
      TaskEventCutnumber                = "00052103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;
  } else if(trainConfig == 8){  // trigger
      TaskEventCutnumber                = "00081103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;

   // Anti isolation variation

  } else if(trainConfig == 10){  // min bias (cuts from PCMEMC 84 + loose iso)
      TaskEventCutnumber                = "00010103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;
      
      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;

      fAntiIsolation[0] = 4;
      fAntiIsolation[1] = 10;
  } else if(trainConfig == 11){  // trigger
      TaskEventCutnumber                = "00052103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;

      fAntiIsolation[0] = 4;
      fAntiIsolation[1] = 10;
  } else if(trainConfig == 12){  // trigger
      TaskEventCutnumber                = "00081103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;

      fAntiIsolation[0] = 4;
      fAntiIsolation[1] = 10;

  } else if(trainConfig == 13){  // min bias (cuts from PCMEMC 84 + loose iso)
      TaskEventCutnumber                = "00010103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;

      fAntiIsolation[0] = 5;
      fAntiIsolation[1] = 20;
  } else if(trainConfig == 14){  // trigger
      TaskEventCutnumber                = "00052103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;

      fAntiIsolation[0] = 5;
      fAntiIsolation[1] = 20;
  } else if(trainConfig == 15){  // trigger
      TaskEventCutnumber                = "00081103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;

      fAntiIsolation[0] = 5;
      fAntiIsolation[1] = 20;


  // pPb 8 TeV
  // ────────────────────────────────────────────────────────────────────────────────
  } else if(trainConfig == 50){  // pPb INT7
      TaskEventCutnumber                = "80010103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
            TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;
  } else if(trainConfig == 51){  // EG2
      TaskEventCutnumber                = "80085103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;
  } else if(trainConfig == 52){  // EG1
      TaskEventCutnumber                = "80083103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;

// Anti isolation variation
  } else if(trainConfig == 60){  // pPb INT7
      TaskEventCutnumber                = "80010103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;

      fAntiIsolation[0] = 4;
      fAntiIsolation[1] = 10;
  } else if(trainConfig == 61){  // EG2
      TaskEventCutnumber                = "80085103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;
      fAntiIsolation[0] = 4;
      fAntiIsolation[1] = 10;
  } else if(trainConfig == 62){  // EG1
      TaskEventCutnumber                = "80083103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;
      fAntiIsolation[0] = 4;
      fAntiIsolation[1] = 10;

 } else if(trainConfig == 63){  // pPb INT7
      TaskEventCutnumber                = "80010103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;

      fAntiIsolation[0] = 5;
      fAntiIsolation[1] = 20;
  } else if(trainConfig == 64){  // EG2
      TaskEventCutnumber                = "80085103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;
      fAntiIsolation[0] = 5;
      fAntiIsolation[1] = 20;
  } else if(trainConfig == 65){  // EG1
      TaskEventCutnumber                = "80083103";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;
      fAntiIsolation[0] = 5;
      fAntiIsolation[1] = 20;

  } else if(trainConfig == 70){  // min bias (cuts from PCMEMC 84 + loose iso)
      TaskEventCutnumber                = "00000000";
      TaskClusterCutnumberEMC           = "1111132060032000000";
      TaskTMCut = TaskClusterCutnumberEMC.Data();
      TaskTMCut.Replace(9,1,"5");
      TaskClusterCutnumberIsolationEMC  = "111113206f022000000";
      TaskClusterCutnumberTaggingEMC    = "111113206f000000000";
      TaskClusterCutnumberPHOS          = "2444411044013300000";
      TaskConvCutnumber                 = "0dm00009f9730000dge0404000";

      minSignalM02 = 0.1;
      maxSignalM02 = 0.5;

      backgroundTrackMatching = kFALSE; // obsolete
      doNeutralIso = kFALSE;
      doChargedIso = kTRUE;
      doTagging = kFALSE;
      doCellIso = kFALSE;
  }
  

  // ================== GetAnalysisManager ===============================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return ;
  }

  // ================== GetInputEventHandler =============================
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

//=========  Set Cutnumber for V0Reader ================================
  TString cutnumberPhoton     = photonCutNumberV0Reader.Data();
  TString cutnumberEvent      = "00000003";
  if(IsHeavyIon==1)
    cutnumberEvent = "10000003";
  else if(IsHeavyIon==2)
    cutnumberEvent = "80000003";

//========= Check V0 Reader in  ANALYSIS manager  =====
  TString V0ReaderName        = Form("V0ReaderV1_%s_%s",cutnumberEvent.Data(),cutnumberPhoton.Data());
  AliV0ReaderV1 *fV0ReaderV1  =  NULL;
  if( !(AliV0ReaderV1*)mgr->GetTask(V0ReaderName.Data()) ){
    cout << "V0Reader: " << V0ReaderName.Data() << " not found!!"<< endl;
    return;
  } else {
    cout << "V0Reader: " << V0ReaderName.Data() << " found!!"<< endl;
  }

  AliAnalysisTaskGammaIsoTree *fQA = new AliAnalysisTaskGammaIsoTree("GammaIsoTree");
  
  fQA->SetV0ReaderName(V0ReaderName);
  fQA->SetIsMC(isMC);
  fQA->SetYCutMC(0.9);
  fQA->SetAntiIsolationE(fAntiIsolation[0],fAntiIsolation[1]);
  fQA->SetRecPtCut(recPtCut);
  fQA->SetGenPtCut(genPtCut);
  
  // fQA->SetSaveClusterCells(doSaveClusterCells);
  // fQA->SetSaveEventProperties(doSaveEventProp);
  // fQA->SetDoAdditionalHistos(makeAdditionalHistos);

  TObjArray *rmaxFacPtHardSetting = settingMaxFacPtHard.Tokenize("_");
  if(rmaxFacPtHardSetting->GetEntries()<1){cout << "ERROR: AddTask_GammaIsoTree during parsing of settingMaxFacPtHard String '" << settingMaxFacPtHard.Data() << "'" << endl; return;}
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

  // add track matcher if do own trackmatching is enabled

  TString TrackMatcherNameSignal = Form("CaloTrackMatcher_Signal_%s_%i",TaskTMCut.Data(),trackMatcherRunningMode);
  TString TrackMatcherNameIsolation = Form("CaloTrackMatcher_Isolation_%s_%i",TaskClusterCutnumberIsolationEMC.Data(),trackMatcherRunningMode);
  TString TrackMatcherNameTagging = Form("CaloTrackMatcher_Tagging_%s_%i",TaskClusterCutnumberTaggingEMC.Data(),trackMatcherRunningMode);
  
  TString clusterTypeStringSignal(TaskClusterCutnumberEMC(0,1));
  Int_t clusterTypeSignal = clusterTypeStringSignal.Atoi();

  TString clusterTypeStringIsolation(TaskClusterCutnumberIsolationEMC(0,1));
  Int_t clusterTypeIsolation = clusterTypeStringIsolation.Atoi();

  TString clusterTypeStringTagging(TaskClusterCutnumberTaggingEMC(0,1));
  Int_t clusterTypeTagging = clusterTypeStringTagging.Atoi();
  if(!doOwnTrackMatching){
    
    // matching for signal clusters
    if(corrTaskSetting.CompareTo("")){
      TrackMatcherNameSignal = TrackMatcherNameSignal+"_"+corrTaskSetting.Data();
      cout << "Using separate track matcher for correction framework setting: " << TrackMatcherNameSignal.Data() << endl;
    }
    if( !(AliCaloTrackMatcher*)mgr->GetTask(TrackMatcherNameSignal.Data()) ){
      AliCaloTrackMatcher* fTrackMatcherSignal = new AliCaloTrackMatcher(TrackMatcherNameSignal.Data(),clusterTypeSignal,trackMatcherRunningMode);
      fTrackMatcherSignal->SetV0ReaderName(V0ReaderName);
      fTrackMatcherSignal->SetCorrectionTaskSetting(corrTaskSetting);
      mgr->AddTask(fTrackMatcherSignal);
      mgr->ConnectInput(fTrackMatcherSignal,0,cinput);
    }
  }

  // Create Event Cuts
  AliConvEventCuts *analysisEventCuts = new AliConvEventCuts();
  analysisEventCuts->SetV0ReaderName(V0ReaderName);
  analysisEventCuts->SetTriggerOverlapRejecion(enableTriggerOverlapRej);
  if(fMinPtHardSet)
    analysisEventCuts->SetMinFacPtHard(minFacPtHard);
  if(fMaxPtHardSet)
    analysisEventCuts->SetMaxFacPtHard(maxFacPtHard);
  if(fSingleMaxPtHardSet)
    analysisEventCuts->SetMaxFacPtHardSingleParticle(maxFacPtHardSingle);
  if(fJetFinderUsage)
      analysisEventCuts->SetUseJetFinderForOutliers(kTRUE);
  if(fUsePtHardFromFile)
    analysisEventCuts->SetUsePtHardBinFromFile(kTRUE);
  if(fUseAddOutlierRej)
    analysisEventCuts->SetUseAdditionalOutlierRejection(kTRUE);
  analysisEventCuts->SetCorrectionTaskSetting(corrTaskSetting);
  if (periodNameV0Reader.CompareTo("") != 0) analysisEventCuts->SetPeriodEnum(periodNameV0Reader);
  analysisEventCuts->InitializeCutsFromCutString(TaskEventCutnumber.Data());
  analysisEventCuts->SetTriggerMimicking(enableTriggerMimicking);
  analysisEventCuts->SetFillCutHistograms("",kFALSE);

  // EMC signal cluster cuts (used to store in tree)
  AliCaloPhotonCuts *analysisClusterCutsEMC = new AliCaloPhotonCuts(isMC,"analysisClusterCutsEMC","analysisClusterCutsEMC");
  analysisClusterCutsEMC->SetV0ReaderName(V0ReaderName);
  analysisClusterCutsEMC->SetCorrectionTaskSetting(corrTaskSetting);
  analysisClusterCutsEMC->SetCaloTrackMatcherName(TrackMatcherNameSignal);
  analysisClusterCutsEMC->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsEMC->InitializeCutsFromCutString(TaskClusterCutnumberEMC.Data());
  analysisClusterCutsEMC->SetFillCutHistograms("");

  // EMC background cluster cuts (used to calculate iso and tagging)
  AliCaloPhotonCuts *analysisClusterCutsIsolationEMC = new AliCaloPhotonCuts(isMC,"analysisClusterCutsIsolationEMC","analysisClusterCutsIsolationEMC");
  analysisClusterCutsIsolationEMC->SetV0ReaderName(V0ReaderName);
  analysisClusterCutsIsolationEMC->SetCorrectionTaskSetting(corrTaskSetting);
  analysisClusterCutsIsolationEMC->SetCaloTrackMatcherName(TrackMatcherNameIsolation);
  analysisClusterCutsIsolationEMC->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsIsolationEMC->InitializeCutsFromCutString(TaskClusterCutnumberIsolationEMC.Data());
  analysisClusterCutsIsolationEMC->SetFillCutHistograms("");

  AliCaloPhotonCuts *analysisClusterCutsTaggingEMC = new AliCaloPhotonCuts(isMC,"analysisClusterCutsTaggingEMC","analysisClusterCutsTaggingEMC");
  analysisClusterCutsTaggingEMC->SetV0ReaderName(V0ReaderName);
  analysisClusterCutsTaggingEMC->SetCorrectionTaskSetting(corrTaskSetting);
  analysisClusterCutsTaggingEMC->SetCaloTrackMatcherName(TrackMatcherNameTagging);
  analysisClusterCutsTaggingEMC->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsTaggingEMC->InitializeCutsFromCutString(TaskClusterCutnumberTaggingEMC.Data());
  analysisClusterCutsTaggingEMC->SetFillCutHistograms("");

  // Used only to handle track matching
  AliCaloPhotonCuts *analysisClusterCutsEMCTrackMatching = new AliCaloPhotonCuts(isMC,"analysisClusterCutsEMCTrackMatching","analysisClusterCutsEMCTrackMatching");
  analysisClusterCutsEMCTrackMatching->SetV0ReaderName(V0ReaderName);
  analysisClusterCutsEMCTrackMatching->SetCorrectionTaskSetting(corrTaskSetting);
  analysisClusterCutsEMCTrackMatching->SetCaloTrackMatcherName(TrackMatcherNameSignal);
  analysisClusterCutsEMCTrackMatching->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsEMCTrackMatching->InitializeCutsFromCutString(TaskTMCut.Data());
  analysisClusterCutsEMCTrackMatching->SetFillCutHistograms("");

  // PHOS cluster cuts
  AliCaloPhotonCuts *analysisClusterCutsPHOS = new AliCaloPhotonCuts(isMC,"analysisClusterCutsPHOS","analysisClusterCutsPHOS");
  analysisClusterCutsPHOS->SetV0ReaderName(V0ReaderName);
  // analysisClusterCutsPHOS->SetCaloTrackMatcherName(TrackMatcherNamePHOS);
  analysisClusterCutsPHOS->SetExtendedMatchAndQA(enableExtMatchAndQA);
  analysisClusterCutsPHOS->InitializeCutsFromCutString(TaskClusterCutnumberPHOS.Data());
  analysisClusterCutsPHOS->SetFillCutHistograms("");
  
  AliConversionPhotonCuts* analysisConvCuts = new AliConversionPhotonCuts();
  analysisConvCuts->SetV0ReaderName(V0ReaderName);
  analysisConvCuts->InitializeCutsFromCutString(TaskConvCutnumber.Data());
  analysisConvCuts->SetFillCutHistograms("");
  
  fQA->SetEventCuts(analysisEventCuts,IsHeavyIon);
  fQA->SetClusterCutsEMC(analysisClusterCutsEMC,IsHeavyIon);
  fQA->SetClusterCutsEMCTrackMatching(analysisClusterCutsEMCTrackMatching,IsHeavyIon);
  fQA->SetClusterCutsIsolationEMC(analysisClusterCutsIsolationEMC,IsHeavyIon);
  fQA->SetClusterCutsTaggingEMC(analysisClusterCutsTaggingEMC,IsHeavyIon);
  fQA->SetClusterCutsPHOS(analysisClusterCutsPHOS,IsHeavyIon);
  fQA->SetConvCuts(analysisConvCuts,IsHeavyIon);
  fQA->SetDoTrackIso(kTRUE);
  fQA->SetDoBackgroundTrackMatching(backgroundTrackMatching);
  fQA->SetDoTrackIso(doChargedIso);
  fQA->SetDoNeutralIso(doNeutralIso);
  fQA->SetDoTagging(doTagging);
  fQA->SetDoCellIso(doCellIso);
  fQA->SetTrackIsoR(trackIsoR);
  fQA->SetNeutralIsoR(neutralIsoR);
  fQA->SetTrackIsoE(trackIsoE);
  fQA->SetNeutralIsoE(neutralIsoE);
  fQA->SetCorrectionTaskSetting(corrTaskSetting);
  fQA->SetSaveConversions(storeConversions);
  fQA->SetSaveEMCClusters(storeEMCalCluster);
  fQA->SetSavePHOSClusters(storePHOSCluster);
  fQA->SetSaveTracks(storeTracks);
  fQA->SetBuffSize(60*1024*1024);
  fQA->SetTrackMatcherRunningMode(trackMatcherRunningMode);
  fQA->SetDoOwnTrackMatching(doOwnTrackMatching);
  fQA->SetUseHistograms(useHistograms);
  fQA->SetMinClsTPC(fMinClsTPC);
  fQA->SetMinClsITS(fMinClsITS);
  fQA->SetChi2PerClsTPC(fChi2PerClsTPC);
  fQA->SetEtaCut(fEtaCut);
  fQA->SetMinPtCut(fPtCut);

  fQA->SetSignalMinM02(minSignalM02);
  fQA->SetSignalMaxM02(maxSignalM02);
  
  mgr->AddTask(fQA);

  mgr->ConnectInput(fQA, 0,  cinput );
  AliAnalysisDataContainer *coutput = NULL;
  AliAnalysisDataContainer *histos = NULL;

  if(corrTaskSetting.CompareTo("")){
    coutput =mgr->CreateContainer( Form("GammaIsoTree_%d_%s",trainConfig,corrTaskSetting.Data()),
                                                              TTree::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              Form("GammaIsoTree_%d.root",trainConfig));
    histos = mgr->CreateContainer( Form("GammaIsoTree_histos_%d_%s",trainConfig,corrTaskSetting.Data()),
                                                              TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              Form("GammaIsoTree_histos_%d.root",trainConfig));
  } else{
    coutput =mgr->CreateContainer( Form("GammaIsoTree_%d",trainConfig),
                                                              TTree::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              Form("GammaIsoTree_%d.root",trainConfig));
    histos = mgr->CreateContainer( Form("GammaIsoTree_histos_%d",trainConfig),
                                                              TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              Form("GammaIsoTree_histos_%d.root",trainConfig));
   
  }
  
  mgr->ConnectOutput (fQA, 1, histos );
  mgr->ConnectOutput (fQA, 2, coutput );

  // mgr->ConnectOutput (fQA, 2, "GammaIsoTreeHistos", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TreeAnalysisHistos", mgr->GetCommonFileName())) );
  return;
}
