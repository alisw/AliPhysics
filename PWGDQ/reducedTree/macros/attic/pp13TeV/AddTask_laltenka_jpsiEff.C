//===============================================================================================================
// addtask for J/psi efficiency
//===============================================================================================================

// MC signals, in order of appearance in AddTask_laltenka_dst_jpsiEff.C
enum MCFilters {
  kJpsiInclusive=0,
  kJpsiNonPrompt,
  kJpsiPrompt,
  kJpsiDecayElectron,
  kJpsiNonPromptDecayElectron,
  kJpsiPromptDecayElectron
};

// trigger modes
enum triggerModes {
  kMB=0,
  kV0HM,
  kEG1DG1,
  kEG2DG2,
  kNTriggerModes=4
};

// settings flags
Bool_t gUseEG1DG1PidMC = kFALSE;
Bool_t gUseEG2DG2PidMC = kFALSE;

// forward declarations
void SetEventCuts(AliReducedAnalysisJpsi2ee* task);
void SetClusterCuts(AliReducedAnalysisJpsi2ee* task);
void SetJpsiElectronTrackCuts(AliReducedAnalysisJpsi2ee* task);
void SetPairCuts(AliReducedAnalysisJpsi2ee* task);
void SetClusterTrackMatcher(AliReducedAnalysisJpsi2ee* task);
void SetMCSignalCutsLeg(AliReducedAnalysisJpsi2ee* task);
void SetMCSignalCutsJPsi(AliReducedAnalysisJpsi2ee* task);
void DefineHistograms(AliReducedAnalysisJpsi2ee* task);

//_______________________________________________________________________________________________________________
AliAnalysisTask* AddTask_laltenka_jpsiEff(Bool_t isAliRoot=kTRUE, Int_t runMode=1, TString prodString="LHC16l") {

  // settings according to production
  //-----------------------------------------------------------------------------------------------------------
  TObjArray*  prodStringArr = prodString.Tokenize("_");
  TString     prod          = prodStringArr->At(0)->GetName();
  
  // MC choice
  Bool_t setMC = kFALSE;
  if (prod.Contains("LHC17h2a") ||
      prod.Contains("LHC17h2b") ||
      prod.Contains("LHC17h2d") ||
      prod.Contains("LHC17h2e") ||
      prod.Contains("LHC17h2f") ||
      prod.Contains("LHC17h2g") ||
      prod.Contains("LHC17h2h2") ||
      prod.Contains("LHC17h2i2") ||
      prod.Contains("LHC17h2j") ||
      prod.Contains("LHC17h2k"))  setMC = kTRUE; // 2016 J/psi inj. MC
  if (prod.Contains("LHC18b1a"))  setMC = kTRUE; // 2017 J/psi inj. MC
  if (prod.Contains("LHC19c6"))   setMC = kTRUE; // 2018 J/psi inj. MC
  if (!setMC) {
    printf("ERROR: AddTask_laltenka_jpsiEff() only defined to run on MC!");
    return 0;
  }

  // J/psi PID choice
  TString triggerString = "";
  if (prodStringArr->GetEntries()>1) {
    triggerString = prodStringArr->At(1)->GetName();
    if (!triggerString.CompareTo("EG1DG1")) gUseEG1DG1PidMC = kTRUE;
    if (!triggerString.CompareTo("EG2DG2")) gUseEG2DG2PidMC = kTRUE;
  }

  // set up analysis task
  //-----------------------------------------------------------------------------------------------------------
  AliReducedAnalysisJpsi2ee* jpsi2eeAnalysis = new AliReducedAnalysisJpsi2ee("Jpsi2eeAnalysis","Jpsi->ee analysis");
  jpsi2eeAnalysis->SetRunOverMC(setMC);
  jpsi2eeAnalysis->SetFillCaloClusterHistograms(kTRUE);
  jpsi2eeAnalysis->SetRunEventMixing(kFALSE);
  jpsi2eeAnalysis->SetRunLikeSignPairing(kFALSE);
  jpsi2eeAnalysis->SetRunPrefilter(kFALSE);

  // set MC pT weights
  //-----------------------------------------------------------------------------------------------------------
  if (setMC) {
    TFile*  fW      = TFile::Open("/Users/lal035/Documents/PWGDQ/data/calibration/MC/ratio_13TeV_two.root");
    TH1F*   weights = (TH1F*)fW->Get("weightsNew");
    if (weights) printf("INFO on AddTask_laltenka_jpsiEff(): Setting MC pT weights!\n");
    jpsi2eeAnalysis->SetMCJpsiPtWeights(weights);
  }

  // set analysis cuts
  //-----------------------------------------------------------------------------------------------------------
  SetEventCuts(jpsi2eeAnalysis);
  SetClusterCuts(jpsi2eeAnalysis);
  SetJpsiElectronTrackCuts(jpsi2eeAnalysis);
  SetPairCuts(jpsi2eeAnalysis);
  SetClusterTrackMatcher(jpsi2eeAnalysis);

  // set MC signal selection
  //-----------------------------------------------------------------------------------------------------------
  if (setMC) {
    SetMCSignalCutsLeg(jpsi2eeAnalysis);
    SetMCSignalCutsJPsi(jpsi2eeAnalysis);
  }

  // initialize analysis task
  //-----------------------------------------------------------------------------------------------------------
  jpsi2eeAnalysis->Init();

  // define histograms histograms
  //-----------------------------------------------------------------------------------------------------------
  AliReducedVarManager::SetDefaultVarNames();
  DefineHistograms(jpsi2eeAnalysis);
  AliReducedVarManager::SetUseVars((Bool_t*)jpsi2eeAnalysis->GetHistogramManager()->GetUsedVars());

  // initialize wrapper AliAnalysisTask
  // (in order to run AliReducedAnalysisJpsi2ee in an aliroot analysis train )
  //-----------------------------------------------------------------------------------------------------------
  AliAnalysisTaskReducedEventProcessor* task = new AliAnalysisTaskReducedEventProcessor("ReducedEventAnalysisManager", runMode);
  task->AddTask(jpsi2eeAnalysis);

  // intercept isAliRoot=kFALSE (nothing to be done yet)
  //-----------------------------------------------------------------------------------------------------------
  if (!isAliRoot) return 0;

  // get analysis manager
  //-----------------------------------------------------------------------------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_laltenka_jpsiEff", "No analysis manager found.");
    return 0;
  }

  // get data container
  //-----------------------------------------------------------------------------------------------------------
  AliAnalysisDataContainer* cReducedEvent = NULL;
  if (runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) {
    printf("INFO on AddTask_laltenka_jpsiEff(): use on the fly events\n");
    cReducedEvent = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ReducedEventDQ");
    if (!cReducedEvent) {
      printf("ERROR: In AddTask_laltenka_jpsiEff(), ReducedEvent exchange container could not be found!\n");
      return 0x0;
    }
  }

  // add task to analysis manager
  //-----------------------------------------------------------------------------------------------------------
  mgr->AddTask(task);

  // connect input data containers
  //-----------------------------------------------------------------------------------------------------------
  if (runMode==AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree)              mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  else if (runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents)  mgr->ConnectInput(task, 0, cReducedEvent);
  else {
    printf("ERROR: In AddTask_laltenka_jpsiEff(), runMode %d not defined!\n", runMode);
    return 0;
  }

  // connect output data containers
  //-----------------------------------------------------------------------------------------------------------
  TString                           outputName = "dstCorrelationsAnalysisHistograms.root";
  if (triggerString.CompareTo(""))  outputName = Form("dstCorrelationsAnalysisHistograms_%s.root", triggerString.Data());
  if (setMC)                        outputName.ReplaceAll(".root", "_MC.root");
  printf("INFO on AddTask_laltenka_jpsiEff(): output container name = %s\n", outputName.Data());
  AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("jpsi2eeHistos", THashList::Class(), AliAnalysisManager::kOutputContainer, outputName);
  mgr->ConnectOutput(task, 1, cOutputHist);
  
  // done
  //-----------------------------------------------------------------------------------------------------------
  return task;
}

//_______________________________________________________________________________________________________________
void SetEventCuts(AliReducedAnalysisJpsi2ee* task) {
  
  // default selection (vertex z and PS)
  AliReducedEventCut* defaultCut = new AliReducedEventCut("DefaultCut", "Default selection");
  defaultCut->AddEventTriggerFilter(AliReducedVarManager::kINT7);
  defaultCut->AddCut(AliReducedVarManager::kVtxZ, -10.0, 10.0);
  if (!task->GetRunOverMC()) defaultCut->AddCut(AliReducedVarManager::kIsPhysicsSelection, 0.1, 2.);

  // SPD pileup cut (https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsPileup)
  // event->IsPileupFromSPD(3,0.8,3.,2.,5.)
  AliReducedEventCut* spdPileupCut = new AliReducedEventCut("SpdPileupCut", "SPD pileup rejection");
  spdPileupCut->AddEventTagFilterBit(9, kTRUE);

  // MV pileup cut (https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsPileup)
  // reject multi-vertexer pileup events, min. weight. dist. 15
  AliReducedEventCut* mvPileupCut = new AliReducedEventCut("MvPileupCut", "MV pileup rejection");
  mvPileupCut->AddEventTagFilterBit(1, kTRUE);

  // combine all cuts
  AliReducedCompositeCut* eventCut = new AliReducedCompositeCut("eventCut", "", kTRUE);
  eventCut->AddCut(defaultCut);
  eventCut->AddCut(spdPileupCut);
  eventCut->AddCut(mvPileupCut);
  task->AddEventCut(eventCut);
}

//_______________________________________________________________________________________________________________
void SetClusterCuts(AliReducedAnalysisJpsi2ee* task) {
  AliReducedCaloClusterCut* clusterCut = new AliReducedCaloClusterCut("clusters","cluster selection");
  clusterCut->AddCut(AliReducedVarManager::kEMCALdetector, AliReducedCaloClusterInfo::kEMCAL, AliReducedCaloClusterInfo::kEMCAL);
  clusterCut->AddCut(AliReducedVarManager::kEMCALm02, 0.01, 0.70);
  task->AddClusterCut(clusterCut);
}

//_______________________________________________________________________________________________________________
void SetJpsiElectronTrackCuts(AliReducedAnalysisJpsi2ee* task) {

  // default cuts to be used in either cut variation
  //_____________________________________________________________________________________________________________
  AliReducedTrackCut* defaultCut = new AliReducedTrackCut("defaultCut", "Electron default selection");
  defaultCut->SetTrackFilterBit(0); // fQualityFlag bit 32 (see AddTask_laltenka_dst_jpsiEff.C)
  defaultCut->AddCut(AliReducedVarManager::kPt,                 1.0,   1.0e+30);
  defaultCut->AddCut(AliReducedVarManager::kEta,               -0.9,   0.9);
  defaultCut->AddCut(AliReducedVarManager::kDcaXY,             -1.0,   1.0);
  defaultCut->AddCut(AliReducedVarManager::kDcaZ,              -3.0,   3.0);
  defaultCut->AddCut(AliReducedVarManager::kTPCncls,           70.0, 161.0);
  defaultCut->AddCut(AliReducedVarManager::kTPCchi2,            0.0,   4.0);
  defaultCut->AddCut(AliReducedVarManager::kITSchi2,            0.0,  30.0);
  defaultCut->AddCut(AliReducedVarManager::kTPCNclusBitsFired,  6.0,   9.0);
  defaultCut->SetRejectKinks();
  defaultCut->SetRequestITSrefit();
  defaultCut->SetRequestTPCrefit();
  defaultCut->SetRequestSPDany();

  // standard electron cut
  //_____________________________________________________________________________________________________________
  AliReducedTrackCut* standardCut = (AliReducedTrackCut*)defaultCut->Clone("standardElectron");
  standardCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron), -1.5, 3.0);
  if (gUseEG1DG1PidMC || gUseEG2DG2PidMC) {
    // EMCal/DCal trigger -> good EMCal runs
    standardCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton),  3.5, 1.e8);
    standardCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion),    3.5, 1.e8, kFALSE,
                        AliReducedVarManager::kEMCALmatchedEnergy, -9999., -9999.);
    standardCut->AddCut(AliReducedVarManager::kEMCALmatchedEOverP,                                                                   0.8,  1.3, kFALSE,
                        AliReducedVarManager::kEMCALmatchedEnergy, -9999., -9999., kTRUE);
  } else {
    // MB/HM trigger
    standardCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton),  3.5, 1.e8);
    standardCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion),    3.5, 1.e8);
  }
  task->AddTrackCut(standardCut);
    
  // PID strict electron cut
  //_____________________________________________________________________________________________________________
  AliReducedTrackCut* PIDstrictCut = (AliReducedTrackCut*)defaultCut->Clone("PIDstrictElectron");
  PIDstrictCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron), -1.0, 2.5);
  if (gUseEG1DG1PidMC || gUseEG2DG2PidMC) {
    // EMCal/DCal trigger -> good EMCal runs
    PIDstrictCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton),  4.0, 1.e8);
    PIDstrictCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion),    4.0, 1.e8, kFALSE,
                         AliReducedVarManager::kEMCALmatchedEnergy, -9999., -9999.);
    PIDstrictCut->AddCut(AliReducedVarManager::kEMCALmatchedEOverP,                                                                   0.9, 1.2, kFALSE,
                         AliReducedVarManager::kEMCALmatchedEnergy, -9999., -9999., kTRUE);
  } else {
    // MB/HM trigger
    PIDstrictCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton),  4.0, 1.e8);
    PIDstrictCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion),    4.0, 1.e8);
  }
  task->AddTrackCut(PIDstrictCut);
    
  // PID open electron cut
  //_____________________________________________________________________________________________________________
  AliReducedTrackCut* PIDopenCut = (AliReducedTrackCut*)defaultCut->Clone("PIDopenElectron");
  PIDopenCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron), -2.0, 3.5);
  if (gUseEG1DG1PidMC || gUseEG2DG2PidMC) {
    // EMCal/DCal trigger -> good EMCal runs
    PIDopenCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton),  3.0, 1.e8);
    PIDopenCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion),    3.0, 1.e8, kFALSE,
                       AliReducedVarManager::kEMCALmatchedEnergy, -9999., -9999.);
    PIDopenCut->AddCut(AliReducedVarManager::kEMCALmatchedEOverP,                                                                   0.7, 1.4, kFALSE,
                       AliReducedVarManager::kEMCALmatchedEnergy, -9999., -9999., kTRUE);
  } else {
    // MB/HM trigger
    PIDopenCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton),  3.0, 1.e8);
    PIDopenCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion),    3.0, 1.e8);
  }
  task->AddTrackCut(PIDopenCut);
}

//_______________________________________________________________________________________________________________
void SetPairCuts(AliReducedAnalysisJpsi2ee* task) {
  
  // default pair cut
  AliReducedTrackCut* defaultPairCut = new AliReducedTrackCut("defaultPairCut", "Pair default selection");
  defaultPairCut->AddCut(AliReducedVarManager::kRap,  -0.9, 0.9);

  // standard pair cut
  if (gUseEG1DG1PidMC) {
    AliReducedTrackCut* standardPairCutLeg1 = (AliReducedTrackCut*)defaultPairCut->Clone("standardPairCutLeg1");
    AliReducedTrackCut* standardPairCutLeg2 = (AliReducedTrackCut*)defaultPairCut->Clone("standardPairCutLeg2");

    standardPairCutLeg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0), 10.0, 1.0e+30);
    standardPairCutLeg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1), 10.0, 1.0e+30);

    AliReducedCompositeCut* standardPairCut = new AliReducedCompositeCut("standardPair", "", kFALSE);
    standardPairCut->AddCut(standardPairCutLeg1);
    standardPairCut->AddCut(standardPairCutLeg2);
    task->AddPairCut(standardPairCut);
  } else if (gUseEG2DG2PidMC) {
    AliReducedTrackCut* standardPairCutLeg1 = (AliReducedTrackCut*)defaultPairCut->Clone("standardPairCutLeg1");
    AliReducedTrackCut* standardPairCutLeg2 = (AliReducedTrackCut*)defaultPairCut->Clone("standardPairCutLeg2");

    standardPairCutLeg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0),  5.0, 1.0e+30);
    standardPairCutLeg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1),  5.0, 1.0e+30);

    AliReducedCompositeCut* standardPairCut = new AliReducedCompositeCut("standardPair", "", kFALSE);
    standardPairCut->AddCut(standardPairCutLeg1);
    standardPairCut->AddCut(standardPairCutLeg2);
    task->AddPairCut(standardPairCut);
  } else {
    AliReducedTrackCut* standardPairCut = (AliReducedTrackCut*)defaultPairCut->Clone("standardPair");
    task->AddPairCut(standardPairCut);
  }
  
  // non-prompt pair cuts
  TFile* cutFile = TFile::Open("/Users/lal035/Documents/PWGDQ/analysis_MC_decayLengthCut/decayLengthResolution_MB_MC_standardElectron_paper.root");

  // default cut with decay-length cut applied - 5% prompt contamination
  AliReducedTrackCut* defaultPairCut5perc = (AliReducedTrackCut*)defaultPairCut->Clone("defaultPairCut5perc");
  defaultPairCut5perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_SS_5percent_fit"), 1.0e+30, kFALSE,
                              AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                              AliReducedVarManager::kPairTypeSPD,  0.,   0.,       kFALSE);
  defaultPairCut5perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_FS_5percent_fit"), 1.0e+30, kFALSE,
                              AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                              AliReducedVarManager::kPairTypeSPD,  1.,   1.,       kFALSE);
  defaultPairCut5perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_FF_5percent_fit"), 1.0e+30, kFALSE,
                              AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                              AliReducedVarManager::kPairTypeSPD,  2.,   2.,       kFALSE);
    
  // default cut with decay-length cut applied - 10% prompt contamination
  AliReducedTrackCut* defaultPairCut10perc = (AliReducedTrackCut*)defaultPairCut->Clone("defaultPairCut10perc");
  defaultPairCut10perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_SS_10percent_fit"), 1.0e+30, kFALSE,
                               AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                               AliReducedVarManager::kPairTypeSPD,  0.,   0.,       kFALSE);
  defaultPairCut10perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_FS_10percent_fit"), 1.0e+30, kFALSE,
                               AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                               AliReducedVarManager::kPairTypeSPD,  1.,   1.,       kFALSE);
  defaultPairCut10perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_FF_10percent_fit"), 1.0e+30, kFALSE,
                               AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                               AliReducedVarManager::kPairTypeSPD,  2.,   2.,       kFALSE);

  // default cut with decay-length cut applied - 20% prompt contamination
  AliReducedTrackCut* defaultPairCut20perc = (AliReducedTrackCut*)defaultPairCut->Clone("defaultPairCut20perc");
  defaultPairCut20perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_SS_20percent_fit"), 1.0e+30, kFALSE,
                               AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                               AliReducedVarManager::kPairTypeSPD,  0.,   0.,       kFALSE);
  defaultPairCut20perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_FS_20percent_fit"), 1.0e+30, kFALSE,
                               AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                               AliReducedVarManager::kPairTypeSPD,  1.,   1.,       kFALSE);
  defaultPairCut20perc->AddCut(AliReducedVarManager::kPseudoProperDecayTime, (TF1*)gDirectory->Get("decayLength_cut_FF_20percent_fit"), 1.0e+30, kFALSE,
                               AliReducedVarManager::kPt,           0.0,  1.0e+30,  kFALSE,
                               AliReducedVarManager::kPairTypeSPD,  2.,   2.,       kFALSE);
    
  // set pair cuts
  if (gUseEG1DG1PidMC) {
    // non-prompt, 5% prompt contamination
    AliReducedTrackCut* pairCutNonPromptLeg1 = (AliReducedTrackCut*)defaultPairCut5perc->Clone("pairCutNonPromptLeg1");
    AliReducedTrackCut* pairCutNonPromptLeg2 = (AliReducedTrackCut*)defaultPairCut5perc->Clone("pairCutNonPromptLeg2");

    pairCutNonPromptLeg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0), 10.0, 1.0e+30);
    pairCutNonPromptLeg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1), 10.0, 1.0e+30);

    AliReducedCompositeCut* pairCutNonPrompt = new AliReducedCompositeCut("nonPromptPair5perc", "", kFALSE);
    pairCutNonPrompt->AddCut(pairCutNonPromptLeg1);
    pairCutNonPrompt->AddCut(pairCutNonPromptLeg2);
    task->AddPairCut(pairCutNonPrompt);

    // non-prompt, 10% prompt contamination
    AliReducedTrackCut* pairCutNonPrompt1Leg1 = (AliReducedTrackCut*)defaultPairCut10perc->Clone("pairCutNonPrompt1Leg1");
    AliReducedTrackCut* pairCutNonPrompt1Leg2 = (AliReducedTrackCut*)defaultPairCut10perc->Clone("pairCutNonPrompt1Leg2");

    pairCutNonPrompt1Leg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0), 10.0, 1.0e+30);
    pairCutNonPrompt1Leg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1), 10.0, 1.0e+30);

    AliReducedCompositeCut* pairCutNonPrompt1 = new AliReducedCompositeCut("nonPromptPair10perc", "", kFALSE);
    pairCutNonPrompt1->AddCut(pairCutNonPrompt1Leg1);
    pairCutNonPrompt1->AddCut(pairCutNonPrompt1Leg2);
    task->AddPairCut(pairCutNonPrompt1);

    // non-prompt, 20% prompt contamination
    AliReducedTrackCut* pairCutNonPrompt2Leg1 = (AliReducedTrackCut*)defaultPairCut20perc->Clone("pairCutNonPrompt2Leg1");
    AliReducedTrackCut* pairCutNonPrompt2Leg2 = (AliReducedTrackCut*)defaultPairCut20perc->Clone("pairCutNonPrompt2Leg2");

    pairCutNonPrompt2Leg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0), 10.0, 1.0e+30);
    pairCutNonPrompt2Leg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1), 10.0, 1.0e+30);

    AliReducedCompositeCut* pairCutNonPrompt2 = new AliReducedCompositeCut("nonPromptPair20perc", "", kFALSE);
    pairCutNonPrompt2->AddCut(pairCutNonPrompt2Leg1);
    pairCutNonPrompt2->AddCut(pairCutNonPrompt2Leg2);
    task->AddPairCut(pairCutNonPrompt2);
  } else if (gUseEG2DG2PidMC) {
    // non-prompt, 5% prompt contamination
    AliReducedTrackCut* pairCutNonPromptLeg1 = (AliReducedTrackCut*)defaultPairCut5perc->Clone("pairCutNonPromptLeg1");
    AliReducedTrackCut* pairCutNonPromptLeg2 = (AliReducedTrackCut*)defaultPairCut5perc->Clone("pairCutNonPromptLeg2");

    pairCutNonPromptLeg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0),  5.0, 1.0e+30);
    pairCutNonPromptLeg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1),  5.0, 1.0e+30);

    AliReducedCompositeCut* pairCutNonPrompt = new AliReducedCompositeCut("nonPromptPair5perc", "", kFALSE);
    pairCutNonPrompt->AddCut(pairCutNonPromptLeg1);
    pairCutNonPrompt->AddCut(pairCutNonPromptLeg2);
    task->AddPairCut(pairCutNonPrompt);

    // non-prompt, 10% prompt contamination
    AliReducedTrackCut* pairCutNonPrompt1Leg1 = (AliReducedTrackCut*)defaultPairCut10perc->Clone("pairCutNonPrompt1Leg1");
    AliReducedTrackCut* pairCutNonPrompt1Leg2 = (AliReducedTrackCut*)defaultPairCut10perc->Clone("pairCutNonPrompt1Leg2");

    pairCutNonPrompt1Leg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0),  5.0, 1.0e+30);
    pairCutNonPrompt1Leg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1),  5.0, 1.0e+30);

    AliReducedCompositeCut* pairCutNonPrompt1 = new AliReducedCompositeCut("nonPromptPair10perc", "", kFALSE);
    pairCutNonPrompt1->AddCut(pairCutNonPrompt1Leg1);
    pairCutNonPrompt1->AddCut(pairCutNonPrompt1Leg2);
    task->AddPairCut(pairCutNonPrompt1);

    // non-prompt, 20% prompt contamination
    AliReducedTrackCut* pairCutNonPrompt2Leg1 = (AliReducedTrackCut*)defaultPairCut20perc->Clone("pairCutNonPrompt2Leg1");
    AliReducedTrackCut* pairCutNonPrompt2Leg2 = (AliReducedTrackCut*)defaultPairCut20perc->Clone("pairCutNonPrompt2Leg2");

    pairCutNonPrompt2Leg1->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+0),  5.0, 1.0e+30);
    pairCutNonPrompt2Leg2->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kPairLegEMCALmatchedEnergy+1),  5.0, 1.0e+30);

    AliReducedCompositeCut* pairCutNonPrompt2 = new AliReducedCompositeCut("nonPromptPair20perc", "", kFALSE);
    pairCutNonPrompt2->AddCut(pairCutNonPrompt2Leg1);
    pairCutNonPrompt2->AddCut(pairCutNonPrompt2Leg2);
    task->AddPairCut(pairCutNonPrompt2);
  } else {
    // non-prompt, 5% prompt contamination
    AliReducedTrackCut* pairCutNonPrompt = (AliReducedTrackCut*)defaultPairCut5perc->Clone("nonPromptPair5perc");
    task->AddPairCut(pairCutNonPrompt);

    // non-prompt, 10% prompt contamination
    AliReducedTrackCut* pairCutNonPrompt1 = (AliReducedTrackCut*)defaultPairCut10perc->Clone("nonPromptPair10perc");
    task->AddPairCut(pairCutNonPrompt1);

    // non-prompt, 20% prompt contamination
    AliReducedTrackCut* pairCutNonPrompt2 = (AliReducedTrackCut*)defaultPairCut20perc->Clone("nonPromptPair20perc");
    task->AddPairCut(pairCutNonPrompt2);
  }
}

//_______________________________________________________________________________________________________________
void SetClusterTrackMatcher(AliReducedAnalysisJpsi2ee* task) {
  AliReducedCaloClusterTrackMatcher* matcher = new AliReducedCaloClusterTrackMatcher();
  matcher->SetMaximumMatchingDistance(0.005);
  task->SetClusterTrackMatcher(matcher);
}

//_______________________________________________________________________________________________________________
void SetMCSignalCutsLeg(AliReducedAnalysisJpsi2ee* task) {
  // electron from inclusive J/psi
  AliReducedTrackCut* trueElectron = new AliReducedTrackCut("TrueElectron", "reconstructed electrons from inclusive Jpsi with MC truth");
  trueElectron->SetMCFilterBit(kJpsiDecayElectron);
  task->AddLegCandidateMCcut(trueElectron);
  
  // electron from non-prompt J/psi
  AliReducedTrackCut* trueElectronNonPrompt = new AliReducedTrackCut("TrueElectronNonPrompt", "reconstructed electrons from non-prompt Jpsi with MC truth");
  trueElectronNonPrompt->SetMCFilterBit(kJpsiNonPromptDecayElectron);
  task->AddLegCandidateMCcut(trueElectronNonPrompt);
  
  // electron from prompt J/psi
  AliReducedTrackCut* trueElectronPrompt = new AliReducedTrackCut("TrueElectronPrompt", "reconstructed electrons from prompt Jpsi with MC truth");
  trueElectronPrompt->SetMCFilterBit(kJpsiPromptDecayElectron);
  task->AddLegCandidateMCcut(trueElectronPrompt);
}

//_______________________________________________________________________________________________________________
void SetMCSignalCutsJPsi(AliReducedAnalysisJpsi2ee* task) {
  // inclusive J/psi
  AliReducedTrackCut* mcTruthJpsi = new AliReducedTrackCut("mcTruthJpsi", "Pure MC truth inclusive J/psi");
  mcTruthJpsi->SetMCFilterBit(kJpsiInclusive);
  mcTruthJpsi->SetApplyReweightMCpt(kFALSE); //reweight is applied only for the prompt (kTRUE will reweight twice)
  mcTruthJpsi->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);

  AliReducedTrackCut* mcTruthJpsiElectron = new AliReducedTrackCut("mcTruthJpsiElectron", "Pure MC truth electron from inclusive J/psi");
  mcTruthJpsiElectron->SetMCFilterBit(kJpsiDecayElectron);
  mcTruthJpsiElectron->AddCut(AliReducedVarManager::kEta,   -0.9, 0.9);
  mcTruthJpsiElectron->AddCut(AliReducedVarManager::kPt,    1.0,  1.0e+30);
  task->AddJpsiMotherMCCut(mcTruthJpsi, mcTruthJpsiElectron);
  
  // non-prompt J/psi
  AliReducedTrackCut* mcTruthJpsiNonPrompt = new AliReducedTrackCut("mcTruthJpsiNonPrompt", "Pure MC truth non-prompt J/psi");
  mcTruthJpsiNonPrompt->SetMCFilterBit(kJpsiNonPrompt);
  mcTruthJpsiNonPrompt->SetApplyReweightMCpt(kFALSE); //no reweight
  mcTruthJpsiNonPrompt->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  
  AliReducedTrackCut* mcTruthJpsiElectronNonPrompt = new AliReducedTrackCut("mcTruthJpsiElectronNonPrompt", "Pure MC truth electron from non-prompt J/psi");
  mcTruthJpsiElectronNonPrompt->SetMCFilterBit(kJpsiNonPromptDecayElectron);
  mcTruthJpsiElectronNonPrompt->AddCut(AliReducedVarManager::kEta,  -0.9, 0.9);
  mcTruthJpsiElectronNonPrompt->AddCut(AliReducedVarManager::kPt,   1.0,  1.0e+30);
  task->AddJpsiMotherMCCut(mcTruthJpsiNonPrompt, mcTruthJpsiElectronNonPrompt);

  // prompt J/psi
  AliReducedTrackCut* mcTruthJpsiPrompt = new AliReducedTrackCut("mcTruthJpsiPrompt", "Pure MC truth prompt J/psi");
  mcTruthJpsiPrompt->SetMCFilterBit(kJpsiPrompt);
  mcTruthJpsiPrompt->SetApplyReweightMCpt(kTRUE); //Apply reweighting only for prompt J/psi
  mcTruthJpsiPrompt->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  
  AliReducedTrackCut* mcTruthJpsiElectronPrompt = new AliReducedTrackCut("mcTruthJpsiElectronPrompt", "Pure MC truth electron from prompt J/psi");
  mcTruthJpsiElectronPrompt->SetMCFilterBit(kJpsiPromptDecayElectron);
  mcTruthJpsiElectronPrompt->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  mcTruthJpsiElectronPrompt->AddCut(AliReducedVarManager::kPt,  1.0,  1.0e+30);
  task->AddJpsiMotherMCCut(mcTruthJpsiPrompt, mcTruthJpsiElectronPrompt);
}

//_______________________________________________________________________________________________________________
void DefineHistograms(AliReducedAnalysisJpsi2ee* task) {
  AliHistogramManager* man = task->GetHistogramManager();

  //-------------------------------------------------------------------------------------------------------------
  // binning
  //-------------------------------------------------------------------------------------------------------------

  // mass binning
  const Int_t kNMassBins = 71;      // [1.4, 4.2]
  Double_t    massBins[kNMassBins];
  for (Int_t i=0; i<kNMassBins; i++) massBins[i] = 1.4+i*0.04;

  // pt binning
  const Int_t kNPtBins          = 11;             // [0., 40.]
  Double_t    ptBins[kNPtBins]  = {0.0, 1.0, 3.0, 5.0, 7.0, 8.0, 12.0, 15.0, 20.0, 30.0, 40.0};
  
  const Int_t kNPtBinsEff = 201;                  // [0., 40.]
  Double_t    ptBinsEff[kNPtBinsEff];
  for (Int_t i=0; i<kNPtBinsEff; i++) ptBinsEff[i] = 0. + i*0.2;

  // phi binning
  const Int_t kNPhiBins = 41;                     // [0., 2*pi]
  Double_t    phiBins[kNPhiBins];
  for (Int_t i=0; i<kNPhiBins; i++) phiBins[i] = i*TMath::Pi()/20;

  // eta binning
  const Int_t kNEtaBins = 31;                     // [-1.5, 1.5]
  Double_t    etaBins[kNEtaBins];
  for (Int_t i=0; i<kNEtaBins; i++) etaBins[i] = -1.5+i*0.1;
  
  // SPD pair type
  const Int_t kNPairTypeSPDBins                   = 4;
  Double_t    pairTypeSPDBins[kNPairTypeSPDBins]  = {-0.5, 0.5, 1.5, 2.5};

  //-------------------------------------------------------------------------------------------------------------
  // histogram classes
  //-------------------------------------------------------------------------------------------------------------
  
  // event quantities
  TString histClasses = "";
  histClasses += "Event_BeforeCuts;";
  histClasses += "Event_AfterCuts;";
  
  // calorimeter quantities
  if (task->GetFillCaloClusterHistograms()) {
    histClasses += "CaloCluster_BeforeCuts;";
    for (Int_t i=0; i<task->GetNClusterCuts(); ++i) {
      TString cutName = task->GetClusterCutName(i);
      histClasses += Form("CaloCluster_%s;", cutName.Data());
    }
  }

  // electron track quantities
  if (task->GetLoopOverTracks()) {
    histClasses += "Track_BeforeCuts;";
    for (Int_t i=0; i<task->GetNTrackCuts(); ++i) {
      TString cutName = task->GetTrackCutName(i);
      histClasses += Form("Track_%s;", cutName.Data());
      if (task->GetRunOverMC()) {
        for (Int_t mcSel=0; mcSel<task->GetNLegCandidateMCcuts(); ++mcSel) {
          histClasses += Form("Track_%s_%s;", cutName.Data(), task->GetLegCandidateMCcutName(mcSel));
        }
      }
    }
  }
  
  // J/psi MC gen. tracks
  if (task->GetRunOverMC()) {
    for (Int_t mcSel=0; mcSel<task->GetNJpsiMotherMCCuts(); ++mcSel) {
      histClasses += Form("%s_PureMCTruth_BeforeSelection;", task->GetJpsiMotherMCcutName(mcSel));
      histClasses += Form("%s_PureMCTruth_AfterSelection;", task->GetJpsiMotherMCcutName(mcSel));
    }
  }

  // pair quantities
  if (task->GetNPairCuts()>1) {
    for (Int_t i=0; i<task->GetNTrackCuts(); ++i) {
      for (Int_t j=0; j<task->GetNPairCuts(); ++j) {
        TString trackCutName  = task->GetTrackCutName(i);
        TString pairCutName   = task->GetPairCutName(j);
        if (!task->GetRunLikeSignPairing()) histClasses += Form("PairSEPM_%s_%s;",
                                                                trackCutName.Data(), pairCutName.Data());
        else                                histClasses += Form("PairSEPP_%s_%s;PairSEPM_%s_%s;PairSEMM_%s_%s;",
                                                                trackCutName.Data(), pairCutName.Data(),
                                                                trackCutName.Data(), pairCutName.Data(),
                                                                trackCutName.Data(), pairCutName.Data());
        if (!task->GetRunOverMC() && !pairCutName.Contains("nonPrompt"))  histClasses += Form("PairMEPP_%s_%s;PairMEPM_%s_%s;PairMEMM_%s_%s;",
                                                                                              trackCutName.Data(), pairCutName.Data(),
                                                                                              trackCutName.Data(), pairCutName.Data(),
                                                                                              trackCutName.Data(), pairCutName.Data());
        if (task->GetRunOverMC()) {
          for (Int_t mcSel=0; mcSel<task->GetNLegCandidateMCcuts(); ++mcSel) histClasses += Form("PairSEPM_%s_%s_%s;", trackCutName.Data(), pairCutName.Data(), task->GetLegCandidateMCcutName(mcSel));
        }
      }
    }
  } else {
    for (Int_t i=0; i<task->GetNTrackCuts(); ++i) {
      TString trackCutName  = task->GetTrackCutName(i);
      if (!task->GetRunLikeSignPairing()) histClasses += Form("PairSEPM_%s;", trackCutName.Data());
      else                                histClasses += Form("PairSEPP_%s;PairSEPM_%s;PairSEMM_%s;", trackCutName.Data(), trackCutName.Data(), trackCutName.Data());
      if (!task->GetRunOverMC())          histClasses += Form("PairMEPP_%s;PairMEPM_%s;PairMEMM_%s;", trackCutName.Data(), trackCutName.Data(), trackCutName.Data());
      if (task->GetRunOverMC()) {
        for (Int_t mcSel=0; mcSel<task->GetNLegCandidateMCcuts(); ++mcSel) histClasses += Form("PairSEPM_%s_%s;", trackCutName.Data(), task->GetLegCandidateMCcutName(mcSel));
      }
    }
  }
  
  // add histograms according to class to histogram manager
  TString classesStr(histClasses);
  TObjArray* arr = classesStr.Tokenize(";");

  //-------------------------------------------------------------------------------------------------------------
  // histograms
  //-------------------------------------------------------------------------------------------------------------

  // loop over histogram classes and add histograms
  printf("INFO on AddTask_laltenka_jpsiEff(): Histgram classes included in histogram manager\n");
  for (Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();

    //-----------------------------------------------------------------------------------------------------------
    // event histograms
    //-----------------------------------------------------------------------------------------------------------
    if(classStr.Contains("Event_")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n",classStr.Data());

      // multiplicity
      man->AddHistogram(classStr.Data(), "SPDnTracklets", "SPD #tracklets in |#eta|<1.0", kFALSE,
                        200, 0.0, 200.,                               AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(), "MultiplicityV0", "", kFALSE,
                        500, 0.0, 5000.,                              AliReducedVarManager::kVZEROTotalMult);

      // vertex histograms
      man->AddHistogram(classStr.Data(), "VtxX", "Vtx X", kFALSE,
                        50, -1.0, 1.0,                                AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(), "VtxY", "Vtx Y", kFALSE,
                        50, -1.0, 1.0,                                AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(), "VtxZ", "Vtx Z", kFALSE,
                        150, -15.0, 15.0,                             AliReducedVarManager::kVtxZ);

      continue;
    }
    
    //-----------------------------------------------------------------------------------------------------------
    // calorimeter cluster histograms
    //-----------------------------------------------------------------------------------------------------------
    if(classStr.Contains("CaloCluster_")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n",classStr.Data());
      
      man->AddHistogram(classStr.Data(), "Energy",          "Cluster energy",                 kFALSE, 500,  0.0,  50.0, AliReducedVarManager::kEMCALclusterEnergy);
      man->AddHistogram(classStr.Data(), "M02",             "Cluster M02",                    kFALSE, 100,  0.0,   1.0, AliReducedVarManager::kEMCALm02);
      man->AddHistogram(classStr.Data(), "M20",             "Cluster M20",                    kFALSE, 100,  0.0,   1.0, AliReducedVarManager::kEMCALm20);
      man->AddHistogram(classStr.Data(), "NCells",          "No. cells per cluster",          kFALSE, 100,  0.0, 100.0, AliReducedVarManager::kEMCALnCells);
      man->AddHistogram(classStr.Data(), "NMatchedTracks",  "No. matched tracks per cluster", kFALSE, 100,  0.0, 100.0, AliReducedVarManager::kEMCALnMatchedTracks);
      
      continue;
    }
    
    //-----------------------------------------------------------------------------------------------------------
    // track and associated track histograms
    //-----------------------------------------------------------------------------------------------------------
    if(classStr.Contains("Track_")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n",classStr.Data());
      
      // kinematics
      man->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution", kFALSE,
                        200, 0., 100.,                                AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE,
                        200, -5.0, 5.0,                               AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAxy_Pt", "DCAxy vs. pT", kFALSE,
                        100, -2.0, 2.0,                               AliReducedVarManager::kDcaXY,
                        100, 0.0, 100.,                               AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE,
                        200, -5.0, 5.0,                               AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(), "DCAz_Pt", "DCAz vs. pT", kFALSE,
                        100, -3.0, 3.0,                               AliReducedVarManager::kDcaZ,
                        100, 0.0, 100.,                               AliReducedVarManager::kPt);

      // PID
      man->AddHistogram(classStr.Data(), "TPCnsigElectron_Pin", "TPC N_{#sigma} electron vs. inner param P", kFALSE,
                        200, 0., 100.,                                AliReducedVarManager::kPin,
                        100, -5.0, 5.0,                               AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(), "EMCMatchedEOverP_Pt", "E/P from the calorimeter matched cluster vs. pT", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kPt,
                        150, 0.0, 1.5,                                AliReducedVarManager::kEMCALmatchedEOverP);
      man->AddHistogram(classStr.Data(), "TPCnsigElectron_EMCMatchedEOverP", "nSigma electron from TPC vs. E/p", kFALSE,
                        150, 0.0, 1.5,                                AliReducedVarManager::kEMCALmatchedEOverP,
                        100, -5.0, 5.0,                               AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);

      continue;
    }
    
    //-----------------------------------------------------------------------------------------------------------
    // MC histograms
    //-----------------------------------------------------------------------------------------------------------
    if(classStr.Contains("MCTruth")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n",classStr.Data());
      
      man->AddHistogram(classStr.Data(), "MassMC",      "MC mass",      kFALSE,  70,  1.4,   4.2, AliReducedVarManager::kMassMC);
      man->AddHistogram(classStr.Data(), "PtMC",        "p_{T} MC",     kFALSE, 1000,  0.0, 50.0, AliReducedVarManager::kPtMC);
      man->AddHistogram(classStr.Data(), "RapidityMC",  "MC rapidity",  kFALSE,  80, -2.0,   2.0, AliReducedVarManager::kRapMC);
      man->AddHistogram(classStr.Data(), "PhiMC",       "MC #varphi",   kFALSE, 100,  0.0,   6.3, AliReducedVarManager::kPhiMC);
      man->AddHistogram(classStr.Data(), "EtaMC",       "MC #eta",      kFALSE, 100, -1.5,   1.5, AliReducedVarManager::kEtaMC);
      man->AddHistogram(classStr.Data(), "PtMC_RapidityMC", "p_{T} MC vs. rapidity MC", kFALSE,
                        80, -2.0,   2.0, AliReducedVarManager::kRapMC,
                        200,  0.0, 100.0, AliReducedVarManager::kPtMC);
      
      // pair histogram for efficiency evaluation
      const Int_t kNVarsEff = 4;
      Int_t varsEff[kNVarsEff] = {AliReducedVarManager::kMassMC, AliReducedVarManager::kPtMC, AliReducedVarManager::kEtaMC, AliReducedVarManager::kPhiMC};
      TArrayD pairHistEffBinLimits[kNVarsEff];
      pairHistEffBinLimits[0] = TArrayD(kNMassBins,   massBins);
      pairHistEffBinLimits[1] = TArrayD(kNPtBinsEff,  ptBinsEff);
      pairHistEffBinLimits[2] = TArrayD(kNEtaBins,    etaBins);
      pairHistEffBinLimits[3] = TArrayD(kNPhiBins,    phiBins);
      man->AddHistogram(classStr.Data(), "PairInvMass_Eff", "Differential pair inv. mass", kNVarsEff, varsEff, pairHistEffBinLimits, 0x0, -1, kTRUE); // kTRUE = THnSparseF

      continue;
    }

    //-----------------------------------------------------------------------------------------------------------
    // pair histograms
    //-----------------------------------------------------------------------------------------------------------
    if (classStr.Contains("PairSE") || classStr.Contains("PairME")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n",classStr.Data());
      
      // pair histogram
      const Int_t kNVars = 2;
      Int_t vars[kNVars] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt};
      TArrayD pairHistBinLimits[kNVars];
      pairHistBinLimits[0] = TArrayD(kNMassBins,          massBins);
      pairHistBinLimits[1] = TArrayD(kNPtBins,            ptBins);
      man->AddHistogram(classStr.Data(), "PairInvMass", "Differential pair inv. mass", kNVars, vars, pairHistBinLimits, 0x0, -1, kFALSE); // kTRUE = THnSparseF
      
      // pair histogram for efficiency evaluation
      if (task->GetRunOverMC()) {
        const Int_t kNVarsEff = 4;
        Int_t varsEff[kNVarsEff] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, AliReducedVarManager::kEta, AliReducedVarManager::kPhi};
        TArrayD pairHistEffBinLimits[kNVarsEff];
        pairHistEffBinLimits[0] = TArrayD(kNMassBins,   massBins);
        pairHistEffBinLimits[1] = TArrayD(kNPtBinsEff,  ptBinsEff);
        pairHistEffBinLimits[2] = TArrayD(kNEtaBins,    etaBins);
        pairHistEffBinLimits[3] = TArrayD(kNPhiBins,    phiBins);
        man->AddHistogram(classStr.Data(), "PairInvMass_Eff", "Differential pair inv. mass", kNVarsEff, varsEff, pairHistEffBinLimits, 0x0, -1, kTRUE); // kTRUE = THnSparseF
      }

      // kinematic quantities
      man->AddHistogram(classStr.Data(), "Pt", "", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE,
                        125, 0.0, 5.0,                                AliReducedVarManager::kMass);
      man->AddHistogram(classStr.Data(), "Pt_Mass", "#it{p}_{T} vs inv. mass", kFALSE,
                        70, 1.4, 4.2,                                 AliReducedVarManager::kMass,
                        400, 0.0, 40.0,                               AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE,
                        80, -2.0,   2.0,                              AliReducedVarManager::kRap);
      man->AddHistogram(classStr.Data(), "Pt_Rapidity", "pT vs. Rapidity", kFALSE,
                        80, -2.0,   2.0,                              AliReducedVarManager::kRap,
                        400, 0.0, 40.0,                               AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE,
                        100, 0.0, 6.3,                                AliReducedVarManager::kPhi);

      // cluster energy for leg
      man->AddHistogram(classStr.Data(), "Leg1EMCMatchedEnergy", "Energy from the calorimeter matched cluster (leg 1)", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kPairLegEMCALmatchedEnergy+0);
      man->AddHistogram(classStr.Data(), "Leg2EMCMatchedEnergy", "Energy from the calorimeter matched cluster (leg 2)", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kPairLegEMCALmatchedEnergy+1);
      man->AddHistogram(classStr.Data(), "Leg1EMCMatchedEnergy_Leg2EMCMatchedEnergy", "Energy from the calorimeter matched cluster (leg 1 vs. leg 2)", kFALSE,
                        200, 0.0, 100.,                               AliReducedVarManager::kPairLegEMCALmatchedEnergy+0,
                        200, 0.0, 100.,                               AliReducedVarManager::kPairLegEMCALmatchedEnergy+1);

      // SPD pair type
      man->AddHistogram(classStr.Data(), "PairTypeSPD", "Pair type SPD", kFALSE,
                        kNPairTypeSPDBins-1, pairTypeSPDBins,         AliReducedVarManager::kPairTypeSPD);
      man->AddHistogram(classStr.Data(), "PairTypeSPD_Pt", "Pair type SPD vs. pt", kFALSE,
                        kNPtBins-1, ptBins,                           AliReducedVarManager::kPt,
                        kNPairTypeSPDBins-1, pairTypeSPDBins,         AliReducedVarManager::kPairTypeSPD);
      
      // pseudo-proper decay length
      man->AddHistogram(classStr.Data(), "DecayLength", "pseudo-proper decay length", kFALSE,
                        800, -2.0, 2.0,                               AliReducedVarManager::kPseudoProperDecayTime);
      man->AddHistogram(classStr.Data(), "DecayLength_Pt", "pseudo-proper decay length vs. pt", kFALSE,
                        80, 0.0, 40.0,                                AliReducedVarManager::kPt,
                        800, -2.0, 2.0,                               AliReducedVarManager::kPseudoProperDecayTime);
      man->AddHistogram(classStr.Data(), "DecayLength_PairTypeSPD", "pseudo-proper decay length vs. pair type SPD", kFALSE,
                        3, -0.5, 2.5,                                 AliReducedVarManager::kPairTypeSPD,
                        800, -2.0, 2.0,                               AliReducedVarManager::kPseudoProperDecayTime);
      man->AddHistogram(classStr.Data(), "DecayLength_Pt_PairTypeSPD", "pseudo-proper decay length vs. pt vs. SPD pair type", kFALSE,
                        800, -2.0, 2.0,                               AliReducedVarManager::kPseudoProperDecayTime,
                        80, 0.0, 40.0,                                AliReducedVarManager::kPt,
                        3, -0.5, 2.5,                                 AliReducedVarManager::kPairTypeSPD);

      continue;
    }
  }
}
