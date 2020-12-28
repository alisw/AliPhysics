//===============================================================================================================
// addtask for pseudo-proper decay length resolution in MC anchored to pp 13TeV data (last updated: 28/08/2019)
//
// NOTE:  - This addtaks needs smeared (improver task) MC trees as input!
//        - The cuts must be identical to AddTask_laltenka_jpsi2ee_correlations!
//
//===============================================================================================================

// incl./prompt/non-prompt flag for weights
enum WeightFlags {
  kWeightJpsiInclusive=0,
  kWeightJpsiNonPrompt,
  kWeightJpsiPrompt,
  kNWeights=3
};

// MC signals, in order of appearance in AddTask_ffionda_dstDCAcorrections_new.C
enum MCFilters {
  kJpsiInclusive=0,
  kJpsiNonPrompt,
  kJpsiPrompt,
  kJpsiRadiative,
  kJpsiNonRadiative,
  kJpsiNonPromptRadiative,    // kJpsiFromBRadiative
  kJpsiNonPromptNonRadiative, // kJpsiFromBNonRadiative
  kJpsiDecayElectron,
  kJpsiNonPromptDecayElectron,
  kJpsiPromptDecayElectron,
  kJpsiRadiativeDecayElectron,
  kJpsiNonRadiativeDecayElectron,
  kJpsiDecayPhoton
};

// forward declaration
void SetEventCuts(AliReducedAnalysisJpsi2ee* task);
void SetJpsiElectronTrackCuts(AliReducedAnalysisJpsi2ee* task);
void SetTrackPrefilterCuts(AliReducedAnalysisJpsi2ee* task);
void SetPairCuts(AliReducedAnalysisJpsi2ee* task);
void SetPairPrefilterCuts(AliReducedAnalysisJpsi2ee* task);
void SetMCSignalCutsLeg(AliReducedAnalysisJpsi2ee* task);
void SetMCSignalCutsJPsi(AliReducedAnalysisJpsi2ee* task, Int_t weightChoice=kWeightJpsiInclusive);
void SetupHistogramManager(AliReducedAnalysisJpsi2ee* task);
void DefineHistograms(AliReducedAnalysisJpsi2ee* task);

//_______________________________________________________________________________________________________________
AliAnalysisTask* AddTask_laltenka_decayLengthResolution(Bool_t isAliRoot=kTRUE, Int_t runMode=1, TString prodString="") {

  //
  // arguments:
  //  - isAliroot = kTRUE for ESD/AOD analysis, = KFALSE for reduced trees
  //  - runMode = 1 for AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents, = 2 for AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree
  //
  printf("INFO on AddTask_laltenka_decayLengthResolution(): (isAliRoot, runMode) :: (%d,%d)\n", isAliRoot, runMode);

  // settings according to production
  //-----------------------------------------------------------------------------------------------------------
  TObjArray*  prodStringArr = prodString.Tokenize("_");
  TString     prod          = prodStringArr->At(0)->GetName();
  
  // weight choice
  Int_t   weightChoice = kWeightJpsiInclusive;
  TString weightString = "";
  if (prodStringArr->GetEntries()>1) {
    weightString = prodStringArr->At(1)->GetName();
    if (!weightString.CompareTo("nonPrompt"))  weightChoice = kWeightJpsiNonPrompt;
    if (!weightString.CompareTo("prompt"))     weightChoice = kWeightJpsiPrompt;
  }
  
  // initialize analysis task
  //-----------------------------------------------------------------------------------------------------------
  AliReducedAnalysisJpsi2ee* jpsiAnalysis = new AliReducedAnalysisJpsi2ee("jpsiAnalysis", "Decay length resolution analysis");
  jpsiAnalysis->Init();
  // ----------------------------------------------------------------------------------
  jpsiAnalysis->SetLoopOverTracks(kTRUE);
  jpsiAnalysis->SetRunEventMixing(kFALSE);
  jpsiAnalysis->SetRunLikeSignPairing(kFALSE);
  jpsiAnalysis->SetRunPairing(kTRUE);
  jpsiAnalysis->SetRunPrefilter(kTRUE);
  jpsiAnalysis->SetStoreJpsiCandidates(kTRUE);
  jpsiAnalysis->SetRunOverMC(kTRUE);

  // set analysis and prefilter cuts
  //-----------------------------------------------------------------------------------------------------------
  SetEventCuts(jpsiAnalysis);
  SetJpsiElectronTrackCuts(jpsiAnalysis);
  SetTrackPrefilterCuts(jpsiAnalysis);
  SetPairCuts(jpsiAnalysis);
  SetPairPrefilterCuts(jpsiAnalysis);

  // set MC signal selection
  //-----------------------------------------------------------------------------------------------------------
  SetMCSignalCutsLeg(jpsiAnalysis);
  SetMCSignalCutsJPsi(jpsiAnalysis, weightChoice);

  // set MC pT weights, depends on choice of weights (prompt vs. non-prompt)
  //-----------------------------------------------------------------------------------------------------------
  TFile*                                  fileWeights = NULL;
  if (weightChoice==kWeightJpsiInclusive) fileWeights = TFile::Open("/Users/lal035/Documents/PWGDQ/data/calibration/MC/ratio_13TeV_two.root");
  if (weightChoice==kWeightJpsiNonPrompt) fileWeights = TFile::Open("/Users/lal035/Documents/PWGDQ/data/calibration/MC/WeightsFONLLMC.root");
  if (weightChoice==kWeightJpsiPrompt)    fileWeights = TFile::Open("/Users/lal035/Documents/PWGDQ/data/calibration/MC/weightsForPromptJpsi.root");
  TH1F*                                   weights     = NULL;
  if (weightChoice==kWeightJpsiInclusive) weights     = (TH1F*)fileWeights->Get("weightsNew");
  if (weightChoice==kWeightJpsiNonPrompt) weights     = (TH1F*)fileWeights->Get("weights");
  if (weightChoice==kWeightJpsiPrompt)    weights     = (TH1F*)fileWeights->Get("weights");
  weights->SetDirectory(0x0);
  fileWeights->Close();
  jpsiAnalysis->SetMCJpsiPtWeights(weights);
  
  // define histograms histograms
  //-----------------------------------------------------------------------------------------------------------
  SetupHistogramManager(jpsiAnalysis);
  
  // initialize wrapper AliAnalysisTask
  // (in order to run AliReducedAnalysisJpsi2ee in an aliroot analysis train )
  //-----------------------------------------------------------------------------------------------------------
  AliAnalysisTaskReducedEventProcessor* task = new AliAnalysisTaskReducedEventProcessor("ReducedEventAnalysisManager", runMode);
  task->AddTask(jpsiAnalysis);

  // intercept isAliRoot=kFALSE (nothing to be done yet)
  //-----------------------------------------------------------------------------------------------------------
  if (!isAliRoot) return 0;
  
  // get analysis manager
  //-----------------------------------------------------------------------------------------------------------
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_laltenka_decayLengthResolution", "No analysis manager found.");
    return 0;
  }

  // get data container
  //-----------------------------------------------------------------------------------------------------------
  AliAnalysisDataContainer* cReducedEvent = NULL;
  if (runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) {
    printf("INFO on AddTask_laltenka_decayLengthResolution(): use on the fly events\n");
    cReducedEvent = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ReducedEventDQ");
    if (!cReducedEvent) {
      printf("ERROR: In AddTask_laltenka_decayLengthResolution(), ReducedEvent exchange container could not be found!\n");
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
    printf("ERROR: In AddTask_laltenka_decayLengthResolution(), runMode %d not defined!\n", runMode);
    return 0;
  }

  // connect output data containers
  //-----------------------------------------------------------------------------------------------------------
  TString outputName = Form("pseudoProperDecayLengthResolutionHistograms_MB_MC_%s.root", weightString.Data());
  printf("INFO on AddTask_laltenka_decayLengthResolution(): output container name = %s\n", outputName.Data());
  AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("histos", THashList::Class(), AliAnalysisManager::kOutputContainer, outputName);
  mgr->ConnectOutput(task, 1, cOutputHist);
  
  // done
  //-----------------------------------------------------------------------------------------------------------
  return task;
}

//_______________________________________________________________________________________________________________
void SetEventCuts(AliReducedAnalysisJpsi2ee* task) {
  AliReducedEventCut* eventCut = new AliReducedEventCut("EventCut", "Event selection");
  eventCut->AddEventTriggerFilter(AliReducedVarManager::kINT7);
  eventCut->AddCut(AliReducedVarManager::kVtxZ, -10.0, 10.0);
  task->AddEventCut(eventCut);
}

//_______________________________________________________________________________________________________________
void SetJpsiElectronTrackCuts(AliReducedAnalysisJpsi2ee* task) {
  // standard cut
  AliReducedTrackCut* standardCut = new AliReducedTrackCut("standardElectron","Electron selection standard");
  standardCut->SetTrackFilterBit(0); // fQualityFlag bit 32 (see AddTask_ffionda_dstDCAcorrections_new.C)
  standardCut->AddCut(AliReducedVarManager::kPt,                  1.0,  30.0);
  standardCut->AddCut(AliReducedVarManager::kEta,                -0.9,   0.9);
  standardCut->AddCut(AliReducedVarManager::kDcaXY,              -1.0,   1.0);
  standardCut->AddCut(AliReducedVarManager::kDcaZ,               -3.0,   3.0);
  standardCut->AddCut(AliReducedVarManager::kTPCncls,            70.0, 160.0);
  standardCut->AddCut(AliReducedVarManager::kTPCchi2,             0.0,   4.0);
  standardCut->AddCut(AliReducedVarManager::kITSchi2,             0.0,  30.0);
  standardCut->AddCut(AliReducedVarManager::kTPCNclusBitsFired,   6.0,   9.0);
  standardCut->SetRejectKinks();
  standardCut->SetRequestITSrefit();
  standardCut->SetRequestTPCrefit();
  standardCut->SetRequestSPDany();
  standardCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron), -1.5, 3.0);
  standardCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton),    3.5, 1.e8);
  standardCut->AddCut(static_cast<AliReducedVarManager::Variables>(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion),      3.5, 1.e8);
  task->AddTrackCut(standardCut);
}

//_______________________________________________________________________________________________________________
void SetTrackPrefilterCuts(AliReducedAnalysisJpsi2ee* task) {
  AliReducedTrackCut* trackPrefilterCut = new AliReducedTrackCut("TrackPrefilterCut", "Track prefilter selection");
  trackPrefilterCut->SetTrackFilterBit(1); // fQualityFlag bit 33 (see AddTask_ffionda_dstDCAcorrections_new.C)
  trackPrefilterCut->AddCut(AliReducedVarManager::kPt, 0.7, 1.0e+30);
  task->AddPrefilterTrackCut(trackPrefilterCut);
}

//_______________________________________________________________________________________________________________
void SetPairCuts(AliReducedAnalysisJpsi2ee* task) {
  AliReducedTrackCut* pairCut = new AliReducedTrackCut("PairCut", "Pair selection");
  pairCut->AddCut(AliReducedVarManager::kPt,    0.0, 1.0e+30);
  pairCut->AddCut(AliReducedVarManager::kRap,  -0.9, 0.9);
  task->AddPairCut(pairCut);
}

//_______________________________________________________________________________________________________________
void SetPairPrefilterCuts(AliReducedAnalysisJpsi2ee* task) {
  AliReducedVarCut* pairPrefilterCut = new AliReducedVarCut("PairPrefilterCut", "Pair prefilter selection");
  pairPrefilterCut->AddCut(AliReducedVarManager::kMass, 0.0, 0.05, kTRUE);
  task->AddPrefilterPairCut(pairPrefilterCut);
}

//_______________________________________________________________________________________________________________
void SetMCSignalCutsLeg(AliReducedAnalysisJpsi2ee* task) {
  // electron from inclusive J/psi
  AliReducedTrackCut* trueElectron = new AliReducedTrackCut("TrueElectron", "reconstructed electrons with MC truth");
  trueElectron->SetMCFilterBit(kJpsiDecayElectron);
  task->AddLegCandidateMCcut(trueElectron);
  
  // electron from non-prompt J/psi
  AliReducedTrackCut* trueElectronNonPrompt = new AliReducedTrackCut("TrueElectronNonPrompt", "reconstructed electrons with MC truth");
  trueElectronNonPrompt->SetMCFilterBit(kJpsiNonPromptDecayElectron);
  task->AddLegCandidateMCcut(trueElectronNonPrompt);
  
  // electron from prompt J/psi
  AliReducedTrackCut* trueElectronPrompt = new AliReducedTrackCut("TrueElectronPrompt", "reconstructed electrons with MC truth");
  trueElectronPrompt->SetMCFilterBit(kJpsiPromptDecayElectron);
  task->AddLegCandidateMCcut(trueElectronPrompt);
}

//_______________________________________________________________________________________________________________
void SetMCSignalCutsJPsi(AliReducedAnalysisJpsi2ee* task, Int_t weightChoice/*=kWeightJpsiInclusive*/) {
  // inclusive J/psi
  AliReducedTrackCut* mcTruthJpsi = new AliReducedTrackCut("mcTruthJpsi", "Pure MC truth J/psi");
  mcTruthJpsi->SetMCFilterBit(kJpsiInclusive);
  mcTruthJpsi->SetApplyReweightMCpt(kFALSE); //reweight is applied only for the prompt or non-prompt (kTRUE will reweight twice)
  mcTruthJpsi->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  AliReducedTrackCut* mcTruthJpsiElectron = new AliReducedTrackCut("mcTruthJpsiElectron", "Pure MC truth electron from J/psi");
  mcTruthJpsiElectron->SetMCFilterBit(kJpsiDecayElectron);
  mcTruthJpsiElectron->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  mcTruthJpsiElectron->AddCut(AliReducedVarManager::kPt, 1.0, 100.0);
  task->AddJpsiMotherMCCut(mcTruthJpsi, mcTruthJpsiElectron);
  
  // non-prompt J/psi
  AliReducedTrackCut* mcTruthJpsiNonPrompt = new AliReducedTrackCut("mcTruthJpsiNonPrompt", "Pure MC truth J/psi");
  mcTruthJpsiNonPrompt->SetMCFilterBit(kJpsiNonPrompt);
  if (weightChoice==kWeightJpsiNonPrompt) mcTruthJpsiNonPrompt->SetApplyReweightMCpt(kTRUE);  // reweight
  else                                    mcTruthJpsiNonPrompt->SetApplyReweightMCpt(kFALSE); // no reweight
  mcTruthJpsiNonPrompt->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  AliReducedTrackCut* mcTruthJpsiElectronNonPrompt = new AliReducedTrackCut("mcTruthJpsiElectronNonPrompt", "Pure MC truth electron from J/psi");
  mcTruthJpsiElectronNonPrompt->SetMCFilterBit(kJpsiNonPromptDecayElectron);
  mcTruthJpsiElectronNonPrompt->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  mcTruthJpsiElectronNonPrompt->AddCut(AliReducedVarManager::kPt, 1.0, 100.0);
  task->AddJpsiMotherMCCut(mcTruthJpsiNonPrompt, mcTruthJpsiElectronNonPrompt);
  
  // prompt J/psi
  AliReducedTrackCut* mcTruthJpsiPrompt = new AliReducedTrackCut("mcTruthJpsiPrompt", "Pure MC truth J/psi");
  mcTruthJpsiPrompt->SetMCFilterBit(kJpsiPrompt);
  if (weightChoice==kWeightJpsiPrompt || weightChoice==kWeightJpsiInclusive)  mcTruthJpsiPrompt->SetApplyReweightMCpt(kTRUE);  // reweight
  else                                                                        mcTruthJpsiPrompt->SetApplyReweightMCpt(kFALSE); // no reweight
  mcTruthJpsiPrompt->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  AliReducedTrackCut* mcTruthJpsiElectronPrompt = new AliReducedTrackCut("mcTruthJpsiElectronPrompt", "Pure MC truth electron from J/psi");
  mcTruthJpsiElectronPrompt->SetMCFilterBit(kJpsiPromptDecayElectron);
  mcTruthJpsiElectronPrompt->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  mcTruthJpsiElectronPrompt->AddCut(AliReducedVarManager::kPt, 1.0, 100.0);
  task->AddJpsiMotherMCCut(mcTruthJpsiPrompt, mcTruthJpsiElectronPrompt);
}

//_______________________________________________________________________________________________________________
void SetupHistogramManager(AliReducedAnalysisJpsi2ee* task) {
  AliReducedVarManager::SetDefaultVarNames();
  DefineHistograms(task);
  AliReducedVarManager::SetUseVars((Bool_t*)task->GetHistogramManager()->GetUsedVars());
}

//_______________________________________________________________________________________________________________
void DefineHistograms(AliReducedAnalysisJpsi2ee* task) {
  AliHistogramManager* man = task->GetHistogramManager();

  //-------------------------------------------------------------------------------------------------------------
  // histogram classes
  //-------------------------------------------------------------------------------------------------------------
  TString histClasses = "";
  for (Int_t i=0; i<task->GetNTrackCuts(); ++i) {
    TString cutName = task->GetTrackCutName(i);
    histClasses += Form("PairSEPM_%s;", cutName.Data());
    if (task->GetRunOverMC()) {
      for(Int_t mcSel=0; mcSel<task->GetNLegCandidateMCcuts(); ++mcSel)
        histClasses += Form("PairSEPM_%s_%s;", cutName.Data(), task->GetLegCandidateMCcutName(mcSel));
    }
  }
  for(Int_t mcSel=0; mcSel<task->GetNJpsiMotherMCCuts(); ++mcSel) {
    histClasses += Form("%s_PureMCTruth_BeforeSelection;", task->GetJpsiMotherMCcutName(mcSel));
    histClasses += Form("%s_PureMCTruth_AfterSelection;", task->GetJpsiMotherMCcutName(mcSel));
  }

  //-------------------------------------------------------------------------------------------------------------
  // histograms
  //-------------------------------------------------------------------------------------------------------------
  
  // pT binning
  const Int_t kNPtBins = 401; // [0., 40.]
  Double_t ptBins[kNPtBins];
  for (Int_t i=0; i<kNPtBins; i++) ptBins[i] = i*0.1;
  
  // pseudo proper decay length binning
  const Int_t kNDecayLengthBins = 1001; // [-1., 1.]
  Double_t decayLengthBins[kNDecayLengthBins];
  for (Int_t i=0; i<kNDecayLengthBins; i++) decayLengthBins[i] = -1.0 + i*0.002; // 20 mum bin width

  // pair type SPD binning
  const Int_t kNPairTypeSPDBins = 4;
  Double_t pairTypeSPDBins[kNPairTypeSPDBins] = {-0.5, 0.5, 1.5, 2.5};
  
  // add histograms according to class to histogram manager
  TString classesStr(histClasses);
  TObjArray* arr = classesStr.Tokenize(";");
  
  // loop over histogram classes and add histograms
  printf("INFO on AddTask_laltenka_decayLengthResolution(): Histgram classes included in histogram manager\n");
  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();

    // MC rec.
    if(classStr.Contains("Pair")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n", classStr.Data());

      man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE,
                        kNPtBins-1,           ptBins[0],          ptBins[kNPtBins-1],                           AliReducedVarManager::kPtMC);
      man->AddHistogram(classStr.Data(), "Pt", "p_{T}", kFALSE,
                        kNPtBins-1,           ptBins[0],          ptBins[kNPtBins-1],                           AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "PairTypeSPD", "SPD pair type", kFALSE,
                        kNPairTypeSPDBins-1,  pairTypeSPDBins[0], pairTypeSPDBins[kNPairTypeSPDBins-1],         AliReducedVarManager::kPairTypeSPD);
      man->AddHistogram(classStr.Data(), "DecayLength", "pseudo-proper decay length", kFALSE,
                        kNDecayLengthBins-1,  decayLengthBins[0], decayLengthBins[kNDecayLengthBins-1],         AliReducedVarManager::kPseudoProperDecayTime);
      man->AddHistogram(classStr.Data(), "DecayLengthErr", "pseudo-proper decay length uncertainty", kFALSE,
                        kNDecayLengthBins-1,  decayLengthBins[0], decayLengthBins[kNDecayLengthBins-1],         AliReducedVarManager::kPseudoProperDecayTimeErr);
      man->AddHistogram(classStr.Data(), "DecayLengthNorm", "normlized pseudo-proper decay length", kFALSE,
                        kNDecayLengthBins-1,  decayLengthBins[0]*10., decayLengthBins[kNDecayLengthBins-1]*10., AliReducedVarManager::kPseudoProperDecayTimeNorm);
      man->AddHistogram(classStr.Data(), "Pt_DecayLength", "pseudo-proper decay length", kFALSE,
                        kNPtBins-1,           ptBins[0],          ptBins[kNPtBins-1],                           AliReducedVarManager::kPt,
                        kNDecayLengthBins-1,  decayLengthBins[0], decayLengthBins[kNDecayLengthBins-1],         AliReducedVarManager::kPseudoProperDecayTime);
      man->AddHistogram(classStr.Data(), "Pt_DecayLengthErr", "pseudo-proper decay length uncertainty vs. pT", kFALSE,
                        kNPtBins-1,           ptBins[0],          ptBins[kNPtBins-1],                           AliReducedVarManager::kPt,
                        kNDecayLengthBins-1,  decayLengthBins[0], decayLengthBins[kNDecayLengthBins-1],         AliReducedVarManager::kPseudoProperDecayTimeErr);
      man->AddHistogram(classStr.Data(), "Pt_DecayLengthNorm", "normalized pseudo-proper decay length vs. pT", kFALSE,
                        kNPtBins-1,           ptBins[0],          ptBins[kNPtBins-1],                           AliReducedVarManager::kPt,
                        kNDecayLengthBins-1,  decayLengthBins[0]*10., decayLengthBins[kNDecayLengthBins-1]*10., AliReducedVarManager::kPseudoProperDecayTimeNorm);
      man->AddHistogram(classStr.Data(), "PtMC_DecayLength", "pseudo-proper decay length", kFALSE,
                        kNPtBins-1,           ptBins[0],          ptBins[kNPtBins-1],                           AliReducedVarManager::kPtMC,
                        kNDecayLengthBins-1,  decayLengthBins[0], decayLengthBins[kNDecayLengthBins-1],         AliReducedVarManager::kPseudoProperDecayTime);
      man->AddHistogram(classStr.Data(), "PtMC_DecayLengthErr", "pseudo-proper decay length uncertainty vs. MC pT", kFALSE,
                        kNPtBins-1,           ptBins[0],          ptBins[kNPtBins-1],                           AliReducedVarManager::kPtMC,
                        kNDecayLengthBins-1,  decayLengthBins[0], decayLengthBins[kNDecayLengthBins-1],         AliReducedVarManager::kPseudoProperDecayTimeErr);
      man->AddHistogram(classStr.Data(), "PtMC_DecayLengthNorm", "normalized pseudo-proper decay length vs. MC pT", kFALSE,
                        kNPtBins-1,           ptBins[0],          ptBins[kNPtBins-1],                           AliReducedVarManager::kPtMC,
                        kNDecayLengthBins-1,  decayLengthBins[0]*10., decayLengthBins[kNDecayLengthBins-1]*10., AliReducedVarManager::kPseudoProperDecayTimeNorm);
      man->AddHistogram(classStr.Data(), "Pt_DecayLength_PairTypeSPD", "pt vs. pseudo-proper decay length vs. SPD pair type", kFALSE,
                        kNPtBins-1,           ptBins[0],          ptBins[kNPtBins-1],                           AliReducedVarManager::kPt,
                        kNDecayLengthBins-1,  decayLengthBins[0], decayLengthBins[kNDecayLengthBins-1],         AliReducedVarManager::kPseudoProperDecayTime,
                        kNPairTypeSPDBins-1,  pairTypeSPDBins[0], pairTypeSPDBins[kNPairTypeSPDBins-1],         AliReducedVarManager::kPairTypeSPD);
      man->AddHistogram(classStr.Data(), "Pt_DecayLengthErr_PairTypeSPD", "pt vs. pseudo-proper decay length uncertainty vs. SPD pair type", kFALSE,
                        kNPtBins-1,           ptBins[0],          ptBins[kNPtBins-1],                           AliReducedVarManager::kPt,
                        kNDecayLengthBins-1,  decayLengthBins[0], decayLengthBins[kNDecayLengthBins-1],         AliReducedVarManager::kPseudoProperDecayTimeErr,
                        kNPairTypeSPDBins-1,  pairTypeSPDBins[0], pairTypeSPDBins[kNPairTypeSPDBins-1],         AliReducedVarManager::kPairTypeSPD);
      man->AddHistogram(classStr.Data(), "Pt_DecayLengthNorm_PairTypeSPD", "pt vs. normalized pseudo-proper decay length vs. SPD pair type", kFALSE,
                        kNPtBins-1,           ptBins[0],          ptBins[kNPtBins-1],                           AliReducedVarManager::kPt,
                        kNDecayLengthBins-1,  decayLengthBins[0]*10., decayLengthBins[kNDecayLengthBins-1]*10., AliReducedVarManager::kPseudoProperDecayTimeNorm,
                        kNPairTypeSPDBins-1,  pairTypeSPDBins[0], pairTypeSPDBins[kNPairTypeSPDBins-1],         AliReducedVarManager::kPairTypeSPD);

      continue;
    }
    
    // MC gen.
    if (classStr.Contains("PureMCTruth")) {
      man->AddHistClass(classStr.Data());
      printf("%s\n", classStr.Data());

      man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE,
                        kNPtBins-1,           ptBins[0],          ptBins[kNPtBins-1],                   AliReducedVarManager::kPtMC);
      man->AddHistogram(classStr.Data(), "DecayLengthMC", "pseudo-proper decay length MC", kFALSE,
                        kNDecayLengthBins-1,  decayLengthBins[0], decayLengthBins[kNDecayLengthBins-1], AliReducedVarManager::kPseudoProperDecayTimeMC);
      man->AddHistogram(classStr.Data(), "PtMC_DecayLengthMC", "pseudo-proper decay length MC", kFALSE,
                        kNPtBins-1,           ptBins[0],          ptBins[kNPtBins-1],                   AliReducedVarManager::kPtMC,
                        kNDecayLengthBins-1,  decayLengthBins[0], decayLengthBins[kNDecayLengthBins-1], AliReducedVarManager::kPseudoProperDecayTimeMC);

      continue;
    }
  }
}
