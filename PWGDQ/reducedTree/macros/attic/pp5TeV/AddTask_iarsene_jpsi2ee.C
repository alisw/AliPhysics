//void Setup(AliReducedAnalysisJpsi2ee* processor, TString prod="LHC10h");
//AliHistogramManager* SetupHistogramManager(TString prod="LHC10h");
//void DefineHistograms(AliHistogramManager* man, TString prod="LHC10h");

enum TrackFilters {
   kBasicCut=0,
   kPrefilterCut,
   kNTrackFilters
};

enum MCFilters {
   kJpsiInclusive=0,
   kJpsiNonPrompt,
   kJpsiPrompt,
   kJpsiRadiative,
   kJpsiNonRadiative,
   kJpsiDecayElectron,
   kJpsiNonPromptDecayElectron,
   kJpsiPromptDecayElectron,
   kJpsiRadiativeDecayElectron,
   kJpsiNonRadiativeDecayElectron,
   kJpsiDecayPhoton
};

enum TrackQualityFilterBits {       // as in AliReducedBaseTrack::fQualityFlags
   kTPCEventPlaneContributor = 0,
   kPhotonConversionV0Leg, 
   kK0sV0Leg,
   kLambdaV0Leg,
   kAntiLambdaV0Leg,
   kKink0PositiveIndex,
   kKink1PositiveIndex,
   kKink2PositiveIndex,
   kPurePhotonConversionV0Leg,
   kPureK0sV0Leg, 
   kPureLambdaV0Leg,
   kAntiLambdaV0Leg, 
   kKink0NegativeIndex, 
   kKink1NegativeIndex,
   kKink2NegativeIndex,
   kAODbits=15,
   kTRDmatch=26
};

// run list string used for run IDs in the VarManager
const Int_t gkNRuns = 72;
TString gRunList = "244340;244343;244351;244355;244359;244364;244377;244411;244416;244418;244421;244453;244456;244480;244481;244482;244483;244484;244531;244540;244542;244617;244618;244619;244626;244627;244628;282008;282016;282021;282025;282030;282031;282050;282051;282078;282098;282099;282118;282119;282120;282122;282123;282125;282126;282127;282146;282147;282189;282206;282224;282227;282229;282230;282247;282302;282303;282304;282305;282306;282307;282309;282312;282313;282314;282340;282341;282342;282343;282365;282366;282367";

const Int_t gkNZBins = 11;
Double_t gZBinLims[gkNZBins] = {
   -10.0, -9.0, -8.0, -7.0, -4.0, 0.0,
   4.0, 7.0, 8.0, 9.0, 10.
};
const Int_t gkNZBinsMC = 2;
Double_t gZBinLimsMC[gkNZBinsMC] = {
   -10.0, 10.
};

/*const Int_t gkNSPDtrkBins = 8;
Double_t gSPDtrkBinLims[gkNSPDtrkBins] = {
   0., 5., 10., 15., 25., 35., 50., 110.
};*/
const Int_t gkNSPDtrkBins = 2;
Double_t gSPDtrkBinLims[gkNSPDtrkBins] = {
   0., 110.
};

// pt binning
const Int_t gkNPtBins = 15;
Double_t gPtLims[gkNPtBins] = {
   0.0, 0.15, 0.30, 0.50, 1.0, 
   1.3, 1.50, 2.00, 3.00, 4.0, 
   5.0, 7.0, 10.,  15.,  20.
};
// pt binning for MC
const Int_t gkNPtBinsMC = 400;
Double_t gPtLimsMC[gkNPtBinsMC+1] = {0.};

// mass binning
const Int_t gkNMassBins = 126;
Double_t gMassBins[gkNMassBins] = {0.};


//__________________________________________________________________________________________
AliAnalysisTask* AddTask_iarsene_jpsi2ee(Bool_t isAliRoot=kTRUE, Int_t runMode=1, TString prod="LHC10h"){    
   //
   // isAliRoot=kTRUE for ESD/AOD analysis in AliROOT, kFALSE for root analysis on reduced trees
   // runMode=1 (AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents)
   //               =2 (AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree)
   //
   //get the current analysis manager

  printf("INFO on AddTask_iarsene_jpsi2ee(): (isAliRoot, runMode) :: (%d,%d) \n", isAliRoot, runMode);

  if(!isAliRoot) return 0x0;
  
  AliReducedAnalysisJpsi2ee* jpsi2eeAnalysis = new AliReducedAnalysisJpsi2ee("Jpsi2eeAnalysis","Jpsi->ee analysis");
  jpsi2eeAnalysis->Init();
  Setup(jpsi2eeAnalysis, prod);
  // initialize an AliAnalysisTask which will wrapp the AliReducedAnalysisJpsi2ee such that it can be run in an aliroot analysis train (e.g. LEGO, local analysis etc)
  AliAnalysisTaskReducedEventProcessor* task = new AliAnalysisTaskReducedEventProcessor("ReducedEventAnalysisManager", runMode);
  task->AddTask(jpsi2eeAnalysis);
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
     Error("AddTask_iarsene_dst", "No analysis manager found.");
     return 0;
  }
     
  AliAnalysisDataContainer* cReducedEvent = NULL;
  if(runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) {
    printf("INFO on AddTask_iarsene_jpsi2ee(): use on the fly events ");
    cReducedEvent = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ReducedEventDQ");
    if(!cReducedEvent) {
      printf("ERROR: In AddTask_iarsene_jpsi2ee(), couldn't find exchange container with ReducedEvent! ");
      printf("             You need to use AddTask_iarsene_dst() in the on-the-fly reduced events mode!");
      return 0x0;
    }
  }
            
  mgr->AddTask(task);
      
  if(runMode==AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree) 
     mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
      
  if(runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) 
     mgr->ConnectInput(task, 0, cReducedEvent);
  
  TString outputFilename = "AnalysisHistograms_jpsi2ee_pp2017.root";
  if(jpsi2eeAnalysis->GetRunOverMC()) outputFilename = "AnalysisHistograms_jpsi2ee_pp2017_MC.root";
  AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("jpsi2eeHistos", THashList::Class(),
                                                               AliAnalysisManager::kOutputContainer, outputFilename.Data());
  mgr->ConnectOutput(task, 1, cOutputHist );
    
  return task;
}


//_________________________________________________________________
void Setup(AliReducedAnalysisJpsi2ee* processor, TString prod /*="LHC10h"*/) {
  //
  // Configure the analysis task
  // Setup histograms, handlers, cuts, etc.
  //
  //processor->SetLoopOverTracks(kFALSE);
  //processor->SetRunEventMixing(kFALSE);
  //processor->SetRunPairing(kFALSE);
  //processor->SetRunOverMC(kTRUE);
  //processor->SetRunLikeSignPairing(kFALSE);
  //processor->SetRunPrefilter(kFALSE);
  
  TFile* wFile = 0x0; 
  if(prod.Contains("LHC15n")) {
     wFile = TFile::Open("/home/iarsene/work/ALICE/pp2017Analysis/weights.root");
     TH1F* weights = (TH1F*)wFile->Get("weights2");
     processor->SetMCJpsiPtWeights(weights);
  }
  
  if(prod.Contains("LHC17p") && !processor->GetRunOverMC()) {
     TFile* corrFile = TFile::Open("/home/iarsene/work/ALICE/pp2017Analysis/nSigmaCalibMap_17p.root");
     AliReducedVarManager::SetTPCelectronCorrectionMaps((TH2F*)corrFile->Get("nSigmaCentroidMap"), (TH2F*)corrFile->Get("nSigmaWidthMap"), AliReducedVarManager::kEta, AliReducedVarManager::kPin);
     corrFile->Close();
  }
  if(prod.Contains("LHC17q") && !processor->GetRunOverMC()) {
     TFile* corrFile = TFile::Open("/home/iarsene/work/ALICE/pp2017Analysis/nSigmaCalibMap_17q.root");
     AliReducedVarManager::SetTPCelectronCorrectionMaps((TH2F*)corrFile->Get("nSigmaCentroidMap"), (TH2F*)corrFile->Get("nSigmaWidthMap"), AliReducedVarManager::kEta, AliReducedVarManager::kPin);
     corrFile->Close();
  }
   
  // Set event cuts
  AliReducedEventCut* evCut1 = new AliReducedEventCut("Centrality","Centrality selection");
  evCut1->AddCut(AliReducedVarManager::kVtxZ, -10.0, 10.0);
  evCut1->AddCut(AliReducedVarManager::kIsPhysicsSelection, 0.1, 2.);   // request physics selection
  evCut1->EnableVertexDistanceCut();
  evCut1->AddEventTriggerFilter(AliReducedVarManager::kINT7);
  
  processor->AddEventCut(evCut1);
  
  SetupCandidateElectronCuts(processor, prod);
  
  
  if(processor->GetRunOverMC()) SetupMCsignals(processor);
  
  SetupHistogramManager(processor, prod);
  SetupMixingHandler(processor);
}


//_________________________________________________________________
void SetupCandidateElectronCuts(AliReducedAnalysisJpsi2ee* task, TString prod) {
   //
   //  setup a standard electron candidate cut
   //   
   
   //SetupTrackingSystematicsCutsSingleLeg(task, kFALSE);
   
   //SetupTrackingSystematicsCuts(task, kTRUE, kFALSE);
   //SetupPIDSystematicsCuts(task, kTRUE, kFALSE, kTRUE);
   //SetupDCASystematicsCuts(task, kTRUE, kFALSE);
   
   //AddStandardCutForComparisons(task, kFALSE);
   //AddStandardCut(task, kTRUE, kFALSE);
   
   //*******************************************************************************
   // To run for systematics (one configuration at a time)
   
   // Standard cut alone
   
   if(prod.Contains("LHC15n")) AddStandardCut(task, kTRUE, kFALSE);                                // standard cut with SPDany and no PID calibration   
   else                       AddStandardCut(task, kTRUE, kTRUE);                                 // standard cut with SPDany and with PID calibration   
   
   //   Tracking systematics -----------------------------------------------------------------------
   /*if(prod.Contains("LHC15n")) SetupTrackingSystematicCuts(task, kTRUE, kFALSE);                     // tracking systematics + standard cut + no PID calibration   
   else                        SetupTrackingSystematicCuts(task, kTRUE, kTRUE);                    // tracking systematics + standard cut + with PID calibration   
   */
   //    Single electron PID systematics -----------------------------------------------------------
   //   Remember to switch off the prefilter
   //if(prod.Contains("LHC15n")) SetupPIDSystematicsCuts(task, kTRUE, kFALSE, kTRUE, kTRUE);        // PID systematics + SPDany + no PID calibration + no standard cut + with conversion electrons 
   //else                        SetupPIDSystematicsCuts(task, kTRUE, kTRUE, kFALSE, kTRUE);        // PID systematics + SPDany + with PID calibration + no standard cut + with conversion electrons 
   
   //   Dielectron pair PID systematics ------------------------------------------------------------
   /*if(prod.Contains("LHC15n")) SetupPIDSystematicsCuts(task, kTRUE, kFALSE, kTRUE, kFALSE);                 // PID systematics + SPDany + with PID calibration + standard cut + with primary electrons
   else                        SetupPIDSystematicsCuts(task, kTRUE, kTRUE, kTRUE, kFALSE);                // PID systematics + SPDany + no PID calibration + standard cut + with primary electrons
   */
   //*******************************************************************************
   
   // set track prefilter cuts
   AliReducedTrackCut* prefTrackCut = new AliReducedTrackCut("prefCut","prefilter selection");
   prefTrackCut->SetTrackFilterBit(kPrefilterCut);
   prefTrackCut->AddCut(AliReducedVarManager::kPt, 0.3,100.0);
   prefTrackCut->SetRequestITSrefit();
   task->AddPrefilterTrackCut(prefTrackCut);  
   
   // set pair prefilter cuts
   AliReducedVarCut* prefPairCut = new AliReducedVarCut("prefCutM50MeV","prefilter pair cuts");
   prefPairCut->AddCut(AliReducedVarManager::kMass, 0.0, 0.05, kTRUE);
   task->AddPrefilterPairCut(prefPairCut);
   
   // Set pair cuts
   AliReducedTrackCut* pairCut1 = new AliReducedTrackCut("Ptpair","Pt pair selection");
   pairCut1->AddCut(AliReducedVarManager::kPt, 0.0,100.0);
   pairCut1->AddCut(AliReducedVarManager::kRap, -0.9,0.9);
   task->AddPairCut(pairCut1);
}

//__________________________________________________________________________________________
void SetupPIDSystematicsCuts(AliReducedAnalysisJpsi2ee* task, Bool_t reqSPDany=kTRUE, Bool_t isCalibrated = kFALSE, Bool_t addStandardCut=kTRUE, Bool_t useWithConversions=kFALSE) {
   //
   // Apply tracking cuts variations on the electron candidates
   //
   Bool_t isMC = task->GetRunOverMC();
   
   // common cuts
   AliReducedTrackCut* commonCuts = new AliReducedTrackCut("commonCuts","");
   if(!useWithConversions) commonCuts->SetTrackFilterBit(kBasicCut);        // eta cut, DCA cuts, ITS and TPC refit already applied
   else {
      if(!isMC) commonCuts->SetTrackQualityFilterBit(kPhotonConversionV0Leg, kFALSE);
      commonCuts->SetRequestITSrefit();
   }
   commonCuts->AddCut(AliReducedVarManager::kPt, 0.0, 1.0, kTRUE);  
   commonCuts->AddCut(AliReducedVarManager::kITSchi2, 0.0, 10.);
   commonCuts->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
   if(!useWithConversions) {
      if(reqSPDany) commonCuts->SetRequestSPDany();
      else          commonCuts->SetRequestSPDfirst();
   }
   commonCuts->SetRejectKinks();
   if(!useWithConversions) {
      commonCuts->AddCut(AliReducedVarManager::kITSnclsShared, 1.9, 6., kTRUE);
      commonCuts->AddCut(AliReducedVarManager::kDcaXY, -0.2, 0.2);
      commonCuts->AddCut(AliReducedVarManager::kDcaZ, -0.4, 0.4);
   }
   if(addStandardCut) {
      AddStandardCut(task, reqSPDany, isCalibrated);
      AddStandardCutForComparisons(task, isCalibrated);
   }
   
   if(useWithConversions) {
      AliReducedCompositeCut* basicElectronSelection = new AliReducedCompositeCut("basicElectronCut", "");
      basicElectronSelection->AddCut(commonCuts);
      task->AddTrackCut(basicElectronSelection);
      
      AliReducedCompositeCut* electronFullBand = new AliReducedCompositeCut("electronFullBand", "", kTRUE);
      electronFullBand->AddCut(basicElectronSelection);
      AliReducedTrackCut* electronFullBandVarCut = new AliReducedTrackCut();
      if(isMC || !isCalibrated) electronFullBandVarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
      else       electronFullBandVarCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
      electronFullBandVarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 2.0, 30000.0);
      electronFullBand->AddCut(electronFullBandVarCut);
      task->AddTrackCut(electronFullBand);
      
      AliReducedCompositeCut* electronLowerHalf = new AliReducedCompositeCut("electronLowerHalf", "", kTRUE);
      electronLowerHalf->AddCut(basicElectronSelection);
      AliReducedTrackCut* electronLowerHalfVarCut = new AliReducedTrackCut();
      if(isMC || !isCalibrated) electronLowerHalfVarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 0.0);
      else       electronLowerHalfVarCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 0.0);
      electronLowerHalf->AddCut(electronLowerHalfVarCut);
      task->AddTrackCut(electronLowerHalf);
      
      AliReducedCompositeCut* electronUpperHalf = new AliReducedCompositeCut("electronUpperHalf", "", kTRUE);
      electronUpperHalf->AddCut(basicElectronSelection);
      AliReducedTrackCut* electronUpperHalfVarCut = new AliReducedTrackCut();
      if(isMC || !isCalibrated) electronUpperHalfVarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, 0.0, 3.0);
      else       electronUpperHalfVarCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 0.0, 3.0);
      electronUpperHalf->AddCut(electronUpperHalfVarCut);
      task->AddTrackCut(electronUpperHalf);
   }
   
   const Int_t kNsystEleVariations = 5;
   Float_t systEleVariations[kNsystEleVariations][2] = {
      {-2.0, 3.0}, {-1.5,3.0}, {-2.5,3.0}, {-2.0, 2.0}, {-2.0, 4.0}
   };
   const Int_t kNsystProtVariations = 3;
   Float_t systProtVariations[kNsystProtVariations] = {3.0, 2.5, 3.5};
   const Int_t kNsystPionVariations = 3;
   Float_t systPionVariations[kNsystPionVariations] = {3.0, 2.5, 3.5};
   
   for(Int_t iele=0; iele<kNsystEleVariations; ++iele) {
      for(Int_t ipro=0; ipro<kNsystProtVariations; ++ipro) {
         for(Int_t ipio=0; ipio<kNsystPionVariations; ++ipio) {
            
            if(addStandardCut && iele==0 && ipro==0 && ipio==0) continue;
                        
            TString cutName="Standard";
            if(iele!=0 || ipro!=0 || ipio!=0)
               cutName = Form("electron%d_prot%d_pion%d", iele, ipro, ipio);
            
            // for the single leg study with conversions we vary just one of the cuts, while keeping the others at their default
            Int_t nZeroes = 0;
            if(iele==0) nZeroes++;
            if(ipro==0) nZeroes++;
            if(ipio==0) nZeroes++;
            if(useWithConversions && nZeroes<2) continue;  
            
            AliReducedCompositeCut* cut = new AliReducedCompositeCut(cutName.Data(), "", kTRUE);
            cut->AddCut(commonCuts);
            AliReducedTrackCut* varCut = new AliReducedTrackCut(Form("%s_varCut", cutName.Data()), "");
            if(task->GetRunOverMC() || !isCalibrated)
               varCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, systEleVariations[iele][0], systEleVariations[iele][1]);
            else
               varCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, systEleVariations[iele][0], systEleVariations[iele][1]);
            
            varCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, systProtVariations[ipro], 30000.0);
            varCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, systPionVariations[ipio], 30000.0);
            cut->AddCut(varCut);
            
            task->AddTrackCut(cut);
         }  // end loop over pion rej settings
      }  // end loop over proton rej settings
   }  // end loop over electron inclusion setting
}

//__________________________________________________________________________________________
void AddStandardCutForComparisons(AliReducedAnalysisJpsi2ee* task, Bool_t isCalibrated = kFALSE) {
   //
   //  Here we add to the task a cut that is commonly agreed between myself and the other analyzers of the 5 TeV pp dataset
   //
   
   Bool_t isMC = task->GetRunOverMC();
   
   // standard cut
   AliReducedTrackCut* standardCuts = new AliReducedTrackCut("StandardCommon","");
   standardCuts->SetTrackFilterBit(kBasicCut);        // eta cut, DCA cuts, ITS and TPC refit already applied
   standardCuts->AddCut(AliReducedVarManager::kPt, 0.0, 1.0, kTRUE);  
   if(isMC || !isCalibrated) standardCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2.0, 3.0);
   else                      standardCuts->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -2.0, 3.0);
   standardCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.0, 30000.0);
   standardCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.0, 30000.0);
   standardCuts->AddCut(AliReducedVarManager::kITSchi2, 0.0, 10.);
   standardCuts->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
   standardCuts->SetRejectKinks();
   standardCuts->SetRequestSPDany();
   
   task->AddTrackCut(standardCuts);
}

//__________________________________________________________________________________________
void AddStandardCut(AliReducedAnalysisJpsi2ee* task, Bool_t useSPDany=kTRUE, Bool_t isCalibrated = kFALSE) {
   //
   //  This is the standard set of cuts
   //
   
   Bool_t isMC = task->GetRunOverMC();
   
   // standard cut
   AliReducedTrackCut* standardCuts = new AliReducedTrackCut("Standard","");
   standardCuts->SetTrackFilterBit(kBasicCut);        // eta cut, DCA cuts, ITS and TPC refit already applied
   standardCuts->AddCut(AliReducedVarManager::kPt, 0.0, 1.0, kTRUE);  
   if(isMC || !isCalibrated) standardCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2.0, 3.0);
   else                      standardCuts->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -2.0, 3.0);
   standardCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.0, 30000.0);
   standardCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.0, 30000.0);
   standardCuts->AddCut(AliReducedVarManager::kITSchi2, 0.0, 10.);
   standardCuts->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
   standardCuts->SetRejectKinks();
   standardCuts->SetRequestSPDany();
   standardCuts->AddCut(AliReducedVarManager::kITSnclsShared, 1.9, 6., kTRUE);
   standardCuts->AddCut(AliReducedVarManager::kDcaXY, -0.2, 0.2);
   standardCuts->AddCut(AliReducedVarManager::kDcaZ, -0.4, 0.4);
   
   task->AddTrackCut(standardCuts);
}


//__________________________________________________________________________________________
void SetupDCASystematicsCuts(AliReducedAnalysisJpsi2ee* task, Bool_t addStandardCut=kTRUE, Bool_t isCalibrated = kFALSE) {
   //
   // Apply DCA cut variations cuts variations on the electron candidates
   //
   Bool_t isMC = task->GetRunOverMC();
   
   // common cuts
   AliReducedTrackCut* commonCuts = new AliReducedTrackCut("commonCuts","");
   commonCuts->SetTrackFilterBit(kBasicCut);        // eta cut, DCA cuts, ITS and TPC refit already applied
   commonCuts->AddCut(AliReducedVarManager::kPt, 0.0, 1.0, kTRUE);  
   if(isMC || !isCalibrated) commonCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2.0, 3.0);
   else                      commonCuts->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -2.0, 3.0);
   commonCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.0, 30000.0);
   commonCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.0, 30000.0);
   commonCuts->AddCut(AliReducedVarManager::kITSchi2, 0.0, 10.);
   commonCuts->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
   commonCuts->SetRejectKinks();
   commonCuts->SetRequestSPDany();
   
   // DCAxy cut variations
   const Int_t kNDCAxyCuts = 6;
   Double_t dcaXYcuts[kNDCAxyCuts] = {1.0, 0.5, 0.3, 0.2, 0.15, 0.1};
   // DCAz cut variations
   const Int_t kNDCAzCuts = 5;
   Double_t dcaZcuts[kNDCAzCuts] = {3.0, 1.0, 0.6, 0.4, 0.2};
   
   if(addStandardCut) AddStandardCutForComparisons(task, isCalibrated);
   
   for(Int_t ixy=0;ixy<kNDCAxyCuts; ++ixy) {
      AliReducedTrackCut* xyCut = new AliReducedTrackCut(Form("DCAxy%.2f", dcaXYcuts[ixy]),"");
      xyCut->AddCut(AliReducedVarManager::kDcaXY, -1.0*dcaXYcuts[ixy], dcaXYcuts[ixy]);
      for(Int_t iz=0;iz<kNDCAzCuts; ++iz) {
         if(ixy==0 && iz==0) continue;
         AliReducedTrackCut* zCut = new AliReducedTrackCut(Form("DCAz%.2f", dcaZcuts[iz]),"");
         zCut->AddCut(AliReducedVarManager::kDcaZ, -1.0*dcaZcuts[iz], dcaZcuts[iz]);
         
         AliReducedCompositeCut* cut = new AliReducedCompositeCut(Form("DCAxy%.0fDCAz%.0f", dcaXYcuts[ixy]*100, dcaZcuts[iz]*100), "");
         cut->AddCut(commonCuts);
         cut->AddCut(xyCut);
         cut->AddCut(zCut);
         task->AddTrackCut(cut);
      }
   }
}

//__________________________________________________________________________________________
void SetupTrackingSystematicCuts(AliReducedAnalysisJpsi2ee* task, Bool_t addStandardCut=kTRUE, Bool_t isCalibrated = kFALSE) {
   //
   // Apply tracking cuts variations on the electron candidates
   //
   Bool_t isMC = task->GetRunOverMC();
   
   // common cuts
   AliReducedTrackCut* commonCuts = new AliReducedTrackCut("commonCuts","");
   commonCuts->SetTrackFilterBit(kBasicCut);        // eta cut, DCA cuts, ITS and TPC refit already applied
   commonCuts->AddCut(AliReducedVarManager::kPt, 0.0, 1.0, kTRUE);  
   if(isMC || !isCalibrated) commonCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2.0, 3.0);
   else                      commonCuts->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -2.0, 3.0);
   commonCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.0, 30000.0);
   commonCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.0, 30000.0);
   commonCuts->AddCut(AliReducedVarManager::kITSchi2, 0.0, 10.);
   commonCuts->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
   commonCuts->SetRejectKinks();
   
   if(addStandardCut) {
      AddStandardCut(task, kTRUE, isCalibrated);
      AddStandardCutForComparisons(task, isCalibrated);
   }
   
   const Int_t kNITSclsReqCuts = 5;      // SPDany, SPDfirst, SPDboth, ITS3cls
   AliReducedTrackCut* itsClsReqCuts[kNITSclsReqCuts];
   itsClsReqCuts[0] = new AliReducedTrackCut("SPDany","");
   itsClsReqCuts[0]->SetRequestSPDany();
   itsClsReqCuts[1] = new AliReducedTrackCut("SPDfirst","");
   itsClsReqCuts[1]->SetRequestSPDfirst();
   itsClsReqCuts[2] = new AliReducedTrackCut("SPDboth","");
   itsClsReqCuts[2]->SetRequestSPDboth();
   itsClsReqCuts[3] = new AliReducedTrackCut("ITS3cls","");
   itsClsReqCuts[3]->AddCut(AliReducedVarManager::kITSncls, 0., 2.1, kTRUE);
   itsClsReqCuts[4] = new AliReducedTrackCut("ITS4cls","");
   itsClsReqCuts[4]->AddCut(AliReducedVarManager::kITSncls, 0., 3.1, kTRUE);
   
   const Int_t kNITSsharedClsCuts = 4;      // ITSshared0, itsShared1, itsShared2, ITSsharedAny
   AliReducedTrackCut* itsSharedClsCuts[kNITSsharedClsCuts];
   itsSharedClsCuts[1] = new AliReducedTrackCut("ITSshared0","");
   itsSharedClsCuts[1]->AddCut(AliReducedVarManager::kITSnclsShared, 0.9, 6., kTRUE);   // no shared cls allowed
   itsSharedClsCuts[0] = new AliReducedTrackCut("ITSshared1","");
   itsSharedClsCuts[0]->AddCut(AliReducedVarManager::kITSnclsShared, 1.9, 6., kTRUE);   // 1 shared cls allowed
   itsSharedClsCuts[2] = new AliReducedTrackCut("ITSshared2","");
   itsSharedClsCuts[2]->AddCut(AliReducedVarManager::kITSnclsShared, 2.9, 6., kTRUE);   // 2 shared cls allowed
   itsSharedClsCuts[3] = new AliReducedTrackCut("ITSshared3","");
   itsSharedClsCuts[3]->AddCut(AliReducedVarManager::kITSnclsShared, 0., 6.);   // no shared cls cut
   
   const Int_t kNDCACuts = 3;
   AliReducedTrackCut* dcaCuts[kNDCACuts];
   dcaCuts[0] = new AliReducedTrackCut("DCAxy20DCAz40","");
   dcaCuts[0]->AddCut(AliReducedVarManager::kDcaXY, -0.2, 0.2);
   dcaCuts[0]->AddCut(AliReducedVarManager::kDcaZ, -0.4, 0.4);
   dcaCuts[1] = new AliReducedTrackCut("DCAxy15DCAz30","");
   dcaCuts[1]->AddCut(AliReducedVarManager::kDcaXY, -0.15, 0.15);
   dcaCuts[1]->AddCut(AliReducedVarManager::kDcaZ, -0.3, 0.3);
   dcaCuts[2] = new AliReducedTrackCut("DCAxy50DCAz100","");
   dcaCuts[2]->AddCut(AliReducedVarManager::kDcaXY, -0.5, 0.5);
   dcaCuts[2]->AddCut(AliReducedVarManager::kDcaZ, -1.0, 1.0);
 
   for(Int_t i=0;i<kNITSclsReqCuts;++i) {
      for(Int_t j=0;j<kNITSsharedClsCuts;++j) {
         for(Int_t k=0;k<kNDCACuts;++k) {
            if(addStandardCut && i==0 && j==0 && k==0) continue;     // this is the standard cut and was already added
            
            AliReducedCompositeCut* cut = new AliReducedCompositeCut(Form("%s_%s_%s", itsClsReqCuts[i]->GetName(), itsSharedClsCuts[j]->GetName(), dcaCuts[k]->GetName()), "");
            cut->AddCut(commonCuts);
            cut->AddCut(itsClsReqCuts[i]);
            cut->AddCut(itsSharedClsCuts[j]);
            cut->AddCut(dcaCuts[k]);
            task->AddTrackCut(cut);
         }
      }
   }
}

//__________________________________________________________________________________________
void SetupTrackingSystematicsCutsSingleLeg(AliReducedAnalysisJpsi2ee* task, Bool_t isCalibrated = kFALSE) {
   //
   // Apply tracking cuts variations on the electron candidates
   //
   Bool_t isMC = task->GetRunOverMC();
   
   // reference cut
   AliReducedTrackCut* referenceCut = new AliReducedTrackCut("referenceCut","");
   referenceCut->SetTrackFilterBit(kBasicCut);
   referenceCut->AddCut(AliReducedVarManager::kPt, 0.0, 1.0, kTRUE);  
   if(isMC || !isCalibrated) referenceCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -1.5, 1.5);
   else                      referenceCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -1.5, 1.5);
   referenceCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.0, 30000.0);
   referenceCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.0, 30000.0);
   task->AddTrackCut(referenceCut);
   
   // kink rejection
   AliReducedCompositeCut* kinkRej = new AliReducedCompositeCut("kinkRej", "");
   kinkRej->AddCut(referenceCut);
   AliReducedTrackCut* kinkRejSpecial = new AliReducedTrackCut("kinkRejSpecial", "");
   kinkRejSpecial->SetRejectKinks();
   kinkRej->AddCut(kinkRejSpecial);
   task->AddTrackCut(kinkRej);
   
   // ITS cluster cuts
   AliReducedTrackCut* spdFirstSpecial = new AliReducedTrackCut("spdFirstSpecial","");
   spdFirstSpecial->SetRequestSPDfirst();
   AliReducedCompositeCut* spdFirst = new AliReducedCompositeCut("spdFirst", "");
   spdFirst->AddCut(referenceCut);
   spdFirst->AddCut(spdFirstSpecial);
   task->AddTrackCut(spdFirst);
   
   AliReducedTrackCut* spdAnySpecial = new AliReducedTrackCut("spdAnySpecial","");
   spdAnySpecial->SetRequestSPDany();
   AliReducedCompositeCut* spdAny = new AliReducedCompositeCut("spdAny", "");
   spdAny->AddCut(referenceCut);
   spdAny->AddCut(spdAnySpecial);
   task->AddTrackCut(spdAny);
   
   AliReducedTrackCut* its3ClsSpecial = new AliReducedTrackCut("its3ClsSpecial","");
   its3ClsSpecial->AddCut(AliReducedVarManager::kITSncls, 0., 2.1, kTRUE);
   AliReducedCompositeCut* its3Cls = new AliReducedCompositeCut("its3Cls", "");
   its3Cls->AddCut(referenceCut);
   its3Cls->AddCut(its3ClsSpecial);
   task->AddTrackCut(its3Cls);
   
   AliReducedTrackCut* its4ClsSpecial = new AliReducedTrackCut("its4ClsSpecial","");
   its4ClsSpecial->AddCut(AliReducedVarManager::kITSncls, 0., 3.1, kTRUE);
   AliReducedCompositeCut* its4Cls = new AliReducedCompositeCut("its4Cls", "");
   its4Cls->AddCut(referenceCut);
   its4Cls->AddCut(its4ClsSpecial);
   task->AddTrackCut(its4Cls);
   
   // ITS shared clusters
   AliReducedTrackCut* itsSharedCluster3Special = new AliReducedTrackCut("itsSharedCluster3Special","");
   itsSharedCluster3Special->AddCut(AliReducedVarManager::kITSnclsShared, 2.9, 6., kTRUE);
   AliReducedCompositeCut* itsSharedCluster3 = new AliReducedCompositeCut("itsSharedCluster3", "");
   itsSharedCluster3->AddCut(referenceCut);
   itsSharedCluster3->AddCut(itsSharedCluster3Special);
   task->AddTrackCut(itsSharedCluster3);
   
   AliReducedTrackCut* itsSharedCluster2Special = new AliReducedTrackCut("itsSharedCluster2Special","");
   itsSharedCluster2Special->AddCut(AliReducedVarManager::kITSnclsShared, 1.9, 6., kTRUE);
   AliReducedCompositeCut* itsSharedCluster2 = new AliReducedCompositeCut("itsSharedCluster2", "");
   itsSharedCluster2->AddCut(referenceCut);
   itsSharedCluster2->AddCut(itsSharedCluster2Special);
   task->AddTrackCut(itsSharedCluster2);
   
   AliReducedTrackCut* itsSharedCluster1Special = new AliReducedTrackCut("itsSharedCluster1Special","");
   itsSharedCluster1Special->AddCut(AliReducedVarManager::kITSnclsShared, 0.9, 6., kTRUE);
   AliReducedCompositeCut* itsSharedCluster1 = new AliReducedCompositeCut("itsSharedCluster1", "");
   itsSharedCluster1->AddCut(referenceCut);
   itsSharedCluster1->AddCut(itsSharedCluster1Special);
   task->AddTrackCut(itsSharedCluster1);
   
   // ITS chi2
   AliReducedTrackCut* itsChi2_30Special = new AliReducedTrackCut("itsChi2_30Special","");
   itsChi2_30Special->AddCut(AliReducedVarManager::kITSchi2, 0.0, 30.);
   AliReducedCompositeCut* itsChi2_30 = new AliReducedCompositeCut("itsChi2_30", "");
   itsChi2_30->AddCut(referenceCut);
   itsChi2_30->AddCut(itsChi2_30Special);
   task->AddTrackCut(itsChi2_30);
   
   AliReducedTrackCut* itsChi2_20Special = new AliReducedTrackCut("itsChi2_20Special","");
   itsChi2_20Special->AddCut(AliReducedVarManager::kITSchi2, 0.0, 20.);
   AliReducedCompositeCut* itsChi2_20 = new AliReducedCompositeCut("itsChi2_20", "");
   itsChi2_20->AddCut(referenceCut);
   itsChi2_20->AddCut(itsChi2_20Special);
   task->AddTrackCut(itsChi2_20);
   
   AliReducedTrackCut* itsChi2_15Special = new AliReducedTrackCut("itsChi2_15Special","");
   itsChi2_15Special->AddCut(AliReducedVarManager::kITSchi2, 0.0, 15.);
   AliReducedCompositeCut* itsChi2_15 = new AliReducedCompositeCut("itsChi2_15", "");
   itsChi2_15->AddCut(referenceCut);
   itsChi2_15->AddCut(itsChi2_15Special);
   task->AddTrackCut(itsChi2_15);
   
   AliReducedTrackCut* itsChi2_10Special = new AliReducedTrackCut("itsChi2_10Special","");
   itsChi2_10Special->AddCut(AliReducedVarManager::kITSchi2, 0.0, 10.);
   AliReducedCompositeCut* itsChi2_10 = new AliReducedCompositeCut("itsChi2_10", "");
   itsChi2_10->AddCut(referenceCut);
   itsChi2_10->AddCut(itsChi2_10Special);
   task->AddTrackCut(itsChi2_10);
   
   // TPC number of clusters
   AliReducedTrackCut* tpcNcls90Special = new AliReducedTrackCut("tpcNcls90Special","");
   tpcNcls90Special->AddCut(AliReducedVarManager::kTPCncls, 90., 160.);
   AliReducedCompositeCut* tpcNcls90 = new AliReducedCompositeCut("tpcNcls90", "");
   tpcNcls90->AddCut(referenceCut);
   tpcNcls90->AddCut(tpcNcls90Special);
   task->AddTrackCut(tpcNcls90);
   
   AliReducedTrackCut* tpcNcls110Special = new AliReducedTrackCut("tpcNcls110Special","");
   tpcNcls110Special->AddCut(AliReducedVarManager::kTPCncls, 110., 160.);
   AliReducedCompositeCut* tpcNcls110 = new AliReducedCompositeCut("tpcNcls110", "");
   tpcNcls110->AddCut(referenceCut);
   tpcNcls110->AddCut(tpcNcls110Special);
   task->AddTrackCut(tpcNcls110);
   
   AliReducedTrackCut* tpcNcls120Special = new AliReducedTrackCut("tpcNcls120Special","");
   tpcNcls120Special->AddCut(AliReducedVarManager::kTPCncls, 120., 160.);
   AliReducedCompositeCut* tpcNcls120 = new AliReducedCompositeCut("tpcNcls120", "");
   tpcNcls120->AddCut(referenceCut);
   tpcNcls120->AddCut(tpcNcls120Special);
   task->AddTrackCut(tpcNcls120);
   
   AliReducedTrackCut* tpcNcls130Special = new AliReducedTrackCut("tpcNcls130Special","");
   tpcNcls130Special->AddCut(AliReducedVarManager::kTPCncls, 130., 160.);
   AliReducedCompositeCut* tpcNcls130 = new AliReducedCompositeCut("tpcNcls130", "");
   tpcNcls130->AddCut(referenceCut);
   tpcNcls130->AddCut(tpcNcls130Special);
   task->AddTrackCut(tpcNcls130);
   
   AliReducedTrackCut* tpcNcls140Special = new AliReducedTrackCut("tpcNcls140Special","");
   tpcNcls140Special->AddCut(AliReducedVarManager::kTPCncls, 140., 160.);
   AliReducedCompositeCut* tpcNcls140 = new AliReducedCompositeCut("tpcNcls140", "");
   tpcNcls140->AddCut(referenceCut);
   tpcNcls140->AddCut(tpcNcls140Special);
   task->AddTrackCut(tpcNcls140);
   
   // TPC number of segments
   AliReducedTrackCut* tpcNsegm4Special = new AliReducedTrackCut("tpcNsegm4Special","");
   tpcNsegm4Special->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 4., 9.);
   AliReducedCompositeCut* tpcNsegm4 = new AliReducedCompositeCut("tpcNsegm4", "");
   tpcNsegm4->AddCut(referenceCut);
   tpcNsegm4->AddCut(tpcNsegm4Special);
   task->AddTrackCut(tpcNsegm4);
   
   AliReducedTrackCut* tpcNsegm5Special = new AliReducedTrackCut("tpcNsegm5Special","");
   tpcNsegm5Special->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 5., 9.);
   AliReducedCompositeCut* tpcNsegm5 = new AliReducedCompositeCut("tpcNsegm5", "");
   tpcNsegm5->AddCut(referenceCut);
   tpcNsegm5->AddCut(tpcNsegm5Special);
   task->AddTrackCut(tpcNsegm5);
   
   AliReducedTrackCut* tpcNsegm6Special = new AliReducedTrackCut("tpcNsegm6Special","");
   tpcNsegm6Special->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 6., 9.);
   AliReducedCompositeCut* tpcNsegm6 = new AliReducedCompositeCut("tpcNsegm6", "");
   tpcNsegm6->AddCut(referenceCut);
   tpcNsegm6->AddCut(tpcNsegm6Special);
   task->AddTrackCut(tpcNsegm6);
   
   AliReducedTrackCut* tpcNsegm7Special = new AliReducedTrackCut("tpcNsegm7Special","");
   tpcNsegm7Special->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 7., 9.);
   AliReducedCompositeCut* tpcNsegm7 = new AliReducedCompositeCut("tpcNsegm7", "");
   tpcNsegm7->AddCut(referenceCut);
   tpcNsegm7->AddCut(tpcNsegm7Special);
   task->AddTrackCut(tpcNsegm7);
   
   AliReducedTrackCut* tpcNsegm8Special = new AliReducedTrackCut("tpcNsegm8Special","");
   tpcNsegm8Special->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 8., 9.);
   AliReducedCompositeCut* tpcNsegm8 = new AliReducedCompositeCut("tpcNsegm8", "");
   tpcNsegm8->AddCut(referenceCut);
   tpcNsegm8->AddCut(tpcNsegm8Special);
   task->AddTrackCut(tpcNsegm8);
   
   // TPC chi2
   AliReducedTrackCut* tpcChi2_4Special = new AliReducedTrackCut("tpcChi2_4Special","");
   tpcChi2_4Special->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
   AliReducedCompositeCut* tpcChi2_4 = new AliReducedCompositeCut("tpcChi2_4", "");
   tpcChi2_4->AddCut(referenceCut);
   tpcChi2_4->AddCut(tpcChi2_4Special);
   task->AddTrackCut(tpcChi2_4);
   
   AliReducedTrackCut* tpcChi2_35Special = new AliReducedTrackCut("tpcChi2_35Special","");
   tpcChi2_35Special->AddCut(AliReducedVarManager::kTPCchi2, 0., 3.5);
   AliReducedCompositeCut* tpcChi2_35 = new AliReducedCompositeCut("tpcChi2_35", "");
   tpcChi2_35->AddCut(referenceCut);
   tpcChi2_35->AddCut(tpcChi2_35Special);
   task->AddTrackCut(tpcChi2_35);
   
   AliReducedTrackCut* tpcChi2_3Special = new AliReducedTrackCut("tpcChi2_3Special","");
   tpcChi2_3Special->AddCut(AliReducedVarManager::kTPCchi2, 0., 3.);
   AliReducedCompositeCut* tpcChi2_3 = new AliReducedCompositeCut("tpcChi2_3", "");
   tpcChi2_3->AddCut(referenceCut);
   tpcChi2_3->AddCut(tpcChi2_3Special);
   task->AddTrackCut(tpcChi2_3);
   
   AliReducedTrackCut* tpcChi2_25Special = new AliReducedTrackCut("tpcChi2_25Special","");
   tpcChi2_25Special->AddCut(AliReducedVarManager::kTPCchi2, 0., 2.5);
   AliReducedCompositeCut* tpcChi2_25 = new AliReducedCompositeCut("tpcChi2_25", "");
   tpcChi2_25->AddCut(referenceCut);
   tpcChi2_25->AddCut(tpcChi2_25Special);
   task->AddTrackCut(tpcChi2_25);
   
   AliReducedTrackCut* tpcChi2_2Special = new AliReducedTrackCut("tpcChi2_2Special","");
   tpcChi2_2Special->AddCut(AliReducedVarManager::kTPCchi2, 0., 2.);
   AliReducedCompositeCut* tpcChi2_2 = new AliReducedCompositeCut("tpcChi2_2", "");
   tpcChi2_2->AddCut(referenceCut);
   tpcChi2_2->AddCut(tpcChi2_2Special);
   task->AddTrackCut(tpcChi2_2);
   
   // TPC golden chi2
   AliReducedTrackCut* goldenChi2_256Special = new AliReducedTrackCut("goldenChi2_256Special","");
   goldenChi2_256Special->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 0., 256.);
   AliReducedCompositeCut* goldenChi2_256 = new AliReducedCompositeCut("goldenChi2_256", "");
   goldenChi2_256->AddCut(referenceCut);
   goldenChi2_256->AddCut(goldenChi2_256Special);
   task->AddTrackCut(goldenChi2_256);
   
   AliReducedTrackCut* goldenChi2_196Special = new AliReducedTrackCut("goldenChi2_196Special","");
   goldenChi2_196Special->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 0., 196.);
   AliReducedCompositeCut* goldenChi2_196 = new AliReducedCompositeCut("goldenChi2_196", "");
   goldenChi2_196->AddCut(referenceCut);
   goldenChi2_196->AddCut(goldenChi2_196Special);
   task->AddTrackCut(goldenChi2_196);
   
   AliReducedTrackCut* goldenChi2_144Special = new AliReducedTrackCut("goldenChi2_144Special","");
   goldenChi2_144Special->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 0., 144.);
   AliReducedCompositeCut* goldenChi2_144 = new AliReducedCompositeCut("goldenChi2_144", "");
   goldenChi2_144->AddCut(referenceCut);
   goldenChi2_144->AddCut(goldenChi2_144Special);
   task->AddTrackCut(goldenChi2_144);
   
   AliReducedTrackCut* goldenChi2_100Special = new AliReducedTrackCut("goldenChi2_100Special","");
   goldenChi2_100Special->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 0., 100.);
   AliReducedCompositeCut* goldenChi2_100 = new AliReducedCompositeCut("goldenChi2_100", "");
   goldenChi2_100->AddCut(referenceCut);
   goldenChi2_100->AddCut(goldenChi2_100Special);
   task->AddTrackCut(goldenChi2_100);
   
   AliReducedTrackCut* goldenChi2_81Special = new AliReducedTrackCut("goldenChi2_81Special","");
   goldenChi2_81Special->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 0., 81.);
   AliReducedCompositeCut* goldenChi2_81 = new AliReducedCompositeCut("goldenChi2_81", "");
   goldenChi2_81->AddCut(referenceCut);
   goldenChi2_81->AddCut(goldenChi2_81Special);
   task->AddTrackCut(goldenChi2_81);
   
   AliReducedTrackCut* goldenChi2_64Special = new AliReducedTrackCut("goldenChi2_64Special","");
   goldenChi2_64Special->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 0., 64.);
   AliReducedCompositeCut* goldenChi2_64 = new AliReducedCompositeCut("goldenChi2_64", "");
   goldenChi2_64->AddCut(referenceCut);
   goldenChi2_64->AddCut(goldenChi2_64Special);
   task->AddTrackCut(goldenChi2_64);
   
   AliReducedTrackCut* goldenChi2_49Special = new AliReducedTrackCut("goldenChi2_49Special","");
   goldenChi2_49Special->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 0., 49.);
   AliReducedCompositeCut* goldenChi2_49 = new AliReducedCompositeCut("goldenChi2_49", "");
   goldenChi2_49->AddCut(referenceCut);
   goldenChi2_49->AddCut(goldenChi2_49Special);
   task->AddTrackCut(goldenChi2_49);
   
   AliReducedTrackCut* goldenChi2_36Special = new AliReducedTrackCut("goldenChi2_36Special","");
   goldenChi2_36Special->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 0., 36.);
   AliReducedCompositeCut* goldenChi2_36 = new AliReducedCompositeCut("goldenChi2_36", "");
   goldenChi2_36->AddCut(referenceCut);
   goldenChi2_36->AddCut(goldenChi2_36Special);
   task->AddTrackCut(goldenChi2_36);
   
   AliReducedTrackCut* goldenChi2_25Special = new AliReducedTrackCut("goldenChi2_25Special","");
   goldenChi2_25Special->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 0., 25.);
   AliReducedCompositeCut* goldenChi2_25 = new AliReducedCompositeCut("goldenChi2_25", "");
   goldenChi2_25->AddCut(referenceCut);
   goldenChi2_25->AddCut(goldenChi2_25Special);
   task->AddTrackCut(goldenChi2_25);
   
   AliReducedTrackCut* goldenChi2_16Special = new AliReducedTrackCut("goldenChi2_16Special","");
   goldenChi2_16Special->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 0., 16.);
   AliReducedCompositeCut* goldenChi2_16 = new AliReducedCompositeCut("goldenChi2_16", "");
   goldenChi2_16->AddCut(referenceCut);
   goldenChi2_16->AddCut(goldenChi2_16Special);
   task->AddTrackCut(goldenChi2_16);
   
   AliReducedTrackCut* goldenChi2_9Special = new AliReducedTrackCut("goldenChi2_9Special","");
   goldenChi2_9Special->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 0., 9.);
   AliReducedCompositeCut* goldenChi2_9 = new AliReducedCompositeCut("goldenChi2_9", "");
   goldenChi2_9->AddCut(referenceCut);
   goldenChi2_9->AddCut(goldenChi2_9Special);
   task->AddTrackCut(goldenChi2_9);

}

//_________________________________________________________________
void SetupMCsignals(AliReducedAnalysisJpsi2ee* task) {
   //
   // setup the needed signals for efficiency calculations
   //
   // MC truth info for reconstructed electron candidates
   AliReducedTrackCut* trueElectron = new AliReducedTrackCut("TrueElectron", "reconstructed electrons with MC truth");
   trueElectron->SetMCFilterBit(kJpsiDecayElectron);
   task->AddLegCandidateMCcut(trueElectron);  
   
   AliReducedTrackCut* trueElectronNonPrompt = new AliReducedTrackCut("TrueElectronNonPrompt", "reconstructed electrons with MC truth");
   trueElectronNonPrompt->SetMCFilterBit(kJpsiNonPromptDecayElectron);
   task->AddLegCandidateMCcut(trueElectronNonPrompt);
   
   AliReducedTrackCut* trueElectronPrompt = new AliReducedTrackCut("TrueElectronPrompt", "reconstructed electrons with MC truth");
   trueElectronPrompt->SetMCFilterBit(kJpsiPromptDecayElectron);
   task->AddLegCandidateMCcut(trueElectronPrompt);    
   
   // MC truth selections(s) for the J/psi mother and electron daughters
   AliReducedTrackCut* mcTruthJpsi = new AliReducedTrackCut("mcTruthJpsi", "Pure MC truth J/psi");
   mcTruthJpsi->SetMCFilterBit(kJpsiInclusive);
   mcTruthJpsi->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
   AliReducedTrackCut* mcTruthJpsiElectron = new AliReducedTrackCut("mcTruthJpsiElectron", "Pure MC truth electron from J/psi");
   mcTruthJpsiElectron->SetMCFilterBit(kJpsiDecayElectron);
   mcTruthJpsiElectron->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
   mcTruthJpsiElectron->AddCut(AliReducedVarManager::kPt, 1.0, 100.0);
   task->AddJpsiMotherMCCut(mcTruthJpsi, mcTruthJpsiElectron);  
   
   AliReducedTrackCut* mcTruthJpsiNonPrompt = new AliReducedTrackCut("mcTruthJpsiNonPrompt", "Pure MC truth J/psi");
   mcTruthJpsiNonPrompt->SetMCFilterBit(kJpsiNonPrompt);
   mcTruthJpsiNonPrompt->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
   AliReducedTrackCut* mcTruthJpsiElectronNonPrompt = new AliReducedTrackCut("mcTruthJpsiElectronNonPrompt", "Pure MC truth electron from J/psi");
   mcTruthJpsiElectronNonPrompt->SetMCFilterBit(kJpsiNonPromptDecayElectron);
   mcTruthJpsiElectronNonPrompt->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
   mcTruthJpsiElectronNonPrompt->AddCut(AliReducedVarManager::kPt, 1.0, 100.0);
   task->AddJpsiMotherMCCut(mcTruthJpsiNonPrompt, mcTruthJpsiElectronNonPrompt);
   
   AliReducedTrackCut* mcTruthJpsiPrompt = new AliReducedTrackCut("mcTruthJpsiPrompt", "Pure MC truth J/psi");
   mcTruthJpsiPrompt->SetMCFilterBit(kJpsiPrompt);
   mcTruthJpsiPrompt->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
   AliReducedTrackCut* mcTruthJpsiElectronPrompt = new AliReducedTrackCut("mcTruthJpsiElectronPrompt", "Pure MC truth electron from J/psi");
   mcTruthJpsiElectronPrompt->SetMCFilterBit(kJpsiPromptDecayElectron);
   mcTruthJpsiElectronPrompt->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
   mcTruthJpsiElectronPrompt->AddCut(AliReducedVarManager::kPt, 1.0, 100.0);
   task->AddJpsiMotherMCCut(mcTruthJpsiPrompt, mcTruthJpsiElectronPrompt);    
}


//_________________________________________________________________
void SetupMixingHandler(AliReducedAnalysisJpsi2ee* task) {
   //
   // setup the mixing handler
   //
   AliMixingHandler* handler = task->GetMixingHandler();
   handler->SetPoolDepth(100);
   handler->SetMixingThreshold(1.0);
   handler->SetDownscaleEvents(1);
   handler->SetDownscaleTracks(1);
   
   
   handler->AddMixingVariable(AliReducedVarManager::kVtxZ, gkNZBins, gZBinLims);
   //handler->AddMixingVariable(AliReducedVarManager::kSPDntracklets, gkNSPDtrkBins, gSPDtrkBinLims);
}


//_________________________________________________________________
void SetupHistogramManager(AliReducedAnalysisJpsi2ee* task, TString prod /*="LHC10h"*/) {
  //
  // setup the histograms manager
  //
  AliReducedVarManager::SetDefaultVarNames();
  
  DefineHistograms(task, prod);
  
  AliReducedVarManager::SetUseVars(task->GetHistogramManager()->GetUsedVars());
}


//_________________________________________________________________
void DefineHistograms(AliReducedAnalysisJpsi2ee* task, TString prod /*="LHC10h"*/) {
  //
  // define histograms
  // NOTE: The DefineHistograms need to be called after the track cuts are defined because the name of the cuts
  //           are used in the histogram lists
   // NOTE: The name of the pair histogram lists for event mixing need to contain "PairMEPP", "PairMEPM", "PairMEMM", in their name, see below.
   //  TODO: make needed changes such that this becomes less prone to mistakes
   
  AliHistogramManager* man = task->GetHistogramManager(); 
   
  TString histClasses = "";
  histClasses += "Event_BeforeCuts;";
  histClasses += "Event_AfterCuts;";
  if(!task->GetRunOverMC()) {
    histClasses += "EventTag_BeforeCuts;";   
    histClasses += "EventTag_AfterCuts;";   
    histClasses += "EventTriggers_BeforeCuts;";
    histClasses += "EventTriggers_AfterCuts;";   
  }
  if(task->GetLoopOverTracks())
    histClasses += "Track_BeforeCuts;";
  if(!task->GetRunOverMC()) {
    histClasses += "TrackStatusFlags_BeforeCuts;";
  }
  if(task->GetLoopOverTracks()) {
    histClasses += "TrackITSclusterMap_BeforeCuts;";
    histClasses += "TrackITSsharedClusterMap_BeforeCuts;";
    histClasses += "TrackTPCclusterMap_BeforeCuts;";
  }
  if(task->GetRunOverMC()) {
     for(Int_t mcSel=0; mcSel<task->GetNJpsiMotherMCCuts(); ++mcSel) {
        histClasses += Form("%s_PureMCTruth_BeforeSelection;", task->GetJpsiMotherMCcutName(mcSel));
        histClasses += Form("%s_PureMCTruth_AfterSelection;", task->GetJpsiMotherMCcutName(mcSel));
     }
  }
  if(task->GetLoopOverTracks()) {
    for(Int_t i=0; i<task->GetNTrackCuts(); ++i) {
      TString cutName = task->GetTrackCutName(i);
      histClasses += Form("Track_%s;", cutName.Data());
      if(!task->GetRunOverMC()) {
        histClasses += Form("TrackStatusFlags_%s;", cutName.Data());
      }
      histClasses += Form("TrackITSclusterMap_%s;", cutName.Data());
      histClasses += Form("TrackITSsharedClusterMap_%s;", cutName.Data());
      histClasses += Form("TrackTPCclusterMap_%s;", cutName.Data());
      if(task->GetRunOverMC()) {
         for(Int_t mcSel=0; mcSel<task->GetNLegCandidateMCcuts(); ++mcSel) {
            histClasses += Form("Track_%s_%s;", cutName.Data(), task->GetLegCandidateMCcutName(mcSel));
            histClasses += Form("TrackStatusFlags_%s_%s;", cutName.Data(), task->GetLegCandidateMCcutName(mcSel));
            histClasses += Form("TrackITSclusterMap_%s_%s;", cutName.Data(), task->GetLegCandidateMCcutName(mcSel));
            histClasses += Form("TrackITSsharedClusterMap_%s_%s;", cutName.Data(), task->GetLegCandidateMCcutName(mcSel));
            histClasses += Form("TrackTPCclusterMap_%s_%s;", cutName.Data(), task->GetLegCandidateMCcutName(mcSel));
         }
      }
      if(task->GetRunPairing()) {
         if(!task->GetRunOverMC()) {
            if(!task->GetRunLikeSignPairing())
               histClasses += Form("PairSEPM_%s;", cutName.Data());
            else
               histClasses += Form("PairSEPP_%s;PairSEPM_%s;PairSEMM_%s;", cutName.Data(), cutName.Data(), cutName.Data());
         }
         if(!task->GetRunOverMC())
            histClasses += Form("PairMEPP_%s;PairMEPM_%s;PairMEMM_%s;", cutName.Data(), cutName.Data(), cutName.Data());
         if(task->GetRunOverMC()) {
            for(Int_t mcSel=0; mcSel<task->GetNLegCandidateMCcuts(); ++mcSel)
               histClasses += Form("PairSEPM_%s_%s;", cutName.Data(), task->GetLegCandidateMCcutName(mcSel));
         }
      }
    }
  }
  
  AliReducedVarManager::SetRunNumbers(gRunList);
  
  const Int_t kNBCBins = 3600;
  Double_t bcHistRange[2] = {-0.5,3599.5};
  
  // set the mass bin limits
  for(Int_t i=0; i<126;++i) gMassBins[i] = 0.0+i*0.04;
  
  // set the pt bin limits for MC
  for(Int_t i=0; i<=gkNPtBinsMC;++i) gPtLimsMC[i] = 0.0+i*0.05;
  
  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");
  
  cout << "Histogram classes included in the Histogram Manager" << endl;
  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    
    if(classStr.Contains("PureMCTruth_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       
  //     man->AddHistogram(classStr.Data(), "MassMC_Pt_CentVZERO", "", kFALSE, gkNMassBins-1, gMassBins, AliReducedVarManager::kMassMC, gkNPtBins-1, gPtLims, AliReducedVarManager::kPtMC, gkNCentBins-1, gCentLims, AliReducedVarManager::kCentVZERO);
       
       man->AddHistogram(classStr.Data(), "MassMC", "MC mass", kFALSE, 200, 0., 5.0, AliReducedVarManager::kMassMC);
       man->AddHistogram(classStr.Data(), "RapidityMC", "MC rapidity", kFALSE, 48, -1.2, 1.2, AliReducedVarManager::kRapMC);
       man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE, gkNPtBinsMC, gPtLimsMC, AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "PtMC_RunID", "p_{T} MC vs Run", kFALSE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 
                                                                           600, 0., 30.0, AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "PtMC_coarse", "p_{T} MC", kFALSE, 20, 0., 20.0, AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "PhiMC", "MC #varphi", kFALSE, 100, 0., 6.3, AliReducedVarManager::kPhiMC);
       man->AddHistogram(classStr.Data(), "EtaMC", "MC #eta", kFALSE, 100, -1.5, 1.5, AliReducedVarManager::kEtaMC);
       man->AddHistogram(classStr.Data(), "PtMC_RapMC", "", kFALSE, 100, -1.2, 1.2, AliReducedVarManager::kRapMC, 100, 0., 15., AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "CosThetaStarCS", "cos(#theta^{*})_{CS}", kFALSE, 88, -1.1, 1.1, AliReducedVarManager::kPairThetaCS);
       man->AddHistogram(classStr.Data(), "CosThetaStarCS_ptMC", "cos(#theta^{*})_{CS} vs MC pt", kFALSE, 88, -1.1, 1.1, AliReducedVarManager::kPairThetaCS, 400, 0., 20., AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "CosThetaStarHE", "cos(#theta^{*})_{HE}", kFALSE, 88, -1.1, 1.1, AliReducedVarManager::kPairThetaHE);
       man->AddHistogram(classStr.Data(), "CosThetaStarHE_ptMC", "cos(#theta^{*})_{HE} vs MC pt", kFALSE, 88, -1.1, 1.1, AliReducedVarManager::kPairThetaHE, 400, 0., 20., AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "CosThetaStarHE_ptMC_coarse", "cos(#theta^{*})_{HE} vs MC pt", kFALSE, 10, -1.0, 1.0, AliReducedVarManager::kPairThetaHE, 20, 0., 20., AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "PhiStarCS", "#varphi^{*}_{CS}", kFALSE, 66, -3.3, 3.3, AliReducedVarManager::kPairPhiCS);
       man->AddHistogram(classStr.Data(), "PhiStarHE", "#varphi^{*}_{HE}", kFALSE, 66, -3.3, 3.3, AliReducedVarManager::kPairPhiHE);
       
       
       continue;
    }
    
    // Event wise histograms
    if(classStr.Contains("EventTag_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       TString tagNames = "";
       tagNames += "AliAnalysisUtils 2013 selection;";
       tagNames += "AliAnalysisUtils MV pileup;";
       tagNames += "AliAnalysisUtils MV pileup, no BC check;";
       tagNames += "AliAnalysisUtils MV pileup, min wght dist 10;";
       tagNames += "AliAnalysisUtils MV pileup, min wght dist 5;";
       tagNames += "IsPileupFromSPD(3,0.6,3.,2.,5.);";
       tagNames += "IsPileupFromSPD(4,0.6,3.,2.,5.);";
       tagNames += "IsPileupFromSPD(5,0.6,3.,2.,5.);";
       tagNames += "IsPileupFromSPD(6,0.6,3.,2.,5.);";
       tagNames += "IsPileupFromSPD(3,0.8,3.,2.,5.);";
       tagNames += "IsPileupFromSPD(4,0.8,3.,2.,5.);";
       tagNames += "IsPileupFromSPD(5,0.8,3.,2.,5.);";
       tagNames += "IsPileupFromSPD(6,0.8,3.,2.,5.);";
       tagNames += "vtx distance selected;";
       tagNames += "unbiased event;";
       man->AddHistogram(classStr.Data(), "EventTags", "Event tags", kFALSE,
                         20, -0.5, 19.5, AliReducedVarManager::kEventTag, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, tagNames.Data());
       man->AddHistogram(classStr.Data(), "EventTags_PhysSelection", "Event tags", kFALSE,
                         20, -0.5, 19.5, AliReducedVarManager::kEventTag, 2,-0.5,1.5,AliReducedVarManager::kIsPhysicsSelection, 0, 0.0, 0.0, AliReducedVarManager::kNothing, tagNames.Data(), "off;on");
       man->AddHistogram(classStr.Data(), "EventTags_SPDtracklets0_Run_prof", "<SPD tracklets> vs run", kTRUE,
                         gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 20, -0.5, 19.5, AliReducedVarManager::kEventTag, 300, 0., 300., AliReducedVarManager::kSPDntracklets, gRunList.Data(), tagNames.Data());
       man->AddHistogram(classStr.Data(), "EventTags_VZEROMult_Run_prof", "<VZERO mult> vs run", kTRUE,
                         gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 20, -0.5, 19.5, AliReducedVarManager::kEventTag, 300, 0., 300., AliReducedVarManager::kVZEROTotalMult, gRunList.Data(), tagNames.Data());
       man->AddHistogram(classStr.Data(), "EventTags_NTotalTracksAnalyzed_Run_prof", "<total tracks analyzed> vs run", kTRUE,
                         gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 20, -0.5, 19.5, AliReducedVarManager::kEventTag, 300, 0., 300., AliReducedVarManager::kNtracksAnalyzed, gRunList.Data(), tagNames.Data());
       man->AddHistogram(classStr.Data(), "EventTags_NTotalPairsAnalyzed_Run_prof", "<total tracks analyzed> vs run", kTRUE,
                         gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 20, -0.5, 19.5, AliReducedVarManager::kEventTag, 300, 0., 300., AliReducedVarManager::kNpairsSelected, gRunList.Data(), tagNames.Data());
       man->AddHistogram(classStr.Data(), "EventTags_VtxZ_Run", "VtxZ vs run", kFALSE,
                         gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 20, -0.5, 19.5, AliReducedVarManager::kEventTag, 300, -15., 15., AliReducedVarManager::kVtxZ, gRunList.Data(), tagNames.Data());
       //man->AddHistogram(classStr.Data(), "EventTags_CentVZERO", "Event tags vs VZERO centrality", kFALSE,
       //                  20, -0.5, 19.5, AliReducedVarManager::kEventTag, 9, 0.0, 90.0, AliReducedVarManager::kCentVZERO, 0, 0.0, 0.0, AliReducedVarManager::kNothing, tagNames.Data());
       continue;
    }
    
    if(classStr.Contains("EventTriggers_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       TString triggerNames = "";
       for(Int_t i=0; i<64; ++i) {triggerNames += AliReducedVarManager::fgkOfflineTriggerNames[i]; triggerNames+=";";}
       
       man->AddHistogram(classStr.Data(), "Triggers2", "", kFALSE,
                         64, -0.5, 63.5, AliReducedVarManager::kOnlineTriggerFired2, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, triggerNames.Data());
       //man->AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "", kFALSE,
       //                  64, -0.5, 63.5, AliReducedVarManager::kOnlineTriggerFired2, 9, 0.0, 90.0, AliReducedVarManager::kCentVZERO, 0, 0.0, 0.0, AliReducedVarManager::kNothing, triggerNames.Data());
       continue;
    }
    
    if(classStr.Contains("Event_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(),"RunNo","Run numbers",kFALSE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.,0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(),"VtxX","Vtx X",kFALSE,500,-1.0,1.0,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxX_vs_Run_prof", "Vtx <X>", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kVtxX, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(),"VtxY","Vtx Y",kFALSE,500,-1.0,1.0,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxY_vs_Run_prof", "Vtx <Y>", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kVtxY, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(),"VtxXY","Vtx X vs vtx Y",kFALSE,100,0.,0.2,AliReducedVarManager::kVtxX, 100,0.3,0.45,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxXZ","Vtx X vs vtx Z",kFALSE, 100,-10.,10.,AliReducedVarManager::kVtxZ, 100,0.0,0.2,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxYZ","Vtx Y vs vtx Z",kFALSE, 100,-10.,10.,AliReducedVarManager::kVtxZ, 100,0.3,0.45,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxZ","Vtx Z",kFALSE,200,-10.,10.,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"VtxZ_vs_Run_prof", "Vtx <Z>", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kVtxZ, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(),"DeltaVtxZtpc","Z - Z_{TPC}",kFALSE,300,-15.,15.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"DeltaVtxztpc_vs_Run_prof", "<Z-Z_{TPC}>", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kDeltaVtxZ, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(),"DeltaVtxZtpc_nTPCout","Z - Z_{TPC} vs nTPCout tracks",kFALSE,
                        100,0.,1200.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout,                        100,-15.,15.,AliReducedVarManager::kDeltaVtxZ);
      
      man->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0", kFALSE, 300, 0., 300., AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(),"SPDntracklets_vs_Run_prof", "SPD #tracklets in |#eta|<1.0", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kSPDntracklets, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(),"SPDfiredChipsL1", "SPD fired chips in layer 1", kFALSE, 300, 0., 300., AliReducedVarManager::kSPDFiredChips+0);
      man->AddHistogram(classStr.Data(),"SPDfiredChipsL2", "SPD fired chips in layer 2", kFALSE, 300, 0., 300., AliReducedVarManager::kSPDFiredChips+1);
      man->AddHistogram(classStr.Data(),"SPDnSigleClusters", "SPD #single clusters", kFALSE, 300, 0., 300., AliReducedVarManager::kSPDnSingleClusters);
      for(Int_t il=0; il<6;++il)
         man->AddHistogram(classStr.Data(),Form("ITSnClustersLayer%d",il+1), Form("Number of clusters in ITS layer %d", il+1), kFALSE, 300, 0., 300., AliReducedVarManager::kITSnClusters+il);
      
      man->AddHistogram(classStr.Data(),"SPDntracklets_SPDfiredChipsL1", "SPD #tracklets in |#eta|<1.0 vs SPD fired chips in L1", kFALSE, 200, 0., 200., AliReducedVarManager::kSPDntracklets, 200, 0., 200., AliReducedVarManager::kSPDFiredChips+0);
      man->AddHistogram(classStr.Data(),"SPDntracklets_SPDfiredChipsL2", "SPD #tracklets in |#eta|<1.0 vs SPD fired chips in L2", kFALSE, 200, 0., 200., AliReducedVarManager::kSPDntracklets, 200, 0., 200., AliReducedVarManager::kSPDFiredChips+1);
      man->AddHistogram(classStr.Data(),"SPDntracklets_SPDnSingleClusters", "SPD #tracklets in |#eta|<1.0 vs SPD single clusters", kFALSE, 200, 0., 200., AliReducedVarManager::kSPDntracklets, 200, 0., 200., AliReducedVarManager::kSPDnSingleClusters);
      man->AddHistogram(classStr.Data(),"SPDntracklets_ITSclsLayer1", "SPD #tracklets in |#eta|<1.0 vs ITS clusters in L1", kFALSE, 200, 0., 200., AliReducedVarManager::kSPDntracklets, 200, 0., 200., AliReducedVarManager::kITSnClusters+0);
      man->AddHistogram(classStr.Data(),"SPDntracklets_ITSclsLayer2", "SPD #tracklets in |#eta|<1.0 vs ITS clusters in L2", kFALSE, 200, 0., 200., AliReducedVarManager::kSPDntracklets, 200, 0., 200., AliReducedVarManager::kITSnClusters+1);
      man->AddHistogram(classStr.Data(),"SPDntracklets_ITSclsLayer5", "SPD #tracklets in |#eta|<1.0 vs ITS clusters in L5", kFALSE, 200, 0., 200., AliReducedVarManager::kSPDntracklets, 200, 0., 200., AliReducedVarManager::kITSnClusters+4);
      man->AddHistogram(classStr.Data(),"SPDntracklets_ITSclsLayer6", "SPD #tracklets in |#eta|<1.0 vs ITS clusters in L6", kFALSE, 200, 0., 200., AliReducedVarManager::kSPDntracklets, 200, 0., 200., AliReducedVarManager::kITSnClusters+5);
      
      
      man->AddHistogram(classStr.Data(),"VZEROmult", "", kFALSE, 200, 0.0, 2000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"VZEROmult_vs_Run_prof", "", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kVZEROTotalMult, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(),"VZEROmult_VtxContributors", "", kFALSE, 
		   200, 0.0, 2000., AliReducedVarManager::kVZEROTotalMult, 200, 0.0, 200., AliReducedVarManager::kNVtxContributors);
      man->AddHistogram(classStr.Data(),"VZEROmult_SPDntracklets", "", kFALSE, 
                        200, 0.0, 200., AliReducedVarManager::kSPDntracklets, 200, 0.0, 2000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsPhysicsSelection, 0,0.0,0.0,AliReducedVarManager::kNothing, 0,0.0,0.0,AliReducedVarManager::kNothing, "off;on");

      man->AddHistogram(classStr.Data(),"NPosTracksAnalyzed","#positrons per event",kFALSE,10,0.,10.,AliReducedVarManager::kNtracksPosAnalyzed);
      man->AddHistogram(classStr.Data(),"NPosTracksAnalyzed_vs_Run_prof", "# positrons per event", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kNtracksPosAnalyzed, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(),"NNegTracksAnalyzed","#electrons per event",kFALSE,10,0.,10.,AliReducedVarManager::kNtracksNegAnalyzed);
      man->AddHistogram(classStr.Data(),"NNegTracksAnalyzed_vs_Run_prof", "# electrons per event", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kNtracksNegAnalyzed, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(),"NTotalTracksAnalyzed","#leg candidates per event",kFALSE,10,0.,10.,AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalTracksAnalyzed_vs_Run_prof", "# leg candidates per event", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kNtracksAnalyzed, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(),"NTotalPairsAnalyzed","#dielectron candidates per event",kFALSE,10,0.,10.,AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"NTotalPairsAnalyzed_vs_Run_prof", "#dielectron candidates per event", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kNpairsSelected, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      
      man->AddHistogram(classStr.Data(),"NTracksTPCout","",kFALSE,200,0.,1200.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout);
      man->AddHistogram(classStr.Data(),"NTracksTPCout_vs_Run_prof", "", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(),"NTracksTPCout_vs_Run", "", kFALSE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 100, 0., 1500., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(),"NTracksTPCout_VZEROmult","",kFALSE,200,0.,1200.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 200, 0., 2000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTPCout_SPDtracklets","",kFALSE,200,0.,1200.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 200, 0., 200., AliReducedVarManager::kSPDntracklets);
      
      
      man->AddHistogram(classStr.Data(),"NTracksAnalyzed_VZEROmult","",kFALSE,10,0.,10.,AliReducedVarManager::kNtracksAnalyzed, 200, 0., 2000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksAnalyzed_VZEROmult_prof","",kTRUE,200, 0., 2000., AliReducedVarManager::kVZEROTotalMult, 10,0.,10.,AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NPairsAnalyzed_VZEROmult_prof","",kTRUE,200, 0., 2000., AliReducedVarManager::kVZEROTotalMult, 10,0.,10.,AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"NTracksAnalyzed_NTracksTPCout_prof","",kTRUE, 30,0.,1500.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 10,0.,10.,AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NPairsAnalyzed_NTracksTPCout_prof","",kTRUE, 30,0.,1500.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 10,0.,10.,AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutVsVZEROmult","",kFALSE,200,0.,20.,AliReducedVarManager::kNTracksTPCoutVsVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutVsVZEROmult_VZEROmult","",kFALSE, 100, 0., 2000., AliReducedVarManager::kVZEROTotalMult, 100,0.,20.,AliReducedVarManager::kNTracksTPCoutVsVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksAnalyzed_RunNo_prof","",kTRUE, gkNRuns, -0.5, -0.5 + gkNRuns, AliReducedVarManager::kRunID, 10, 0., 10., AliReducedVarManager::kNtracksAnalyzed, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(),"NPairsAnalyzed_RunNo_prof","",kTRUE, gkNRuns, -0.5, -0.5 + gkNRuns, AliReducedVarManager::kRunID, 10, 0., 10., AliReducedVarManager::kNpairsSelected, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      continue;
    }  // end if className contains "Event"    
    
    // Track histograms
    if(classStr.Contains("TrackITSclusterMap_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       man->AddHistogram(classStr.Data(), "ITSlayerHit", "Hits in the ITS layers", kFALSE,
                         6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
       if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(), "ITSlayerHit_Phi", "Hits in the ITS layers vs #varphi", kFALSE,
                           180, 0.0, 6.29, AliReducedVarManager::kPhi, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
         man->AddHistogram(classStr.Data(), "ITSlayerHit_Eta", "Hits in the ITS layers vs #eta", kFALSE,
                           100, -1.0, 1.0, AliReducedVarManager::kEta, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
       }
       continue;
    }  // end of ITSclusterMap histogram definitions
    
    if(classStr.Contains("TrackITSsharedClusterMap_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       man->AddHistogram(classStr.Data(), "ITSlayerShared", "Shared cls in the ITS layers", kFALSE,
                         6, 0.5, 6.5, AliReducedVarManager::kITSlayerShared);
       if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(), "ITSlayerShared_Phi", "Shared cls in the ITS layers vs #varphi", kFALSE,
                           180, 0.0, 6.29, AliReducedVarManager::kPhi, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerShared);
         man->AddHistogram(classStr.Data(), "ITSlayerShared_Eta", "Shared cls in the ITS layers vs #eta", kFALSE,
                           100, -1.0, 1.0, AliReducedVarManager::kEta, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerShared);
       }
       continue;
    }  // end of ITSclusterMap histogram definitions
    
    if(classStr.Contains("TrackTPCclusterMap_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       man->AddHistogram(classStr.Data(), "TPCclusterMap", "TPC cluster map", kFALSE,
                         8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
       if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(), "TPCclusterMap_Phi", "TPC cluster map vs #varphi", kFALSE,
                           180, 0.0, 6.29, AliReducedVarManager::kPhi, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
         man->AddHistogram(classStr.Data(), "TPCclusterMap_Eta", "TPC cluster map vs #eta", kFALSE,
                           100, -1.0, 1.0, AliReducedVarManager::kEta, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
         man->AddHistogram(classStr.Data(), "TPCclusterMap_Pt", "TPC cluster map vs p_{T}", kFALSE,
                           100, 0.0, 10.0, AliReducedVarManager::kPt, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
       }
       continue;
    }  // end of TPCclusterMap histogram definitions
    
    TString trkStatusNames = "";
    for(Int_t iflag=0; iflag<AliReducedVarManager::kNTrackingStatus; ++iflag) {
       trkStatusNames += AliReducedVarManager::fgkTrackingStatusNames[iflag];
       trkStatusNames += ";";
    }
    if(classStr.Contains("TrackStatusFlags_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       man->AddHistogram(classStr.Data(), "TrackingFlags", "Tracking flags;;", kFALSE,
                         AliReducedVarManager::kNTrackingStatus, -0.5, AliReducedVarManager::kNTrackingStatus-0.5, AliReducedVarManager::kTrackingFlag, 
                         0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, trkStatusNames.Data());
       continue;
    }
    
    if(classStr.Contains("Track_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      
      //man->AddHistogram(classStr.Data(), "P_Pin", "P vs Pin", kFALSE, 100, 0.0, 10.0, AliReducedVarManager::kP, 100, 0.0, 10.0, AliReducedVarManager::kPin);
      man->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution", kFALSE, 200, 0.0, 10.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "P", "p distribution", kFALSE, 100, 0.0, 10.0, AliReducedVarManager::kP);
      man->AddHistogram(classStr.Data(), "Pin", "p_{IN} distribution", kFALSE, 200, 0.0, 10.0, AliReducedVarManager::kPin);
      man->AddHistogram(classStr.Data(), "Eta", "", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEta);
      if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(), "Eta_Phi", "", kFALSE, 40, -1.0, 1.0, AliReducedVarManager::kEta, 36, 0., 6.29, AliReducedVarManager::kPhi);
         man->AddHistogram(classStr.Data(), "Eta_P", "", kFALSE, 40, -1.0, 1.0, AliReducedVarManager::kEta, 50, 0., 5.0, AliReducedVarManager::kP);
      }
      man->AddHistogram(classStr.Data(), "Phi", "", kFALSE, 36, 0.0, 6.29, AliReducedVarManager::kPhi);
      if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(), "DCAxy_Pt", "DCAxy", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kDcaXY, 50, 0.0, 10.0, AliReducedVarManager::kPt);
         man->AddHistogram(classStr.Data(), "DCAxy_DCAz", "DCAxy_DCAz", kFALSE, 200, -1.0, 1.0, AliReducedVarManager::kDcaXY, 200, -3.0, 3.0, AliReducedVarManager::kDcaZ);
      }
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE, 1000, -2.5, 2.5, AliReducedVarManager::kDcaXY);
      if(classStr.Contains("Standard")) 
         man->AddHistogram(classStr.Data(), "DCAxy_vs_Run_prof", "", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kDcaXY, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE, 1000, -5.0, 5.0, AliReducedVarManager::kDcaZ);
      if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(), "DCAz_vs_Run_prof", "", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kDcaZ, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
         man->AddHistogram(classStr.Data(), "DCAz_Pt", "DCAz", kFALSE, 100, -3.0, 3.0, AliReducedVarManager::kDcaZ, 50, 0.0, 10.0, AliReducedVarManager::kPt);
         man->AddHistogram(classStr.Data(), "DCAxy_goldenChi2", "DCAxy vs golden chi2", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kDcaXY, 200, 0.0, 1000.0, AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
         man->AddHistogram(classStr.Data(), "DCAz_goldenChi2", "DCAz vs golden chi2", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kDcaZ, 200, 0.0, 1000.0, AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
      }
      
      man->AddHistogram(classStr.Data(),"ITSncls", "ITS nclusters", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSncls);
      if(classStr.Contains("Standard")) 
         man->AddHistogram(classStr.Data(), "ITSncls_vs_Run_prof", "", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kITSncls, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      man->AddHistogram(classStr.Data(),"ITSnclsShared", "ITS nclusters shared", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSnclsShared);
      if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(), "ITSnclsShared_vs_Run_prof", "", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kITSnclsShared, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
         man->AddHistogram(classStr.Data(),"Eta_Phi_ITSncls_prof","ITS <nclusters> vs (#eta,#phi)",kTRUE,
                   36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kITSncls);
         man->AddHistogram(classStr.Data(),"Eta_Phi_ITSnclsShared_prof","ITS <nclusters-shared> vs (#eta,#phi)",kTRUE,
                        36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kITSnclsShared);
      }
      man->AddHistogram(classStr.Data(),"ITSchi2", "ITS #chi^{2}", kFALSE, 200,0.0,50.0, AliReducedVarManager::kITSchi2);
      //man->AddHistogram(classStr.Data(),"ITSchi2_RunID_prof","",kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 10, 0., 10., AliReducedVarManager::kITSchi2);
      if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(),"ITSchi2_ITSncls", "ITS #chi^{2} vs ITS clusters", kFALSE, 100,0.0,20.0, AliReducedVarManager::kITSchi2,7,-0.5,6.5, AliReducedVarManager::kITSncls);
         /*man->AddHistogram(classStr.Data(),"Eta_Phi_ITSchi2_prof","ITS <#chi^{2}> vs (#eta,#phi)",kTRUE,
                        36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 100, 0.0, 1000, AliReducedVarManager::kITSchi2);        */
         man->AddHistogram(classStr.Data(),"ITSchi2_NTracksTPCout_prof","",kTRUE,30,0.,1500.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 10, 0., 10., AliReducedVarManager::kITSchi2);
      }
      man->AddHistogram(classStr.Data(),"Chi2TPCConstrainedVsGlobal","",kFALSE,1000,0.,200.,AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
      if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(),"Chi2TPCConstrainedVsGlobal_NTracksTPCout_prof","",kTRUE,30,0.,1500.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout,1000,0.,200.,AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
         man->AddHistogram(classStr.Data(), "Chi2TPCConstrainedVsGlobal_vs_Run_prof", "", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
      }
      
      man->AddHistogram(classStr.Data(),"TPCncls","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCncls);
      if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(),"TPCncls_RunID_prof","",kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 10, 0., 10., AliReducedVarManager::kTPCncls, 0, 0.,0.,AliReducedVarManager::kNothing, gRunList.Data());
         man->AddHistogram(classStr.Data(),"TPCncls_NTracksTPCout_prof","",kTRUE,30,0.,1500.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 10, 0., 10., AliReducedVarManager::kTPCncls);
      }
      man->AddHistogram(classStr.Data(),"TPCcrossedRows","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCcrossedRows);
      man->AddHistogram(classStr.Data(),"TPCfindableClusters","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsF);
      man->AddHistogram(classStr.Data(),"TPCcrossedRowsOverFindableClusters","", kFALSE, 200,0.0,1.5,AliReducedVarManager::kTPCcrossedRowsOverFindableClusters);
      man->AddHistogram(classStr.Data(),"TPCnclsShared","",kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsShared);
      if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(),"TPCnclsShared_NTracksTPCout_prof","",kTRUE, 30,0.,1500.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout,160,-0.5,159.5,AliReducedVarManager::kTPCnclsShared);
      }
      man->AddHistogram(classStr.Data(),"TPCnclsSharedRatio","",kFALSE, 200,0.,1.,AliReducedVarManager::kTPCnclsSharedRatio);        
      
      man->AddHistogram(classStr.Data(),"TPCchi2","",kFALSE, 200,0.0,10.0,AliReducedVarManager::kTPCchi2);
      if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(),"TPCchi2_NTracksTPCout_prof","",kTRUE, 30,0.,1500.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout,10,0.0,10.,AliReducedVarManager::kTPCchi2);
         man->AddHistogram(classStr.Data(), "TPCchi2_vs_Run_prof", "", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kTPCchi2, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
         man->AddHistogram(classStr.Data(),"TPCchi2_Pt","",kFALSE, 200,0.0,10.0,AliReducedVarManager::kTPCchi2, 100, 0., 10.0, AliReducedVarManager::kPt);
         man->AddHistogram(classStr.Data(),"TPCchi2_TPCncls", "TPC #chi^{2} vs TPC clusters", kFALSE, 100,0.0,10.0, AliReducedVarManager::kTPCchi2,160,0.,160., AliReducedVarManager::kTPCncls);
         /*man->AddHistogram(classStr.Data(),"TPCchi2_CentVZERO", "TPC #chi^{2} vs CentVZERO", kFALSE, 100,0.0,10.0, AliReducedVarManager::kTPCchi2,20,0.,100., AliReducedVarManager::kCentVZERO);*/
         man->AddHistogram(classStr.Data(),"Eta_TPCchi2","TPC #chi^{2} vs #eta",kFALSE,
                        36, -0.9, 0.9, AliReducedVarManager::kEta, 100, 0.0, 10., AliReducedVarManager::kTPCchi2);
      }
      man->AddHistogram(classStr.Data(),"TPCsegments","",kFALSE, 9,0.0,9.0,AliReducedVarManager::kTPCNclusBitsFired);
      
      /*man->AddHistogram(classStr.Data(),"Eta_Phi_TPCncls_prof","TPC <nclusters> vs (#eta,#phi)",kTRUE,
       *             36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCncls);    
       m an->Ad*dHistogram(classStr.Data(),"Eta_Phi_TPCcrossedRows_prof","TPC <n crossed rows> vs (#eta,#phi)",kTRUE,
       36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCcrossedRows);
       man->AddHistogram(classStr.Data(),"Eta_Phi_TPCnclsF_prof","TPC <nclusters findable> vs (#eta,#phi)",kTRUE,
       36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCnclsF);*/
      if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(),"Eta_Phi_TPCchi2_prof","TPC <#chi^{2}> vs (#eta,#phi)",kTRUE,
                           36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCchi2);
      }
      
      /*man->AddHistogram(classStr.Data(),"TPCchi2_NTracksTPCout_VZEROmult_prof","",kTRUE,20,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20, 0., 30000., AliReducedVarManager::kVZEROTotalMult, 10, 0., 10., AliReducedVarManager::kTPCchi2);*/
      /*man->AddHistogram(classStr.Data(),"Eta_Phi_TPCchi2","TPC #chi^{2} vs (#eta,#phi)",kFALSE,
       *                  36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 100, 0., 10.0, AliReducedVarManager::kTPCchi2); */
      man->AddHistogram(classStr.Data(),"TPCsignal_Pin","TPC dE/dx vs. inner param P",kFALSE,
                        100,0.0,10.0,AliReducedVarManager::kPin,100,49.5,150.5,AliReducedVarManager::kTPCsignal);
      if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(),"TPCsignal_NTracksTPCout_prof","",kTRUE, 30,0.,1500.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout,160,49.5,150.5,AliReducedVarManager::kTPCsignal);
      }
      man->AddHistogram(classStr.Data(),"TPCsignalN","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCsignalN);
      /*man->AddHistogram(classStr.Data(),"Eta_Phi_TPCsignalN_prof","TPC <nclusters pid> vs (#eta,#phi)",kTRUE,
       *             36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCsignalN);    */
      man->AddHistogram(classStr.Data(),"TPCnsigElectron_Pin","TPC N_{#sigma} electron vs. inner param P",kFALSE,
                        100,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      if(classStr.Contains("Standard")) {
         man->AddHistogram(classStr.Data(), "TPCnsigElectron_vs_Run_prof", "", kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 300, 0., 300., AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, 0, 0., 0., AliReducedVarManager::kNothing, gRunList.Data());
         man->AddHistogram(classStr.Data(),"TPCnsigElectron_NTracksTPCout_prof","",kTRUE, 30,0.,1500.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      }
      man->AddHistogram(classStr.Data(),"TPCnsigElectron_Eta","TPC N_{#sigma} electron vs. #eta",kFALSE,
                        36,-0.9,0.9,AliReducedVarManager::kEta,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin","TPC N_{#sigma} electron corrected vs. inner param P",kFALSE,
                        100,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Eta","TPC N_{#sigma} electron corrected vs. #eta",kFALSE,
                        36,-0.9,0.9,AliReducedVarManager::kEta,100,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron);
      
      /*man->AddHistogram(classStr.Data(),"TPCchi2_vs_TPCnsigElectron_Pin_prof","<TPC chi2> vs TPC N_{#sigma} electron and inner param P",kTRUE,
                        20,0.0,5.0,AliReducedVarManager::kPin,30,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, 10, 0., 10., AliReducedVarManager::kTPCchi2);*/
      /*man->AddHistogram(classStr.Data(),"TPCchi2_vs_TPCnsigElectron_Pin","TPC chi2 vs TPC N_{#sigma} electron and inner param P",kFALSE,
                        20,0.0,5.0,AliReducedVarManager::kPin,30,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, 80, 0., 8., AliReducedVarManager::kTPCchi2);*/
      /*man->AddHistogram(classStr.Data(),"TPCchi2_vs_TPCnsigElectronCorrected_Pin_prof","<TPC chi2> vs TPC N_{#sigma} electron corrected and inner param P",kTRUE,
                        20,0.0,5.0,AliReducedVarManager::kPin,30,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 10, 0., 10., AliReducedVarManager::kTPCchi2);*/
      if(classStr.Contains("Standard")) {  
         man->AddHistogram(classStr.Data(),"TOFbeta_P","TOF #beta vs P",kFALSE,
                 100,0.0,10.0,AliReducedVarManager::kP, 110,0.0,1.1,AliReducedVarManager::kTOFbeta);
         man->AddHistogram(classStr.Data(),"TOFnSigEle_P","TOF nSigEle vs P",kFALSE,
                        100,0.0,10.0,AliReducedVarManager::kP, 100,-5.0,5.0,AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron);
         man->AddHistogram(classStr.Data(),"TOFnSigProt_P","TOF nSigProton vs P",kFALSE,
                        100,0.0,10.0,AliReducedVarManager::kP, 100,-5.0,5.0,AliReducedVarManager::kTOFnSig+AliReducedVarManager::kProton);
         man->AddHistogram(classStr.Data(),"TOFmismatchProbab","TOF mismatch probability",kFALSE,
                        2,0.0,2.0,AliReducedVarManager::kTOFmismatchProbability);
      }
      
      //man->AddHistogram(classStr.Data(),"TPCchi2_ITSchi2","",kFALSE, 100,0.0,10.0,AliReducedVarManager::kTPCchi2, 100,0.0,40.0,AliReducedVarManager::kITSchi2);
      //man->AddHistogram(classStr.Data(),"ITSchi2_Chi2TPCconstrainedVsGlobal","",kFALSE, 500,0.0,1000.0,AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 100,0.0,40.0,AliReducedVarManager::kITSchi2);
      //man->AddHistogram(classStr.Data(),"TPCchi2_Chi2TPCconstrainedVsGlobal","",kFALSE, 500,0.0,1000.0,AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 100,0.0,10.0,AliReducedVarManager::kTPCchi2);
      
      
        if(classStr.Contains("TrueElectron")) {
          man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE, 150, 0., 15.0, AliReducedVarManager::kPtMC);
          //man->AddHistogram(classStr.Data(), "PtRec_PtMC", "p_{T} MC vs p_{T} reconstructed", kFALSE, 100, 0., 10.0, AliReducedVarManager::kPtMC, 100, 0., 10.0, AliReducedVarManager::kPt);
          man->AddHistogram(classStr.Data(), "PhiMC", "#varphi MC", kFALSE, 180, 0., 6.3, AliReducedVarManager::kPhiMC);
          //man->AddHistogram(classStr.Data(), "PhiRec_PhiMC", "#varphi MC vs #varphi reconstructed", kFALSE, 180, 0., 6.3, AliReducedVarManager::kPhiMC, 180, 0., 6.3, AliReducedVarManager::kPhi);
          man->AddHistogram(classStr.Data(), "EtaMC", "#eta MC", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEtaMC);
          //man->AddHistogram(classStr.Data(), "EtaRec_EtaMC", "#eta MC vs #eta reconstructed", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEtaMC, 100, -1.0, 1.0, AliReducedVarManager::kEta);          
          //man->AddHistogram(classStr.Data(), "PDGcode0", "PDG code of the track", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC);
          //man->AddHistogram(classStr.Data(), "PDGcode1", "PDG code of the track's mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+1);
          //man->AddHistogram(classStr.Data(), "PDGcode2", "PDG code of the track's grand-mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+2);
          //man->AddHistogram(classStr.Data(), "PDGcode3", "PDG code of the track's grand-grand mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+3);
        }
      continue;
    }  // end if "TrackQA"
        
     Int_t vars[3] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, AliReducedVarManager::kVtxZ};
    //Int_t vars[3] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, AliReducedVarManager::kSPDntracklets};
    //Int_t vars[4] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, AliReducedVarManager::kVtxZ, AliReducedVarManager::kSPDntracklets};
       
    TArrayD pairHistBinLimits[3];
    pairHistBinLimits[0] = TArrayD(gkNMassBins, gMassBins);
    if(task->GetRunOverMC()) {
       pairHistBinLimits[1] = TArrayD(gkNPtBinsMC+1, gPtLimsMC);
       pairHistBinLimits[2] = TArrayD(gkNZBinsMC, gZBinLimsMC);
    }
    else {
      pairHistBinLimits[1] = TArrayD(gkNPtBins, gPtLims);
      pairHistBinLimits[2] = TArrayD(gkNZBins, gZBinLims);
    }
    //pairHistBinLimits[3] = TArrayD(gkNSPDtrkBins, gSPDtrkBinLims);
    
    // Histograms for pairs
    if(classStr.Contains("Pair")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "PairInvMass", "Differential pair inv. mass", 3, vars, pairHistBinLimits);
      if(classStr.Contains("PairSE") || classStr.Contains("PairPrefilterSE")) {
        man->AddHistogram(classStr.Data(), "PairType", "Pair type", kFALSE, 4, -0.5, 3.5, AliReducedVarManager::kPairType);
	man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE, gkNMassBins-1, gMassBins, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "MassWide", "Invariant mass", kFALSE, 300, 0.0, 12.0, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPt);
        man->AddHistogram(classStr.Data(), "Pt_coarse", "", kFALSE, 20, 0.0, 20.0, AliReducedVarManager::kPt);
        //man->AddHistogram(classStr.Data(), "Leg1Pt_Leg2Pt", "", kFALSE, 100, 0.0, 10.0, AliReducedVarManager::kPairLegPt, 
        //                  100, 0.0, 10.0, AliReducedVarManager::kPairLegPt+1);
        //man->AddHistogram(classStr.Data(), "Leg1Pt_Leg2Pt_Mass", "", kFALSE, 100, 0.0, 10.0, AliReducedVarManager::kPairLegPt, 
        //                  100, 0.0, 10.0, AliReducedVarManager::kPairLegPt+1, 25, 2.5, 3.5, AliReducedVarManager::kMass);
        
        //man->AddHistogram(classStr.Data(), "LegPtSum", "", kFALSE, 100, 0.0, 10.0, AliReducedVarManager::kPairLegPtSum);
        //man->AddHistogram(classStr.Data(), "LegPtSum_PairPt", "", kFALSE, 100, 0.0, 10.0, AliReducedVarManager::kPt, 
        //                                                                  100, 0.0, 10.0, AliReducedVarManager::kPairLegPtSum);
        //man->AddHistogram(classStr.Data(), "LegPtSum_Mass", "", kFALSE, 100, 0.0, 10.0, AliReducedVarManager::kPairLegPtSum, 
        //                  100, 0.0, 10.0, AliReducedVarManager::kMass);
        
	//man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kRap);
	//man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 315, 0.0, 6.3, AliReducedVarManager::kPhi);
      }   // end if "QA"
      
      if(classStr.Contains("TrueElectron")) {
         man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE, 200, 0., 10.0, AliReducedVarManager::kPtMC);
         man->AddHistogram(classStr.Data(), "PtMC_coarse", "p_{T} MC", kFALSE, 20, 0., 20.0, AliReducedVarManager::kPtMC);
         //man->AddHistogram(classStr.Data(), "MassMC_Mass", "Invariant mass, MC vs reconstructed", kFALSE, 150, 2.0, 3.5, AliReducedVarManager::kMass,
          //  150, 2.0, 3.5, AliReducedVarManager::kMassMC);
         //man->AddHistogram(classStr.Data(), "PtMC_Pt", "pair pT, MC vs reconstructed", kFALSE, 150, 0.0, 15., AliReducedVarManager::kPt,
          //                 150, 0.0, 15.0, AliReducedVarManager::kPtMC);
         /*man->AddHistogram(classStr.Data(), "PtMC_MassMC", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPtMC,
                           30, 2.7, 3.3, AliReducedVarManager::kMassMC);
         man->AddHistogram(classStr.Data(), "Pt_Mass", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPt,
                           30, 2.7, 3.3, AliReducedVarManager::kMass);
         man->AddHistogram(classStr.Data(), "PtMC_Mass", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPtMC,
                           30, 2.7, 3.3, AliReducedVarManager::kMass);
         man->AddHistogram(classStr.Data(), "Pt_MassMC", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPt,
                           30, 2.7, 3.3, AliReducedVarManager::kMassMC);  */
         //man->AddHistogram(classStr.Data(), "PtMC_Pt_lowPtZoom", "pair pT, MC vs reconstructed", kFALSE, 100, 0.0, 0.5, AliReducedVarManager::kPt,
          //                 100, 0.0, 0.5, AliReducedVarManager::kPtMC);
          //man->AddHistogram(classStr.Data(), "PtMC_Pt_Mass_coarse", "pair pT, MC vs reconstructed, vs mass", kFALSE, 25, 0.0, 0.5, AliReducedVarManager::kPt,
          //                25, 0.0, 0.5, AliReducedVarManager::kPtMC, 15, 2.72, 3.32, AliReducedVarManager::kMass);
          if(classStr.Contains("Standard")) {
            man->AddHistogram(classStr.Data(), "PtMC_Pt_Mass", "pair pT, MC vs reconstructed, vs mass", kFALSE, 200, 0.0, 20.0, AliReducedVarManager::kPt,
                              200, 0.0, 20.0, AliReducedVarManager::kPtMC, 15, 2.72, 3.32, AliReducedVarManager::kMass);
            man->AddHistogram(classStr.Data(), "CosThetaStarCS", "cos(#theta^{*})_{CS}", kFALSE, 88, -1.1, 1.1, AliReducedVarManager::kPairThetaCS);
            man->AddHistogram(classStr.Data(), "CosThetaStarCS_pt_mass", "cos(#theta^{*})_{CS} vs reconstructed pt vs mass", kFALSE, 88, -1.1, 1.1, AliReducedVarManager::kPairThetaCS, 400, 0., 20., AliReducedVarManager::kPt, 15, 2.72, 3.32, AliReducedVarManager::kMass);
            man->AddHistogram(classStr.Data(), "CosThetaStarHE", "cos(#theta^{*})_{HE}", kFALSE, 88, -1.1, 1.1, AliReducedVarManager::kPairThetaHE);
            man->AddHistogram(classStr.Data(), "CosThetaStarHE_pt_mass", "cos(#theta^{*})_{HE} vs reconstructed pt vs mass", kFALSE, 88, -1.1, 1.1, AliReducedVarManager::kPairThetaHE, 400, 0., 20., AliReducedVarManager::kPt, 15, 2.72, 3.32, AliReducedVarManager::kMass);
            
          }
      }
      
      continue;
    }   // end if for Pair classes of histograms
  }  // end loop over histogram classes
}
