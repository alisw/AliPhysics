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
const Int_t gkNRuns = 2;
TString gRunList = "280234;280235";

// centrality binning used for event mixing and main pair histograms
const Int_t gkNCentBins = 15;
Double_t gCentLims[gkNCentBins] = {
   0.0, 2.5, 5.0, 7.5, 10., 
   15.0, 20., 25., 30., 40., 
   50.0,  60.0, 70.0, 80., 90.0
};

// pt binning
const Int_t gkNPtBins = 12;
Double_t gPtLims[gkNPtBins] = {
   0.0, 0.50, 1.0, 1.5, 2.0, 
   2.5, 3.0, 4.0, 5.0, 7.0,
   10., 20.0
};

// vertex-z binning
const Int_t gkNVtxZBins = 7;
Double_t gVtxZLims[gkNVtxZBins] = {
   -10.0, -7.0, -3.0, 0.0, 3.0, 7.0, 10.0
};

// mass binning
const Int_t gkNMassBins = 91;
Double_t gMassBins[gkNMassBins] = {0.};


//__________________________________________________________________________________________
AliAnalysisTask* AddTask_iarsene_jpsi2ee_XeXe(Bool_t isAliRoot=kTRUE, Int_t runMode=1, TString prod="LHC10h"){    
   //
   // isAliRoot=kTRUE for ESD/AOD analysis in AliROOT, kFALSE for root analysis on reduced trees
   // runMode=1 (AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents)
   //               =2 (AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree)
   //
   //get the current analysis manager

  printf("INFO on AddTask_iarsene_jpsi2ee(): (isAliRoot, runMode) :: (%d,%d) \n", isAliRoot, runMode);

  AliReducedAnalysisJpsi2ee* jpsi2eeAnalysis = new AliReducedAnalysisJpsi2ee("Jpsi2eeAnalysis","Jpsi->ee analysis");
  jpsi2eeAnalysis->Init();
  //jpsi2eeAnalysis->SetLoopOverTracks(kFALSE);
  //jpsi2eeAnalysis->SetRunEventMixing(kFALSE);
  //jpsi2eeAnalysis->SetRunPairing(kFALSE);
  jpsi2eeAnalysis->SetRunOverMC(kTRUE);
  //jpsi2eeAnalysis->SetRunLikeSignPairing(kFALSE);
  //jpsi2eeAnalysis->SetRunPrefilter(kFALSE);
  Setup(jpsi2eeAnalysis, prod);
  // initialize an AliAnalysisTask which will wrapp the AliReducedAnalysisJpsi2ee such that it can be run in an aliroot analysis train (e.g. LEGO, local analysis etc)
  AliAnalysisTaskReducedEventProcessor* task = new AliAnalysisTaskReducedEventProcessor("ReducedEventAnalysisManager", runMode);
  task->AddTask(jpsi2eeAnalysis);
  
  if(isAliRoot){
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
  
     TString outputFilename = "AnalysisHistograms_jpsi2ee_XeXe.root";
     if(jpsi2eeAnalysis->GetRunOverMC()) outputFilename = "AnalysisHistograms_jpsi2ee_XeXe_MC.root";
     AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("jpsi2eeHistos", THashList::Class(),
                                                                  AliAnalysisManager::kOutputContainer, outputFilename.Data());
     mgr->ConnectOutput(task, 1, cOutputHist );
  }
  else {
    // nothing at the moment   
  }
  
  return task;
}


//_________________________________________________________________
void Setup(AliReducedAnalysisJpsi2ee* processor, TString prod /*="LHC10h"*/) {
  //
  // Configure the analysis task
  // Setup histograms, handlers, cuts, etc.
  //
   Bool_t isMC = processor->GetRunOverMC();
   
  if(!processor->GetRunOverMC()) {
      //TFile* corrFile = TFile::Open("/home/iarsene/work/ALICE/XeXe2017/tpcPostcalibrationMaps_cent0_30.root");
      TFile* corrFile = TFile::Open("/home/iarsene/work/ALICE/XeXe2017/tpcPostcalibrationMaps_cent30_90.root");
      AliReducedVarManager::SetTPCelectronCorrectionMaps((TH2F*)corrFile->Get("Mean"), (TH2F*)corrFile->Get("Width"), 
                                                         AliReducedVarManager::kEta, AliReducedVarManager::kCentVZERO);
      corrFile->Close();   
  }
   
  // Set event cuts
  AliReducedEventCut* evCut1 = new AliReducedEventCut("Centrality","Centrality selection");
  evCut1->AddCut(AliReducedVarManager::kCentVZERO, 0., 90.);
  evCut1->AddCut(AliReducedVarManager::kVtxZ, -10.0, 10.0);
  if(!isMC)  evCut1->AddCut(AliReducedVarManager::kIsPhysicsSelection, 0.1, 2.);   // request physics selection
  evCut1->EnableVertexDistanceCut();
  TF1* cutCorrTPCoutVZEROmult = new TF1("cutCorrTPCoutVZEROmult", "[0]+[1]*x+[2]*x*x", 0., 1.e+5);
  cutCorrTPCoutVZEROmult->SetParameters(-1000., 2.8, 1.2e-5);
  if(prod.Contains("LHC17n") && !processor->GetRunOverMC()) 
     evCut1->AddCut(AliReducedVarManager::kVZEROTotalMult, cutCorrTPCoutVZEROmult, 99999., kFALSE, AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 0., 99998.);
    
  processor->AddEventCut(evCut1);
  
  // define some generic track cuts 
  AliReducedTrackCut* spdAny = new AliReducedTrackCut("spdAny","");
  spdAny->SetRequestSPDany();
  AliReducedTrackCut* spdFirst = new AliReducedTrackCut("spdFirst","");
  spdFirst->SetRequestSPDfirst();
  AliReducedTrackCut* conversionElectron = new AliReducedTrackCut("conversionElectron","");
  conversionElectron->SetTrackQualityFilterBit(kPhotonConversionV0Leg, kFALSE);
  
  SetupTrackingSystematicsCuts(processor, "prot36_pion34", kFALSE, 3.0, 3.6, 3.4);
  //SetupTrackingSystematicsCuts(processor, "prot30_pion30", kTRUE, 3.0, 3.0, 3.0);
  SetupTrackingSystematicsCuts(processor, "prot40_pion36", kFALSE, 3.0, 4.0, 3.6);
  
  //SetupPIDSystematicsCuts(processor, kTRUE);   // kFALSE -> use SPDfirst
  
  //SetupPIDSystematicsCutsSingleLeg(processor, kTRUE);       // calibrated electron band
  //SetupPIDSystematicsCutsSingleLeg(processor, kFALSE);      // uncalibrated
  
  
  // set track prefilter cuts
  AliReducedTrackCut* prefTrackCut1 = new AliReducedTrackCut("prefCutPt09","prefilter P selection");
  prefTrackCut1->SetTrackFilterBit(kPrefilterCut);
  prefTrackCut1->AddCut(AliReducedVarManager::kPt, 0.8,100.0);
  //prefTrackCut1->SetRequestTPCrefit();
  prefTrackCut1->SetRequestITSrefit();
  processor->AddPrefilterTrackCut(prefTrackCut1);  
  
  // MC truth info for reconstructed electron candidates
  AliReducedTrackCut* trueElectron = new AliReducedTrackCut("TrueElectron", "reconstructed electrons with MC truth");
  trueElectron->SetMCFilterBit(kJpsiDecayElectron);
  processor->AddLegCandidateMCcut(trueElectron);  
  
  AliReducedTrackCut* trueElectronNonPrompt = new AliReducedTrackCut("TrueElectronNonPrompt", "reconstructed electrons with MC truth");
  trueElectronNonPrompt->SetMCFilterBit(kJpsiNonPromptDecayElectron);
  processor->AddLegCandidateMCcut(trueElectronNonPrompt);

  AliReducedTrackCut* trueElectronPrompt = new AliReducedTrackCut("TrueElectronPrompt", "reconstructed electrons with MC truth");
  trueElectronPrompt->SetMCFilterBit(kJpsiPromptDecayElectron);
  processor->AddLegCandidateMCcut(trueElectronPrompt);    
  
  // MC truth selections(s) for the J/psi mother and electron daughters
  AliReducedTrackCut* mcTruthJpsi = new AliReducedTrackCut("mcTruthJpsi", "Pure MC truth J/psi");
  mcTruthJpsi->SetMCFilterBit(kJpsiInclusive);
  mcTruthJpsi->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  AliReducedTrackCut* mcTruthJpsiElectron = new AliReducedTrackCut("mcTruthJpsiElectron", "Pure MC truth electron from J/psi");
  mcTruthJpsiElectron->SetMCFilterBit(kJpsiDecayElectron);
  mcTruthJpsiElectron->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  mcTruthJpsiElectron->AddCut(AliReducedVarManager::kPt, 1.0, 100.0);
  processor->AddJpsiMotherMCCut(mcTruthJpsi, mcTruthJpsiElectron);  
  
  AliReducedTrackCut* mcTruthJpsiNonPrompt = new AliReducedTrackCut("mcTruthJpsiNonPrompt", "Pure MC truth J/psi");
  mcTruthJpsiNonPrompt->SetMCFilterBit(kJpsiNonPrompt);
  mcTruthJpsiNonPrompt->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  AliReducedTrackCut* mcTruthJpsiElectronNonPrompt = new AliReducedTrackCut("mcTruthJpsiElectronNonPrompt", "Pure MC truth electron from J/psi");
  mcTruthJpsiElectronNonPrompt->SetMCFilterBit(kJpsiNonPromptDecayElectron);
  mcTruthJpsiElectronNonPrompt->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  mcTruthJpsiElectronNonPrompt->AddCut(AliReducedVarManager::kPt, 1.0, 100.0);
  processor->AddJpsiMotherMCCut(mcTruthJpsiNonPrompt, mcTruthJpsiElectronNonPrompt);

  AliReducedTrackCut* mcTruthJpsiPrompt = new AliReducedTrackCut("mcTruthJpsiPrompt", "Pure MC truth J/psi");
  mcTruthJpsiPrompt->SetMCFilterBit(kJpsiPrompt);
  mcTruthJpsiPrompt->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  AliReducedTrackCut* mcTruthJpsiElectronPrompt = new AliReducedTrackCut("mcTruthJpsiElectronPrompt", "Pure MC truth electron from J/psi");
  mcTruthJpsiElectronPrompt->SetMCFilterBit(kJpsiPromptDecayElectron);
  mcTruthJpsiElectronPrompt->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  mcTruthJpsiElectronPrompt->AddCut(AliReducedVarManager::kPt, 1.0, 100.0);
  processor->AddJpsiMotherMCCut(mcTruthJpsiPrompt, mcTruthJpsiElectronPrompt);    
  
  // set pair prefilter cuts
  AliReducedVarCut* prefPairCut = new AliReducedVarCut("prefCutM50MeV","prefilter pair cuts");
  prefPairCut->AddCut(AliReducedVarManager::kMass, 0.0, 0.05, kTRUE);
  processor->AddPrefilterPairCut(prefPairCut);
  
  // Set pair cuts
  AliReducedTrackCut* pairCut1 = new AliReducedTrackCut("Ptpair","Pt pair selection");
  pairCut1->AddCut(AliReducedVarManager::kPt, 0.0,100.0);
  pairCut1->AddCut(AliReducedVarManager::kRap, -0.9,0.9);
  processor->AddPairCut(pairCut1);
    
  SetupHistogramManager(processor, prod);
  SetupMixingHandler(processor);
}


//__________________________________________________________________________________________
void SetupTrackingSystematicsCuts(AliReducedAnalysisJpsi2ee* processor, const Char_t* settingName, Bool_t useSPDany, Double_t eleCut, Double_t protCut, Double_t pionCut) {
   //
   // add tracking systematic cuts
   //
   Bool_t isMC = processor->GetRunOverMC();
   
   AliReducedTrackCut* commonVarCuts = new AliReducedTrackCut();
   commonVarCuts->SetTrackFilterBit(kBasicCut);
   commonVarCuts->AddCut(AliReducedVarManager::kPt, 1.0, 100.0);  
   commonVarCuts->SetRejectKinks();
   commonVarCuts->AddCut(AliReducedVarManager::kTPCncls, 70., 161.);   
   if(isMC)   commonVarCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -1.0*eleCut, eleCut);
   else      commonVarCuts->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -1.0*eleCut, eleCut);
   commonVarCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, protCut, 30000.0);
   commonVarCuts->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, pionCut, 30000.0);
   
   AliReducedTrackCut* spdAny = new AliReducedTrackCut("spdAny","");
   spdAny->SetRequestSPDany();
   AliReducedTrackCut* spdFirst = new AliReducedTrackCut("spdFirst","");
   spdFirst->SetRequestSPDfirst();
   
   // add track cuts to the jpsi2ee task
   AliReducedCompositeCut* spdAnyCut = new AliReducedCompositeCut(Form("%s_spdAny", settingName), "", kTRUE);
   spdAnyCut->AddCut(commonVarCuts);
   spdAnyCut->AddCut(spdAny);
   AliReducedVarCut* itsExtraCuts = new AliReducedTrackCut();
   itsExtraCuts->AddCut(AliReducedVarManager::kITSchi2, 0., 15.);
   //itsExtraCuts->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
   itsExtraCuts->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 6., 9.);
   spdAnyCut->AddCut(itsExtraCuts);
   processor->AddTrackCut(spdAnyCut);   
   
   AliReducedCompositeCut* spdFirstCut = new AliReducedCompositeCut(Form("%s_spdFirst", settingName), "", kTRUE);
   spdFirstCut->AddCut(commonVarCuts);
   spdFirstCut->AddCut(itsExtraCuts);
   spdFirstCut->AddCut(spdFirst);
   processor->AddTrackCut(spdFirstCut);
   
   AliReducedCompositeCut* itsCls3 = new AliReducedCompositeCut(Form("%s_itsCls3", settingName), "", kTRUE);
   itsCls3->AddCut(commonVarCuts);
   itsCls3->AddCut(itsExtraCuts);
   AliReducedTrackCut* itsCls3ExtraCut = new AliReducedTrackCut();
   itsCls3ExtraCut->AddCut(AliReducedVarManager::kITSncls, 2.9, 6.1);
   itsCls3->AddCut(itsCls3ExtraCut);
   processor->AddTrackCut(itsCls3);
   
   AliReducedCompositeCut* itsCls4 = new AliReducedCompositeCut(Form("%s_itsCls4", settingName), "", kTRUE);
   itsCls4->AddCut(commonVarCuts);
   itsCls4->AddCut(itsExtraCuts);
   AliReducedTrackCut* itsCls4ExtraCut = new AliReducedTrackCut();
   itsCls4ExtraCut->AddCut(AliReducedVarManager::kITSncls, 3.9, 6.1);
   itsCls4->AddCut(itsCls4ExtraCut);
   processor->AddTrackCut(itsCls4);
   
   AliReducedCompositeCut* itsCls5 = new AliReducedCompositeCut(Form("%s_itsCls5", settingName), "", kTRUE);
   itsCls5->AddCut(commonVarCuts);
   itsCls5->AddCut(itsExtraCuts);
   AliReducedTrackCut* itsCls5ExtraCut = new AliReducedTrackCut();
   itsCls5ExtraCut->AddCut(AliReducedVarManager::kITSncls, 4.9, 6.1);
   itsCls5->AddCut(itsCls5ExtraCut);
   //processor->AddTrackCut(itsCls5);
   
   AliReducedCompositeCut* itsChi2_11 = new AliReducedCompositeCut(Form("%s_itsChi2_11", settingName), "", kTRUE);
   itsChi2_11->AddCut(commonVarCuts);
   AliReducedTrackCut* itsChi2CommonExtraCuts = new AliReducedTrackCut();
   if(useSPDany)  itsChi2CommonExtraCuts->SetRequestSPDany();
   else                 itsChi2CommonExtraCuts->SetRequestSPDfirst();
   itsChi2CommonExtraCuts->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
   itsChi2CommonExtraCuts->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 6., 9.);
   itsChi2_11->AddCut(itsChi2CommonExtraCuts);
   AliReducedTrackCut* itsChi2_11_cut = new AliReducedTrackCut();
   itsChi2_11_cut->AddCut(AliReducedVarManager::kITSchi2, 0., 11.);
   itsChi2_11->AddCut(itsChi2_11_cut);
   processor->AddTrackCut(itsChi2_11);
   
   AliReducedCompositeCut* itsChi2_9 = new AliReducedCompositeCut(Form("%s_itsChi2_9", settingName), "", kTRUE);
   itsChi2_9->AddCut(commonVarCuts);
   itsChi2_9->AddCut(itsChi2CommonExtraCuts);
   AliReducedTrackCut* itsChi2_9_cut = new AliReducedTrackCut();
   itsChi2_9_cut->AddCut(AliReducedVarManager::kITSchi2, 0., 9.);
   itsChi2_9->AddCut(itsChi2_9_cut);
   processor->AddTrackCut(itsChi2_9);
   
   AliReducedCompositeCut* itsChi2_7 = new AliReducedCompositeCut(Form("%s_itsChi2_7", settingName), "", kTRUE);
   itsChi2_7->AddCut(commonVarCuts);
   itsChi2_7->AddCut(itsChi2CommonExtraCuts);
   AliReducedTrackCut* itsChi2_7_cut = new AliReducedTrackCut();
   itsChi2_7_cut->AddCut(AliReducedVarManager::kITSchi2, 0., 7.);
   itsChi2_7->AddCut(itsChi2_7_cut);
   processor->AddTrackCut(itsChi2_7);
   
   AliReducedCompositeCut* itsChi2_13 = new AliReducedCompositeCut(Form("%s_itsChi2_13", settingName), "", kTRUE);
   itsChi2_13->AddCut(commonVarCuts);
   itsChi2_13->AddCut(itsChi2CommonExtraCuts);
   AliReducedTrackCut* itsChi2_13_cut = new AliReducedTrackCut();
   itsChi2_13_cut->AddCut(AliReducedVarManager::kITSchi2, 0., 13.);
   itsChi2_13->AddCut(itsChi2_13_cut);
   processor->AddTrackCut(itsChi2_13);
   
   AliReducedCompositeCut* itsChi2_17 = new AliReducedCompositeCut(Form("%s_itsChi2_17", settingName), "", kTRUE);
   itsChi2_17->AddCut(commonVarCuts);
   itsChi2_17->AddCut(itsChi2CommonExtraCuts);
   AliReducedTrackCut* itsChi2_17_cut = new AliReducedTrackCut();
   itsChi2_17_cut->AddCut(AliReducedVarManager::kITSchi2, 0., 17.);
   itsChi2_17->AddCut(itsChi2_17_cut);
   processor->AddTrackCut(itsChi2_17);
   
   AliReducedCompositeCut* itsChi2_19 = new AliReducedCompositeCut(Form("%s_itsChi2_19", settingName), "", kTRUE);
   itsChi2_19->AddCut(commonVarCuts);
   itsChi2_19->AddCut(itsChi2CommonExtraCuts);
   AliReducedTrackCut* itsChi2_19_cut = new AliReducedTrackCut();
   itsChi2_19_cut->AddCut(AliReducedVarManager::kITSchi2, 0., 19.);
   itsChi2_19->AddCut(itsChi2_19_cut);
   processor->AddTrackCut(itsChi2_19);
   
   AliReducedCompositeCut* tpcChi2_60 = new AliReducedCompositeCut(Form("%s_tpcChi2_60", settingName), "", kTRUE);
   tpcChi2_60->AddCut(commonVarCuts);
   AliReducedTrackCut* tpcChi2CommonExtraCuts = new AliReducedTrackCut();
   if(useSPDany) tpcChi2CommonExtraCuts->SetRequestSPDany();
   else                tpcChi2CommonExtraCuts->SetRequestSPDfirst();
   tpcChi2CommonExtraCuts->AddCut(AliReducedVarManager::kITSchi2, 0., 15.);
   tpcChi2CommonExtraCuts->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 6., 9.);
   tpcChi2_60->AddCut(tpcChi2CommonExtraCuts);
   AliReducedTrackCut* tpcChi2_60_cut = new AliReducedTrackCut();
   tpcChi2_60_cut->AddCut(AliReducedVarManager::kTPCchi2, 0., 6.);
   tpcChi2_60->AddCut(tpcChi2_60_cut);
   processor->AddTrackCut(tpcChi2_60);
   
   AliReducedCompositeCut* tpcChi2_50 = new AliReducedCompositeCut(Form("%s_tpcChi2_50", settingName), "", kTRUE);
   tpcChi2_50->AddCut(commonVarCuts);
   tpcChi2_50->AddCut(tpcChi2CommonExtraCuts);
   AliReducedTrackCut* tpcChi2_50_cut = new AliReducedTrackCut();
   tpcChi2_50_cut->AddCut(AliReducedVarManager::kTPCchi2, 0., 5.0);
   tpcChi2_50->AddCut(tpcChi2_50_cut);
   processor->AddTrackCut(tpcChi2_50);
   
   AliReducedCompositeCut* tpcChi2_35 = new AliReducedCompositeCut(Form("%s_tpcChi2_35", settingName), "", kTRUE);
   tpcChi2_35->AddCut(commonVarCuts);
   tpcChi2_35->AddCut(tpcChi2CommonExtraCuts);
   AliReducedTrackCut* tpcChi2_35_cut = new AliReducedTrackCut();
   tpcChi2_35_cut->AddCut(AliReducedVarManager::kTPCchi2, 0., 3.5);
   tpcChi2_35->AddCut(tpcChi2_35_cut);
   processor->AddTrackCut(tpcChi2_35);
   
   AliReducedCompositeCut* tpcChi2_30 = new AliReducedCompositeCut(Form("%s_tpcChi2_30", settingName), "", kTRUE);
   tpcChi2_30->AddCut(commonVarCuts);
   tpcChi2_30->AddCut(tpcChi2CommonExtraCuts);
   AliReducedTrackCut* tpcChi2_30_cut = new AliReducedTrackCut();
   tpcChi2_30_cut->AddCut(AliReducedVarManager::kTPCchi2, 0., 3.0);
   tpcChi2_30->AddCut(tpcChi2_30_cut);
   processor->AddTrackCut(tpcChi2_30);
   
   AliReducedCompositeCut* tpcChi2_25 = new AliReducedCompositeCut(Form("%s_tpcChi2_25", settingName), "", kTRUE);
   tpcChi2_25->AddCut(commonVarCuts);
   tpcChi2_25->AddCut(tpcChi2CommonExtraCuts);
   AliReducedTrackCut* tpcChi2_25_cut = new AliReducedTrackCut();
   tpcChi2_25_cut->AddCut(AliReducedVarManager::kTPCchi2, 0., 2.5);
   tpcChi2_25->AddCut(tpcChi2_25_cut);
   processor->AddTrackCut(tpcChi2_25);
   
   AliReducedCompositeCut* tpcSegm4 = new AliReducedCompositeCut(Form("%s_tpcSegm4", settingName), "", kTRUE);
   tpcSegm4->AddCut(commonVarCuts);
   AliReducedTrackCut* tpcSegmCommonExtraCuts = new AliReducedTrackCut();
   if(useSPDany) tpcSegmCommonExtraCuts->SetRequestSPDany();
   else                tpcSegmCommonExtraCuts->SetRequestSPDfirst();
   tpcSegmCommonExtraCuts->AddCut(AliReducedVarManager::kITSchi2, 0., 15.);
   tpcSegmCommonExtraCuts->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
   tpcSegm4->AddCut(tpcSegmCommonExtraCuts);
   AliReducedTrackCut* tpcSegm4_cut = new AliReducedTrackCut();
   tpcSegm4_cut->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 4., 9.);
   tpcSegm4->AddCut(tpcSegm4_cut);
   processor->AddTrackCut(tpcSegm4);
   
   AliReducedCompositeCut* tpcSegm5 = new AliReducedCompositeCut(Form("%s_tpcSegm5", settingName), "", kTRUE);
   tpcSegm5->AddCut(commonVarCuts);
   tpcSegm5->AddCut(tpcSegmCommonExtraCuts);
   AliReducedTrackCut* tpcSegm5_cut = new AliReducedTrackCut();
   tpcSegm5_cut->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 5., 9.);
   tpcSegm5->AddCut(tpcSegm5_cut);
   processor->AddTrackCut(tpcSegm5);
   
   AliReducedCompositeCut* tpcSegm7 = new AliReducedCompositeCut(Form("%s_tpcSegm7", settingName), "", kTRUE);
   tpcSegm7->AddCut(commonVarCuts);
   tpcSegm7->AddCut(tpcSegmCommonExtraCuts);
   AliReducedTrackCut* tpcSegm7_cut = new AliReducedTrackCut();
   tpcSegm7_cut->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 7., 9.);
   tpcSegm7->AddCut(tpcSegm7_cut);
   processor->AddTrackCut(tpcSegm7);
   
   AliReducedCompositeCut* tpcCls90 = new AliReducedCompositeCut(Form("%s_tpcCls90", settingName), "", kTRUE);
   tpcCls90->AddCut(commonVarCuts);
   AliReducedTrackCut* tpcClsCommonExtraCuts = new AliReducedTrackCut();
   if(useSPDany) tpcClsCommonExtraCuts->SetRequestSPDany();
   else                tpcClsCommonExtraCuts->SetRequestSPDfirst();
   tpcClsCommonExtraCuts->AddCut(AliReducedVarManager::kITSchi2, 0., 15.);
   tpcClsCommonExtraCuts->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
   tpcClsCommonExtraCuts->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 6., 9.);
   tpcCls90->AddCut(tpcClsCommonExtraCuts);
   AliReducedTrackCut* tpcCls90_cut = new AliReducedTrackCut();
   tpcCls90_cut->AddCut(AliReducedVarManager::kTPCncls, 90.,160.0);
   tpcCls90->AddCut(tpcCls90_cut);
   processor->AddTrackCut(tpcCls90);
   
   AliReducedCompositeCut* tpcCls100 = new AliReducedCompositeCut(Form("%s_tpcCls100", settingName), "", kTRUE);
   tpcCls100->AddCut(commonVarCuts);
   tpcCls100->AddCut(tpcClsCommonExtraCuts);
   AliReducedTrackCut* tpcCls100_cut = new AliReducedTrackCut();
   tpcCls100_cut->AddCut(AliReducedVarManager::kTPCncls, 100.,160.0);
   tpcCls100->AddCut(tpcCls100_cut);
   processor->AddTrackCut(tpcCls100);
   
   AliReducedCompositeCut* tpcCls110 = new AliReducedCompositeCut(Form("%s_tpcCls110", settingName), "", kTRUE);
   tpcCls110->AddCut(commonVarCuts);
   tpcCls110->AddCut(tpcClsCommonExtraCuts);
   AliReducedTrackCut* tpcCls110_cut = new AliReducedTrackCut();
   tpcCls110_cut->AddCut(AliReducedVarManager::kTPCncls, 110.,160.0);
   tpcCls110->AddCut(tpcCls110_cut);
   processor->AddTrackCut(tpcCls110);
   
   AliReducedCompositeCut* tpcCls115 = new AliReducedCompositeCut(Form("%s_tpcCls115", settingName), "", kTRUE);
   tpcCls115->AddCut(commonVarCuts);
   tpcCls115->AddCut(tpcClsCommonExtraCuts);
   AliReducedTrackCut* tpcCls115_cut = new AliReducedTrackCut();
   tpcCls115_cut->AddCut(AliReducedVarManager::kTPCncls, 115.,160.0);
   tpcCls115->AddCut(tpcCls115_cut);
   processor->AddTrackCut(tpcCls115);
   
   AliReducedCompositeCut* tpcCls120 = new AliReducedCompositeCut(Form("%s_tpcCls120", settingName), "", kTRUE);
   tpcCls120->AddCut(commonVarCuts);
   tpcCls120->AddCut(tpcClsCommonExtraCuts);
   AliReducedTrackCut* tpcCls120_cut = new AliReducedTrackCut();
   tpcCls120_cut->AddCut(AliReducedVarManager::kTPCncls, 120.,160.0);
   tpcCls120->AddCut(tpcCls120_cut);
   processor->AddTrackCut(tpcCls120);
   
   AliReducedCompositeCut* tpcCls125 = new AliReducedCompositeCut(Form("%s_tpcCls125", settingName), "", kTRUE);
   tpcCls125->AddCut(commonVarCuts);
   tpcCls125->AddCut(tpcClsCommonExtraCuts);
   AliReducedTrackCut* tpcCls125_cut = new AliReducedTrackCut();
   tpcCls125_cut->AddCut(AliReducedVarManager::kTPCncls, 125.,160.0);
   tpcCls125->AddCut(tpcCls125_cut);
   //processor->AddTrackCut(tpcCls125);
   
   AliReducedCompositeCut* tpcCls130 = new AliReducedCompositeCut(Form("%s_tpcCls130", settingName), "", kTRUE);
   tpcCls130->AddCut(commonVarCuts);
   tpcCls130->AddCut(tpcClsCommonExtraCuts);
   AliReducedTrackCut* tpcCls130_cut = new AliReducedTrackCut();
   tpcCls130_cut->AddCut(AliReducedVarManager::kTPCncls, 130.,160.0);
   tpcCls130->AddCut(tpcCls130_cut);
   //processor->AddTrackCut(tpcCls130);
}


//__________________________________________________________________________________________
void SetupPIDSystematicsCuts(AliReducedAnalysisJpsi2ee* task, Bool_t useSPDany) {
   //
   // add pid systematics cuts
   //
   AliReducedTrackCut* commonVarCuts = new AliReducedTrackCut("commonVarCuts","");
   commonVarCuts->SetTrackFilterBit(kBasicCut);
   commonVarCuts->AddCut(AliReducedVarManager::kPt, 1.0, 100.0);  
   commonVarCuts->SetRejectKinks();
   commonVarCuts->AddCut(AliReducedVarManager::kTPCncls, 70., 161.);
   commonVarCuts->AddCut(AliReducedVarManager::kITSchi2, 0., 15.);
   commonVarCuts->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
   commonVarCuts->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 6., 9.);
   if(useSPDany) commonVarCuts->SetRequestSPDany();
   else                commonVarCuts->SetRequestSPDfirst();
   
   Double_t eleCuts[3] = {3.0, 2.8, 2.6};
   Double_t protCuts[4] = {3.6, 3.8, 4.0, 4.2};
   Double_t pionCuts[5] = {3.2, 3.4, 3.6, 3.8, 4.0};
   for(Int_t iele=0; iele<3; ++iele) {
      for(Int_t ipro=0; ipro<4; ++ipro) {
         for(Int_t ipio=0; ipio<5; ++ipio) {
            
            AliReducedCompositeCut* cut = new AliReducedCompositeCut(Form("electron%.0f_prot%.0f_pion%.0f", 
                                                                          eleCuts[iele]*10, protCuts[ipro]*10, pionCuts[ipio]*10), "", kTRUE);
            cut->AddCut(commonVarCuts);
            AliReducedTrackCut* varCut = new AliReducedTrackCut(Form("electron%.0f_prot%.0f_pion%.0f_varCut", 
                                                                     eleCuts[iele]*10, protCuts[ipro]*10, pionCuts[ipio]*10), "");
            if(task->GetRunOverMC())
               varCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -1.0*eleCuts[iele], eleCuts[iele]);
            else
               varCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -1.0*eleCuts[iele], eleCuts[iele]);
            
            varCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, protCuts[ipro], 30000.0);
            varCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, pionCuts[ipio], 30000.0);
            cut->AddCut(varCut);
            
            task->AddTrackCut(cut);
         }  // end loop over pion rej settings
      }  // end loop over proton rej settings
   }  // end loop over electron inclusion setting
}

//__________________________________________________________________________________________
void SetupPIDSystematicsCutsSingleLeg(AliReducedAnalysisJpsi2ee* processor, Bool_t calibrated) {
   //
   //
   //
   Bool_t isMC = processor->GetRunOverMC();
   
   AliReducedTrackCut* conversionElectron = new AliReducedTrackCut("conversionElectron","");
   conversionElectron->SetTrackQualityFilterBit(kPhotonConversionV0Leg, kFALSE);
   
   AliReducedTrackCut* commonVarCuts = new AliReducedTrackCut("commonVarCuts","");
   commonVarCuts->SetTrackFilterBit(kBasicCut);
   //commonVarCuts->AddCut(AliReducedVarManager::kPt, 1.0, 100.0);  
   //commonVarCuts->SetRejectKinks();
   commonVarCuts->AddCut(AliReducedVarManager::kTPCncls, 70., 161.);
   //commonVarCuts->AddCut(AliReducedVarManager::kITSchi2, 0., 15.);
   //commonVarCuts->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
   commonVarCuts->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 6., 9.);
   //if(useSPDany) commonVarCuts->SetRequestSPDany();
   //else                commonVarCuts->SetRequestSPDfirst();
   
   AliReducedCompositeCut* basicElectronSelection = new AliReducedCompositeCut("basicElectronCut", "");
   if(!isMC) basicElectronSelection->AddCut(conversionElectron);
   basicElectronSelection->AddCut(commonVarCuts);
   processor->AddTrackCut(basicElectronSelection);
   
   AliReducedCompositeCut* electron30 = new AliReducedCompositeCut(Form("electron30%s", (!calibrated ? "_uncalibrated" : "")), "", kTRUE);
   electron30->AddCut(basicElectronSelection);
   AliReducedTrackCut* electron30VarCut = new AliReducedTrackCut();
   if(isMC || !calibrated)  electron30VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
   else       electron30VarCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
   electron30VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 4.0, 30000.0);
   electron30VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.6, 30000.0);
   electron30->AddCut(electron30VarCut);
   processor->AddTrackCut(electron30);
   
   AliReducedCompositeCut* electron28 = new AliReducedCompositeCut(Form("electron28%s", (!calibrated ? "_uncalibrated" : "")), "", kTRUE);
   electron28->AddCut(basicElectronSelection);
   AliReducedTrackCut* electron28VarCut = new AliReducedTrackCut();
   if(isMC || !calibrated)   electron28VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2.8, 2.8);
   else        electron28VarCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -2.8, 2.8);
   electron28VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 4.0, 30000.0);
   electron28VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.6, 30000.0);
   electron28->AddCut(electron28VarCut);
   processor->AddTrackCut(electron28);
   
   AliReducedCompositeCut* electron26 = new AliReducedCompositeCut(Form("electron26%s", (!calibrated ? "_uncalibrated" : "")), "", kTRUE);
   electron26->AddCut(basicElectronSelection);
   AliReducedTrackCut* electron26VarCut = new AliReducedTrackCut();
   if(isMC || !calibrated)   electron26VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2.6, 2.6);
   else        electron26VarCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -2.6, 2.6);
   electron26VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 4.0, 30000.0);
   electron26VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.6, 30000.0);
   electron26->AddCut(electron26VarCut);
   processor->AddTrackCut(electron26);
   
   AliReducedCompositeCut* proton42 = new AliReducedCompositeCut(Form("proton42%s", (!calibrated ? "_uncalibrated" : "")), "", kTRUE);
   proton42->AddCut(basicElectronSelection);
   AliReducedTrackCut* proton42VarCut = new AliReducedTrackCut();
   if(isMC || !calibrated)   proton42VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
   else        proton42VarCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
   proton42VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 4.2, 30000.0);
   proton42VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.6, 30000.0);
   proton42->AddCut(proton42VarCut);
   processor->AddTrackCut(proton42);
   
   AliReducedCompositeCut* proton40 = new AliReducedCompositeCut(Form("proton40%s", (!calibrated ? "_uncalibrated" : "")), "", kTRUE);
   proton40->AddCut(basicElectronSelection);
   AliReducedTrackCut* proton40VarCut = new AliReducedTrackCut();
   if(isMC || !calibrated)   proton40VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
   else        proton40VarCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
   proton40VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 4.0, 30000.0);
   proton40VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.6, 30000.0);
   proton40->AddCut(proton40VarCut);
   processor->AddTrackCut(proton40);
   
   AliReducedCompositeCut* proton38 = new AliReducedCompositeCut(Form("proton38%s", (!calibrated ? "_uncalibrated" : "")), "", kTRUE);
   proton38->AddCut(basicElectronSelection);
   AliReducedTrackCut* proton38VarCut = new AliReducedTrackCut();
   if(isMC || !calibrated)   proton38VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
   else        proton38VarCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
   proton38VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.8, 30000.0);
   proton38VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.6, 30000.0);
   proton38->AddCut(proton38VarCut);
   processor->AddTrackCut(proton38);
   
   AliReducedCompositeCut* proton36 = new AliReducedCompositeCut(Form("proton36%s", (!calibrated ? "_uncalibrated" : "")), "", kTRUE);
   proton36->AddCut(basicElectronSelection);
   AliReducedTrackCut* proton36VarCut = new AliReducedTrackCut();
   if(isMC || !calibrated)   proton36VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
   else        proton36VarCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
   proton36VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.6, 30000.0);
   proton36VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.6, 30000.0);
   proton36->AddCut(proton36VarCut);
   processor->AddTrackCut(proton36);
   
   AliReducedCompositeCut* pion40 = new AliReducedCompositeCut(Form("pion40%s", (!calibrated ? "_uncalibrated" : "")), "", kTRUE);
   pion40->AddCut(basicElectronSelection);
   AliReducedTrackCut* pion40VarCut = new AliReducedTrackCut();
   if(isMC || !calibrated)   pion40VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
   else        pion40VarCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
   pion40VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 4.0, 30000.0);
   pion40VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 4.0, 30000.0);
   pion40->AddCut(pion40VarCut);
   processor->AddTrackCut(pion40);
   
   AliReducedCompositeCut* pion38 = new AliReducedCompositeCut(Form("pion38%s", (!calibrated ? "_uncalibrated" : "")), "", kTRUE);
   pion38->AddCut(basicElectronSelection);
   AliReducedTrackCut* pion38VarCut = new AliReducedTrackCut();
   if(isMC || !calibrated)   pion38VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
   else        pion38VarCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
   pion38VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 4.0, 30000.0);
   pion38VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.8, 30000.0);
   pion38->AddCut(pion38VarCut);
   processor->AddTrackCut(pion38);
   
   AliReducedCompositeCut* pion36 = new AliReducedCompositeCut(Form("pion36%s", (!calibrated ? "_uncalibrated" : "")), "", kTRUE);
   pion36->AddCut(basicElectronSelection);
   AliReducedTrackCut* pion36VarCut = new AliReducedTrackCut();
   if(isMC || !calibrated)   pion36VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
   else        pion36VarCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
   pion36VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 4.0, 30000.0);
   pion36VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.6, 30000.0);
   pion36->AddCut(pion36VarCut);
   processor->AddTrackCut(pion36);
   
   AliReducedCompositeCut* pion34 = new AliReducedCompositeCut(Form("pion34%s", (!calibrated ? "_uncalibrated" : "")), "", kTRUE);
   pion34->AddCut(basicElectronSelection);
   AliReducedTrackCut* pion34VarCut = new AliReducedTrackCut();
   if(isMC || !calibrated)   pion34VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
   else        pion34VarCut->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 3.0);
   pion34VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 4.0, 30000.0);
   pion34VarCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.4, 30000.0);
   pion34->AddCut(pion34VarCut);
   processor->AddTrackCut(pion34);
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
   
   handler->AddMixingVariable(AliReducedVarManager::kCentVZERO, gkNCentBins, gCentLims);
   //handler->AddMixingVariable(AliReducedVarManager::kVtxZ, gkNVtxZBins, gVtxZLims);
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
      histClasses += Form("TrackTPCclusterMap_%s;", cutName.Data());
      if(task->GetRunOverMC()) {
         for(Int_t mcSel=0; mcSel<task->GetNLegCandidateMCcuts(); ++mcSel) {
            histClasses += Form("Track_%s_%s;", cutName.Data(), task->GetLegCandidateMCcutName(mcSel));
            histClasses += Form("TrackStatusFlags_%s_%s;", cutName.Data(), task->GetLegCandidateMCcutName(mcSel));
            histClasses += Form("TrackITSclusterMap_%s_%s;", cutName.Data(), task->GetLegCandidateMCcutName(mcSel));
            histClasses += Form("TrackTPCclusterMap_%s_%s;", cutName.Data(), task->GetLegCandidateMCcutName(mcSel));
         }
      }
      if(!task->GetRunLikeSignPairing())
         histClasses += Form("PairSEPM_%s;", cutName.Data());
      else
         histClasses += Form("PairSEPP_%s;PairSEPM_%s;PairSEMM_%s;", cutName.Data(), cutName.Data(), cutName.Data());
      if(!task->GetRunOverMC())
        histClasses += Form("PairMEPP_%s;PairMEPM_%s;PairMEMM_%s;", cutName.Data(), cutName.Data(), cutName.Data());
      if(task->GetRunOverMC()) {
         for(Int_t mcSel=0; mcSel<task->GetNLegCandidateMCcuts(); ++mcSel)
            histClasses += Form("PairSEPM_%s_%s;", cutName.Data(), task->GetLegCandidateMCcutName(mcSel));
      }
    }
  }
  
  AliReducedVarManager::SetRunNumbers(gRunList);
  
  const Int_t kNBCBins = 3600;
  Double_t bcHistRange[2] = {-0.5,3599.5};
  
  // set the mass bin limits
  for(Int_t i=0; i<gkNMassBins;++i) gMassBins[i] = 1.0+i*0.04;
  
  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");
  
  cout << "Histogram classes included in the Histogram Manager" << endl;
  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    
    if(classStr.Contains("PureMCTruth_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       
       man->AddHistogram(classStr.Data(), "MassMC_Pt_CentVZERO", "", kFALSE, gkNMassBins-1, gMassBins, AliReducedVarManager::kMassMC, gkNPtBins-1, gPtLims, AliReducedVarManager::kPtMC, gkNCentBins-1, gCentLims, AliReducedVarManager::kCentVZERO);
       
       man->AddHistogram(classStr.Data(), "MassMC", "MC mass", kFALSE, 200, 0., 5.0, AliReducedVarManager::kMassMC);
       man->AddHistogram(classStr.Data(), "RapidityMC", "MC rapidity", kFALSE, 48, -1.2, 1.2, AliReducedVarManager::kRapMC);
       man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE, 1000, 0., 10.0, AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "PtMC_coarse", "p_{T} MC", kFALSE, 20, 0., 20.0, AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "PhiMC", "MC #varphi", kFALSE, 100, 0., 6.3, AliReducedVarManager::kPhiMC);
       man->AddHistogram(classStr.Data(), "EtaMC", "MC #eta", kFALSE, 100, -1.5, 1.5, AliReducedVarManager::kEtaMC);
       man->AddHistogram(classStr.Data(), "PtMC_RapMC", "", kFALSE, 100, -1.2, 1.2, AliReducedVarManager::kRapMC, 100, 0., 15., AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "CosThetaStarCS", "cos(#theta^{*})_{CS}", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaCS);
       man->AddHistogram(classStr.Data(), "CosThetaStarCS_ptMC", "cos(#theta^{*})_{CS} vs MC pt", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaCS, 50, 0., 1., AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "CosThetaStarHE", "cos(#theta^{*})_{HE}", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaHE);
       man->AddHistogram(classStr.Data(), "CosThetaStarHE_ptMC", "cos(#theta^{*})_{HE} vs MC pt", kFALSE, 22, -1.1, 1.1, AliReducedVarManager::kPairThetaHE, 100, 0., 1., AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "CosThetaStarHE_ptMC_coarse", "cos(#theta^{*})_{HE} vs MC pt", kFALSE, 10, -1.0, 1.0, AliReducedVarManager::kPairThetaHE, 20, 0., 20., AliReducedVarManager::kPtMC);
       man->AddHistogram(classStr.Data(), "PhiStarCS", "#varphi^{*}_{CS}", kFALSE, 22, -3.3, 3.3, AliReducedVarManager::kPairPhiCS);
       man->AddHistogram(classStr.Data(), "PhiStarHE", "#varphi^{*}_{HE}", kFALSE, 22, -3.3, 3.3, AliReducedVarManager::kPairPhiHE);
       
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
       man->AddHistogram(classStr.Data(), "EventTags", "Event tags", kFALSE,
                         20, -0.5, 19.5, AliReducedVarManager::kEventTag, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, tagNames.Data());
       man->AddHistogram(classStr.Data(), "EventTags_CentVZERO", "Event tags vs VZERO centrality", kFALSE,
                         20, -0.5, 19.5, AliReducedVarManager::kEventTag, 9, 0.0, 90.0, AliReducedVarManager::kCentVZERO, 0, 0.0, 0.0, AliReducedVarManager::kNothing, tagNames.Data());
       continue;
    }
    
    if(classStr.Contains("EventTriggers_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       TString triggerNames = "";
       for(Int_t i=0; i<64; ++i) {triggerNames += AliReducedVarManager::fgkOfflineTriggerNames[i]; triggerNames+=";";}
       
       man->AddHistogram(classStr.Data(), "Triggers2", "", kFALSE,
                         64, -0.5, 63.5, AliReducedVarManager::kOnlineTriggerFired2, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, triggerNames.Data());
       man->AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "", kFALSE,
                         64, -0.5, 63.5, AliReducedVarManager::kOnlineTriggerFired2, 9, 0.0, 90.0, AliReducedVarManager::kCentVZERO, 0, 0.0, 0.0, AliReducedVarManager::kNothing, triggerNames.Data());
       continue;
    }
    
    if(classStr.Contains("Event_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(),"RunNo","Run numbers",kFALSE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID);
      man->AddHistogram(classStr.Data(),"VtxX","Vtx X",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxY","Vtx Y",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxZ","Vtx Z",kFALSE,300,-15.,15.,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentVZERO);
      man->AddHistogram(classStr.Data(),"CentVZERO_VtxZ","Centrality(VZERO) vs vtxZ",kFALSE, 50, 0.0, 100.0, AliReducedVarManager::kCentVZERO,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentSPD);
      man->AddHistogram(classStr.Data(),"CentSPD_VtxZ","Centrality(SPD) vs vtxZ",kFALSE, 50, 0.0, 100.0, AliReducedVarManager::kCentSPD,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentQuality","Centrality quality",kFALSE, 500, 0., 500., AliReducedVarManager::kCentQuality);
      man->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0", kFALSE, 200, 0., 5000., AliReducedVarManager::kSPDntracklets);

      man->AddHistogram(classStr.Data(),"VZEROmult", "", kFALSE, 500, 0.0, 50000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"VZEROmult_VtxContributors", "", kFALSE, 
		   200, 0.0, 40000., AliReducedVarManager::kVZEROTotalMult, 200, 0.0, 10000., AliReducedVarManager::kNVtxContributors);
      man->AddHistogram(classStr.Data(),"VZEROmult_SPDntracklets", "", kFALSE, 
                        200, 0.0, 5000., AliReducedVarManager::kSPDntracklets, 200, 0.0, 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsPhysicsSelection, 0,0.0,0.0,AliReducedVarManager::kNothing, 0,0.0,0.0,AliReducedVarManager::kNothing, "off;on");

      man->AddHistogram(classStr.Data(),"NPosTracksAnalyzed","#positrons per event",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksPosAnalyzed);
      man->AddHistogram(classStr.Data(),"NNegTracksAnalyzed","#electrons per event",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksNegAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalTracksAnalyzed","#leg candidates per event",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalPairsAnalyzed","#dielectron candidates per event",kFALSE,100,0.,100.,AliReducedVarManager::kNpairsSelected);

      man->AddHistogram(classStr.Data(),"NTracksTPCout","",kFALSE,1800,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout);
      man->AddHistogram(classStr.Data(),"NTracksTPCout_VZEROmult","",kFALSE,100,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 100, 0., 30000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTotalTracksAnalyzed_NTracksTPCout_VZEROmult_prof","",kTRUE,20,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20, 0., 30000., AliReducedVarManager::kVZEROTotalMult, 10, 0., 10., AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalPairsAnalyzed_NTracksTPCout_VZEROmult_prof","",kTRUE,20,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20, 0., 30000., AliReducedVarManager::kVZEROTotalMult, 10, 0., 10., AliReducedVarManager::kNpairsSelected);
      
      man->AddHistogram(classStr.Data(),"NTracksAnalyzed_VZEROmult","",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksAnalyzed, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutVsVZEROmult","",kFALSE,1000,0.,10.,AliReducedVarManager::kNTracksTPCoutVsVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutVsVZEROmult_CentVZERO","",kFALSE, 20, 0., 100., AliReducedVarManager::kCentVZERO, 1000,0.,10.,AliReducedVarManager::kNTracksTPCoutVsVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksAnalyzed_RunNo_prof","",kTRUE, gkNRuns, -0.5, -0.5 + gkNRuns, AliReducedVarManager::kRunID, 10, 0., 10., AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NPairsAnalyzed_RunNo_prof","",kTRUE, gkNRuns, -0.5, -0.5 + gkNRuns, AliReducedVarManager::kRunID, 10, 0., 10., AliReducedVarManager::kNpairsSelected);
      continue;
    }  // end if className contains "Event"    
    
    // Track histograms
    if(classStr.Contains("TrackITSclusterMap_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       man->AddHistogram(classStr.Data(), "ITSlayerHit", "Hits in the ITS layers", kFALSE,
                         6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
       //man->AddHistogram(classStr.Data(), "ITSlayerHit_Phi", "Hits in the ITS layers vs #varphi", kFALSE,
         //                180, 0.0, 6.29, AliReducedVarManager::kPhi, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
       //man->AddHistogram(classStr.Data(), "ITSlayerHit_Eta", "Hits in the ITS layers vs #eta", kFALSE,
         //                100, -1.0, 1.0, AliReducedVarManager::kEta, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
       continue;
    }  // end of ITSclusterMap histogram definitions
    
    if(classStr.Contains("TrackTPCclusterMap_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       man->AddHistogram(classStr.Data(), "TPCclusterMap", "TPC cluster map", kFALSE,
                         8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
       //man->AddHistogram(classStr.Data(), "TPCclusterMap_Phi", "TPC cluster map vs #varphi", kFALSE,
        //                 180, 0.0, 6.29, AliReducedVarManager::kPhi, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
       //man->AddHistogram(classStr.Data(), "TPCclusterMap_Eta", "TPC cluster map vs #eta", kFALSE,
         //                100, -1.0, 1.0, AliReducedVarManager::kEta, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
       //man->AddHistogram(classStr.Data(), "TPCclusterMap_Pt", "TPC cluster map vs p_{T}", kFALSE,
         //                100, 0.0, 10.0, AliReducedVarManager::kPt, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
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
      man->AddHistogram(classStr.Data(), "Eta_Phi", "", kFALSE, 40, -1.0, 1.0, AliReducedVarManager::kEta, 36, 0., 6.29, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "Eta_P", "", kFALSE, 40, -1.0, 1.0, AliReducedVarManager::kEta, 50, 0., 5.0, AliReducedVarManager::kP);
      man->AddHistogram(classStr.Data(), "Phi", "", kFALSE, 36, 0.0, 6.29, AliReducedVarManager::kPhi);
      //man->AddHistogram(classStr.Data(), "DCAxy_Pt", "DCAxy", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 50, 0.0, 5.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE, 100, -2.5, 2.5, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kDcaZ);
      //man->AddHistogram(classStr.Data(), "DCAz_Pt", "DCAz", kFALSE, 100, -3.0, 3.0, AliReducedVarManager::kDcaZ, 50, 0.0, 5.0, AliReducedVarManager::kPt);
      
        man->AddHistogram(classStr.Data(),"ITSncls", "ITS nclusters", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"ITSnclsShared", "ITS nclusters shared", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSnclsShared);
        /*man->AddHistogram(classStr.Data(),"Eta_Phi_ITSncls_prof","ITS <nclusters> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"Eta_Phi_ITSnclsShared_prof","ITS <nclusters-shared> vs (#eta,#phi)",kTRUE,
                          36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kITSnclsShared);*/
	man->AddHistogram(classStr.Data(),"ITSchi2", "ITS #chi^{2}", kFALSE, 200,0.0,50.0, AliReducedVarManager::kITSchi2);
        //man->AddHistogram(classStr.Data(),"ITSchi2_RunID_prof","",kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 10, 0., 10., AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"ITSchi2_ITSncls", "ITS #chi^{2} vs ITS clusters", kFALSE, 100,0.0,20.0, AliReducedVarManager::kITSchi2,7,-0.5,6.5, AliReducedVarManager::kITSncls);
	man->AddHistogram(classStr.Data(),"Eta_Phi_ITSchi2_prof","ITS <#chi^{2}> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 100, 0.0, 1000, AliReducedVarManager::kITSchi2);        
        man->AddHistogram(classStr.Data(),"ITSchi2_NTracksTPCout_VZEROmult_prof","",kTRUE,20,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20, 0., 30000., AliReducedVarManager::kVZEROTotalMult, 10, 0., 10., AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"Chi2TPCConstrainedVsGlobal","",kFALSE,1000,0.,200.,AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
        
        man->AddHistogram(classStr.Data(),"TPCncls","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCncls);
        //man->AddHistogram(classStr.Data(),"TPCncls_RunID_prof","",kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 10, 0., 10., AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCncls_NTracksTPCout_VZEROmult_prof","",kTRUE,20,0.,10000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20, 0., 30000., AliReducedVarManager::kVZEROTotalMult, 10, 0., 10., AliReducedVarManager::kTPCncls);
	man->AddHistogram(classStr.Data(),"TPCcrossedRows","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCcrossedRows);
        man->AddHistogram(classStr.Data(),"TPCfindableClusters","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsF);
        man->AddHistogram(classStr.Data(),"TPCcrossedRowsOverFindableClusters","", kFALSE, 200,0.0,1.5,AliReducedVarManager::kTPCcrossedRowsOverFindableClusters);
	man->AddHistogram(classStr.Data(),"TPCnclsShared","",kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsShared);
        man->AddHistogram(classStr.Data(),"TPCnclsSharedRatio","",kFALSE, 200,0.,1.,AliReducedVarManager::kTPCnclsSharedRatio);        
        
        man->AddHistogram(classStr.Data(),"TPCchi2","",kFALSE, 200,0.0,10.0,AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_Pt","",kFALSE, 200,0.0,10.0,AliReducedVarManager::kTPCchi2, 100, 0., 10.0, AliReducedVarManager::kPt);
        man->AddHistogram(classStr.Data(),"TPCchi2_TPCncls", "TPC #chi^{2} vs TPC clusters", kFALSE, 100,0.0,10.0, AliReducedVarManager::kTPCchi2,160,0.,160., AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCchi2_CentVZERO", "TPC #chi^{2} vs CentVZERO", kFALSE, 100,0.0,10.0, AliReducedVarManager::kTPCchi2,20,0.,100., AliReducedVarManager::kCentVZERO);
        man->AddHistogram(classStr.Data(),"Eta_TPCchi2","TPC #chi^{2} vs #eta",kFALSE,
                          36, -0.9, 0.9, AliReducedVarManager::kEta, 100, 0.0, 10., AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCsegments","",kFALSE, 9,0.0,9.0,AliReducedVarManager::kTPCNclusBitsFired);
        
        /*man->AddHistogram(classStr.Data(),"Eta_Phi_TPCncls_prof","TPC <nclusters> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCncls);    
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCcrossedRows_prof","TPC <n crossed rows> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCcrossedRows);
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCnclsF_prof","TPC <nclusters findable> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCnclsF);*/
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCchi2_prof","TPC <#chi^{2}> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCchi2);
        
        man->AddHistogram(classStr.Data(),"TPCchi2_NTracksTPCout_VZEROmult_prof","",kTRUE,20,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20, 0., 30000., AliReducedVarManager::kVZEROTotalMult, 10, 0., 10., AliReducedVarManager::kTPCchi2);
        /*man->AddHistogram(classStr.Data(),"Eta_Phi_TPCchi2","TPC #chi^{2} vs (#eta,#phi)",kFALSE,
                          36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 100, 0., 10.0, AliReducedVarManager::kTPCchi2); */
        man->AddHistogram(classStr.Data(),"TPCsignal_Pin","TPC dE/dx vs. inner param P",kFALSE,
                     100,0.0,10.0,AliReducedVarManager::kPin,100,49.5,150.5,AliReducedVarManager::kTPCsignal);
	man->AddHistogram(classStr.Data(),"TPCsignalN","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCsignalN);
	/*man->AddHistogram(classStr.Data(),"Eta_Phi_TPCsignalN_prof","TPC <nclusters pid> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCsignalN);    */
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Pin","TPC N_{#sigma} electron vs. inner param P",kFALSE,
                     100,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Eta","TPC N_{#sigma} electron vs. #eta",kFALSE,
                          36,-0.9,0.9,AliReducedVarManager::kEta,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin","TPC N_{#sigma} electron corrected vs. inner param P",kFALSE,
                          100,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Eta","TPC N_{#sigma} electron corrected vs. #eta",kFALSE,
                          36,-0.9,0.9,AliReducedVarManager::kEta,100,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron);
        
        man->AddHistogram(classStr.Data(),"TPCchi2_vs_TPCnsigElectron_Pin_prof","<TPC chi2> vs TPC N_{#sigma} electron and inner param P",kTRUE,
                          20,0.0,5.0,AliReducedVarManager::kPin,30,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, 10, 0., 10., AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_vs_TPCnsigElectron_Pin","TPC chi2 vs TPC N_{#sigma} electron and inner param P",kFALSE,
                          20,0.0,5.0,AliReducedVarManager::kPin,30,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, 80, 0., 8., AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_vs_TPCnsigElectronCorrected_Pin_prof","<TPC chi2> vs TPC N_{#sigma} electron corrected and inner param P",kTRUE,
                          20,0.0,5.0,AliReducedVarManager::kPin,30,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 10, 0., 10., AliReducedVarManager::kTPCchi2);
        
        //man->AddHistogram(classStr.Data(),"TOFbeta_P","TOF #beta vs P",kFALSE,
          //           100,0.0,10.0,AliReducedVarManager::kP, 110,0.0,1.1,AliReducedVarManager::kTOFbeta);
        
        //man->AddHistogram(classStr.Data(),"TPCchi2_ITSchi2","",kFALSE, 100,0.0,10.0,AliReducedVarManager::kTPCchi2, 100,0.0,40.0,AliReducedVarManager::kITSchi2);
        //man->AddHistogram(classStr.Data(),"ITSchi2_Chi2TPCconstrainedVsGlobal","",kFALSE, 500,0.0,1000.0,AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 100,0.0,40.0,AliReducedVarManager::kITSchi2);
        //man->AddHistogram(classStr.Data(),"TPCchi2_Chi2TPCconstrainedVsGlobal","",kFALSE, 500,0.0,1000.0,AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 100,0.0,10.0,AliReducedVarManager::kTPCchi2);
        
        if(classStr.Contains("MCTruth")) {
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
        
    //Int_t vars[4] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, AliReducedVarManager::kCentVZERO, AliReducedVarManager::kVtxZ};
    Int_t vars[3] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, AliReducedVarManager::kCentVZERO};   
    
    TArrayD pairHistBinLimits[3];
    pairHistBinLimits[0] = TArrayD(gkNMassBins, gMassBins);
    pairHistBinLimits[1] = TArrayD(gkNPtBins, gPtLims);
    pairHistBinLimits[2] = TArrayD(gkNCentBins, gCentLims);
    //pairHistBinLimits[3] = TArrayD(gkNVtxZBins, gVtxZLims);
    
    // Histograms for pairs
    if(classStr.Contains("Pair")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "PairInvMass", "Differential pair inv. mass", 3, vars, pairHistBinLimits);
      if(classStr.Contains("PairSE") || classStr.Contains("PairPrefilterSE")) {
        man->AddHistogram(classStr.Data(), "PairType", "Pair type", kFALSE, 4, -0.5, 3.5, AliReducedVarManager::kPairType);
	man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE, gkNMassBins-1, gMassBins, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPt);
        man->AddHistogram(classStr.Data(), "Pt_coarse", "", kFALSE, 20, 0.0, 20.0, AliReducedVarManager::kPt);
	man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kRap);
	man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 315, 0.0, 6.3, AliReducedVarManager::kPhi);
      }   // end if "QA"
      
      if(classStr.Contains("MCTruth")) {
         man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE, 200, 0., 10.0, AliReducedVarManager::kPtMC);
         man->AddHistogram(classStr.Data(), "PtMC_coarse", "p_{T} MC", kFALSE, 20, 0., 20.0, AliReducedVarManager::kPtMC);
         //man->AddHistogram(classStr.Data(), "MassMC_Mass", "Invariant mass, MC vs reconstructed", kFALSE, 150, 2.0, 3.5, AliReducedVarManager::kMass,
          //  150, 2.0, 3.5, AliReducedVarManager::kMassMC);
         //man->AddHistogram(classStr.Data(), "PtMC_Pt", "pair pT, MC vs reconstructed", kFALSE, 150, 0.0, 15., AliReducedVarManager::kPt,
          //                 150, 0.0, 15.0, AliReducedVarManager::kPtMC);
         man->AddHistogram(classStr.Data(), "PtMC_MassMC", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPtMC,
                           30, 2.7, 3.3, AliReducedVarManager::kMassMC);
         man->AddHistogram(classStr.Data(), "Pt_Mass", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPt,
                           30, 2.7, 3.3, AliReducedVarManager::kMass);
         man->AddHistogram(classStr.Data(), "PtMC_Mass", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPtMC,
                           30, 2.7, 3.3, AliReducedVarManager::kMass);
         man->AddHistogram(classStr.Data(), "Pt_MassMC", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPt,
                           30, 2.7, 3.3, AliReducedVarManager::kMassMC);
         //man->AddHistogram(classStr.Data(), "PtMC_Pt_lowPtZoom", "pair pT, MC vs reconstructed", kFALSE, 100, 0.0, 0.5, AliReducedVarManager::kPt,
          //                 100, 0.0, 0.5, AliReducedVarManager::kPtMC);
         //man->AddHistogram(classStr.Data(), "PtMC_Pt_Mass", "pair pT, MC vs reconstructed, vs mass", kFALSE, 100, 0.0, 0.5, AliReducedVarManager::kPt,
          //                 100, 0.0, 0.5, AliReducedVarManager::kPtMC, 15, 2.72, 3.32, AliReducedVarManager::kMass);
         man->AddHistogram(classStr.Data(), "PtMC_Pt_Mass_coarse", "pair pT, MC vs reconstructed, vs mass", kFALSE, 25, 0.0, 0.5, AliReducedVarManager::kPt,
                           25, 0.0, 0.5, AliReducedVarManager::kPtMC, 15, 2.72, 3.32, AliReducedVarManager::kMass);
         
         man->AddHistogram(classStr.Data(), "Mass_Pt_CentVZERO", "", kFALSE, gkNMassBins-1, gMassBins, AliReducedVarManager::kMass, gkNPtBins-1, gPtLims, AliReducedVarManager::kPt, gkNCentBins-1, gCentLims, AliReducedVarManager::kCentVZERO);
      }
      
      continue;
    }   // end if for Pair classes of histograms
  }  // end loop over histogram classes
}
