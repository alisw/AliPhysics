//void Setup(AliReducedAnalysisJpsi2ee* processor, TString prod="LHC10h");
//AliHistogramManager* SetupHistogramManager(TString prod="LHC10h");
//void DefineHistograms(AliHistogramManager* man, TString prod="LHC10h");

enum TrackFilters {
   kBasicCut=0,
   kSPDany,
   kSPDfirst,
   kTPC100,
   kTPC120,
   kTPCchi2_3,
   kITSchi2_10,
   kEleInclusion_35,
   kEleInclusion_30,
   kProtRejection_30,
   kProtRejection_35,
   kProtRejection_40,
   kPionRejection_30,
   kPionRejection_35,
   kPionRejection_40,
   kNTrackFilters
};

enum MCFilters {
   kJpsiInclusive=0,
   kJpsiNonPrompt,
   kJpsiPrompt,
   kJpsiRadiative,
   kJpsiNonRadiative,
   kJpsiDecayElectron,
   kJpsiDecayPhoton
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
const Int_t gkNPtBins = 23;
Double_t gPtLims[gkNPtBins] = {
   0.0, 0.25, 0.50, 0.75, 1.0, 
   1.5, 2.0, 2.5, 3.0, 3.5, 
   4.0, 4.5, 5.0, 6.0, 7.0,
   8.0, 9.0, 10., 12., 15., 
   20., 25., 30. 
};

// mass binning
const Int_t gkNMassBins = 126;
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
  
     AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("jpsi2eeHistos", THashList::Class(),
                                                                  AliAnalysisManager::kOutputContainer, "AnalysisHistograms_jpsi2ee_XeXe.root");
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
  
  // Set event cuts
  AliReducedEventCut* evCut1 = new AliReducedEventCut("Centrality","Centrality selection");
  evCut1->AddCut(AliReducedVarManager::kCentVZERO, 0., 90.);
  evCut1->AddCut(AliReducedVarManager::kVtxZ, -10.0, 10.0);
  evCut1->AddCut(AliReducedVarManager::kIsPhysicsSelection, 0.1, 2.);   // request physics selection
  evCut1->EnableVertexDistanceCut();
  
  processor->AddEventCut(evCut1);
  
  /*
  // Set track cuts
  AliReducedTrackCut* standardCut = new AliReducedTrackCut("standard","");
  standardCut->AddCut(AliReducedVarManager::kP, 1.0,30.0);
  standardCut->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  standardCut->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  standardCut->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  standardCut->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 30000.0);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.0, 30000.0);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2.0, 3.0, kFALSE, AliReducedVarManager::kEta, -0.9, -0.75);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2.6, 3.0, kFALSE, AliReducedVarManager::kEta, 0.75, 0.9);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2.4, 3.0, kFALSE, AliReducedVarManager::kEta, -0.4, 0.4);
  standardCut->SetRejectKinks();
  standardCut->SetRequestITSrefit();
  standardCut->SetRequestTPCrefit();
  //standardCut->SetRequestTOFout();
  standardCut->SetRequestSPDany();
  standardCut->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  standardCut->AddCut(AliReducedVarManager::kITSchi2, 0., 10.);
  standardCut->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 6., 9.);
  standardCut->AddCut(AliReducedVarManager::kTPCnclsSharedRatio, 0.3, 2., kTRUE);
  standardCut->AddCut(AliReducedVarManager::kTPCcrossedRowsOverFindableClusters, 0.8, 2.);
  //standardCut->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 36., 1.0e+10, kTRUE);
  processor->AddTrackCut(standardCut); 
  */
  
  AliReducedTrackCut* standardCut = new AliReducedTrackCut("standardCut","");
  standardCut->AddCut(AliReducedVarManager::kPt, 1.0, 100.);
  standardCut->SetRejectKinks();
  standardCut->SetTrackFilterBit(kSPDany);
  standardCut->SetTrackFilterBit(kEleInclusion_30);
  standardCut->SetTrackFilterBit(kProtRejection_30);
  standardCut->SetTrackFilterBit(kPionRejection_30);  
  processor->AddTrackCut(standardCut);
  
  AliReducedTrackCut* spdFirstCut = new AliReducedTrackCut("spdFirstCut","");
  spdFirstCut->AddCut(AliReducedVarManager::kPt, 1.0, 100.);
  spdFirstCut->SetRejectKinks();
  spdFirstCut->SetTrackFilterBit(kSPDfirst);
  spdFirstCut->SetTrackFilterBit(kEleInclusion_30);
  spdFirstCut->SetTrackFilterBit(kProtRejection_30);
  spdFirstCut->SetTrackFilterBit(kPionRejection_30);  
  processor->AddTrackCut(spdFirstCut);
  
  AliReducedTrackCut* protRej30pionRej35 = new AliReducedTrackCut("protRej30pionRej35","");
  protRej30pionRej35->AddCut(AliReducedVarManager::kPt, 1.0, 100.);
  protRej30pionRej35->SetRejectKinks();
  protRej30pionRej35->SetTrackFilterBit(kSPDany);
  protRej30pionRej35->SetTrackFilterBit(kEleInclusion_30);
  protRej30pionRej35->SetTrackFilterBit(kProtRejection_30);
  protRej30pionRej35->SetTrackFilterBit(kPionRejection_35);
  processor->AddTrackCut(protRej30pionRej35);
  
  AliReducedTrackCut* protRej30pionRej40 = new AliReducedTrackCut("protRej30pionRej40","");
  protRej30pionRej40->AddCut(AliReducedVarManager::kPt, 1.0, 100.);
  protRej30pionRej40->SetRejectKinks();
  protRej30pionRej40->SetTrackFilterBit(kSPDany);
  protRej30pionRej40->SetTrackFilterBit(kEleInclusion_30);
  protRej30pionRej40->SetTrackFilterBit(kProtRejection_30);
  protRej30pionRej40->SetTrackFilterBit(kPionRejection_40);
  processor->AddTrackCut(protRej30pionRej40);
  
  AliReducedTrackCut* protRej35pionRej30 = new AliReducedTrackCut("protRej35pionRej30","");
  protRej35pionRej30->AddCut(AliReducedVarManager::kPt, 1.0, 100.);
  protRej35pionRej30->SetRejectKinks();
  protRej35pionRej30->SetTrackFilterBit(kSPDany);
  protRej35pionRej30->SetTrackFilterBit(kEleInclusion_30);
  protRej35pionRej30->SetTrackFilterBit(kProtRejection_35);
  protRej35pionRej30->SetTrackFilterBit(kPionRejection_30);
  processor->AddTrackCut(protRej35pionRej30);
  
  AliReducedTrackCut* protRej35pionRej35 = new AliReducedTrackCut("protRej35pionRej35","");
  protRej35pionRej35->AddCut(AliReducedVarManager::kPt, 1.0, 100.);
  protRej35pionRej35->SetRejectKinks();
  protRej35pionRej35->SetTrackFilterBit(kSPDany);
  protRej35pionRej35->SetTrackFilterBit(kEleInclusion_30);
  protRej35pionRej35->SetTrackFilterBit(kProtRejection_35);
  protRej35pionRej35->SetTrackFilterBit(kPionRejection_35);
  processor->AddTrackCut(protRej35pionRej35);
  
  AliReducedTrackCut* protRej35pionRej40 = new AliReducedTrackCut("protRej35pionRej40","");
  protRej35pionRej40->AddCut(AliReducedVarManager::kPt, 1.0, 100.);
  protRej35pionRej40->SetRejectKinks();
  protRej35pionRej40->SetTrackFilterBit(kSPDany);
  protRej35pionRej40->SetTrackFilterBit(kEleInclusion_30);
  protRej35pionRej40->SetTrackFilterBit(kProtRejection_35);
  protRej35pionRej40->SetTrackFilterBit(kPionRejection_40);
  processor->AddTrackCut(protRej35pionRej40);
  
  AliReducedTrackCut* protRej40pionRej30 = new AliReducedTrackCut("protRej40pionRej30","");
  protRej40pionRej30->AddCut(AliReducedVarManager::kPt, 1.0, 100.);
  protRej40pionRej30->SetRejectKinks();
  protRej40pionRej30->SetTrackFilterBit(kSPDany);
  protRej40pionRej30->SetTrackFilterBit(kEleInclusion_30);
  protRej40pionRej30->SetTrackFilterBit(kProtRejection_40);
  protRej40pionRej30->SetTrackFilterBit(kPionRejection_30);
  processor->AddTrackCut(protRej40pionRej30);
  
  AliReducedTrackCut* protRej40pionRej35 = new AliReducedTrackCut("protRej40pionRej35","");
  protRej40pionRej35->AddCut(AliReducedVarManager::kPt, 1.0, 100.);
  protRej40pionRej35->SetRejectKinks();
  protRej40pionRej35->SetTrackFilterBit(kSPDany);
  protRej40pionRej35->SetTrackFilterBit(kEleInclusion_30);
  protRej40pionRej35->SetTrackFilterBit(kProtRejection_40);
  protRej40pionRej35->SetTrackFilterBit(kPionRejection_35);
  processor->AddTrackCut(protRej40pionRej35);
  
  AliReducedTrackCut* protRej40pionRej40 = new AliReducedTrackCut("protRej40pionRej40","");
  protRej40pionRej40->AddCut(AliReducedVarManager::kPt, 1.0, 100.);
  protRej40pionRej40->SetRejectKinks();
  protRej40pionRej40->SetTrackFilterBit(kSPDany);
  protRej40pionRej40->SetTrackFilterBit(kEleInclusion_30);
  protRej40pionRej40->SetTrackFilterBit(kProtRejection_40);
  protRej40pionRej40->SetTrackFilterBit(kPionRejection_40);
  processor->AddTrackCut(protRej40pionRej40);
  
  // set track prefilter cuts
  AliReducedTrackCut* prefTrackCut1 = new AliReducedTrackCut("prefCutPt09","prefilter P selection");
  prefTrackCut1->AddCut(AliReducedVarManager::kPt, 0.8,100.0);
  prefTrackCut1->SetRequestTPCrefit();
  processor->AddPrefilterTrackCut(prefTrackCut1);  
  
  // MC truth info for reconstructed electron candidates
  AliReducedTrackCut* trueElectron = new AliReducedTrackCut("TrueElectron", "reconstructed electrons with MC truth");
  trueElectron->SetMCFilterBit(kJpsiDecayElectron);
  processor->AddLegCandidateMCcut(trueElectron);  
  
  // MC truth selections(s) for the J/psi mother and electron daughters
  AliReducedTrackCut* mcTruthJpsi = new AliReducedTrackCut("mcTruthJpsi", "Pure MC truth J/psi");
  mcTruthJpsi->SetMCFilterBit(kJpsiInclusive);
  mcTruthJpsi->AddCut(AliReducedVarManager::kRap, -0.9, 0.9);
  AliReducedTrackCut* mcTruthJpsiElectron = new AliReducedTrackCut("mcTruthJpsiElectron", "Pure MC truth electron from J/psi");
  mcTruthJpsiElectron->SetMCFilterBit(kJpsiDecayElectron);
  mcTruthJpsiElectron->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  mcTruthJpsiElectron->AddCut(AliReducedVarManager::kPt, 1.0, 100.0);
  processor->AddJpsiMotherMCCut(mcTruthJpsi, mcTruthJpsiElectron);  
  
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


//_________________________________________________________________
void SetupMixingHandler(AliReducedAnalysisJpsi2ee* task) {
   //
   // setup the mixing handler
   //
   AliMixingHandler* handler = task->GetMixingHandler();
   handler->SetPoolDepth(50);
   handler->SetMixingThreshold(1.0);
   handler->SetDownscaleEvents(1);
   handler->SetDownscaleTracks(1);
   
   handler->AddMixingVariable(AliReducedVarManager::kCentVZERO, gkNCentBins, gCentLims);
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
  for(Int_t i=0; i<gkNMassBins;++i) gMassBins[i] = 0.0+i*0.04;
  
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
      man->AddHistogram(classStr.Data(),"NTracksTPCout_VZEROmult","",kFALSE,900,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      
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
       man->AddHistogram(classStr.Data(), "ITSlayerHit_Phi", "Hits in the ITS layers vs #varphi", kFALSE,
                         180, 0.0, 6.29, AliReducedVarManager::kPhi, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
       man->AddHistogram(classStr.Data(), "ITSlayerHit_Eta", "Hits in the ITS layers vs #eta", kFALSE,
                         100, -1.0, 1.0, AliReducedVarManager::kEta, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
       continue;
    }  // end of ITSclusterMap histogram definitions
    
    if(classStr.Contains("TrackTPCclusterMap_")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       man->AddHistogram(classStr.Data(), "TPCclusterMap", "TPC cluster map", kFALSE,
                         8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
       man->AddHistogram(classStr.Data(), "TPCclusterMap_Phi", "TPC cluster map vs #varphi", kFALSE,
                         180, 0.0, 6.29, AliReducedVarManager::kPhi, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
       man->AddHistogram(classStr.Data(), "TPCclusterMap_Eta", "TPC cluster map vs #eta", kFALSE,
                         100, -1.0, 1.0, AliReducedVarManager::kEta, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
       man->AddHistogram(classStr.Data(), "TPCclusterMap_Pt", "TPC cluster map vs p_{T}", kFALSE,
                         100, 0.0, 10.0, AliReducedVarManager::kPt, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
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
      man->AddHistogram(classStr.Data(), "P_Pin", "P vs Pin", kFALSE, 200, 0.0, 10.0, AliReducedVarManager::kP, 200, 0.0, 10.0, AliReducedVarManager::kPin);
      man->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution", kFALSE, 1000, 0.0, 50.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Eta", "", kFALSE, 1000, -1.5, 1.5, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Eta_Phi", "", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEta, 100, 0., 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "Phi", "", kFALSE, 1000, 0.0, 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "DCAxy_Pt", "DCAxy", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 50, 0.0, 5.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE, 200, -5.0, 5.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE, 200, -5.0, 5.0, AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(), "DCAz_Pt", "DCAz", kFALSE, 100, -3.0, 3.0, AliReducedVarManager::kDcaZ, 50, 0.0, 5.0, AliReducedVarManager::kPt);
      
        man->AddHistogram(classStr.Data(),"ITSncls", "ITS nclusters", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"ITSnclsShared", "ITS nclusters shared", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSnclsShared);
        man->AddHistogram(classStr.Data(),"Eta_Phi_ITSncls_prof","ITS <nclusters> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"Eta_Phi_ITSnclsShared_prof","ITS <nclusters-shared> vs (#eta,#phi)",kTRUE,
                          36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kITSnclsShared);
	man->AddHistogram(classStr.Data(),"ITSchi2", "ITS #chi^{2}", kFALSE, 200,0.0,50.0, AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"ITSchi2_RunID_prof","",kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 10, 0., 10., AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"ITSchi2_ITSncls", "ITS #chi^{2} vs ITS clusters", kFALSE, 100,0.0,100.0, AliReducedVarManager::kITSchi2,7,-0.5,6.5, AliReducedVarManager::kITSncls);
	man->AddHistogram(classStr.Data(),"Eta_Phi_ITSchi2_prof","ITS <#chi^{2}> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 100, 0.0, 1000, AliReducedVarManager::kITSchi2);        
        man->AddHistogram(classStr.Data(),"Chi2TPCConstrainedVsGlobal","",kFALSE,200,0.,200.,AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
        
        man->AddHistogram(classStr.Data(),"TPCncls","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCncls_RunID_prof","",kTRUE, gkNRuns, -0.5, -0.5+gkNRuns, AliReducedVarManager::kRunID, 10, 0., 10., AliReducedVarManager::kTPCncls);
	man->AddHistogram(classStr.Data(),"TPCcrossedRows","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCcrossedRows);
        man->AddHistogram(classStr.Data(),"TPCfindableClusters","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsF);
        man->AddHistogram(classStr.Data(),"TPCcrossedRowsOverFindableClusters","", kFALSE, 200,0.0,1.5,AliReducedVarManager::kTPCcrossedRowsOverFindableClusters);
	man->AddHistogram(classStr.Data(),"TPCnclsShared","",kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsShared);
        man->AddHistogram(classStr.Data(),"TPCnclsSharedRatio","",kFALSE, 200,0.,1.,AliReducedVarManager::kTPCnclsSharedRatio);        
        
        man->AddHistogram(classStr.Data(),"TPCchi2","",kFALSE, 200,0.0,10.0,AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_Pt","",kFALSE, 200,0.0,10.0,AliReducedVarManager::kTPCchi2, 100, 0., 10.0, AliReducedVarManager::kPt);
        man->AddHistogram(classStr.Data(),"TPCchi2_TPCncls", "TPC #chi^{2} vs TPC clusters", kFALSE, 100,0.0,10.0, AliReducedVarManager::kTPCchi2,160,0.,160., AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCsegments","",kFALSE, 9,0.0,9.0,AliReducedVarManager::kTPCNclusBitsFired);

        man->AddHistogram(classStr.Data(),"Eta_Phi_TPCncls_prof","TPC <nclusters> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCncls);    
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCcrossedRows_prof","TPC <n crossed rows> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCcrossedRows);
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCnclsF_prof","TPC <nclusters findable> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCnclsF);
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCchi2_prof","TPC <#chi^{2}> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"Eta_Phi_TPCchi2","TPC #chi^{2} vs (#eta,#phi)",kFALSE,
                          36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 100, 0., 10.0, AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCsignal_Pin","TPC dE/dx vs. inner param P",kFALSE,
                     300,0.0,30.0,AliReducedVarManager::kPin,200,-0.5,199.5,AliReducedVarManager::kTPCsignal);
	man->AddHistogram(classStr.Data(),"TPCsignalN","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCsignalN);
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCsignalN_prof","TPC <nclusters pid> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCsignalN);    
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Pin","TPC N_{#sigma} electron vs. inner param P",kFALSE,
                     100,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Eta","TPC N_{#sigma} electron vs. #eta",kFALSE,
                          36,-0.9,0.9,AliReducedVarManager::kEta,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TOFbeta_P","TOF #beta vs P",kFALSE,
                     200,0.0,20.0,AliReducedVarManager::kP, 220,0.0,1.1,AliReducedVarManager::kTOFbeta);
        
        man->AddHistogram(classStr.Data(),"TPCchi2_ITSchi2","",kFALSE, 100,0.0,10.0,AliReducedVarManager::kTPCchi2, 100,0.0,40.0,AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"ITSchi2_Chi2TPCconstrainedVsGlobal","",kFALSE, 500,0.0,1000.0,AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 100,0.0,40.0,AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_Chi2TPCconstrainedVsGlobal","",kFALSE, 500,0.0,1000.0,AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 100,0.0,10.0,AliReducedVarManager::kTPCchi2);
        
        if(classStr.Contains("MCTruth")) {
          man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE, 150, 0., 15.0, AliReducedVarManager::kPtMC);
          man->AddHistogram(classStr.Data(), "PtRec_PtMC", "p_{T} MC vs p_{T} reconstructed", kFALSE, 150, 0., 15.0, AliReducedVarManager::kPtMC, 150, 0., 15.0, AliReducedVarManager::kPt);
          man->AddHistogram(classStr.Data(), "PhiMC", "#varphi MC", kFALSE, 180, 0., 6.3, AliReducedVarManager::kPhiMC);
          man->AddHistogram(classStr.Data(), "PhiRec_PhiMC", "#varphi MC vs #varphi reconstructed", kFALSE, 180, 0., 6.3, AliReducedVarManager::kPhiMC, 180, 0., 6.3, AliReducedVarManager::kPhi);
          man->AddHistogram(classStr.Data(), "EtaMC", "#eta MC", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEtaMC);
          man->AddHistogram(classStr.Data(), "EtaRec_EtaMC", "#eta MC vs #eta reconstructed", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEtaMC, 100, -1.0, 1.0, AliReducedVarManager::kEta);          
          man->AddHistogram(classStr.Data(), "PDGcode0", "PDG code of the track", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC);
          man->AddHistogram(classStr.Data(), "PDGcode1", "PDG code of the track's mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+1);
          man->AddHistogram(classStr.Data(), "PDGcode2", "PDG code of the track's grand-mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+2);
          man->AddHistogram(classStr.Data(), "PDGcode3", "PDG code of the track's grand-grand mother", kFALSE, 12000, -6000., 6000.0, AliReducedVarManager::kPdgMC+3);
        }
      continue;
    }  // end if "TrackQA"
        
     Int_t vars[3] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, AliReducedVarManager::kCentVZERO};
       
    TArrayD pairHistBinLimits[3];
    pairHistBinLimits[0] = TArrayD(gkNMassBins, gMassBins);
    pairHistBinLimits[1] = TArrayD(gkNPtBins, gPtLims);
    pairHistBinLimits[2] = TArrayD(gkNCentBins, gCentLims);
    
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
         man->AddHistogram(classStr.Data(), "PtMC", "p_{T} MC", kFALSE, 1000, 0., 10.0, AliReducedVarManager::kPtMC);
         man->AddHistogram(classStr.Data(), "PtMC_coarse", "p_{T} MC", kFALSE, 20, 0., 20.0, AliReducedVarManager::kPtMC);
         man->AddHistogram(classStr.Data(), "MassMC_Mass", "Invariant mass, MC vs reconstructed", kFALSE, 150, 2.0, 3.5, AliReducedVarManager::kMass,
            150, 2.0, 3.5, AliReducedVarManager::kMassMC);
         man->AddHistogram(classStr.Data(), "PtMC_Pt", "pair pT, MC vs reconstructed", kFALSE, 150, 0.0, 15., AliReducedVarManager::kPt,
                           150, 0.0, 15.0, AliReducedVarManager::kPtMC);
         man->AddHistogram(classStr.Data(), "PtMC_MassMC", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPtMC,
                           30, 2.7, 3.3, AliReducedVarManager::kMassMC);
         man->AddHistogram(classStr.Data(), "Pt_Mass", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPt,
                           30, 2.7, 3.3, AliReducedVarManager::kMass);
         man->AddHistogram(classStr.Data(), "PtMC_Mass", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPtMC,
                           30, 2.7, 3.3, AliReducedVarManager::kMass);
         man->AddHistogram(classStr.Data(), "Pt_MassMC", "", kFALSE, 100, 0.0, 1., AliReducedVarManager::kPt,
                           30, 2.7, 3.3, AliReducedVarManager::kMassMC);
         man->AddHistogram(classStr.Data(), "PtMC_Pt_lowPtZoom", "pair pT, MC vs reconstructed", kFALSE, 100, 0.0, 0.5, AliReducedVarManager::kPt,
                           100, 0.0, 0.5, AliReducedVarManager::kPtMC);
         man->AddHistogram(classStr.Data(), "PtMC_Pt_Mass", "pair pT, MC vs reconstructed, vs mass", kFALSE, 100, 0.0, 0.5, AliReducedVarManager::kPt,
                           100, 0.0, 0.5, AliReducedVarManager::kPtMC, 15, 2.72, 3.32, AliReducedVarManager::kMass);
         man->AddHistogram(classStr.Data(), "PtMC_Pt_Mass_coarse", "pair pT, MC vs reconstructed, vs mass", kFALSE, 25, 0.0, 0.5, AliReducedVarManager::kPt,
                           25, 0.0, 0.5, AliReducedVarManager::kPtMC, 15, 2.72, 3.32, AliReducedVarManager::kMass);
         
         man->AddHistogram(classStr.Data(), "Mass_Pt_CentVZERO", "", kFALSE, gkNMassBins-1, gMassBins, AliReducedVarManager::kMass, gkNPtBins-1, gPtLims, AliReducedVarManager::kPt, gkNCentBins-1, gCentLims, AliReducedVarManager::kCentVZERO);
      }
      
      continue;
    }   // end if for Pair classes of histograms
  }  // end loop over histogram classes
}
