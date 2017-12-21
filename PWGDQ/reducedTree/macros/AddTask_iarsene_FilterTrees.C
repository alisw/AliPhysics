//void Setup(AliReducedAnalysisJpsi2ee* processor, TString prod="LHC10h");
//AliHistogramManager* SetupHistogramManager(TString prod="LHC10h");
//void DefineHistograms(AliHistogramManager* man, TString prod="LHC10h");

//__________________________________________________________________________________________
AliAnalysisTask* AddTask_iarsene_FilterTrees(Bool_t isAliRoot=kTRUE, Int_t runMode=1, TString prod="LHC10h"){    
   //
   // isAliRoot=kTRUE for ESD/AOD analysis in AliROOT, kFALSE for root analysis on reduced trees
   // runMode=1 (AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents)
   //               =2 (AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree)
   //
   //get the current analysis manager

  printf("INFO on AddTask_iarsene_FilterTrees(): (isAliRoot, runMode) :: (%d,%d) \n", isAliRoot, runMode);

  AliReducedAnalysisFilterTrees* filterTask = new AliReducedAnalysisFilterTrees("FilterTrees","filter DST trees");
  filterTask->Init();
  filterTask->SetFilteredTreeWritingOption(AliReducedAnalysisTaskSE::kFullEventsWithFullTracks);
  filterTask->SetWriteFilteredTracks(kFALSE);
  filterTask->SetWriteFilteredPairs(kFALSE);
  filterTask->SetBuildCandidatePairs(AliReducedPairInfo::kJpsiToEE);
  //filterTask->SetBuildCandidatePairs(AliReducedPairInfo::kADzeroToKplusPiminus);
  filterTask->SetBuildCandidateLikePairs(kTRUE);
  filterTask->SetRunCandidatePrefilter(kTRUE);
  filterTask->SetRunCandidatePrefilterOnSameCharge(kFALSE);
  Setup(filterTask, prod);
  // initialize an AliAnalysisTask which will wrapp the AliReducedAnalysisFilterTrees such that it can be run in an aliroot analysis train (e.g. LEGO, local analysis etc)
  AliAnalysisTaskReducedEventProcessor* task = new AliAnalysisTaskReducedEventProcessor("ReducedEventAnalysisManager", runMode, kTRUE);
  task->AddTask(filterTask);
  
  if(isAliRoot){
     AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
     if (!mgr) {
        Error("AddTask_iarsene_dst", "No analysis manager found.");
        return 0;
     }
     
     AliAnalysisDataContainer* cReducedEvent = NULL;
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) {
       printf("INFO on AddTask_iarsene_FilterTrees(): use on the fly events ");
       cReducedEvent = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ReducedEventDQ");
       if(!cReducedEvent) {
         printf("ERROR: In AddTask_iarsene_FilterTrees(), couldn't find exchange container with ReducedEvent! ");
         printf("             You need to use AddTask_iarsene_dst() in the on-the-fly reduced events mode!");
         return 0x0;
       }
     }
            
     mgr->AddTask(task);
      
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree) 
        mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
      
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) 
       mgr->ConnectInput(task, 0, cReducedEvent);
  
     AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("filterQAhistos", THashList::Class(),
                                                                  AliAnalysisManager::kOutputContainer, "dstAnalysisHistograms.root");
     mgr->ConnectOutput(task, 1, cOutputHist );
     
     AliAnalysisDataContainer *cOutputTree = mgr->CreateContainer("filteredDstTree", TTree::Class(),
                                                                  AliAnalysisManager::kOutputContainer, "dstTree.root");
     mgr->ConnectOutput(task, 2, cOutputTree );
  }
  else {
    // nothing at the moment   
  }
  
  return task;
}


//_________________________________________________________________
void Setup(AliReducedAnalysisFilterTrees* processor, TString prod /*="LHC10h"*/) {
  //
  // Configure the analysis task
  // Setup histograms, handlers, cuts, etc.
  // 
  
  // Set event cuts
  AliReducedEventCut* evCut1 = new AliReducedEventCut("Centrality","Centrality selection");
  if(prod.Contains("LHC15o")) evCut1->AddCut(AliReducedVarManager::kCentVZERO, 50., 90.);
  evCut1->AddCut(AliReducedVarManager::kVtxZ, -5.0, 6.5);
  evCut1->AddCut(AliReducedVarManager::kIsPhysicsSelection, 0.1, 2.);   // request physics selection
  //evCut1->AddCut(AliReducedVarManager::kDeltaVtxZ, -0.2, 0.2);   // select based on the difference between the vtxZ and the tPC vtxZ
  //evCut1->AddCut(AliReducedVarManager::kIRIntClosestIntMap+1, 0.99, 5000., kTRUE);   // exclude out-of-bunch pileup
  //evCut1->EnableVertexDistanceCut();
  //evCut1->AddCut(AliReducedVarManager::kTZEROpileup, -0.1, 0.1);
  TF1* cutCorrTPCoutVZEROmult = new TF1("cutCorrTPCoutVZEROmult", "[0]+[1]*x+[2]*x*x", 0., 1.e+5);
  cutCorrTPCoutVZEROmult->SetParameters(-2200., 2.5, 1.2e-5);
  if(prod.Contains("LHC15o") && !processor->GetRunOverMC()) 
      evCut1->AddCut(AliReducedVarManager::kVZEROTotalMult, cutCorrTPCoutVZEROmult, 99999., kFALSE, AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 0., 99998.);
  
  processor->AddEventCut(evCut1);
  
  // Set track cuts
  AliReducedTrackCut* standardCut = new AliReducedTrackCut("standard","");
  standardCut->AddCut(AliReducedVarManager::kP, 1.2,30.0);
  standardCut->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  standardCut->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  standardCut->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  standardCut->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.0, 30000.0);
  standardCut->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.0, 30000.0);
  standardCut->SetRejectKinks();
  standardCut->SetRequestITSrefit();
  standardCut->SetRequestTPCrefit();
  //standardCut->SetRequestTOFout();
  standardCut->SetRequestSPDany();
  TF1* chi2Cut=new TF1("chi2Cut","[0]+[1]*x",0.,15000.);
  chi2Cut->SetParameters(1.9, 1.1e-4);
  //standardCut->AddCut(AliReducedVarManager::kTPCchi2, 0., chi2Cut, kFALSE, AliReducedVarManager::kNTracksPerTrackingStatus        +AliReducedVarManager::kTPCout, 0., 99998.);
  standardCut->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  standardCut->AddCut(AliReducedVarManager::kITSchi2, 0., 10.);
  standardCut->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 6., 9.);
  standardCut->AddCut(AliReducedVarManager::kTPCnclsSharedRatio, 0.3, 2., kTRUE);
  //standardCut->AddCut(AliReducedVarManager::kTPCcrossedRowsOverFindableClusters, 0.8, 2.);
  //standardCut->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 36., 1.0e+10, kTRUE);
  processor->AddTrackCut(standardCut); 
  
  AliReducedTrackCut* cut1 = new AliReducedTrackCut("cut1","");
  cut1->AddCut(AliReducedVarManager::kPt, 0.5,30.0);
  cut1->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  cut1->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  cut1->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  cut1->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  cut1->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3.0, 3.0);
  //cut1->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 1.0, 30000.0);
  //cut1->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 1.0, 30000.0);
  //cut1->SetRejectKinks();
  cut1->SetRequestITSrefit();
  cut1->SetRequestTPCrefit();
  //cut1->SetRequestTOFout();
  //cut1->SetRequestSPDany();
  cut1->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  cut1->AddCut(AliReducedVarManager::kITSchi2, 0., 10.);
  cut1->AddCut(AliReducedVarManager::kTPCNclusBitsFired, 6., 9.);
  cut1->AddCut(AliReducedVarManager::kTPCnclsSharedRatio, 0.3, 2., kTRUE);
  //cut1->AddCut(AliReducedVarManager::kTPCcrossedRowsOverFindableClusters, 0.8, 2.);
  //cut1->AddCut(AliReducedVarManager::kChi2TPCConstrainedVsGlobal, 36., 1.0e+10, kTRUE);
  processor->AddTrackCut(cut1); 
  
  
  // TODO: allow multiple track cuts
  
  // set track prefilter cuts  TODO: add prefilter for creating jpsi candidates
/*  AliReducedTrackCut* prefTrackCut1 = new AliReducedTrackCut("prefCutPt09","prefilter P selection");
  prefTrackCut1->AddCut(AliReducedVarManager::kP, 0.7,100.0);
  prefTrackCut1->SetRequestTPCrefit();
  processor->AddPrefilterTrackCut(prefTrackCut1);  
  
  // set pair prefilter cuts
  AliReducedVarCut* prefPairCut = new AliReducedVarCut("prefCutM50MeV","prefilter pair cuts");
  prefPairCut->AddCut(AliReducedVarManager::kMass, 0.0, 0.05, kTRUE);
  processor->AddPrefilterPairCut(prefPairCut);
*/
  
  // Set pair cuts
  AliReducedTrackCut* pairCut1 = new AliReducedTrackCut("Ptpair1","Pt pair selection");
  pairCut1->AddCut(AliReducedVarManager::kPt, 0.0,1.0);
  pairCut1->AddCut(AliReducedVarManager::kRap, -0.9,0.9);
  processor->AddPairCut(pairCut1);
  AliReducedTrackCut* pairCut2 = new AliReducedTrackCut("Ptpair2","Pt pair selection");
  pairCut2->AddCut(AliReducedVarManager::kPt, 1.0,2.0);
  pairCut2->AddCut(AliReducedVarManager::kRap, -0.9,0.9);
  processor->AddPairCut(pairCut2);
  AliReducedTrackCut* pairCut3 = new AliReducedTrackCut("Ptpair3","Pt pair selection");
  pairCut3->AddCut(AliReducedVarManager::kPt, 2.0,3.0);
  pairCut3->AddCut(AliReducedVarManager::kRap, -0.9,0.9);
  processor->AddPairCut(pairCut3);
  
  // set candidate leg cuts
  // Jpsi candidate legs ---------------------------------------------------------------------
  processor->AddCandidateLeg1Cut(cut1);
  processor->AddCandidateLeg1Cut(standardCut);
  
  AliReducedTrackCut* prefTrackCut1 = new AliReducedTrackCut("prefCutPt07","prefilter P selection");
  prefTrackCut1->AddCut(AliReducedVarManager::kP, 0.7,100.0);
  prefTrackCut1->SetRequestTPCrefit();
  processor->AddCandidateLeg1PrefilterCut(prefTrackCut1);  
  
  AliReducedVarCut* prefPairCut = new AliReducedVarCut("prefCutM50MeV","prefilter pair cuts");
  prefPairCut->AddCut(AliReducedVarManager::kMass, 0.0, 0.05, kTRUE);
  processor->AddCandidateLeg1PairPrefilterCut(prefPairCut);
  
  // D+ -> K+ pi-
  AliReducedTrackCut* kplus = new AliReducedTrackCut("kplus","");
  kplus->AddCut(AliReducedVarManager::kCharge, 0.9, 1.1);
  kplus->AddCut(AliReducedVarManager::kPt, 0.4,10.0);
  kplus->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  kplus->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  kplus->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  kplus->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  kplus->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kKaon, -3.0, 3.0);
  //kplus->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, -3.0, 3.0, kTRUE);
  //kplus->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, -3.0, 3.0, kTRUE);
  kplus->SetRequestITSrefit();
  kplus->SetRequestTPCrefit();
  kplus->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  //processor->AddCandidateLeg1Cut(kplus);
  
  AliReducedTrackCut* piminus = new AliReducedTrackCut("piminus","");
  piminus->AddCut(AliReducedVarManager::kCharge, -1.1, -0.9);
  piminus->AddCut(AliReducedVarManager::kPt, 0.4,10.0);
  piminus->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
  piminus->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  piminus->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  piminus->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  piminus->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, -3.0, 3.0);
  //piminus->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, -3.0, 3.0, kTRUE);
  //piminus->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kKaon, -3.0, 3.0, kTRUE);
  piminus->SetRequestITSrefit();
  piminus->SetRequestTPCrefit();
  piminus->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.0);
  //processor->AddCandidateLeg2Cut(piminus); 
  
  AliReducedTrackCut* prefTrackCutKplus = new AliReducedTrackCut("prefCutKaon","prefilter P selection");
  prefTrackCutKplus->AddCut(AliReducedVarManager::kP, 0.5,100.0);
  prefTrackCutKplus->SetRequestTPCrefit();
  //processor->AddCandidateLeg1PrefilterCut(prefTrackCutKplus);  
  
  AliReducedTrackCut* prefTrackCutPiminus = new AliReducedTrackCut("prefCutPion","prefilter P selection");
  prefTrackCutPiminus->AddCut(AliReducedVarManager::kP, 0.5,100.0);
  prefTrackCutPiminus->SetRequestTPCrefit();
  //processor->AddCandidateLeg2PrefilterCut(prefTrackCutPiminus);  
  
  AliReducedVarCut* prefPairCutKaon = new AliReducedVarCut("prefCutM50MeVKaon","prefilter pair cuts");
  prefPairCutKaon->AddCut(AliReducedVarManager::kMass, 0.0, 0.1, kTRUE);
  //processor->AddCandidateLeg1PairPrefilterCut(prefPairCutKaon);
  AliReducedVarCut* prefPairCutPion = new AliReducedVarCut("prefCutM50MeVPion","prefilter pair cuts");
  prefPairCutPion->AddCut(AliReducedVarManager::kMass, 0.0, 0.1, kTRUE);
  //processor->AddCandidateLeg2PairPrefilterCut(prefPairCutPion);
    
  
  AliReducedTrackCut* pairCut1 = new AliReducedTrackCut("Ptpair","Pt pair selection");
  pairCut1->AddCut(AliReducedVarManager::kPt, 0.0,100.0);
  pairCut1->AddCut(AliReducedVarManager::kRap, -0.9,0.9);
  processor->AddCandidatePairCut(pairCut1);
  
  
  //SetupHistogramManager(processor->GetHistogramManager(), prod);
  SetupHistogramManager(processor, prod);
  //SetupMixingHandler(processor);
}


//_________________________________________________________________
void SetupHistogramManager(AliReducedAnalysisFilterTrees* task, TString prod /*="LHC10h"*/) {
  //
  // setup the histograms manager
  //
  AliReducedVarManager::SetDefaultVarNames();
  
  DefineHistograms(task, prod);
  
  AliReducedVarManager::SetUseVars(task->GetHistogramManager()->GetUsedVars());
}


//_________________________________________________________________
void DefineHistograms(AliReducedAnalysisFilterTrees* task, TString prod /*="LHC10h"*/) {
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
  
  //histClasses += "EventTag_BeforeCuts;";   
  //histClasses += "EventTag_AfterCuts;";   
  histClasses += "EventTriggers_BeforeCuts;";
  histClasses += "EventTriggers_AfterCuts;";   
  
  if(task->GetWriteFilteredTracks()) {
  histClasses += "Track_BeforeCuts;";
     for(Int_t i=0; i<task->GetNTrackCuts(); ++i) {
         TString cutName = task->GetTrackCutName(i);
         histClasses += Form("Track_%s;", cutName.Data());
     }
  }
  
  if(task->GetWriteFilteredPairs()) {
     histClasses += "Pair_BeforeCuts;";
     histClasses += "PairQualityFlags_BeforeCuts;";
     for(Int_t i=0; i<task->GetNPairCuts(); ++i) {
        TString cutName = task->GetPairCutName(i);
        histClasses += Form("Pair_%s;", cutName.Data());
        histClasses += Form("PairQualityFlags_%s;", cutName.Data());
     }
  }
  
  if(task->GetBuildCandidatePairs()) {
     for(Int_t i=0; i<task->GetNCandidateLegCuts(); ++i) {
        histClasses += Form("Track_LEG1_BeforePrefilter_%s;", task->GetCandidateLegCutName(i,1));
        histClasses += Form("Track_LEG2_BeforePrefilter_%s;", task->GetCandidateLegCutName(i,2));
        if(task->GetRunCandidatePrefilter()) {
           histClasses += Form("Track_LEG1_AfterPrefilter_%s;", task->GetCandidateLegCutName(i,1));
           histClasses += Form("Track_LEG2_AfterPrefilter_%s;", task->GetCandidateLegCutName(i,2));
        }
        histClasses += Form("Pair_Candidate12_%s%s;", 
                         task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""));
        if(task->GetBuildCandidateLikePairs()) {
           histClasses += Form("Pair_Candidate11_%s%s;", 
                               task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""));
           histClasses += Form("Pair_Candidate22_%s%s;", 
                               task->GetCandidateLegCutName(i,1), (task->IsAsymmetricDecayChannel() ? Form("_%s", task->GetCandidateLegCutName(i,2)) : ""));
        }
     }
     if(task->GetRunCandidatePrefilter()) {
        histClasses += "Track_LEG1_PrefilterTrack;";
        histClasses += "Track_LEG2_PrefilterTrack;";
     }
  }
  
  Int_t runNBins = 0;
  Double_t runHistRange[2] = {0.0,0.0};
  
  // Pb-Pb from 2010 run range is default: LHC10h
  runNBins = 2500;
  runHistRange[0] = 137100.;
  runHistRange[1] = 139600.;
  
  // Pb-Pb of 2011
  if(prod.Contains("LHC11h")) {
    runNBins = 2700;
    runHistRange[0] = 167900.;
    runHistRange[1] = 170600.;
  }
  
  // Pb-Pb of 2015
  if(prod.Contains("LHC15o")) {
     runNBins = 2100;
     runHistRange[0] = 244900.;
     runHistRange[1] = 247000.;
  }
  
  // p-Pb of 2013
  if(prod.Contains("LHC13b") || prod.Contains("LHC13c")) {
    runNBins = 400;
    runHistRange[0] = 195300.;
    runHistRange[1] = 195700.;
  }
  
  // pp at 13 TeV
  if(prod.Contains("LHC16l")) {
     runNBins = 1140;
     runHistRange[0] = 258880.;
     runHistRange[1] = 260020.;
  }
  
  // p-Pb at 8.16 TeV
  if(prod.Contains("LHC16r")) {
     runNBins = 1000;
     runHistRange[0] = 265400.;
     runHistRange[1] = 266400.;
  }
  
  const Int_t kNBCBins = 3600;
  Double_t bcHistRange[2] = {-0.5,3599.5};
  
  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");
  
  cout << "Histogram classes included in the Histogram Manager" << endl;
  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    
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
      man->AddHistogram(classStr.Data(),"RunNo","Run numbers",kFALSE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo);
      man->AddHistogram(classStr.Data(),"TimeFromSOR","Events vs time",kFALSE, 450, 0., 450., AliReducedVarManager::kTimeRelativeSOR);
      man->AddHistogram(classStr.Data(),"VtxX","Vtx X",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxX_TimeFromSOR_prof","<Vtx X> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxY","Vtx Y",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxY_TimeFromSOR_prof","<Vtx Y> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxZ","Vtx Z",kFALSE,300,-15.,15.,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"VtxZ_TimeFromSOR_prof","<Vtx Z> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-15.0,15.0,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentVZERO);
      man->AddHistogram(classStr.Data(),"CentVZERO_VtxZ","Centrality(VZERO) vs vtxZ",kFALSE, 50, 0.0, 100.0, AliReducedVarManager::kCentVZERO,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentVZERO_TimeFromSOR","Centrality(VZERO) vs time from SOR",kFALSE,
                        90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 50, 0.0, 100.0, AliReducedVarManager::kCentVZERO);
      man->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentSPD);
      man->AddHistogram(classStr.Data(),"CentSPD_VtxZ","Centrality(SPD) vs vtxZ",kFALSE, 50, 0.0, 100.0, AliReducedVarManager::kCentSPD,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentTPC);
      man->AddHistogram(classStr.Data(),"CentTPC_VtxZ","Centrality(TPC) vs vtxZ",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentTPC,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentQuality","Centrality quality",kFALSE, 500, 0., 500., AliReducedVarManager::kCentQuality);
      man->AddHistogram(classStr.Data(),"IRIntClosestInt1Map","Closest out of bunch interactions Int1",kFALSE,7000,-3500.,3500.,AliReducedVarManager::kIRIntClosestIntMap);
      man->AddHistogram(classStr.Data(),"IRIntClosestInt2Map","Closest out of bunch interactions Int2",kFALSE,7000,-3500.,3500.,AliReducedVarManager::kIRIntClosestIntMap+1);
      man->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event",kFALSE,500,0.,30000.,AliReducedVarManager::kNtracksTotal);
      man->AddHistogram(classStr.Data(),"NTracksTotal_BeamIntensity0_prof","Number of total tracks per event",kTRUE,100, 3.0e+12, 9.0e+12, AliReducedVarManager::kBeamIntensity0, 500,0.,20000.,AliReducedVarManager::kNtracksTotal);
      man->AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event",kFALSE,300,0.,300.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"NTracksSelected_TimeFromSOR","Averaged number of selected tracks per event vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"NTracksSelected_CentVZERO_TimeFromSOR","Averaged number of selected tracks per event per centrality vs time from SOR",kTRUE, 20, 0.0, 100.0, AliReducedVarManager::kCentVZERO, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"DeltaVtxZ","Z_{global}-Z_{TPC}",kFALSE,300,-1.,1.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"DeltaVtxZ_TimeFromSOR_prof","Z_{global}-Z_{TPC} vs time from SOR",kTRUE,90, 0.0, 450.,  AliReducedVarManager::kTimeRelativeSOR, 300,-1.,1.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"NVtxTPCContributors","",kFALSE,500,0.,10000.,AliReducedVarManager::kNVtxTPCContributors);
      man->AddHistogram(classStr.Data(),"NVtxTPCContributors_BeamIntensity0","",kTRUE, 100, 3.0e+12, 9.0e+12, AliReducedVarManager::kBeamIntensity0, 500,0.,10000.,AliReducedVarManager::kNVtxTPCContributors);
      man->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0", kFALSE, 200, 0., 5000., AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(),"SPDntracklets_TimeFromSOR_prof", "SPD <#tracklets> in |#eta|<1.0 vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, 0., 5000., AliReducedVarManager::kSPDntracklets);
      for(Int_t il=0; il<2; ++il) {
        man->AddHistogram(classStr.Data(), Form("SPDfiredChips_layer%d",il+1), Form("SPD fired chips in layer %d",il+1), 
			  kFALSE, 200, 0., 1000., AliReducedVarManager::kSPDFiredChips+il);
        man->AddHistogram(classStr.Data(), Form("SPDfiredChips_layer%d_TimeFromSOR_prof",il+1), Form("SPD <#fired chips> in layer %d vs time from SOR",il+1), 
                          kTRUE, 90, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR, 200, 0., 1000., AliReducedVarManager::kSPDFiredChips+il);
      }
      for(Int_t il=0; il<6; ++il) {
        man->AddHistogram(classStr.Data(), Form("ITSclusters_layer%d",il+1), Form("ITS clusters in layer %d",il+1), 
			  kFALSE, 300, 0., 15000., AliReducedVarManager::kITSnClusters+il);
        man->AddHistogram(classStr.Data(), Form("ITSclusters_layer%d_TimeFromSOR_prof",il+1), Form("ITS <clusters> in layer %d vs time from SOR",il+1), 
                          kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300, 0., 15000., AliReducedVarManager::kITSnClusters+il);
      }
      man->AddHistogram(classStr.Data(), "SPDnSingleClusters", "SPD single clusters", 
			kFALSE, 200, 0., 10000., AliReducedVarManager::kSPDnSingleClusters);	
      man->AddHistogram(classStr.Data(), "SPDnSingleClusters_TimeFromSOR_prof", "SPD single <#clusters> vs time from SOR", 
                        kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, 0., 10000., AliReducedVarManager::kSPDnSingleClusters);	
      man->AddHistogram(classStr.Data(), "SPDnSingleClusters_BeamInt0_prof", "SPD single <#clusters> vs beam intensity", 
                        kTRUE, 100, 3.0e+12, 9.0e+12., AliReducedVarManager::kBeamIntensity0, 200, 0., 10000., AliReducedVarManager::kSPDnSingleClusters);	
      man->AddHistogram(classStr.Data(),"VZEROmult", "", kFALSE, 500, 0.0, 50000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"VZEROmult_VtxContributors", "", kFALSE, 
		   200, 0.0, 40000., AliReducedVarManager::kVZEROTotalMult, 200, 0.0, 10000., AliReducedVarManager::kNVtxContributors);
      man->AddHistogram(classStr.Data(),"VZEROmult_SPDntracklets", "", kFALSE, 
                        200, 0.0, 5000., AliReducedVarManager::kSPDntracklets, 200, 0.0, 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsPhysicsSelection, 0,0.0,0.0,AliReducedVarManager::kNothing, 0,0.0,0.0,AliReducedVarManager::kNothing, "off;on");
      man->AddHistogram(classStr.Data(),"IsPhysicsSelection_TimeFromSOR_prof","Fraction of physics selection selected",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 2,-0.5,1.5,AliReducedVarManager::kIsPhysicsSelection);
      man->AddHistogram(classStr.Data(),"NPosTracksAnalyzed","#positrons per event",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksPosAnalyzed);
      man->AddHistogram(classStr.Data(),"NNegTracksAnalyzed","#electrons per event",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksNegAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalTracksAnalyzed","#leg candidates per event",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalPairsAnalyzed","#dielectron candidates per event",kFALSE,100,0.,100.,AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"NPosTracksAnalyzed_TimeFromSOR_prof","<#positrons> per event vs time from SOR",kTRUE,90, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksPosAnalyzed);
      man->AddHistogram(classStr.Data(),"NNegTracksAnalyzed_TimeFromSOR_prof","<#electrons> per event vs time from SOR",kTRUE,90, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR,100,0.,100.,AliReducedVarManager::kNtracksNegAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalTracksAnalyzed_TimeFromSOR_prof","<#leg candidates> per event vs time from SOR",kTRUE,90, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR,100,0.,100.,AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalTracksAnalyzed_BeamInt0_prof","<#leg candidates> per event vs beam intensity",kTRUE,100, 3.0e+12, 9.0e+12, AliReducedVarManager::kBeamIntensity0,100,0.,100.,AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalPairsAnalyzed_TimeFromSOR_prof","<#dielectron> candidates per event vs time from SOR",kTRUE,90, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR,100,0.,100.,AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"NTotalPairsAnalyzed_BeamInt0_prof","<#dielectron> candidates vs beam intensity",kTRUE,100, 3.0e+12, 9.0e+12, AliReducedVarManager::kBeamIntensity0,100,0.,100.,AliReducedVarManager::kNpairsSelected);
      
      /*for(Int_t i=0;i<36;++i)
         man->AddHistogram(classStr.Data(),Form("NTotalTracksAnalyzedInPhiBins_phiSec%d_%s_TimeFromSOR_prof", i%18, (i<18 ? "negEta" : "posEta")), Form("<#leg candidates> per event in #varphi sector %d and %s vs time from SOR", i%18, (i<18 ? "negative #eta" : "positive #eta")),kTRUE,90, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR,100,0.,100.,AliReducedVarManager::kNtracksAnalyzedInPhiBins+i);
      */
      
      man->AddHistogram(classStr.Data(),"TZEROpileup", "TZERO pileup", kFALSE, 2, -0.5, 1.5, AliReducedVarManager::kTZEROpileup);
      man->AddHistogram(classStr.Data(),"TZEROpileup_TimeFromSOR_prof", "TZERO pileup vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 2, -0.5, 1.5, AliReducedVarManager::kTZEROpileup);
      man->AddHistogram(classStr.Data(),"TZEROsatellite", "TZERO satellite", kFALSE, 2, -0.5, 1.5, AliReducedVarManager::kTZEROsatellite);
      man->AddHistogram(classStr.Data(),"TZEROsatellite_TimeFromSOR_prof", "TZERO satellite vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 2, -0.5, 1.5, AliReducedVarManager::kTZEROsatellite);
      
      man->AddHistogram(classStr.Data(),"NTracksITSin","",kFALSE,300,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kITSin);
      man->AddHistogram(classStr.Data(),"NTracksITSout","",kFALSE,300,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kITSout);
      man->AddHistogram(classStr.Data(),"NTracksTPCin","",kFALSE,300,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCin);
      man->AddHistogram(classStr.Data(),"NTracksTPCout","",kFALSE,1800,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout);
      man->AddHistogram(classStr.Data(),"NTracksTRDout","",kFALSE,300,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTRDout);
      man->AddHistogram(classStr.Data(),"NTracksTRDin","",kFALSE,300,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTRDin);
      man->AddHistogram(classStr.Data(),"NTracksTPCout_VZEROmult","",kFALSE,900,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTPCrefit_VZEROmult","",kFALSE,900,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCrefit, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutFromPileup_VZEROmult","",kFALSE,1200,-10000.,30000.,AliReducedVarManager::kNTracksTPCoutFromPileup, 400, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutFromPileup_VZEROmult_nAnalyzedPairs_prof","",kTRUE,40,-10000.,30000.,AliReducedVarManager::kNTracksTPCoutFromPileup, 40, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 10, 0., 100., AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutFromPileup_TimeFromSOR_nAnalyzedPairs_prof","",kTRUE, 45, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 40,-10000.,30000.,AliReducedVarManager::kNTracksTPCoutFromPileup, 10, 0., 100., AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"NTracksTPCout_VZEROmult_nAnalyzedPairs_prof","",kTRUE,60,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 40, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 10, 0., 100., AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutFromPileup","",kFALSE,2000,-10000.,30000.,AliReducedVarManager::kNTracksTPCoutFromPileup);
      man->AddHistogram(classStr.Data(),"NTracksTRDout_VZEROmult","",kFALSE,300,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTRDout, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTRDrefit_VZEROmult","",kFALSE,300,0.,15000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTRDrefit, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksITSout_VZEROmult","",kFALSE,300,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kITSout, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksITSrefit_VZEROmult","",kFALSE,300,0.,30000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kITSrefit, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTOFout_VZEROmult","",kFALSE,300,0.,3000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTOFout, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTOFrefit_VZEROmult","",kFALSE,300,0.,3000.,AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTOFrefit, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTotal_VZEROmult","",kFALSE,300,0.,30000.,AliReducedVarManager::kNtracksTotal, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksSelected_VZEROmult","",kFALSE,300,0.,300.,AliReducedVarManager::kNtracksSelected, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksAnalyzed_VZEROmult","",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksAnalyzed, 300, 0., 40000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutVsVZEROmult","",kFALSE,1000,0.,10.,AliReducedVarManager::kNTracksTPCoutVsVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksTPCoutVsVZEROmult_CentVZERO","",kFALSE, 20, 0., 100., AliReducedVarManager::kCentVZERO, 1000,0.,10.,AliReducedVarManager::kNTracksTPCoutVsVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"NTracksAnalyzed_RunNo_prof","",kTRUE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 10, 0., 10., AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NPairsAnalyzed_RunNo_prof","",kTRUE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 10, 0., 10., AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"EventAverageTPCchi2","Average TPC chi2 per track",kFALSE,100,0.,5.,AliReducedVarManager::kEvAverageTPCchi2);
      man->AddHistogram(classStr.Data(),"EventAverageTPCchi2_run_prof","Average TPC chi2 per track vs run",kTRUE,runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 100,0.,5.,AliReducedVarManager::kEvAverageTPCchi2);
      man->AddHistogram(classStr.Data(),"EventAverageTPCchi2_NTPCout","Average TPC chi2 per track vs n TPCout",kFALSE,100,0.,5.,AliReducedVarManager::kEvAverageTPCchi2, 50, 0.0, 15000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout);
      man->AddHistogram(classStr.Data(),"EventAverageTPCchi2_NTPCout_prof","Average TPC chi2 per track vs n TPCout",kTRUE,50, 0.0, 15000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 100,0.,5.,AliReducedVarManager::kEvAverageTPCchi2);
      
      Bool_t isCalibrated = kTRUE;
      man->AddHistogram(classStr.Data(), "VZEROA_NEmptyChannels_VtxCent_prof", "No. VZERO-A empty channels per event vs. centrality SPD and vertex Z", kTRUE,
                   24, -12.0, +12.0, AliReducedVarManager::kVtxZ, 20, 0.0, 100., AliReducedVarManager::kCentSPD, 32, 0.0, 32.0, AliReducedVarManager::kVZEROAemptyChannels);
      man->AddHistogram(classStr.Data(), "VZEROC_NEmptyChannels_VtxCent_prof", "No. VZERO-C empty channels per event vs. centrality SPD and vertex Z", kTRUE,
                   24, -12.0, +12.0, AliReducedVarManager::kVtxZ, 20, 0.0, 100., AliReducedVarManager::kCentSPD, 32, 0.0, 32.0, AliReducedVarManager::kVZEROCemptyChannels);
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
       //man->AddHistogram(classStr.Data(), "ITSlayerHit_DCAxy", "Hits in the ITS layers vs DCA_{xy}", kFALSE,
       //                  1000, -0.5, 0.5, AliReducedVarManager::kDcaXY, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
       //man->AddHistogram(classStr.Data(), "ITSlayerHit_DCAz", "Hits in the ITS layers vs DCA_{z}", kFALSE,
       //                  1800, -1.0, 1.0, AliReducedVarManager::kDcaZ, 6, 0.5, 6.5, AliReducedVarManager::kITSlayerHit);
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
       /*man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_TrackingFlags","Corrected TPC N_{#sigma} electron vs. inner param P vs tracking flags;;",kFALSE,
                         50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, AliReducedVarManager::kNTrackingStatus, -0.5, AliReducedVarManager::kNTrackingStatus-0.5, AliReducedVarManager::kTrackingFlag, "", "", trkStatusNames.Data());
       man->AddHistogram(classStr.Data(),"TPCncls_Eta_NTracksTPCoutFromPileup_TrackingFlags_prof","",kTRUE, 20,-1.0,1.0,AliReducedVarManager::kEta,
                         43,-1500., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, AliReducedVarManager::kNTrackingStatus, -0.5, AliReducedVarManager::kNTrackingStatus-0.5, AliReducedVarManager::kTrackingFlag, "", "", trkStatusNames.Data(),  AliReducedVarManager::kTPCncls);*/
       continue;
    }
    
    if(classStr.Contains("Track_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "P_Pin", "P vs Pin", kFALSE, 200, 0.0, 10.0, AliReducedVarManager::kP, 200, 0.0, 10.0, AliReducedVarManager::kPin);
      man->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution", kFALSE, 1000, 0.0, 50.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Pt_TimeFromSOR", "<p_{T}> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, 0.0, 50.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Eta", "", kFALSE, 1000, -1.5, 1.5, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Eta_Phi", "", kFALSE, 100, -1.0, 1.0, AliReducedVarManager::kEta, 100, 0., 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "Eta_TimeFromSOR", "<#eta> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, -1.5, 1.5, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Phi", "", kFALSE, 1000, 0.0, 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "Phi_TimeFromSOR", "<#varphi> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, 0.0, 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "DCAxy_Pt", "DCAxy", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 50, 0.0, 5.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAxy_TPCchi2", "DCAxy vs TPC chi2", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 50, 0.0, 5.0, AliReducedVarManager::kTPCchi2);
      man->AddHistogram(classStr.Data(), "DCAz_TPCchi2", "DCAz vs TPC chi2", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kDcaZ, 50, 0.0, 5.0, AliReducedVarManager::kTPCchi2);
      man->AddHistogram(classStr.Data(), "DCAxy_ITSchi2", "DCAxy vs ITS chi2", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 100, 0.0, 25.0, AliReducedVarManager::kITSchi2);
      man->AddHistogram(classStr.Data(), "DCAz_ITSchi2", "DCAz vs ITS chi2", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kDcaZ, 100, 0.0, 25.0, AliReducedVarManager::kITSchi2);
      man->AddHistogram(classStr.Data(), "DCAxy_goldenChi2", "DCAxy vs golden chi2", kFALSE, 100, -2.0, 2.0, AliReducedVarManager::kDcaXY, 100, 0.0, 100.0, AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
      man->AddHistogram(classStr.Data(), "DCAz_goldenChi2", "DCAz vs golden chi2", kFALSE, 100, -5.0, 5.0, AliReducedVarManager::kDcaZ, 100, 0.0, 100.0, AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE, 200, -5.0, 5.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAxy_TimeFromSOR", "<DCAxy> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, -5.0, 5.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE, 200, -5.0, 5.0, AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(), "DCAz_Pt", "DCAz", kFALSE, 100, -3.0, 3.0, AliReducedVarManager::kDcaZ, 50, 0.0, 5.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "DCAz_TimeFromSOR", "<DCAz> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, -5.0, 5.0, AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(),"DCAxy_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(),"DCAz_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kDcaZ);
      
        man->AddHistogram(classStr.Data(),"ITSncls", "ITS nclusters", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"ITSnclsShared", "ITS nclusters shared", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSnclsShared);
        man->AddHistogram(classStr.Data(),"ITSncls_TimeFromSOR", "<ITS nclusters> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 7,-0.5,6.5,AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"ITSnclsShared_TimeFromSOR", "<ITS nclusters shared> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 7,-0.5,6.5,AliReducedVarManager::kITSnclsShared);
        man->AddHistogram(classStr.Data(),"Eta_Phi_ITSncls_prof","ITS <nclusters> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"Eta_Phi_ITSnclsShared_prof","ITS <nclusters-shared> vs (#eta,#phi)",kTRUE,
                          36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kITSnclsShared);
	man->AddHistogram(classStr.Data(),"ITSchi2", "ITS #chi^{2}", kFALSE, 200,0.0,50.0, AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"ITSchi2_RunNo_prof","",kTRUE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 10, 0., 10., AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"ITSchi2_ITSncls", "ITS #chi^{2} vs ITS clusters", kFALSE, 100,0.0,100.0, AliReducedVarManager::kITSchi2,7,-0.5,6.5, AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"ITSchi2_TimeFromSOR", "<ITS #chi^{2}> vs time from SOR", kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200,0.0,20.0, AliReducedVarManager::kITSchi2);
	man->AddHistogram(classStr.Data(),"Eta_Phi_ITSchi2_prof","ITS <#chi^{2}> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 100, 0.0, 1000, AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"ITSncls_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"ITSchi2_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"ITSnSharedCls_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kITSnclsShared);
        
        man->AddHistogram(classStr.Data(),"Chi2TPCConstrainedVsGlobal","",kFALSE,200,0.,200.,AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
        man->AddHistogram(classStr.Data(),"Chi2TPCConstrainedVsGlobal_TimeFromSOR_prof","<Chi2 TPC constrained vs global> vs time from SOR",kTRUE,90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
        man->AddHistogram(classStr.Data(),"Chi2TPCConstrainedVsGlobal_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
        man->AddHistogram(classStr.Data(),"Chi2TPCConstrainedVsGlobal_Eta_MultVZERO_NTracksTPCout_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 20,0., 30000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
        
        man->AddHistogram(classStr.Data(),"TPCncls","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCncls_RunNo_prof","",kTRUE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 10, 0., 10., AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCncls_TimeFromSOR","<TPC #cls> vs time from SOR",kTRUE,90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 160,-0.5,159.5,AliReducedVarManager::kTPCncls);
	man->AddHistogram(classStr.Data(),"TPCcrossedRows","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCcrossedRows);
        man->AddHistogram(classStr.Data(),"TPCfindableClusters","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsF);
        //man->AddHistogram(classStr.Data(),"TPCfindableClusters_TPCcrossedRows","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsF, 160, -0.5, 159.5, AliReducedVarManager::kTPCcrossedRows);
        man->AddHistogram(classStr.Data(),"TPCcrossedRowsOverFindableClusters","", kFALSE, 200,0.0,1.5,AliReducedVarManager::kTPCcrossedRowsOverFindableClusters);
        man->AddHistogram(classStr.Data(),"TPCcrossedRowsOverFindableClusters_TimeFromSOR","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,1.,AliReducedVarManager::kTPCcrossedRowsOverFindableClusters);
        //man->AddHistogram(classStr.Data(),"TPCcrossedRowsOverFindableClusters_TPCnclsSharedRatio","", kFALSE, 200,0.0,1.5,AliReducedVarManager::kTPCcrossedRowsOverFindableClusters, 200, 0.0, 1.0, AliReducedVarManager::kTPCnclsSharedRatio);
	man->AddHistogram(classStr.Data(),"TPCnclsShared","",kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsShared);
        man->AddHistogram(classStr.Data(),"TPCnclsSharedRatio","",kFALSE, 200,0.,1.,AliReducedVarManager::kTPCnclsSharedRatio);
        //man->AddHistogram(classStr.Data(),"TPCnclsShared_TPCncls","",kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsShared, 160, -0.5, 159.5, AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCnclsShared_TimeFromSOR","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsShared);
        man->AddHistogram(classStr.Data(),"TPCnclsSharedRatio_TimeFromSOR","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,1.,AliReducedVarManager::kTPCnclsSharedRatio);
        man->AddHistogram(classStr.Data(),"TPCnclsShared_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCnclsShared);
        man->AddHistogram(classStr.Data(),"TPCnclsShared_Eta_MultVZERO_NTracksTPCout_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 20,0., 30000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCnclsShared);
        man->AddHistogram(classStr.Data(),"TPCncls_Eta_NTracksTPCoutFromPileup_prof","",kTRUE, 20,-1.0,1.0,AliReducedVarManager::kEta,
          43,-1500., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 10, 0., 160., AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCncls_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCncls_Eta_MultVZERO_NTracksTPCout_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 20,0., 30000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCncls);
        
        
        man->AddHistogram(classStr.Data(),"TPCchi2","",kFALSE, 200,0.0,10.0,AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_nTPCout","",kFALSE, 100, 0.0, 15000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 200,0.0,10.0,AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_Pt","",kFALSE, 200,0.0,10.0,AliReducedVarManager::kTPCchi2, 100, 0., 10.0, AliReducedVarManager::kPt);
        man->AddHistogram(classStr.Data(),"TPCchi2_RunNo_prof","",kTRUE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 10, 0., 10., AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_TPCncls", "TPC #chi^{2} vs TPC clusters", kFALSE, 100,0.0,10.0, AliReducedVarManager::kTPCchi2,160,0.,160., AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCsegments","",kFALSE, 9,0.0,9.0,AliReducedVarManager::kTPCNclusBitsFired);
        man->AddHistogram(classStr.Data(),"TPCsegments_TimeFromSOR","",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 8,0.0,8.0,AliReducedVarManager::kTPCNclusBitsFired);
        man->AddHistogram(classStr.Data(),"TPCsegments_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCNclusBitsFired);
        man->AddHistogram(classStr.Data(),"TPCchi2_TimeFromSOR","TPC <#chi^{2}> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200,0.0,10.0,AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_Eta_MultVZERO_NTracksTPCout_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 30,0., 30000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCchi2);
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
        man->AddHistogram(classStr.Data(),"TPCsignalN_TimeFromSOR","TPC <#cls pid> vs time from SOR",kTRUE, 90, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 160,-0.5,159.5,AliReducedVarManager::kTPCsignalN);
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCsignalN_prof","TPC <nclusters pid> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCsignalN);    
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_TimeFromSOR","TPC N_{#sigma} electron vs. time from SOR",kFALSE,
                          90, 0., 450.,AliReducedVarManager::kTimeRelativeSOR,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_TimeFromSOR_prof","<TPC N_{#sigma} electron> vs. time from SOR",kTRUE,
                          90, 0., 450.,AliReducedVarManager::kTimeRelativeSOR,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Pin","TPC N_{#sigma} electron vs. inner param P",kFALSE,
                     100,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        /*man->AddHistogram(classStr.Data(),"TPCnsigElectron_Pin_TimeFromSOR","TPC N_{#sigma} electron vs. inner param P vs time from SOR",kFALSE,
                          100,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, 15, 0., 450., AliReducedVarManager::kTimeRelativeSOR);*/
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Eta","TPC N_{#sigma} electron vs. #eta",kFALSE,
                          36,-0.9,0.9,AliReducedVarManager::kEta,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        /*man->AddHistogram(classStr.Data(),"TPCnsigElectron_Eta_TimeFromSOR","TPC N_{#sigma} electron vs. #eta vs time from SOR",kFALSE,
                          36,-0.9,0.9,AliReducedVarManager::kEta,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, 15, 0., 450., AliReducedVarManager::kTimeRelativeSOR);*/
        man->AddHistogram(classStr.Data(),"TPCnSigEle_Eta_MultVZERO_NTracksTPCoutFromPileup_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 22,-2000., 20000., AliReducedVarManager::kNTracksTPCoutFromPileup, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnSigEle_Eta_MultVZERO_NTracksTPCout_prof","",kTRUE, 10, 0., 40000., AliReducedVarManager::kVZEROTotalMult, 20,0., 30000., AliReducedVarManager::kNTracksPerTrackingStatus+AliReducedVarManager::kTPCout, 20,-1.0,1.0, AliReducedVarManager::kEta, "","","", AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin","Corrected TPC N_{#sigma} electron vs. inner param P",kFALSE,
                          100,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron);
        
        /*
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_DCAxy","DCAxy vs corrected TPC N_{#sigma} electron vs. inner param P",kFALSE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 80, -2., 2., AliReducedVarManager::kDcaXY);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_DCAz","DCAz vs corrected TPC N_{#sigma} electron vs. inner param P",kFALSE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 120, -3., 3., AliReducedVarManager::kDcaZ);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_ITSchi2","ITS chi2 vs corrected TPC N_{#sigma} electron vs. inner param P",kFALSE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 50, 0., 50., AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_ITSchi2_prof","ITS <chi2> vs corrected TPC N_{#sigma} electron vs. inner param P",kTRUE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 50, 0., 50., AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_TPCchi2","TPC chi2 vs corrected TPC N_{#sigma} electron vs. inner param P",kFALSE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 50, 0., 5., AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_TPCchi2_prof","TPC <chi2> vs corrected TPC N_{#sigma} electron vs. inner param P",kTRUE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 50, 0., 5., AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_Chi2TPCConstrainedVsGlobal","TPC constrained vs global chi2 vs corrected TPC N_{#sigma} electron vs. inner param P",kFALSE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 50, 0., 200., AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_Chi2TPCConstrainedVsGlobal_prof","TPC constrained vs global <chi2> vs corrected TPC N_{#sigma} electron vs. inner param P",kTRUE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 50, 0., 200., AliReducedVarManager::kChi2TPCConstrainedVsGlobal);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_NclsITS_prof","<#ITScls> vs corrected TPC N_{#sigma} electron vs. inner param P",kTRUE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 7, 0., 7., AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_NclsTPC_prof","<#TPCcls> vs corrected TPC N_{#sigma} electron vs. inner param P",kTRUE,
                          50,0.0,5.0,AliReducedVarManager::kPin,50,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 80, 0., 160., AliReducedVarManager::kTPCncls);
        */
        
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Eta","TPC N_{#sigma} electron vs. #eta",kFALSE,
                          36,-0.9,0.9,AliReducedVarManager::kEta,100,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_TimeFromSOR_prof","<corrected TPC N_{#sigma} electron> vs. time from SOR",kTRUE,
                          90, 0., 450.,AliReducedVarManager::kTimeRelativeSOR,80,-4.0,4.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron);
        
        /*man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin_MassUsedForTracking_prof","<mass used for tracking> vs corrected TPC N_{#sigma} electron vs. inner param P",kTRUE,
                          100,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 40, 0., 4., AliReducedVarManager::kMassUsedForTracking);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Eta_MassUsedForTracking_prof","<mass used for tracking> vs corrected TPC N_{#sigma} electron vs. eta",kTRUE,
                          36,-0.9,0.9,AliReducedVarManager::kEta,100,-5.0,5.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, 40, 0., 4., AliReducedVarManager::kMassUsedForTracking);*/
	man->AddHistogram(classStr.Data(),"TPCnsigElectron_Run","TPC N_{#sigma} electron vs. run",kTRUE,
                     runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      
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
    
    if(classStr.Contains("PairQualityFlags")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       TString pairQualityFlagNames = " ;K^{0}_{S}#rightarrow#pi^{+}#pi^{-};#Lambda#rightarrow p#pi^{-};#bar{#Lambda}#rightarrow #bar{p}#pi^{+};#gamma#rightarrow e^{+}e^{-};";
       man->AddHistogram(classStr.Data(), "PairQualityFlags", "Pair quality flags;;", kFALSE,
                         32, -0.5, 31.5, AliReducedVarManager::kPairQualityFlag, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, pairQualityFlagNames.Data());
       continue;
    }
    
    Double_t massBinWidth = 0.001;     // *GeV/c^2
    Double_t massRange[2] = {0.0,5.0};
    Int_t nMassBins = TMath::Nint((massRange[1]-massRange[0])/massBinWidth);
    
    TString candidateNames = "#gamma#rightarrow e^{+}e^{-};K^{0}_{S}#rightarrow#pi^{+}#pi^{-};";
    candidateNames += "#Lambda#rightarrow p#pi^{-};#bar{#Lambda}#rightarrow #bar{p}#pi^{+};";
    
    // Histograms for pairs
    if(classStr.Contains("Pair")) {
       man->AddHistClass(classStr.Data());
       cout << classStr.Data() << endl;
       man->AddHistogram(classStr.Data(), "CandidateId", "Candidate id", kFALSE,
                         5, -0.5, 4.5, AliReducedVarManager::kCandidateId, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, candidateNames.Data());
       man->AddHistogram(classStr.Data(), "PairType", "Pair type", kFALSE, 4, -0.5, 3.5, AliReducedVarManager::kPairType);
       man->AddHistogram(classStr.Data(), "PairChi2", "Pair #chi^{2}", kFALSE, 200, 0.0, 50, AliReducedVarManager::kPairChisquare);	
       man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE, nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMass);
       man->AddHistogram(classStr.Data(), "Mass_V0K0s", "Invariant mass, K^{0}_{s} assumption", kFALSE,
                         nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0);
       man->AddHistogram(classStr.Data(), "Mass_V0Lambda", "Invariant mass, #Lambda^{0} assumption", kFALSE,
                         nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0+1);
       man->AddHistogram(classStr.Data(), "Mass_V0ALambda", "Invariant mass, #bar{#Lambda^{0}} assumption", kFALSE,
                         nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0+2);
       man->AddHistogram(classStr.Data(), "Mass_V0Gamma", "Invariant mass, #gamma conversion assumption", kFALSE,
                         nMassBins, massRange[0], massRange[1], AliReducedVarManager::kMassV0+3);
       man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPt);
       man->AddHistogram(classStr.Data(), "Px", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPx);
       man->AddHistogram(classStr.Data(), "Py", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPy);
       man->AddHistogram(classStr.Data(), "Pz", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPz);
       man->AddHistogram(classStr.Data(), "P", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kP);
       man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kRap);
       man->AddHistogram(classStr.Data(), "Eta", "Pseudo-rapidity #eta", kFALSE, 240, -2.0, 2.0, AliReducedVarManager::kEta);
       man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 315, 0.0, 6.3, AliReducedVarManager::kPhi);
       man->AddHistogram(classStr.Data(), "Theta", "#theta distribution", kFALSE, 1000, 0.0, 3.2, AliReducedVarManager::kTheta);
       man->AddHistogram(classStr.Data(), "LxyOrR", "L_{xy}/Decay Radius", kFALSE, 1000, -10.0, 10.0, AliReducedVarManager::kPairLxy);
       man->AddHistogram(classStr.Data(), "OpeningAngle", "Opening angle", kFALSE, 1000, 0.0, 3.2, AliReducedVarManager::kPairOpeningAngle);
       man->AddHistogram(classStr.Data(), "PointingAngle", "Pointing angle", kFALSE, 1000, 0.0, 3.2, AliReducedVarManager::kPairPointingAngle);
       
       continue;
    }   // end if for Pair classes of histograms    
  }  // end loop over histogram classes
}
