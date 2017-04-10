void Setup(AliReducedAnalysisTaskSE* processor, TString prod="LHC10h");
AliHistogramManager* SetupHistogramManager(TString prod="LHC10h");
void DefineHistograms(AliHistogramManager* man, TString prod="LHC10h");

//__________________________________________________________________________________________
AliAnalysisTask* AddTask_iarsene_testTask(Bool_t isAliRoot=kTRUE, Int_t runMode=1, TString prod="LHC10h"){    
   //
   // isAliRoot=kTRUE for ESD/AOD analysis in AliROOT, kFALSE for root analysis on reduced trees
   // runMode=1 (AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents)
   //               =2 (AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree)
   //
   //get the current analysis manager

  printf("INFO on AddTask_iarsene_testTask(): (isAliRoot, runMode) :: (%d,%d)", isAliRoot, runMode);

  AliReducedAnalysisTest* testAnalysis = new AliReducedAnalysisTest("TestAnalysis","Test analysis");
  testAnalysis->Init();
  Setup(testAnalysis, prod);
  // initialize an AliAnalysisTask which will wrapp the AliReducedAnalysisTest such that it can be run in an aliroot analysis train (e.g. LEGO, local analysis etc)
  AliAnalysisTaskReducedEventProcessor* task = new AliAnalysisTaskReducedEventProcessor("ReducedEventAnalysisManager", runMode);
  task->AddTask(testAnalysis);
  
  if(isAliRoot){
     AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
     if (!mgr) {
        Error("AddTask_iarsene_dst", "No analysis manager found.");
        return 0;
     }
     
     AliAnalysisDataContainer* cReducedEvent = NULL;
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) {
       printf("INFO on AddTask_iarsene_testTask(): use on the fly events ");
       cReducedEvent = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ReducedEventDQ");
       if(!cReducedEvent) {
         printf("ERROR: In AddTask_iarsene_testTask(), couldn't find exchange container with ReducedEvent! ");
         printf("             You need to use AddTask_iarsene_dst() in the on-the-fly reduced events mode!");
         return 0x0;
       }
     }
            
     mgr->AddTask(task);
      
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree) 
        mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
      
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) 
       mgr->ConnectInput(task, 0, cReducedEvent);
  
     AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("testHistos", THashList::Class(),
                                                                  AliAnalysisManager::kOutputContainer, "dstAnalysisHistograms.root");
     mgr->ConnectOutput(task, 1, cOutputHist );
  }
  
  return task;
}


//_________________________________________________________________
void Setup(AliReducedAnalysisTest* processor, TString prod /*="LHC10h"*/) {
  //
  // Configure the event processor
  // Setup histograms, handlers, cuts, etc.
  //
   if(prod.Contains("LHC15o")) {
      TFile* grpFile = TFile::Open("/home/iarsene/work/ALICE/treeAnalysis/newdst/development/test/grpData_LHC15o_allRuns.root");
      AliReducedVarManager::SetLHCDataInfo((TH1F*)grpFile->Get("lumiTotal"), (TH1F*)grpFile->Get("intTotal0"), 
                                           (TH1F*)grpFile->Get("intTotal1"), (TH1I*)grpFile->Get("fillNumber"));
      AliReducedVarManager::SetGRPDataInfo((TH1I*)grpFile->Get("dipolePolarity"), (TH1I*)grpFile->Get("l3Polarity"), 
                                           (TH1I*)grpFile->Get("timeStart"), (TH1I*)grpFile->Get("timeStop"));
      grpFile->Close();
   }
   
  SetupHistogramManager(processor->GetHistogramManager(), prod);
  
  // Set event cuts
  AliReducedEventCut* evCut1 = new AliReducedEventCut("Centrality","Centrality selection");
  evCut1->AddCut(AliReducedVarManager::kCentVZERO, 0., 90.);
  AliReducedEventCut* evCut2 = new AliReducedEventCut("VertexZ","Vertex selection");
  evCut2->AddCut(AliReducedVarManager::kVtxZ, -25.0, 25.0);
  //processor->AddEventCut(evCut1);
  processor->AddEventCut(evCut2);
  
  // Set track cuts
  AliReducedTrackCut* trackCut1 = new AliReducedTrackCut("Pt","Pt selection");
  trackCut1->AddCut(AliReducedVarManager::kPt, 0.,100.0);
  trackCut1->AddCut(AliReducedVarManager::kEta, -1.5,1.5);
  processor->AddTrackCut(trackCut1);  
  
  // Set pair cuts
  AliReducedTrackCut* pairCut1 = new AliReducedTrackCut("Ptpair","Pt pair selection");
  pairCut1->AddCut(AliReducedVarManager::kPt, 0.0,100.0);
  pairCut1->AddCut(AliReducedVarManager::kEta, -1.5,1.5);
  processor->AddPairCut(pairCut1);
}


//_________________________________________________________________
void SetupHistogramManager(AliHistogramManager* man, TString prod /*="LHC10h"*/) {
  //
  // setup the histograms manager
  //
  AliReducedVarManager::SetDefaultVarNames();
  
  DefineHistograms(man, prod);
  
  AliReducedVarManager::SetUseVars(man->GetUsedVars());
}


//_________________________________________________________________
void DefineHistograms(AliHistogramManager* man, TString prod /*="LHC10h"*/) {
  //
  // define histograms
  //
  TString histClasses = "";
  histClasses += "Event_NoCuts;";        //ok
  histClasses += "Event_AfterCuts;";     //ok
  histClasses += "EvtTags;L0TriggerInput;L1TriggerInput;L2TriggerInput;";   //ok
  histClasses += "OnlineTriggers_vs_L0TrigInputs;OnlineTriggers_vs_L1TrigInputs;OnlineTriggers_vs_L2TrigInputs;";  //ok
  histClasses += "ITSclusterMap;TPCclusterMap;";   //ok
  histClasses += "CaloClusters;";   //ok
  histClasses += "OnlineTriggers_NoCuts;";  //ok
  histClasses += "OnlineTriggers_AfterCuts;";    //ok
  histClasses += "TrackQA_GammaLeg;";         //ok
  histClasses += "TrackQA_PureGammaLeg;";
  histClasses += "TrackQA_K0sLeg;";
  histClasses += "TrackQA_PureK0sLeg;";
  histClasses += "TrackQA_LambdaPosLeg;";
  histClasses += "TrackQA_PureLambdaPosLeg;";
  histClasses += "TrackQA_LambdaNegLeg;";
  histClasses += "TrackQA_PureLambdaNegLeg;";
  histClasses += "TrackQA_ALambdaPosLeg;";
  histClasses += "TrackQA_PureALambdaPosLeg;";
  histClasses += "TrackQA_ALambdaNegLeg;";
  histClasses += "TrackQA_PureALambdaNegLeg;";
  histClasses += "TrackQA_AllTracks;";            //ok
  histClasses += "TrackingFlags;TrackQualityFlags;PairQualityFlags_Offline;PairQualityFlags_OnTheFly;";   //ok
  histClasses += "PairQA_OfflineGamma;";           // ok
  histClasses += "PairQA_OfflinePureGamma;";
  histClasses += "PairQA_OnTheFlyGamma;";
  histClasses += "PairQA_OnTheFlyPureGamma;";
  histClasses += "PairQA_OfflineK0s;";
  histClasses += "PairQA_OfflinePureK0s;";
  histClasses += "PairQA_OnTheFlyK0s;";
  histClasses += "PairQA_OnTheFlyPureK0s;";
  histClasses += "PairQA_OfflineLambda;";
  histClasses += "PairQA_OfflinePureLambda;";
  histClasses += "PairQA_OnTheFlyLambda;";
  histClasses += "PairQA_OnTheFlyPureLambda;";
  histClasses += "PairQA_OfflineALambda;";
  histClasses += "PairQA_OfflinePureALambda;";
  histClasses += "PairQA_OnTheFlyALambda;";
  histClasses += "PairQA_OnTheFlyPureALambda;";
    
  Int_t runNBins = 0;
  Double_t runHistRange[2] = {0.0,0.0};
  
  // Pb-Pb from 2010 run range is default: LHC10h
  runNBins = 2500;
  runHistRange[0] = 137100.;
  runHistRange[1] = 139600.;
  
  // Pb-Pb of 2011
  if(!prod.CompareTo("LHC11h")) {
    runNBins = 2700;
    runHistRange[0] = 167900.;
    runHistRange[1] = 170600.;
  }
  
  // Pb-Pb of 2015
  if(!prod.CompareTo("LHC15o")) {
     runNBins = 2100;
     runHistRange[0] = 244900.;
     runHistRange[1] = 247000.;
  }
  
  // p-Pb of 2013
  if(!prod.CompareTo("LHC13b") || !prod.CompareTo("LHC13c")) {
    runNBins = 400;
    runHistRange[0] = 195300.;
    runHistRange[1] = 195700.;
  }
  
  const Int_t kNBCBins = 3600;
  Double_t bcHistRange[2] = {-0.5,3599.5};
  
  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");
  
  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    
    // Event wise histograms
    if(classStr.Contains("Event")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(),"RunNo","Run numbers",kFALSE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo);
      man->AddHistogram(classStr.Data(),"VtxX","Vtx X",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxY","Vtx Y",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxZ","Vtx Z",kFALSE,300,-15.,15.,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"NVtxContributors","",kFALSE,500,0.,10000.,AliReducedVarManager::kNVtxContributors);      
      man->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentVZERO);
      man->AddHistogram(classStr.Data(),"CentVZEROA","Centrality(VZERO-A)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentVZEROA);
      man->AddHistogram(classStr.Data(),"CentVZEROC","Centrality(VZERO-C)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentVZEROC);
      man->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentSPD);
      man->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentTPC);
      man->AddHistogram(classStr.Data(),"CentZDC","Centrality(ZDC)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentZDC);
      man->AddHistogram(classStr.Data(),"CentZNA","Centrality(ZNA)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentZNA);
      man->AddHistogram(classStr.Data(),"CentQuality","Centrality quality",kFALSE, 100, -50.5, 49.5, AliReducedVarManager::kCentQuality);
      
      TString multEstimators = "OnlineV0M;OnlineV0A;OnlineV0C;ADM;ADA;ADC;SPDClusters;SPDTracklets;RefMult05;RefMult08";
      TObjArray* arrEstim=multEstimators.Tokenize(";");
      
      for(Int_t i=0;i<10;++i) {
         TString estimStr = arrEstim->At(i)->GetName();
         man->AddHistogram(classStr.Data(), Form("MultEstimator_%s", estimStr.Data()), Form("Multiplicity estimator %s", estimStr.Data()), kFALSE, 100, 0.0, 1000, AliReducedVarManager::kMultEstimatorOnlineV0M+i);
         man->AddHistogram(classStr.Data(), Form("MultEstimatorPercentile_%s", estimStr.Data()), Form("Multiplicity estimator percentile %s", estimStr.Data()), kFALSE, 100, 0.0, 100, AliReducedVarManager::kMultEstimatorPercentileOnlineV0M+i);
      }
      
      man->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event",kFALSE,500,0.,20000.,AliReducedVarManager::kNtracksTotal);
      man->AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event",kFALSE,200,0.,200.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"NV0sTotal","Number of V0 candidates per event",kFALSE,200,0.,30000.,AliReducedVarManager::kNV0total);
      man->AddHistogram(classStr.Data(),"NV0sSelected","Number of selected V0 candidates per event",kFALSE,200,0.,2000.,AliReducedVarManager::kNV0selected);
      man->AddHistogram(classStr.Data(),"EventNumberInESDFile","Event number in ESD file",kFALSE, 1000, 0.0, 1000.0, AliReducedVarManager::kEventNumberInFile);
      man->AddHistogram(classStr.Data(),"BC","Bunch crossing",kFALSE,3500,0.,3500.,AliReducedVarManager::kBC);
      man->AddHistogram(classStr.Data(),"TimeFromSOR","Events vs time",kFALSE, 450, 0., 450., AliReducedVarManager::kTimeRelativeSOR);
      man->AddHistogram(classStr.Data(),"EventType","Event type",kFALSE,100,0.,100.,AliReducedVarManager::kEventType);
      man->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsPhysicsSelection, 0,0.0,0.0,AliReducedVarManager::kNothing, 0,0.0,0.0,AliReducedVarManager::kNothing, "off;on");
      man->AddHistogram(classStr.Data(),"IsSPDPileup","Event has pileup (SPD)",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsSPDPileup, 0,0.0,0.0,AliReducedVarManager::kNothing, 0,0.0,0.0,AliReducedVarManager::kNothing, "no;yes");
      man->AddHistogram(classStr.Data(),"IsSPDPileupMultBins","Event has pileup (SPD multiplicity bins)",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsSPDPileupMultBins, 0,0.0,0.0,AliReducedVarManager::kNothing, 0,0.0,0.0,AliReducedVarManager::kNothing, "no;yes");
      man->AddHistogram(classStr.Data(),"IRIntClosestInt1Map","Closest out of bunch interactions Int1",kFALSE,7000,-3500.,3500.,AliReducedVarManager::kIRIntClosestIntMap);
      man->AddHistogram(classStr.Data(),"IRIntClosestInt2Map","Closest out of bunch interactions Int2",kFALSE,7000,-3500.,3500.,AliReducedVarManager::kIRIntClosestIntMap+1);
      man->AddHistogram(classStr.Data(),"VtxXtpc","Vtx X (TPC)",kFALSE,300,-1.,1.,AliReducedVarManager::kVtxXtpc);
      man->AddHistogram(classStr.Data(),"VtxYtpc","Vtx Y (TPC)",kFALSE,300,-1.,1.,AliReducedVarManager::kVtxYtpc);
      man->AddHistogram(classStr.Data(),"VtxZtpc","Vtx Z (TPC)",kFALSE,300,-15.,15.,AliReducedVarManager::kVtxZtpc);
      man->AddHistogram(classStr.Data(),"DeltaVtxZ","Z_{global}-Z_{TPC}",kFALSE,300,-1.,1.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"NVtxTPCContributors","",kFALSE,500,0.,10000.,AliReducedVarManager::kNVtxTPCContributors);
      man->AddHistogram(classStr.Data(),"VtxXspd","Vtx X (SPD)",kFALSE,300,-1.,1.,AliReducedVarManager::kVtxXspd);
      man->AddHistogram(classStr.Data(),"VtxYspd","Vtx Y (SPD)",kFALSE,300,-1.,1.,AliReducedVarManager::kVtxYspd);
      man->AddHistogram(classStr.Data(),"VtxZspd","Vtx Z (SPD)",kFALSE,300,-15.,15.,AliReducedVarManager::kVtxZspd);
      man->AddHistogram(classStr.Data(),"DeltaVtxZspd","Z_{global}-Z_{SPD}",kFALSE,300,-1.,1.,AliReducedVarManager::kDeltaVtxZspd);
      man->AddHistogram(classStr.Data(),"NVtxSPDContributors","",kFALSE,500,0.,10000.,AliReducedVarManager::kNVtxSPDContributors);
      man->AddHistogram(classStr.Data(),"NSPDpileups","No. SPD pileups",kFALSE,10,0.,10.,AliReducedVarManager::kNSPDpileups);
      man->AddHistogram(classStr.Data(),"NTrackPileups","No. Track pileups",kFALSE,10,0.,10.,AliReducedVarManager::kNTrackPileups);
      man->AddHistogram(classStr.Data(),"NTPCclusters","Number of TPC clusters per event",kFALSE,200,0.,100000.,AliReducedVarManager::kNTPCclusters);
      man->AddHistogram(classStr.Data(),"NPMDtracks","No. PMD tracks",kFALSE,200,0.,200.,AliReducedVarManager::kNPMDtracks);
      man->AddHistogram(classStr.Data(),"NTRDtracks","No. TRD tracks",kFALSE,200,0.,200.,AliReducedVarManager::kNTRDtracks);
      man->AddHistogram(classStr.Data(),"NTRDtracklets","No. TRD tracklets",kFALSE,500,0.,50000.,AliReducedVarManager::kNTRDtracklets);
      man->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0", kFALSE, 200, 0., 5000., AliReducedVarManager::kSPDntracklets);
      for(Int_t i=0; i<32; ++i)
         man->AddHistogram(classStr.Data(),Form("SPDntracklets_%d",i), "", kFALSE, 100, 0., 300., AliReducedVarManager::kSPDntrackletsEta+i);
      for(Int_t il=0; il<2; ++il)
         man->AddHistogram(classStr.Data(), Form("SPDfiredChips_layer%d",il+1), Form("SPD fired chips in layer %d",il+1), 
                           kFALSE, 200, 0., 600., AliReducedVarManager::kSPDFiredChips+il);
      for(Int_t il=0; il<6; ++il)
         man->AddHistogram(classStr.Data(), Form("ITSclusters_layer%d",il+1), Form("ITS clusters in layer %d",il+1), 
                          kFALSE, 200, 0., 10000., AliReducedVarManager::kITSnClusters+il);
      man->AddHistogram(classStr.Data(), "SPDnSingleClusters", "SPD single clusters", 
                          kFALSE, 200, 0., 6000., AliReducedVarManager::kSPDnSingleClusters);	
      for(Int_t i=0; i<64; ++i)
         man->AddHistogram(classStr.Data(), Form("VZEROmult_ch%d",i), "", kFALSE, 200, 0.0, 1000.0, AliReducedVarManager::kVZEROChannelMult+i);
      man->AddHistogram(classStr.Data(),"VZEROAmult", "", kFALSE, 300, 0.0, 20000., AliReducedVarManager::kVZEROATotalMult);
      man->AddHistogram(classStr.Data(),"VZEROCmult", "", kFALSE, 300, 0.0, 20000., AliReducedVarManager::kVZEROCTotalMult);
      man->AddHistogram(classStr.Data(),"VZEROmult", "", kFALSE, 300, 0.0, 30000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"VZEROAemptyChannels", "", kFALSE, 64, 0.0, 64., AliReducedVarManager::kVZEROAemptyChannels);
      man->AddHistogram(classStr.Data(),"VZEROCemptyChannels", "", kFALSE, 64, 0.0, 64., AliReducedVarManager::kVZEROCemptyChannels);
      for(Int_t i=0; i<10; ++i)
         man->AddHistogram(classStr.Data(), Form("ZDCNenergy_ch%d",i), "", kFALSE, 200, 0.0, (i==0 || i==5 ? 100000.0 : 20000.), AliReducedVarManager::kZDCnEnergyCh+i);
      for(Int_t i=0; i<10; ++i)
         man->AddHistogram(classStr.Data(), Form("ZDCPenergy_ch%d",i), "", kFALSE, 200, 0.0, (i==0 || i==5 ? 40000.0 : 20000.), AliReducedVarManager::kZDCpEnergyCh+i);
      for(Int_t i=0; i<26; ++i)
         man->AddHistogram(classStr.Data(), Form("TZEROamp_ch%d",i), "", kFALSE, 200, 0.0, 100.0, AliReducedVarManager::kTZEROAmplitudeCh+i);
      for(Int_t i=0; i<3; ++i)
         man->AddHistogram(classStr.Data(), Form("TZERO_TOF%d",i), "", kFALSE, 200, -1000.0, 1000.0, AliReducedVarManager::kTZEROTOF+i);
      for(Int_t i=0; i<3; ++i)
         man->AddHistogram(classStr.Data(), Form("TZERO_TOFbest%d",i), "", kFALSE, 200, -1000.0, 1000.0, AliReducedVarManager::kTZEROTOFbest+i);
      man->AddHistogram(classStr.Data(),"TZEROzVtx", "", kFALSE, 300, -15.0, 15., AliReducedVarManager::kTZEROzVtx);
      man->AddHistogram(classStr.Data(),"TZEROstartTime", "", kFALSE, 1000, -1000.0, 1000., AliReducedVarManager::kTZEROstartTime);
      man->AddHistogram(classStr.Data(),"TZEROpileup", "TZERO pileup", kFALSE, 2, -0.5, 1.5, AliReducedVarManager::kTZEROpileup);
      man->AddHistogram(classStr.Data(),"TZEROsatellite", "TZERO satellite", kFALSE, 2, -0.5, 1.5, AliReducedVarManager::kTZEROsatellite);

      continue;
    }  // end if className contains "Event"    
    
    if(classStr.Contains("EvtTags")) {
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
	           64, -0.5, 63.5, AliReducedVarManager::kEventTag, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, tagNames.Data());
      continue;
    }
    if(classStr.Contains("L0TriggerInput")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      TString trigInputs = "";
      for(Int_t i=0; i<32; ++i) {trigInputs += Form("L0 Input %d", i); trigInputs+=";";}
      man->AddHistogram(classStr.Data(), "L0TriggerInputs", "L0 trigger inputs", kFALSE,
	           32, -0.5, 31.5, AliReducedVarManager::kL0TriggerInput, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, trigInputs.Data());
      continue;
    }
    if(classStr.Contains("L1TriggerInput")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      TString trigInputs = "";
      for(Int_t i=0; i<32; ++i) {trigInputs += Form("L1 Input %d", i); trigInputs+=";";}
      man->AddHistogram(classStr.Data(), "L1TriggerInputs", "L1 trigger inputs", kFALSE,
	           32, -0.5, 31.5, AliReducedVarManager::kL1TriggerInput, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, trigInputs.Data());
      continue;
    }
    if(classStr.Contains("L2TriggerInput")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      TString trigInputs = "";
      for(Int_t i=0; i<16; ++i) {trigInputs += Form("L2 Input %d", i); trigInputs+=";";}
      man->AddHistogram(classStr.Data(), "L2TriggerInputs", "L2 trigger inputs", kFALSE,
	           16, -0.5, 15.5, AliReducedVarManager::kL2TriggerInput, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, trigInputs.Data());
      continue;
    }
    
    // Offline trigger histograms
    if(classStr.Contains("OnlineTriggers")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      TString triggerNames = "";
      for(Int_t i=0; i<64; ++i) {triggerNames += AliReducedVarManager::fgkOfflineTriggerNames[i]; triggerNames+=";";}
      
      man->AddHistogram(classStr.Data(), "Triggers", "", kFALSE,
	           64, -0.5, 63.5, AliReducedVarManager::kOnlineTrigger, 2, -0.5, 1.5, AliReducedVarManager::kOnlineTriggerFired, 0, 0.0, 0.0, AliReducedVarManager::kNothing, triggerNames.Data(), "off;on");
      man->AddHistogram(classStr.Data(), "Triggers2", "", kFALSE,
	           64, -0.5, 63.5, AliReducedVarManager::kOnlineTriggerFired2, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, triggerNames.Data());
      man->AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "", kFALSE,
	           64, -0.5, 63.5, AliReducedVarManager::kOnlineTriggerFired2, 20, 0.0, 100.0, AliReducedVarManager::kCentVZERO, 0, 0.0, 0.0, AliReducedVarManager::kNothing, triggerNames.Data());
      continue;
    }
    
    if(classStr.Contains("OnlineTriggers_vs_L0TrigInputs")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      TString triggerNames = "";
      for(Int_t i=0; i<64; ++i) {triggerNames += AliReducedVarManager::fgkOfflineTriggerNames[i]; triggerNames+=";";}
      TString trigInputs = "";
      for(Int_t i=0; i<32; ++i) {trigInputs += Form("L0 Input %d", i); trigInputs+=";";}
      man->AddHistogram(classStr.Data(), "OnlineTriggers_L0TriggerInputs", "", kFALSE,
	           64, -0.5, 63.5, AliReducedVarManager::kOnlineTriggerFired2, 32, -0.5, 31.5, AliReducedVarManager::kL0TriggerInput, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 
		   triggerNames.Data(), trigInputs);
    }
    if(classStr.Contains("OnlineTriggers_vs_L1TrigInputs")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      TString triggerNames = "";
      for(Int_t i=0; i<64; ++i) {triggerNames += AliReducedVarManager::fgkOfflineTriggerNames[i]; triggerNames+=";";}
      TString trigInputs = "";
      for(Int_t i=0; i<32; ++i) {trigInputs += Form("L1 Input %d", i); trigInputs+=";";}
      man->AddHistogram(classStr.Data(), "OnlineTriggers_L1TriggerInputs", "", kFALSE,
	           64, -0.5, 63.5, AliReducedVarManager::kOnlineTriggerFired2, 32, -0.5, 31.5, AliReducedVarManager::kL1TriggerInput, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 
		   triggerNames.Data(), trigInputs);
    }
    if(classStr.Contains("OnlineTriggers_vs_L2TrigInputs")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      TString triggerNames = "";
      for(Int_t i=0; i<64; ++i) {triggerNames += AliReducedVarManager::fgkOfflineTriggerNames[i]; triggerNames+=";";}
      TString trigInputs = "";
      for(Int_t i=0; i<16; ++i) {trigInputs += Form("L2 Input %d", i); trigInputs+=";";}
      man->AddHistogram(classStr.Data(), "OnlineTriggers_L2TriggerInputs", "", kFALSE,
	           64, -0.5, 63.5, AliReducedVarManager::kOfflineTriggerFired2, 16, -0.5, 15.5, AliReducedVarManager::kL2TriggerInput, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 
		   triggerNames.Data(), trigInputs);
    }
        
    if(classStr.Contains("ITSclusterMap")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "ITSlayerHit", "Hits in the ITS layers", kFALSE,
                   13, -6.5, 6.5, AliReducedVarManager::kITSlayerHit);
      man->AddHistogram(classStr.Data(), "ITSlayerHit_Phi", "Hits in the ITS layers vs #varphi", kFALSE,
                   100, 0.0, 6.29, AliReducedVarManager::kPhi, 13, -6.5, 6.5, AliReducedVarManager::kITSlayerHit);
      man->AddHistogram(classStr.Data(), "ITSlayerHit_Eta", "Hits in the ITS layers vs #eta", kFALSE,
                   100, -1.0, 1.0, AliReducedVarManager::kEta, 13, -6.5, 6.5, AliReducedVarManager::kITSlayerHit);
      man->AddHistogram(classStr.Data(), "ITSlayerHit_DCAxy", "Hits in the ITS layers vs DCA_{xy}", kFALSE,
                   1000, -0.5, 0.5, AliReducedVarManager::kDcaXY, 13, -6.5, 6.5, AliReducedVarManager::kITSlayerHit);
      man->AddHistogram(classStr.Data(), "ITSlayerHit_DCAz", "Hits in the ITS layers vs DCA_{z}", kFALSE,
                   1800, -1.0, 1.0, AliReducedVarManager::kDcaZ, 13, -6.5, 6.5, AliReducedVarManager::kITSlayerHit);
      continue;
    }  // end of ITSclusterMap histogram definitions
    
    if(classStr.Contains("TPCclusterMap")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "TPCclusterMap", "TPC cluster map", kFALSE,
                   8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
      man->AddHistogram(classStr.Data(), "TPCclusterMap_Phi", "TPC cluster map vs #varphi", kFALSE,
                   360, 0.0, 6.29, AliReducedVarManager::kPhi, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
      man->AddHistogram(classStr.Data(), "TPCclusterMap_Eta", "TPC cluster map vs #eta", kFALSE,
                   100, -1.0, 1.0, AliReducedVarManager::kEta, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
      man->AddHistogram(classStr.Data(), "TPCclusterMap_Pt", "TPC cluster map vs p_{T}", kFALSE,
                   100, 0.0, 10.0, AliReducedVarManager::kPt, 8, -0.5, 7.5, AliReducedVarManager::kTPCclusBitFired);
      continue;
    }  // end of TPCclusterMap histogram definitions
    
    // Calorimeter cluster histograms
    if(classStr.Contains("CaloClusters")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "Energy", "Cluster energy", kFALSE,
	           2, 0.5, 2.5, AliReducedVarManager::kEMCALdetector, 500, 0.0, 20.0, AliReducedVarManager::kEMCALclusterEnergy, 0, 0.0, 0.0, AliReducedVarManager::kNothing, "EMCAL;PHOS");
      man->AddHistogram(classStr.Data(), "Dx", "Cluster dx", kFALSE,
	           2, 0.5, 2.5, AliReducedVarManager::kEMCALdetector, 600, -300.0, 300.0, AliReducedVarManager::kEMCALclusterDx, 0, 0.0, 0.0, AliReducedVarManager::kNothing, "EMCAL;PHOS");
      man->AddHistogram(classStr.Data(), "Dz", "Cluster dz", kFALSE,
	           2, 0.5, 2.5, AliReducedVarManager::kEMCALdetector, 600, -300.0, 300.0, AliReducedVarManager::kEMCALclusterDz, 0, 0.0, 0.0, AliReducedVarManager::kNothing, "EMCAL;PHOS");
      man->AddHistogram(classStr.Data(), "M20", "Cluster M20", kFALSE,
	           2, 0.5, 2.5, AliReducedVarManager::kEMCALdetector, 400, -2.0, 200.0, AliReducedVarManager::kEMCALm20, 0, 0.0, 0.0, AliReducedVarManager::kNothing, "EMCAL;PHOS");
      man->AddHistogram(classStr.Data(), "M02", "Cluster M02", kFALSE,
	           2, 0.5, 2.5, AliReducedVarManager::kEMCALdetector, 400, -2.0, 200.0, AliReducedVarManager::kEMCALm02, 0, 0.0, 0.0, AliReducedVarManager::kNothing, "EMCAL;PHOS");
      man->AddHistogram(classStr.Data(), "Dispersion", "Cluster dispersion", kFALSE,
	           2, 0.5, 2.5, AliReducedVarManager::kEMCALdetector, 400, 0.0, 50.0, AliReducedVarManager::kEMCALdispersion, 0, 0.0, 0.0, AliReducedVarManager::kNothing, "EMCAL;PHOS");
      man->AddHistogram(classStr.Data(), "M20_M02", "", kFALSE, 400, -2.0, 200., AliReducedVarManager::kEMCALm02, 400, -2.0, 200., AliReducedVarManager::kEMCALm20);
      
      man->AddHistogram(classStr.Data(), "CaloClusterEnergy_Run", "Cluster <energy> vs run", kTRUE,
                   runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 2, 0.5, 2.5, AliReducedVarManager::kEMCALdetector, 200, 0.0, 50.0, AliReducedVarManager::kEMCALclusterEnergy, "", "EMCAL;PHOS");
      man->AddHistogram(classStr.Data(), "CaloClusterDx_Run", "Cluster <dx> vs run", kTRUE,
                   runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 2, 0.5, 2.5, AliReducedVarManager::kEMCALdetector, 200, -300.0, 300.0, AliReducedVarManager::kEMCALclusterDx, "", "EMCAL;PHOS");
      man->AddHistogram(classStr.Data(), "CaloClusterDz_Run", "Cluster <dz> vs run", kTRUE,
                   runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 2, 0.5, 2.5, AliReducedVarManager::kEMCALdetector, 200, -300.0, 300.0, AliReducedVarManager::kEMCALclusterDz, "", "EMCAL;PHOS");
    }  // end if "CaloClusters" histograms
    
    TString trkStatusNames = "";
    for(Int_t iflag=0; iflag<AliReducedVarManager::kNTrackingStatus; ++iflag) {
      trkStatusNames += AliReducedVarManager::fgkTrackingStatusNames[iflag];
      trkStatusNames += ";";
    }
    TString trackQualityFlagNames = "EP;#gamma;K^{0}_{S};#Lambda;#bar{#Lambda};kink0;kink1;kink2;";
    trackQualityFlagNames += "pure #gamma;pure K^{0}_{S};pure #Lambda;pure #bar{#Lambda};-kink0;-kink1;-kink2;";
    //trackQualityFlagNames += "bayes e>0.5;bayes #pi>0.5;bayes K>0.5;bayes p>0.5;bayes>0.7;bayes>0.8;bayes>0.9";
    trackQualityFlagNames += "AOD filter bit 0; AOD filter bit 1; AOD filter bit 2; AOD filter bit 3; AOD filter bit 4;";
    trackQualityFlagNames += "AOD filter bit 5; AOD filter bit 6; AOD filter bit 7; AOD filter bit 8; AOD filter bit 9;";
      
    
    if(classStr.Contains("TrackingFlags")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "TrackingFlags", "Tracking flags;;", kFALSE,
	                AliReducedVarManager::kNTrackingStatus, -0.5, AliReducedVarManager::kNTrackingStatus-0.5, AliReducedVarManager::kTrackingFlag, 
			0, 0.0, 0.0, AliReducedVarManager::kNothing, 0, 0.0, 0.0, AliReducedVarManager::kNothing, trkStatusNames.Data());
      continue;
    }
    
    if(classStr.Contains("TrackQualityFlags")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "TrackQualityFlags", "Track quality flags;;", kFALSE,
	                64, -0.5, 63.5, AliReducedVarManager::kTrackQualityFlag, 0, 0.0, 0.0, AliReducedVarManager::kNothing, 
			0, 0.0, 0.0, AliReducedVarManager::kNothing, trackQualityFlagNames.Data());
      continue;
    }
    
    // Track histograms
    if(classStr.Contains("TrackQA")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution", kFALSE, 1000, 0.0, 50.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Eta", "", kFALSE, 1000, -1.5, 1.5, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Phi", "", kFALSE, 1000, 0.0, 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "Charge", "", kFALSE, 2, -1.0, 1.0, AliReducedVarManager::kCharge);
      man->AddHistogram(classStr.Data(), "PtTPC", "p_{T} (TPC) distribution", kFALSE, 1000, 0.0, 50.0, AliReducedVarManager::kPtTPC);
      man->AddHistogram(classStr.Data(), "EtaTPC", "#eta (TPC)", kFALSE, 1000, -1.5, 1.5, AliReducedVarManager::kEtaTPC);
      man->AddHistogram(classStr.Data(), "PhiTPC", "#varphi (TPC)", kFALSE, 1000, 0.0, 6.3, AliReducedVarManager::kPhiTPC);
      man->AddHistogram(classStr.Data(), "Pinner", "p inner param", kFALSE, 1000, 0.0, 50.0, AliReducedVarManager::kPin);
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE, 1000, -10.0, 10.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE, 1000, -10.0, 10.0, AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(), "TrackLength", "Track length", kFALSE, 1000, 0.0, 1000.0, AliReducedVarManager::kTrackLength);
      man->AddHistogram(classStr.Data(),"MassUsedForTracking","",kFALSE,40,0.,4.,AliReducedVarManager::kMassUsedForTracking);
      man->AddHistogram(classStr.Data(),"ITSncls", "ITS nclusters", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSncls);
      man->AddHistogram(classStr.Data(),"Eta_Phi_ITSncls_prof","ITS <nclusters> vs (#eta,#phi)",kTRUE,
                        192, -1.2, 1.2, AliReducedVarManager::kEta, 126, 0.0, 6.3, AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kITSncls);
      man->AddHistogram(classStr.Data(),"ITSsignal", "ITS dE/dx", kFALSE, 400,0.0,1000.0, AliReducedVarManager::kITSsignal);
      man->AddHistogram(classStr.Data(),"Eta_Phi_ITSsignal_prof","ITS <dE/dx> vs (#eta,#phi)",kTRUE,
                        192, -1.2, 1.2, AliReducedVarManager::kEta, 126, 0.0, 6.3, AliReducedVarManager::kPhi, 100, 0.0, 1000, AliReducedVarManager::kITSsignal);
      man->AddHistogram(classStr.Data(),"ITSchi2", "ITS #chi^{2}", kFALSE, 200,0.0,20.0, AliReducedVarManager::kITSchi2);
      man->AddHistogram(classStr.Data(),"Eta_Phi_ITSchi2_prof","ITS <#chi^{2}> vs (#eta,#phi)",kTRUE,
                        192, -1.2, 1.2, AliReducedVarManager::kEta, 126, 0.0, 6.3, AliReducedVarManager::kPhi, 100, 0.0, 1000, AliReducedVarManager::kITSchi2);
      man->AddHistogram(classStr.Data(),"ITSnclsShared", "ITS nclusters shared", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSnclsShared);
      man->AddHistogram(classStr.Data(),"TPCncls","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCncls);
      man->AddHistogram(classStr.Data(),"TPCcrossedRows","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCcrossedRows);
      man->AddHistogram(classStr.Data(),"TPCnclsF","",kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsF);
      man->AddHistogram(classStr.Data(),"TPCnclsShared","",kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsShared);
      man->AddHistogram(classStr.Data(),"TPCchi2","",kFALSE, 200,0.0,10.0,AliReducedVarManager::kTPCchi2);
      man->AddHistogram(classStr.Data(),"Eta_Phi_TPCncls_prof","TPC <nclusters> vs (#eta,#phi)",kTRUE,
                        192, -1.2, 1.2, AliReducedVarManager::kEta, 126, 0.0, 6.3, AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCncls);    
      man->AddHistogram(classStr.Data(),"Eta_Phi_TPCcrossedRows_prof","TPC <n crossed rows> vs (#eta,#phi)",kTRUE,
                        192, -1.2, 1.2, AliReducedVarManager::kEta, 126, 0.0, 6.3, AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCcrossedRows);
      man->AddHistogram(classStr.Data(),"Eta_Phi_TPCnclsF_prof","TPC <nclusters findable> vs (#eta,#phi)",kTRUE,
                        192, -1.2, 1.2, AliReducedVarManager::kEta, 126, 0.0, 6.3, AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCnclsF);
      man->AddHistogram(classStr.Data(),"Eta_Phi_TPCchi2_prof","TPC <#chi^{2}> vs (#eta,#phi)",kTRUE,
                        192, -1.2, 1.2, AliReducedVarManager::kEta, 126, 0.0, 6.3, AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCchi2);
      man->AddHistogram(classStr.Data(),"TPCsignal_Pin","TPC dE/dx vs. inner param P",kFALSE,
                        400,0.0,4.0,AliReducedVarManager::kPin,500,-0.5,499.5,AliReducedVarManager::kTPCsignal);
      man->AddHistogram(classStr.Data(),"TPCsignalN","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCsignalN);
      man->AddHistogram(classStr.Data(),"Eta_Phi_TPCsignalN_prof","TPC <nclusters pid> vs (#eta,#phi)",kTRUE,
                        192, -1.2, 1.2, AliReducedVarManager::kEta, 126, 0.0, 6.3, AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCsignalN);    
      man->AddHistogram(classStr.Data(),"TPCnsigElectron_Pin","TPC N_{#sigma} electron vs. inner param P",kFALSE,
                        200,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(),"TPCnsigElectron_Eta","TPC N_{#sigma} electron vs. #eta",kFALSE,
                        100,-1.0,1.0,AliReducedVarManager::kEta,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(),"TPCnsigPion_Pin","TPC N_{#sigma} pion vs. inner param P",kFALSE,
                        200,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion);
      man->AddHistogram(classStr.Data(),"TPCnsigKaon_Pin","TPC N_{#sigma} kaon vs. inner param P",kFALSE,
                        200,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kKaon);
      man->AddHistogram(classStr.Data(),"TPCnsigProton_Pin","TPC N_{#sigma} proton vs. inner param P",kFALSE,
                        200,0.0,10.0,AliReducedVarManager::kPin,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton);
      man->AddHistogram(classStr.Data(),"TPCnsigElectron_Run","TPC N_{#sigma} electron vs. run",kTRUE,
                        runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(),"TPCnsigPion_Run","TPC N_{#sigma} pion vs. run",kTRUE,
                        runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion);
      man->AddHistogram(classStr.Data(),"TPCnsigKaon_Run","TPC N_{#sigma} kaon vs. run",kTRUE,
                        runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kKaon);
      man->AddHistogram(classStr.Data(),"TPCnsigProton_Run","TPC N_{#sigma} proton vs. run",kTRUE,
                        runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton);
      man->AddHistogram(classStr.Data(),"TOFdeltaBC","TOF #delta BC", kFALSE, 7000, -3500.,3500.,AliReducedVarManager::kTOFdeltaBC);
      man->AddHistogram(classStr.Data(),"TOFtime","TOF time", kFALSE, 500, 0., 1.0e+5,AliReducedVarManager::kTOFtime);
      man->AddHistogram(classStr.Data(),"TOFdx","TOF dx", kFALSE, 500, -5., 5., AliReducedVarManager::kTOFdx);
      man->AddHistogram(classStr.Data(),"TOFdz","TOF dz", kFALSE, 500, -5., 5., AliReducedVarManager::kTOFdz);
      man->AddHistogram(classStr.Data(),"TOFmismatchProbab","TOF mismatch probability", kFALSE, 120, -0.1, 1.1, AliReducedVarManager::kTOFmismatchProbability);
      man->AddHistogram(classStr.Data(),"TOFchi2","TOF #chi^{2}", kFALSE, 500, -1., 10., AliReducedVarManager::kTOFchi2);
      man->AddHistogram(classStr.Data(),"TOFbeta_P","TOF #beta vs P",kFALSE,
                        200,0.0,20.0,AliReducedVarManager::kP, 220,0.0,1.1,AliReducedVarManager::kTOFbeta);
      man->AddHistogram(classStr.Data(),"TOFnsigElectron_P","TOF N_{#sigma} electron vs. P",kFALSE,
                        200,0.0,20.0,AliReducedVarManager::kP, 100,-5.0,5.0,AliReducedVarManager::kTOFnSig+AliReducedVarManager::kElectron);
      man->AddHistogram(classStr.Data(),"TOFnsigPion_P","TOF N_{#sigma} pion vs. P",kFALSE,
                        200,0.0,20.0,AliReducedVarManager::kP, 100,-5.0,5.0,AliReducedVarManager::kTOFnSig+AliReducedVarManager::kPion);
      man->AddHistogram(classStr.Data(),"TOFnsigKaon_P","TOF N_{#sigma} kaon vs. P",kFALSE,
                        200,0.0,20.0,AliReducedVarManager::kP, 100,-5.0,5.0,AliReducedVarManager::kTOFnSig+AliReducedVarManager::kKaon);
      man->AddHistogram(classStr.Data(),"TOFnsigProton_P","TOF N_{#sigma} proton vs. P",kFALSE,
                        200,0.0,20.0,AliReducedVarManager::kP, 100,-5.0,5.0,AliReducedVarManager::kTOFnSig+AliReducedVarManager::kProton);
      man->AddHistogram(classStr.Data(),"TRDntracklets","TRD ntracklets",kFALSE,
                        7,-0.5,6.5,AliReducedVarManager::kTRDntracklets);
      man->AddHistogram(classStr.Data(),"TRDntrackletsPID","TRD ntracklets PID",kFALSE,
                        7,-0.5,6.5,AliReducedVarManager::kTRDntrackletsPID);
      man->AddHistogram(classStr.Data(),"TRDprobabElectronLQ1D","TRD electron probability LQ1D",kFALSE,
                        500,0.0,1.0,AliReducedVarManager::kTRDpidProbabilitiesLQ1D);
      man->AddHistogram(classStr.Data(),"TRDprobabPionLQ1D","TRD pion probability LQ1D",kFALSE,
                        500,0.0,1.0,AliReducedVarManager::kTRDpidProbabilitiesLQ1D+1);
      man->AddHistogram(classStr.Data(),"TRDprobabElectronLQ2D","TRD electron probability LQ2D",kFALSE,
                        500,0.0,1.0,AliReducedVarManager::kTRDpidProbabilitiesLQ2D);
      man->AddHistogram(classStr.Data(),"TRDprobabPionLQ2D","TRD pion probability LQ2D",kFALSE,
                        500,0.0,1.0,AliReducedVarManager::kTRDpidProbabilitiesLQ2D+1);
      man->AddHistogram(classStr.Data(),"Eta_Phi_TRDntracklets_prof","TRD <ntracklets> vs (#eta,#phi)",kTRUE,
                        192, -1.2, 1.2, AliReducedVarManager::kEta, 126, 0.0, 6.3, AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kTRDntracklets);   
      man->AddHistogram(classStr.Data(),"Eta_Phi_TRDntrackletsPID_prof","TRD <ntracklets PID> vs (#eta,#phi)",kTRUE,
                        192, -1.2, 1.2, AliReducedVarManager::kEta, 126, 0.0, 6.3, AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kTRDntrackletsPID);
      man->AddHistogram(classStr.Data(), "MatchedCaloClusterId", "EMCAL/PHOS matched cluster id", kFALSE,
                        5000, 0.0, 5000.0, AliReducedVarManager::kEMCALmatchedClusterId);
      man->AddHistogram(classStr.Data(), "MatchedEnergy", "Energy from the calorimeter matched cluster", kFALSE,
                        200, 0.0, 50.0, AliReducedVarManager::kEMCALmatchedEnergy);
      man->AddHistogram(classStr.Data(), "MatchedEOverP", "E/P from the calorimeter matched cluster", kFALSE,
                        300, 0.0, 1.5, AliReducedVarManager::kEMCALmatchedEOverP);
      man->AddHistogram(classStr.Data(), "MatchedEOverP_P", "E/P from the calorimeter matched cluster vs P", kFALSE,
                        10, 0.0, 10.0, AliReducedVarManager::kPin, 300, 0.0, 1.5, AliReducedVarManager::kEMCALmatchedEOverP);
      man->AddHistogram(classStr.Data(),"MatchedEnergy_Run","Calorimeter matched energy vs run",kTRUE,
                        runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 100, 0.0, 50.0, AliReducedVarManager::kEMCALmatchedEnergy);
      man->AddHistogram(classStr.Data(),"EMCALmatchedEOverP_Run","Calorimeter matched E/P vs run",kTRUE,
                        runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo, 100, 0.0, 2.0, AliReducedVarManager::kEMCALmatchedEOverP);
      man->AddHistogram(classStr.Data(),"EMCALmatchedEnergy_P","Calorimeter matched energy vs momentum",kFALSE,
                        100, 0.0, 10.0, AliReducedVarManager::kP, 100, 0.0, 10.0, AliReducedVarManager::kEMCALmatchedEnergy);
      man->AddHistogram(classStr.Data(),"EMCALmatchedEnergy_CentVZERO","Calorimeter matched energy vs centrality VZERO",kTRUE,
                        20, 0.0, 100.0, AliReducedVarManager::kCentVZERO, 100, 0.0, 10.0, AliReducedVarManager::kEMCALmatchedEnergy);
      man->AddHistogram(classStr.Data(),"Eta_Phi_MatchedEnergy_prof","EMCAL matched cluster <energy> vs (#eta,#phi)",kTRUE,
                        192, -1.2, 1.2, AliReducedVarManager::kEta, 126, 0.0, 6.3, AliReducedVarManager::kPhi, 100, 0.0, 10.0, AliReducedVarManager::kEMCALmatchedEnergy);
      
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
      if(classStr.Contains("QA")) {
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
      }   // end if "QA"
      
      continue;
    }   // end if for Pair classes of histograms
  }  // end loop over histogram classes
}
