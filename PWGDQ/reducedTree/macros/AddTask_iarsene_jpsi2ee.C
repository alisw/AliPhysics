//void Setup(AliReducedAnalysisJpsi2ee* processor, TString prod="LHC10h");
//AliHistogramManager* SetupHistogramManager(TString prod="LHC10h");
//void DefineHistograms(AliHistogramManager* man, TString prod="LHC10h");

//__________________________________________________________________________________________
AliAnalysisTask* AddTask_iarsene_jpsi2ee(Bool_t isAliRoot=kTRUE, Int_t runMode=1, TString prod="LHC10h"){    
   //
   // isAliRoot=kTRUE for ESD/AOD analysis in AliROOT, kFALSE for root analysis on reduced trees
   // runMode=1 (AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents)
   //               =2 (AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree)
   //
   //get the current analysis manager

  printf("INFO on AddTask_iarsene_jpsi2ee(): (isAliRoot, runMode) :: (%d,%d) \n", isAliRoot, runMode);

  AliReducedAnalysisJpsi2ee* jpsi2eeAnalysis = new AliReducedAnalysisJpsi2ee("Jpsi2eeAnalysis","Jpsi->ee analysis");
  jpsi2eeAnalysis->Init();
  Setup(jpsi2eeAnalysis, prod);
  jpsi2eeAnalysis->SetRunEventMixing(kTRUE);
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
                                                                  AliAnalysisManager::kOutputContainer, "dstAnalysisHistograms.root");
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
   
  TFile* corrFile = TFile::Open("/home/iarsene/work/ALICE/treeAnalysis/newdst/development/test/CorrMaps_CentralityNew.root");
  AliReducedVarManager::SetTPCelectronCorrectionMaps((TH2F*)corrFile->Get("CorrMapCentroid"), (TH2F*)corrFile->Get("CorrMapWidth"), AliReducedVarManager::kEta, AliReducedVarManager::kCentVZERO);
  corrFile->Close();
  
  TFile* grpFile = TFile::Open("/home/iarsene/work/ALICE/treeAnalysis/newdst/development/test/grpData_LHC15o.root");
  AliReducedVarManager::SetLHCDataInfo((TH1F*)grpFile->Get("lumiTotal"), (TH1F*)grpFile->Get("intTotal0"), 
                                       (TH1F*)grpFile->Get("intTotal1"), (TH1I*)grpFile->Get("fillNumber"));
  AliReducedVarManager::SetGRPDataInfo((TH1I*)grpFile->Get("dipolePolarity"), (TH1I*)grpFile->Get("l3Polarity"), 
                                       (TH1I*)grpFile->Get("timeStart"), (TH1I*)grpFile->Get("timeStop"));
  grpFile->Close(); 
  
  // Set event cuts
  AliReducedEventCut* evCut1 = new AliReducedEventCut("Centrality","Centrality selection");
  evCut1->AddCut(AliReducedVarManager::kCentVZERO, 70., 90.);
  evCut1->AddCut(AliReducedVarManager::kVtxZ, -10.0, 10.0);
  evCut1->AddCut(AliReducedVarManager::kIsPhysicsSelection, 0.1, 2.);   // request physics selection
  processor->AddEventCut(evCut1);
  
  // Set track cuts
  AliReducedTrackCut* trackCut1 = new AliReducedTrackCut("HighPt","High pt cuts");
  trackCut1->AddCut(AliReducedVarManager::kP, 3.0,100.0);
  trackCut1->AddCut(AliReducedVarManager::kDcaXY, -0.2,0.2);
  trackCut1->AddCut(AliReducedVarManager::kDcaZ, -1.0,1.0);
  trackCut1->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  //trackCut1->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut1->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 2.0, 1.0e+30);
  trackCut1->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -1.5, 4.0);
  trackCut1->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -1.5, 4.0, kFALSE, AliReducedVarManager::kEta, -0.9, -0.6);
  trackCut1->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -1.5, 4.0, kFALSE, AliReducedVarManager::kEta, 0.6, 0.9);
  trackCut1->SetRejectKinks();
  //trackCut1->SetRequestITSrefit();
  trackCut1->SetRequestTPCrefit();
  //trackCut1->SetRequestSPDany();
  trackCut1->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  processor->AddTrackCut(trackCut1);  
  AliReducedTrackCut* trackCut2 = new AliReducedTrackCut("P10","P>1.0 GeV/c, standard cuts");
  trackCut2->AddCut(AliReducedVarManager::kP, 1.0,100.0);
  trackCut2->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  trackCut2->AddCut(AliReducedVarManager::kDcaZ, -1.0,1.0);
  trackCut2->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut2->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut2->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
  trackCut2->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 4.0);
  trackCut2->SetRejectKinks();
  trackCut2->SetRequestITSrefit();
  trackCut2->SetRequestTPCrefit();
  //trackCut2->SetRequestSPDany();
  //processor->AddTrackCut(trackCut2);  
  AliReducedTrackCut* trackCut3 = new AliReducedTrackCut("P10_weakHadrRej","P>1.0 GeV/c, weak hadr rejection");
  trackCut3->AddCut(AliReducedVarManager::kP, 1.0,100.0);
  trackCut3->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  trackCut3->AddCut(AliReducedVarManager::kDcaZ, -1.0,1.0);
  trackCut3->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut3->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.0, 1.0e+30);
  trackCut3->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.0, 1.0e+30);
  trackCut3->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 4.0);
  trackCut3->SetRejectKinks();
  trackCut3->SetRequestITSrefit();
  trackCut3->SetRequestTPCrefit();
  //trackCut3->SetRequestSPDany();
  //processor->AddTrackCut(trackCut3);  
  AliReducedTrackCut* trackCut4 = new AliReducedTrackCut("P10_weakerHadrRej","P>1.0 GeV/c, weaker hadr rejection");
  trackCut4->AddCut(AliReducedVarManager::kP, 1.0,100.0);
  trackCut4->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  trackCut4->AddCut(AliReducedVarManager::kDcaZ, -1.0,1.0);
  trackCut4->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut4->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 2.5, 1.0e+30);
  trackCut4->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.0, 1.0e+30);
  trackCut4->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 4.0);
  trackCut4->SetRejectKinks();
  trackCut4->SetRequestITSrefit();
  trackCut4->SetRequestTPCrefit();
  //trackCut4->SetRequestSPDany();
  //processor->AddTrackCut(trackCut4);  
  AliReducedTrackCut* trackCut5 = new AliReducedTrackCut("P10_weakestHadrRej","P>1.0 GeV/c, weakest hadr rejection");
  trackCut5->AddCut(AliReducedVarManager::kP, 1.0,100.0);
  trackCut5->AddCut(AliReducedVarManager::kDcaXY, -1.0,1.0);
  trackCut5->AddCut(AliReducedVarManager::kDcaZ, -1.0,1.0);
  trackCut5->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut5->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 2.0, 1.0e+30);
  trackCut5->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.0, 1.0e+30);
  trackCut5->AddCut(AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron, -3.0, 4.0);
  trackCut5->SetRejectKinks();
  trackCut5->SetRequestITSrefit();
  trackCut5->SetRequestTPCrefit();
  //trackCut5->SetRequestSPDany();
  //processor->AddTrackCut(trackCut5);  
  
  // set track prefilter cuts
  AliReducedTrackCut* prefTrackCut1 = new AliReducedTrackCut("prefCutPt09","prefilter P selection");
  prefTrackCut1->AddCut(AliReducedVarManager::kP, 0.9,100.0);
  prefTrackCut1->SetRequestTPCrefit();
  processor->AddPrefilterTrackCut(prefTrackCut1);  
  
  // set pair prefilter cuts
  AliReducedVarCut* prefPairCut = new AliReducedVarCut("prefCutM50MeV","prefilter pair cuts");
  prefPairCut->AddCut(AliReducedVarManager::kMass, 0.0, 0.05, kTRUE);
  processor->AddPrefilterPairCut(prefPairCut);
  
  // Set pair cuts
  AliReducedTrackCut* pairCut1 = new AliReducedTrackCut("Ptpair","Pt pair selection");
  pairCut1->AddCut(AliReducedVarManager::kPt, 0.0,100.0);
  //pairCut1->AddCut(AliReducedVarManager::kEta, -1.5,1.5);
  processor->AddPairCut(pairCut1);
  
  //SetupHistogramManager(processor->GetHistogramManager(), prod);
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
   handler->SetEventVariables(AliReducedVarManager::kCentVZERO,AliReducedVarManager::kVtxZ,AliReducedVarManager::kVZERORP+6+1);
   Float_t centLims[11] = {0.0, 5.0, 10.0, 20.0, 30.0, 40., 50., 60., 70.,80.,90.};
   //Float_t zLims[11] = {-10.,-8.,-6.,-4.,-2.,0.,2.,4.,6.,8.,10.};
   Float_t zLims[2] = {-10.,10.};
   Float_t epLims[2] = {-10000.,10000.};
   handler->SetCentralityLimits(11,centLims);
   handler->SetEventVertexLimits(2,zLims);
   handler->SetEventPlaneLimits(2,epLims);
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
  for(Int_t i=0; i<task->GetNTrackCuts(); ++i) {
    TString cutName = task->GetTrackCutName(i);
    histClasses += Form("Track_%s;", cutName.Data());
    //histClasses += Form("Track_NoPrefilter_%s;", cutName.Data());
    //histClasses += Form("PairPrefilterSEPP_%s;", cutName.Data());
    //histClasses += Form("PairPrefilterSEPM_%s;", cutName.Data());
    //histClasses += Form("PairPrefilterSEMM_%s;", cutName.Data());
    //histClasses += Form("PairSEPP_%s;", cutName.Data());
    //histClasses += Form("PairSEPM_%s;", cutName.Data());
    //histClasses += Form("PairSEMM_%s;", cutName.Data());
    histClasses += Form("PairSEPP_%s;PairSEPM_%s;PairSEMM_%s;", cutName.Data(), cutName.Data(), cutName.Data());
    //histClasses += Form("PairMEPP_%s;", cutName.Data());
    //histClasses += Form("PairMEPM_%s;", cutName.Data());
    //histClasses += Form("PairMEMM_%s;", cutName.Data());
    histClasses += Form("PairMEPP_%s;PairMEPM_%s;PairMEMM_%s;", cutName.Data(), cutName.Data(), cutName.Data());
  }
  
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
  
  // pp at 13 TeV
  if(!prod.CompareTo("LHC16l")) {
     runNBins = 1140;
     runHistRange[0] = 258880.;
     runHistRange[1] = 260020.;
  }
  
  const Int_t kNBCBins = 3600;
  Double_t bcHistRange[2] = {-0.5,3599.5};
  
  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");
  
  cout << "Histogram classes included in the Histogram Manager" << endl;
  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    
    // Event wise histograms
    if(classStr.Contains("Event_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(),"RunNo","Run numbers",kFALSE, runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo);
      man->AddHistogram(classStr.Data(),"VtxX","Vtx X",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxX_TimeFromSOR_prof","<Vtx X> vs time from SOR",kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxY","Vtx Y",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxY_TimeFromSOR_prof","<Vtx Y> vs time from SOR",kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxZ","Vtx Z",kFALSE,300,-15.,15.,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"VtxZ_TimeFromSOR_prof","<Vtx Z> vs time from SOR",kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300,-15.0,15.0,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentVZERO);
      man->AddHistogram(classStr.Data(),"CentVZERO_VtxZ","Centrality(VZERO) vs vtxZ",kFALSE, 50, 0.0, 100.0, AliReducedVarManager::kCentVZERO,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentVZERO_TimeFromSOR","Centrality(VZERO) vs time from SOR",kFALSE,
                        225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 50, 0.0, 100.0, AliReducedVarManager::kCentVZERO);
      man->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentSPD);
      man->AddHistogram(classStr.Data(),"CentSPD_VtxZ","Centrality(SPD) vs vtxZ",kFALSE, 50, 0.0, 100.0, AliReducedVarManager::kCentSPD,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC)",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentTPC);
      man->AddHistogram(classStr.Data(),"CentTPC_VtxZ","Centrality(TPC) vs vtxZ",kFALSE, 100, 0.0, 100.0, AliReducedVarManager::kCentTPC,
         50, -10., 10., AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"CentQuality","Centrality quality",kFALSE, 500, 0., 500., AliReducedVarManager::kCentQuality);
      man->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event",kFALSE,500,0.,20000.,AliReducedVarManager::kNtracksTotal);
      man->AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"NTracksSelected_TimeFromSOR","Averaged number of selected tracks per event vs time from SOR",kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"NTracksSelected_CentVZERO_TimeFromSOR","Averaged number of selected tracks per event per centrality vs time from SOR",kTRUE, 20, 0.0, 100.0, AliReducedVarManager::kCentVZERO, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksSelected);
      man->AddHistogram(classStr.Data(),"DeltaVtxZ","Z_{global}-Z_{TPC}",kFALSE,300,-1.,1.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"DeltaVtxZ_TimeFromSOR_prof","Z_{global}-Z_{TPC} vs time from SOR",kTRUE,225, 0.0, 450.,  AliReducedVarManager::kTimeRelativeSOR, 300,-1.,1.,AliReducedVarManager::kDeltaVtxZ);
      man->AddHistogram(classStr.Data(),"NVtxTPCContributors","",kFALSE,500,0.,10000.,AliReducedVarManager::kNVtxTPCContributors);
      man->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0", kFALSE, 200, 0., 5000., AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(),"SPDntracklets_TimeFromSOR_prof", "SPD <#tracklets> in |#eta|<1.0 vs time from SOR", kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, 0., 5000., AliReducedVarManager::kSPDntracklets);
      for(Int_t il=0; il<2; ++il) {
        man->AddHistogram(classStr.Data(), Form("SPDfiredChips_layer%d",il+1), Form("SPD fired chips in layer %d",il+1), 
			  kFALSE, 200, 0., 1000., AliReducedVarManager::kSPDFiredChips+il);
        man->AddHistogram(classStr.Data(), Form("SPDfiredChips_layer%d_TimeFromSOR_prof",il+1), Form("SPD <#fired chips> in layer %d vs time from SOR",il+1), 
                          kTRUE, 225, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR, 200, 0., 1000., AliReducedVarManager::kSPDFiredChips+il);
      }
      for(Int_t il=0; il<6; ++il) {
        man->AddHistogram(classStr.Data(), Form("ITSclusters_layer%d",il+1), Form("ITS clusters in layer %d",il+1), 
			  kFALSE, 300, 0., 15000., AliReducedVarManager::kITSnClusters+il);
        man->AddHistogram(classStr.Data(), Form("ITSclusters_layer%d_TimeFromSOR_prof",il+1), Form("ITS <clusters> in layer %d vs time from SOR",il+1), 
                          kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 300, 0., 15000., AliReducedVarManager::kITSnClusters+il);
      }
      man->AddHistogram(classStr.Data(), "SPDnSingleClusters", "SPD single clusters", 
			kFALSE, 200, 0., 10000., AliReducedVarManager::kSPDnSingleClusters);	
      man->AddHistogram(classStr.Data(), "SPDnSingleClusters_TimeFromSOR_prof", "SPD single <#clusters> vs time from SOR", 
                        kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, 0., 10000., AliReducedVarManager::kSPDnSingleClusters);	
      man->AddHistogram(classStr.Data(),"VZEROmult", "", kFALSE, 300, 0.0, 30000., AliReducedVarManager::kVZEROTotalMult);
      man->AddHistogram(classStr.Data(),"VZEROmult_VtxContributors", "", kFALSE, 
		   200, 0.0, 10000., AliReducedVarManager::kVZEROTotalMult, 200, 0.0, 1000., AliReducedVarManager::kNVtxContributors);
      man->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsPhysicsSelection, 0,0.0,0.0,AliReducedVarManager::kNothing, 0,0.0,0.0,AliReducedVarManager::kNothing, "off;on");
      man->AddHistogram(classStr.Data(),"IsPhysicsSelection_TimeFromSOR_prof","Fraction of physics selection selected",kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 2,-0.5,1.5,AliReducedVarManager::kIsPhysicsSelection);
      man->AddHistogram(classStr.Data(),"NPosTracksAnalyzed","#positrons per event",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksPosAnalyzed);
      man->AddHistogram(classStr.Data(),"NNegTracksAnalyzed","#electrons per event",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksNegAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalTracksAnalyzed","#leg candidates per event",kFALSE,100,0.,100.,AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalPairsAnalyzed","#dielectron candidates per event",kFALSE,100,0.,100.,AliReducedVarManager::kNpairsSelected);
      man->AddHistogram(classStr.Data(),"NPosTracksAnalyzed_TimeFromSOR_prof","<#positrons> per event vs time from SOR",kTRUE,225, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR, 100,0.,100.,AliReducedVarManager::kNtracksPosAnalyzed);
      man->AddHistogram(classStr.Data(),"NNegTracksAnalyzed_TimeFromSOR_prof","<#electrons> per event vs time from SOR",kTRUE,225, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR,100,0.,100.,AliReducedVarManager::kNtracksNegAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalTracksAnalyzed_TimeFromSOR_prof","<#leg candidates> per event vs time from SOR",kTRUE,225, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR,100,0.,100.,AliReducedVarManager::kNtracksAnalyzed);
      man->AddHistogram(classStr.Data(),"NTotalPairsAnalyzed_TimeFromSOR_prof","<#dielectron> candidates per event vs time from SOR",kTRUE,225, 0.0, 450., AliReducedVarManager::kTimeRelativeSOR,100,0.,100.,AliReducedVarManager::kNpairsSelected);
      continue;
    }  // end if className contains "Event"    
    
    // Track histograms
    if(classStr.Contains("Track_")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution", kFALSE, 1000, 0.0, 50.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Pt_TimeFromSOR", "<p_{T}> vs time from SOR", kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, 0.0, 50.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Eta", "", kFALSE, 1000, -1.5, 1.5, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Eta_TimeFromSOR", "<#eta> vs time from SOR", kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, -1.5, 1.5, AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Phi", "", kFALSE, 1000, 0.0, 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "Phi_TimeFromSOR", "<#varphi> vs time from SOR", kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 1000, 0.0, 6.3, AliReducedVarManager::kPhi);
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE, 200, -5.0, 5.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAxy_TimeFromSOR", "<DCAxy> vs time from SOR", kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, -5.0, 5.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE, 200, -5.0, 5.0, AliReducedVarManager::kDcaZ);
      man->AddHistogram(classStr.Data(), "DCAz_TimeFromSOR", "<DCAz> vs time from SOR", kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200, -5.0, 5.0, AliReducedVarManager::kDcaZ);
      
        man->AddHistogram(classStr.Data(),"ITSncls", "ITS nclusters", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"ITSncls_TimeFromSOR", "<ITS nclusters> vs time from SOR", kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 7,-0.5,6.5,AliReducedVarManager::kITSncls);
        man->AddHistogram(classStr.Data(),"Eta_Phi_ITSncls_prof","ITS <nclusters> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kITSncls);
	man->AddHistogram(classStr.Data(),"ITSchi2", "ITS #chi^{2}", kFALSE, 200,0.0,20.0, AliReducedVarManager::kITSchi2);
        man->AddHistogram(classStr.Data(),"ITSchi2_TimeFromSOR", "<ITS #chi^{2}> vs time from SOR", kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200,0.0,20.0, AliReducedVarManager::kITSchi2);
	man->AddHistogram(classStr.Data(),"Eta_Phi_ITSchi2_prof","ITS <#chi^{2}> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 100, 0.0, 1000, AliReducedVarManager::kITSchi2);
      
        man->AddHistogram(classStr.Data(),"TPCncls","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCncls);
        man->AddHistogram(classStr.Data(),"TPCncls_TimeFromSOR","<TPC #cls> vs time from SOR",kTRUE,225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 160,-0.5,159.5,AliReducedVarManager::kTPCncls);
	man->AddHistogram(classStr.Data(),"TPCcrossedRows","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCcrossedRows);
	man->AddHistogram(classStr.Data(),"TPCnclsShared","",kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsShared);
	man->AddHistogram(classStr.Data(),"TPCchi2","",kFALSE, 200,0.0,10.0,AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCchi2_TimeFromSOR","TPC <#chi^{2}> vs time from SOR",kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 200,0.0,10.0,AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"Eta_Phi_TPCncls_prof","TPC <nclusters> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCncls);    
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCcrossedRows_prof","TPC <n crossed rows> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCcrossedRows);
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCnclsF_prof","TPC <nclusters findable> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCnclsF);
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCchi2_prof","TPC <#chi^{2}> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCsignal_Pin","TPC dE/dx vs. inner param P",kFALSE,
                     300,0.0,30.0,AliReducedVarManager::kPin,200,-0.5,199.5,AliReducedVarManager::kTPCsignal);
	man->AddHistogram(classStr.Data(),"TPCsignalN","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCsignalN);
        man->AddHistogram(classStr.Data(),"TPCsignalN_TimeFromSOR","TPC <#cls pid> vs time from SOR",kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 160,-0.5,159.5,AliReducedVarManager::kTPCsignalN);
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCsignalN_prof","TPC <nclusters pid> vs (#eta,#phi)",kTRUE,
                     36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCsignalN);    
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_TimeFromSOR","TPC N_{#sigma} electron vs. time from SOR",kFALSE,
                          225, 0., 450.,AliReducedVarManager::kTimeRelativeSOR,80,-4.0,4.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_TimeFromSOR_prof","<TPC N_{#sigma} electron> vs. time from SOR",kTRUE,
                          225, 0., 450.,AliReducedVarManager::kTimeRelativeSOR,80,-4.0,4.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Pin","TPC N_{#sigma} electron vs. inner param P",kFALSE,
                     100,0.0,10.0,AliReducedVarManager::kPin,80,-4.0,4.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Pin_TimeFromSOR","TPC N_{#sigma} electron vs. inner param P vs time from SOR",kFALSE,
                          100,0.0,10.0,AliReducedVarManager::kPin,80,-4.0,4.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, 15, 0., 450., AliReducedVarManager::kTimeRelativeSOR);
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Eta","TPC N_{#sigma} electron vs. #eta",kFALSE,
                          36,-0.9,0.9,AliReducedVarManager::kEta,80,-4.0,4.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Eta_TimeFromSOR","TPC N_{#sigma} electron vs. #eta vs time from SOR",kFALSE,
                          36,-0.9,0.9,AliReducedVarManager::kEta,80,-4.0,4.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, 15, 0., 450., AliReducedVarManager::kTimeRelativeSOR);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Pin","Corrected TPC N_{#sigma} electron vs. inner param P",kFALSE,
                          100,0.0,10.0,AliReducedVarManager::kPin,80,-4.0,4.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_Eta","TPC N_{#sigma} electron vs. #eta",kFALSE,
                          36,-0.9,0.9,AliReducedVarManager::kEta,80,-4.0,4.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron);
        man->AddHistogram(classStr.Data(),"TPCnsigElectronCorrected_TimeFromSOR_prof","<corrected TPC N_{#sigma} electron> vs. time from SOR",kTRUE,
                          225, 0., 450.,AliReducedVarManager::kTimeRelativeSOR,80,-4.0,4.0,AliReducedVarManager::kTPCnSigCorrected+AliReducedVarManager::kElectron);
	man->AddHistogram(classStr.Data(),"TPCnsigElectron_Run","TPC N_{#sigma} electron vs. run",kTRUE,
                     runNBins, runHistRange[0], runHistRange[1], AliReducedVarManager::kRunNo,100,-5.0,5.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
      
        man->AddHistogram(classStr.Data(),"TOFbeta_P","TOF #beta vs P",kFALSE,
                     200,0.0,20.0,AliReducedVarManager::kP, 220,0.0,1.1,AliReducedVarManager::kTOFbeta);
      
      continue;
    }  // end if "TrackQA"
        
    const Int_t kNMassBins = 150;
    Double_t massBins[kNMassBins];
    for(Int_t i=0; i<kNMassBins; ++i) massBins[i] = 0.0 + i*0.1; 
    
    //const Int_t kNPtBins = 55;
    const Int_t kNPtBins = 10;
    Double_t ptBins[kNPtBins] = {0., 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20., 30., 50.};
    //for(Int_t i=0; i<=50; ++i) ptBins[i] = 0.0+ i*0.01;
    //ptBins[51] = 1.5; ptBins[52] = 4.5; ptBins[53] = 10.0; ptBins[54] = 20.0;
    
    
    const Int_t kNCentBins = 11;
    //Double_t centBins[kNCentBins] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.,20.,25.,30.,35.,40.,50.,60.,70.,90.};
    Double_t centBins[kNCentBins] = {0.0, 5.0, 10.0, 20.0, 30.0, 40., 50., 60., 70.,80.,90.};
    
    const Int_t kNVtxBins = 2;
    Double_t vtxBins[kNVtxBins] = {-10.,10.};
    
    const Int_t kNEPbins = 2;
    Double_t epBins[kNEPbins] = {-0.5*TMath::Pi()*10000.,0.5*TMath::Pi()*10000.};
    
    Int_t vars[5] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, 
       AliReducedVarManager::kCentVZERO, AliReducedVarManager::kVtxZ, AliReducedVarManager::kVZERORP+6+1};
       
    //TArrayD pairHistBinLimits[5] = {TArrayD(kNMassBins,massBins), TArrayD(kNPtBins,ptBins), TArrayD(kNCentBins,centBins), TArrayD(kNVtxBins,vtxBins), TArrayD(kNEPbins,epBins)};
    TArrayD pairHistBinLimits[5];
    pairHistBinLimits[0] = TArrayD(kNMassBins,massBins);
    pairHistBinLimits[1] = TArrayD(kNPtBins,ptBins);
    pairHistBinLimits[2] = TArrayD(kNCentBins,centBins);
    pairHistBinLimits[3] = TArrayD(kNVtxBins,vtxBins);
    pairHistBinLimits[4] = TArrayD(kNEPbins,epBins);
       
    // Histograms for pairs
    if(classStr.Contains("Pair")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "PairInvMass", "Differential pair inv. mass", 5, vars, pairHistBinLimits);
      if(classStr.Contains("PairSE") || classStr.Contains("PairPrefilterSE")) {
        man->AddHistogram(classStr.Data(), "PairType", "Pair type", kFALSE, 4, -0.5, 3.5, AliReducedVarManager::kPairType);
	man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE, kNMassBins, 0.0, 15.0, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "Mass_TimeFromSOR", "Invariant mass vs time from SOR", kFALSE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "Mass_TimeFromSOR_prof", "<Invariant mass> vs time from SOR", kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "Mass_Cent_TimeFromSOR_prof", "<Invariant mass> vs centrality and time from SOR", kTRUE, 225, 0., 450., AliReducedVarManager::kTimeRelativeSOR, 18, 0.0, 90., AliReducedVarManager::kCentVZERO, kNMassBins, 0.0, 5.0, AliReducedVarManager::kMass);
        man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 1000, 0.0, 10.0, AliReducedVarManager::kPt);
	man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kRap);
	man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 315, 0.0, 6.3, AliReducedVarManager::kPhi);
      }   // end if "QA"
      
      continue;
    }   // end if for Pair classes of histograms
  }  // end loop over histogram classes
}
