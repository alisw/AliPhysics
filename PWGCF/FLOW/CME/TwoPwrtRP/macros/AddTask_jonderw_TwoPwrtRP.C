void Setup(AliAnalysisTaskTwoPwrtRP* processor, Int_t period=-1);
AliHistogramManager* SetupHistogramManager(Int_t period=-1);
void DefineHistograms(AliHistogramManager* man, Int_t period=-1, AliAnalysisTaskSE* analysis);
//void GetFilterTPCFlowCuts(AliCMEAnalysisCuts *trackCuts);
//void GetFilterGlobalFlowCuts(AliCMEAnalysisCuts *trackCuts);
void GetFilterTrackCuts(AliCMEAnalysisCuts *trackCuts);
void GetFilterGlobalTrackCuts2010(AliCMEAnalysisCuts *trackCuts);
void GetEventCuts(AliCMEAnalysisCuts *eventCuts);


//__________________________________________________________________________________________
AliAnalysisTaskSE* AddTask_jonderw_TwoPwrtRP(Bool_t isAliRoot=kTRUE){    // called with isAliRoot=kTRUE for ESD/AOD analysis in AliROOT, and kFALSE for root analysis on reduced trees
  //get the current analysis manager


  if(isAliRoot){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
      Error("AddTask_iarsene_dst", "No analysis manager found.");
      return 0;
    }
  }

  const Bool_t TPCstandalone=kFALSE;

  // REDUCED TASKS

  AliAnalysisTaskTwoPwrtRP* task = new AliAnalysisTaskTwoPwrtRP("ChargeCorrelations");
  //task->SelectCollisionCandidates(AliVEvent::kMB);
  task->SetTPCstandaloneTracks(TPCstandalone);

  AliCMEAnalysisCuts* eventCuts=new AliCMEAnalysisCuts();
  AliCMEAnalysisCuts* trackCuts=new AliCMEAnalysisCuts();
  if(task->IsTracksTPCstandalone()) GetFilterTrackCuts(trackCuts);
  else GetFilterGlobalTrackCuts2010(trackCuts);

  GetEventCuts(eventCuts);
  task->SetEventCuts(eventCuts);

  AliCMEAnalysisCuts* TPCtrackCuts=new AliCMEAnalysisCuts();
  AliCMEAnalysisCuts* GlobaltrackCuts=new AliCMEAnalysisCuts();
  //GetFilterTPCFlowCuts(TPCtrackCuts);
  //GetFilterGlobalFlowCuts(GlobaltrackCuts);
  task->SetTPCtrackcuts(TPCtrackCuts);
  task->SetGlobaltrackcuts(GlobaltrackCuts);
  //task->SetRemoveOutliers(kTRUE);

  //Double_t Ctbins[][2]  = {{0.0, 2}, {100.0, 1}};
  //Double_t  Ptbins[][2]  = {{0.2, 2}, {5.0, 1}};
  //Double_t Etabins[][2] = {{-0.8, 2}, {0.8, 1}};
  //Double_t DEtabins[][2] = {{0.0, 2}, {0.8, 1}};
  //Double_t  Ptbins[][2]  = {{0.2, 5}, {0.5, 3}, {2.0, 15}, {3.0, 4}, {5.0, 8}};

  //Double_t Ctbins[][2]  =  {{ 0.0, 3}, {50.0, 50}, {100.0,25}};
  Double_t Ctbins[][2]  =  {{ 0.0, 3}, {10.0, 2}, {80.0,7}};
  Double_t  Ptbins[][2]  = {{ 0.0, 4}, {3.0, 30}, {5.0, 8}, {10.0, 1}};
  Double_t MPtbins[][2]  = {{ 0.0, 4}, {3.0, 30}, {5.0, 8}, {10.0, 1}};
  Double_t DPtbins[][2]  = {{ 0.0, 4}, {3.0, 30}, {5.0, 8}, {10.0, 1}};
  Double_t Etabins[][2]  = {{-0.8, 2}, {0.8, 16}};
  Double_t DEtabins[][2] = {{ 0.0, 2}, {1.6, 16}};


  TAxis axis[7];
  axis[0] =  AliQnCorrectionsAxes::MakeAxis(Ctbins);
  axis[1] =  AliQnCorrectionsAxes::MakeAxis(Ptbins);
  axis[2] =  AliQnCorrectionsAxes::MakeAxis(MPtbins);
  axis[3] =  AliQnCorrectionsAxes::MakeAxis(DPtbins);
  axis[4] =  AliQnCorrectionsAxes::MakeAxis(Etabins);
  axis[5] =  AliQnCorrectionsAxes::MakeAxis(Etabins);
  axis[6] =  AliQnCorrectionsAxes::MakeAxis(DEtabins);

  task->SetAxes(&axis[0], &axis[1], &axis[2], &axis[3], &axis[4], &axis[5], &axis[6]);
  //task->SetAxes(&axis[0], &axis[1], &axis[2], 0x0, &axis[4], 0x0, 0x0);

  TArrayD binLimits[3];
  Char_t* eventplanes[4] = {"VZEROA","VZEROC","FMDA","FMDC"};

  AliChargeOnePwrtRP* oneP[10];

  vec[0] = new AliChargeOnePwrtRP();
  vec[0]->SetVars(0,1,4);
  vec[0]->SetTrackCuts(trackCuts);
  vec[0]->SetRangeHarmonics(1,4);
  vec[0]->SetTHn("PRLcuts","TPCcuts","c","1Pt","1Eta",eventplanes);
  task->Set1PwrtRP(vec[0]);

  //vec[1] = new AliChargeOnePwrtRP();
  //vec[1]->SetVars(0,1,4);
  //vec[1]->SetRangeHarmonics(1,4);
  //vec[1]->SetTHn("PRLcuts","TPCcuts","p","1Pt","1Eta",eventplanes);
  //task->Set1PwrtRP(vec[1],AddProtonCuts(trackCuts));

  //vec[2] = new AliChargeOnePwrtRP();
  //vec[2]->SetVars(0,1,4);
  //vec[2]->SetRangeHarmonics(2,4);
  //vec[2]->SetTHn("PRLcuts","TPCcuts","K","1Pt","1Eta",eventplanes);
  //task->Set1PwrtRP(vec[2],AddKaonCuts(trackCuts));

  //vec[3] = new AliChargeOnePwrtRP();
  //vec[3]->SetVars(0,1,4);
  //vec[3]->SetRangeHarmonics(2,4);
  //vec[3]->SetTHn("PRLcuts","TPCcuts","p","1Pt","1Eta",eventplanes);
  //task->Set1PwrtRP(vec[3],trackCuts);

  //AliChargeTwoPwrtRP* twoP = new AliChargeTwoPwrtRP(2, axis);
  //twoP->SetVars(0,1);
  //twoP->SetTHn("PRLcuts","TPCcuts","TPCcuts","c","c","1Pt");
  ////twoP->SetTHn("PRLcuts","TPCcuts","TPCcuts","c","c","1Pt");
  ////task->Set2PwrtRP(twoP, trackCuts, trackCuts);

  AliChargeTwoPwrtRP* twoP[5];

  //twoP[0] = new AliChargeTwoPwrtRP();
  //twoP[0]->SetVars(0,1,4);
  //twoP[0]->SetTHn("PRLcuts","TPCcuts","TPCcuts","c","c","1Pt","1Eta",eventplanes);
  //twoP[0]->SetTrackCuts(trackCuts, trackCuts);
  //task->Set2PwrtRP(twoP[0], trackCuts, trackCuts);

  twoP[0] = new AliChargeTwoPwrtRP();
  twoP[0]->SetVars(0,3,6);
  twoP[0]->SetTHn("PRLcuts","TPCcuts","TPCcuts","c","c","dPt","dEta",eventplanes);
  twoP[0]->SetTrackCuts(trackCuts, trackCuts);
  task->Set2PwrtRP(twoP[0]);

  twoP[1] = new AliChargeTwoPwrtRP();
  twoP[1]->SetVars(0,2,5);
  twoP[1]->SetTHn("PRLcuts","TPCcuts","TPCcuts","c","c","mPt","mEta",eventplanes);
  twoP[1]->SetTrackCuts(trackCuts, trackCuts);
  task->Set2PwrtRP(twoP[1]);


  //task->Setup();
  Setup(task, 2010);



  if(isAliRoot){
  //AliAnalysisDataContainer* cReducedEvent = mgr->GetContainers()->FindObject("ReducedEventDQ");
  //if(!cReducedEvent) {cout<<"Couldn't find exchange container with ReducedEvent, run AddTask_iarsene_dst.C"<<endl; return;}
  //AliAnalysisTaskReducedEventProcessor* redTask = new AliAnalysisTaskReducedEventProcessor("ReducedEventAnalysisManager");
  //redTask->AddTask(task);

  mgr->AddTask(task);

  AliAnalysisDataContainer *cOutputHist =
    mgr->CreateContainer("ChargeCorrelations",
                         THashList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "ChargeCorrelations.root");

  mgr->ConnectInput( task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cOutputHist );

  }
  return task;
}



//_____________________________________________________________________
AliCMEAnalysisCuts* AddPionCuts(AliCMEAnalysisCuts *trackCuts){
  AliCMEAnalysisCuts* cuts=new AliCMEAnalysisCuts();
  cuts->SetName("Pions");
  cuts->CopyCuts(trackCuts);
  cuts->AddCut(AliCMEVarManager::kBayes+1,0.89,1.0);
  return cuts;
}
//_____________________________________________________________________
AliCMEAnalysisCuts* AddKaonCuts(AliCMEAnalysisCuts *trackCuts){
  AliCMEAnalysisCuts* cuts=new AliCMEAnalysisCuts();
  cuts->SetName("Kaons");
  cuts->CopyCuts(trackCuts);
  cuts->AddCut(AliCMEVarManager::kBayes+2,0.89,1.0);
  return cuts;
}

//_____________________________________________________________________
AliCMEAnalysisCuts* AddProtonCuts(AliCMEAnalysisCuts *trackCuts){
  AliCMEAnalysisCuts* cuts=new AliCMEAnalysisCuts();
  cuts->SetName("Protons");
  cuts->CopyCuts(trackCuts);
  cuts->AddCut(AliCMEVarManager::kBayes+3,0.89,1.0);
  return cuts;
}


//_____________________________________________________________________
void GetFilterGlobalTrackCuts2010(AliCMEAnalysisCuts *trackCuts) {
  trackCuts->SetName("Global");

  trackCuts->AddCut(AliCMEVarManager::kEta,-0.8,0.8);
  trackCuts->AddCut(AliCMEVarManager::kPt,0.2,5.0);
  trackCuts->AddCut(AliCMEVarManager::kImpactPar2D,0.0,1.0);
  //trackCuts->AddCut(AliCMEVarManager::kImpactParXY,-0.3,0.3);  // hardcoded in analysis task
  //trackCuts->AddCut(AliCMEVarManager::kImpactParZ,-0.3,0.3);   // hardcoded in analysis task
  trackCuts->AddCut(AliCMEVarManager::kNclsTPC,69.5,161.0);
  trackCuts->AddCut(AliCMEVarManager::kTPCsignal,10.0,50000.0);
  trackCuts->AddCut(AliCMEVarManager::kTPCchi2Cl,0.1,4.0);
  trackCuts->AddFlag(AliCMEVarManager::kTrackingStatus+AliCMEVarManager::kTPCrefit  , 0.5,1.5);
  trackCuts->AddFlag(AliCMEVarManager::kTrackingStatus+AliCMEVarManager::kITSrefit  , 0.5,1.5);
  trackCuts->AddFlag(AliCMEVarManager::kTrackingStatus+AliCMEVarManager::kKinkIndex0,-0.5,0.5);
}


//_____________________________________________________________________
void GetFilterTrackCuts(AliCMEAnalysisCuts *trackCuts) {
  trackCuts->SetName("TrackCuts");
  trackCuts->AddCut(AliCMEVarManager::kEta,-0.8,0.8);
  trackCuts->AddCut(AliCMEVarManager::kPt,0.2,5.0);
  trackCuts->AddCut(AliCMEVarManager::kImpactPar2D,0.0,1.0);
  //trackCuts->AddCut(AliCMEVarManager::kImpactParXY,-3.0,3.0);
  //trackCuts->AddCut(AliCMEVarManager::kImpactParZ,-3.0,3.0);
  trackCuts->AddCut(AliCMEVarManager::kNclsTPC,69.5,161.0);
  trackCuts->AddCut(AliCMEVarManager::kTPCsignal,10.0,50000.0);
  trackCuts->AddCut(AliCMEVarManager::kTPCchi2Cl,0.2,4.0);
  //trackCuts->AddCut(AliCMEVarManager::kEMCALmatchedEnergy,0.5,1.5);
  //trackCuts->AddFlag(AliCMEVarManager::kTrackingStatus,AliCMEVarManager::kTPCrefit,kTRUE);
  trackCuts->AddFlag(AliCMEVarManager::kKinkIndex0,-0.5,0.5);
}



//_____________________________________________________________________
void GetEventCuts(AliCMEAnalysisCuts *eventCuts) {
  eventCuts->SetName("Event");
  eventCuts->AddCut(AliCMEVarManager::kCentrality,0.0,80.00);
  //eventCuts->AddCut(AliCMEVarManager::kCentVZEROmTPC,-5.0,5.0);
  eventCuts->AddCut(AliCMEVarManager::kZvPrim,-10.0,10.0);
  //eventCuts->AddCut(AliCMEVarManager::kCentQuality,-0.5,0.5);
  //eventCuts->AddCut(AliCMEVarManager::kOnlineTriggersFired+1,0.5,1.5); // kOnlineTriggersFired+i = Log2(AliVEvent::EOfflineTriggerTypes)
}




//_________________________________________________________________
void Setup(AliAnalysisTaskTwoPwrtRP* processor, Int_t period /*=-1*/) {
  //
  // Configure the event processor
  // Setup histograms, handlers, cuts, etc.
  //
  AliHistogramManager* hMan = new AliHistogramManager("HistogramManager", AliCMEVarManager::kNMaxValues);
  //hMan->SetUseDefaultVariableNames(kTRUE);
  //hMan->SetDefaultVarNames(AliCMEVarManager::fgkParticleNames,AliCMEVarManager::fgVariableUnits);
  processor->SetHistogramManager(hMan);
  SetupHistogramManager(hMan, period, processor);

}


//_________________________________________________________________
void SetupHistogramManager(AliHistogramManager* man, Int_t period /*=-1*/, AliAnalysisTaskSE* processor) {
  //
  // setup the histograms manager
  //
  //AliCMEVarManager::SetDefaultVarNames();

  DefineHistograms(man, period, processor);

  //AliCMEVarManager::SetUseVars(man->GetUsedVars());
}


//_________________________________________________________________
void DefineHistograms(AliHistogramManager* man, Int_t period /*=-1*/, AliAnalysisTaskSE* processor) {
  //
  // define histograms
  //
  TString histClasses = "";
  histClasses += "Event_NoCuts;";        //ok
  histClasses += "Event_AfterCuts;";     //ok
  histClasses += "Event_NoEventPlane;";     //ok

  AliAnalysisTaskTwoPwrtRP* analysis=(AliAnalysisTaskTwoPwrtRP*) processor;

  for(Int_t ic=0; ic<analysis->GetNOnePwrtRP(); ic++){
    //histClasses += "TrackQA_"+analysis->GetOnePwrtRP_TrackCuts(ic)->GetCutsName()+"_TPC_TOF_ITS;";
    histClasses += Form("TrackQA_%s_TPC_TOF_ITS;",analysis->Get1PwrtRP(ic)->GetTrackCuts()->GetName());
    histClasses += "TrackQA_Profile_TPC_TOF_ITS;";
  }

  //histClasses += "TrackQA_Pions_TPC_TOF_ITS;";     //ok
  //histClasses += "TrackQA_Kaons_TPC_TOF_ITS;";     //ok
  //histClasses += "TrackQA_Protons_TPC_TOF_ITS;";     //ok
  histClasses += "EvtTags;L0TriggerInput;L1TriggerInput;L2TriggerInput;";   //ok
  histClasses += "OnlineTriggers_vs_L0TrigInputs;OnlineTriggers_vs_L1TrigInputs;OnlineTriggers_vs_L2TrigInputs;";  //ok
  histClasses += "ITSclusterMap;TPCclusterMap;";   //ok
  histClasses += "CaloClusters;";   //ok
  histClasses += "OnlineTriggers_NoCuts;";  //ok
  histClasses += "OnlineTriggers_AfterCuts;";    //ok

  Int_t runNBins = 0;
  Double_t runHistRange[2] = {0.0,0.0};

  // Pb-Pb from 2011 run range is default
  runNBins = 2700;
  runHistRange[0] = 167900.;
  runHistRange[1] = 170600.;

  if(period==2010) {
    runNBins = 2501;
    runHistRange[0] = 137099.5;
    runHistRange[1] = 139600.5;
  }
  if(period==2011) { // ??
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
      man->AddHistogram(classStr.Data(),"RunNo","Run numbers",kFALSE, runNBins, runHistRange[0], runHistRange[1], AliCMEVarManager::kRunNumber);
      man->AddHistogram(classStr.Data(),"VtxX","Vtx X",kFALSE,300,-0.4,0.4,AliCMEVarManager::kXvPrim);
      man->AddHistogram(classStr.Data(),"VtxY","Vtx Y",kFALSE,300,-0.4,0.4,AliCMEVarManager::kYvPrim);
      man->AddHistogram(classStr.Data(),"VtxZ","Vtx Z",kFALSE,300,-15.,15.,AliCMEVarManager::kZvPrim);
      man->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO)",kFALSE, 100, 0.0, 100.0, AliCMEVarManager::kCentrality);
      man->AddHistogram(classStr.Data(),"CentVZEROA","Centrality(VZERO-A)",kFALSE, 100, 0.0, 100.0, AliCMEVarManager::kCentralityV0A);
      man->AddHistogram(classStr.Data(),"CentVZEROC","Centrality(VZERO-C)",kFALSE, 100, 0.0, 100.0, AliCMEVarManager::kCentralityV0C);
      man->AddHistogram(classStr.Data(),"CentVZEROmTPC","Centrality(VZERO)-Centrality(TPC)",kFALSE, 100, 0.0, 100.0, AliCMEVarManager::kCentrality, 200, -100.0, 100.0, AliCMEVarManager::kCentralityV0mSPD);
      man->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD)",kFALSE, 100, 0.0, 100.0, AliCMEVarManager::kCentralitySPD);
//      //man->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC)",kFALSE, 100, 0.0, 100.0, AliCMEVarManager::kCentralityTPC);
//      //man->AddHistogram(classStr.Data(),"CentZDC","Centrality(ZDC)",kFALSE, 100, 0.0, 100.0, AliCMEVarManager::kCentralityZDC);
//      man->AddHistogram(classStr.Data(),"CentZNA","Centrality(ZNA)",kFALSE, 100, 0.0, 100.0, AliCMEVarManager::kCentralityZNA);
//      //man->AddHistogram(classStr.Data(),"CentQuality","Centrality quality",kFALSE, 100, -50.5, 49.5, AliCMEVarManager::kCentQuality);
//      man->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event",kFALSE,500,0.,20000.,AliCMEVarManager::kNtracksTotal);
//      man->AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event",kFALSE,500,0.,20000.,AliCMEVarManager::kNtracksSelected);
//      man->AddHistogram(classStr.Data(),"NTracksTPCvsGlobal","Number of selected TPC tracks vs Global tracks per event",kFALSE,2500,0.,5000.,AliCMEVarManager::kNtracksSelected,2500,0.,5000.,AliCMEVarManager::kNtracksTotal);
//      man->AddHistogram(classStr.Data(),"EventNumberInESDFile","Event number in ESD file",kFALSE, 1000, 0.0, 1000.0, AliCMEVarManager::kEventNumberInFile);
//      man->AddHistogram(classStr.Data(),"BC","Bunch crossing",kFALSE,3500,0.,3500.,AliCMEVarManager::kBC);
//      man->AddHistogram(classStr.Data(),"EventType","Event type",kFALSE,100,0.,100.,AliCMEVarManager::kEventType);
//      man->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag",kFALSE,
//		   2,-0.5,1.5,AliCMEVarManager::kIsPhysicsSelection, 0,0.0,0.0,AliCMEVarManager::kNothing, 0,0.0,0.0,AliCMEVarManager::kNothing, "off;on");
//      man->AddHistogram(classStr.Data(),"IsSPDPileup","Event has pileup (SPD)",kFALSE,
//		   2,-0.5,1.5,AliCMEVarManager::kIsSPDPileup, 0,0.0,0.0,AliCMEVarManager::kNothing, 0,0.0,0.0,AliCMEVarManager::kNothing, "no;yes");
//      man->AddHistogram(classStr.Data(),"IsSPDPileupMultBins","Event has pileup (SPD multiplicity bins)",kFALSE,
//		   2,-0.5,1.5,AliCMEVarManager::kIsSPDPileupMultBins, 0,0.0,0.0,AliCMEVarManager::kNothing, 0,0.0,0.0,AliCMEVarManager::kNothing, "no;yes");
//      man->AddHistogram(classStr.Data(),"IRIntClosestInt1Map","Closest out of bunch interactions Int1",kFALSE,7000,-3500.,3500.,AliCMEVarManager::kIRIntClosestIntMap);
//      man->AddHistogram(classStr.Data(),"IRIntClosestInt2Map","Closest out of bunch interactions Int2",kFALSE,7000,-3500.,3500.,AliCMEVarManager::kIRIntClosestIntMap+1);
//      man->AddHistogram(classStr.Data(),"NVtxContributors","",kFALSE,500,0.,10000.,AliCMEVarManager::kNVtxContributors);
//      man->AddHistogram(classStr.Data(),"VtxXtpc","Vtx X (TPC)",kFALSE,300,-1.,1.,AliCMEVarManager::kVtxXtpc);
//      man->AddHistogram(classStr.Data(),"VtxYtpc","Vtx Y (TPC)",kFALSE,300,-1.,1.,AliCMEVarManager::kVtxYtpc);
//      man->AddHistogram(classStr.Data(),"VtxZtpc","Vtx Z (TPC)",kFALSE,300,-15.,15.,AliCMEVarManager::kVtxZtpc);
//      man->AddHistogram(classStr.Data(),"DeltaVtxZ","Z_{global}-Z_{TPC}",kFALSE,300,-1.,1.,AliCMEVarManager::kDeltaVtxZ);
//      man->AddHistogram(classStr.Data(),"NVtxTPCContributors","",kFALSE,500,0.,10000.,AliCMEVarManager::kNVtxTPCContributors);
//      man->AddHistogram(classStr.Data(),"NSPDpileups","No. SPD pileups",kFALSE,10,0.,10.,AliCMEVarManager::kNSPDpileups);
//      man->AddHistogram(classStr.Data(),"NTrackPileups","No. Track pileups",kFALSE,10,0.,10.,AliCMEVarManager::kNTrackPileups);
//      man->AddHistogram(classStr.Data(),"NPMDtracks","No. PMD tracks",kFALSE,200,0.,200.,AliCMEVarManager::kNPMDtracks);
//      man->AddHistogram(classStr.Data(),"NTRDtracks","No. TRD tracks",kFALSE,200,0.,200.,AliCMEVarManager::kNTRDtracks);
//      man->AddHistogram(classStr.Data(),"NTRDtracklets","No. TRD tracklets",kFALSE,500,0.,50000.,AliCMEVarManager::kNTRDtracklets);
//      man->AddHistogram(classStr.Data(),"NV0sTotal","Number of V0 candidates per event",kFALSE,200,0.,30000.,AliCMEVarManager::kNV0total);
//      man->AddHistogram(classStr.Data(),"NV0sSelected","Number of selected V0 candidates per event",kFALSE,200,0.,2000.,AliCMEVarManager::kNV0selected);
//      man->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0", kFALSE, 200, 0., 5000., AliCMEVarManager::kSPDntracklets);
//      for(Int_t i=0; i<32; ++i)
//	man->AddHistogram(classStr.Data(),Form("SPDntracklets_%d",i), "", kFALSE, 100, 0., 300., AliCMEVarManager::kSPDntrackletsEta+i);
//      for(Int_t il=0; il<2; ++il)
//        man->AddHistogram(classStr.Data(), Form("SPDfiredChips_layer%d",il+1), Form("SPD fired chips in layer %d",il+1),
//			  kFALSE, 200, 0., 600., AliCMEVarManager::kSPDFiredChips+il);
//      for(Int_t il=0; il<6; ++il)
//        man->AddHistogram(classStr.Data(), Form("ITSclusters_layer%d",il+1), Form("ITS clusters in layer %d",il+1),
//			  kFALSE, 200, 0., 10000., AliCMEVarManager::kITSnClusters+il);
//      man->AddHistogram(classStr.Data(), "SPDnSingleClusters", "SPD single clusters",
//			kFALSE, 200, 0., 6000., AliCMEVarManager::kSPDnSingleClusters);
//      for(Int_t i=0; i<64; ++i)
//        man->AddHistogram(classStr.Data(), Form("VZEROmult_ch%d",i), "", kFALSE, 200, 0.0, 1000.0, AliCMEVarManager::kVZEROChannelMult+i);
//      man->AddHistogram(classStr.Data(),"VZEROAmult", "", kFALSE, 300, 0.0, 20000., AliCMEVarManager::kVZEROATotalMult);
//      man->AddHistogram(classStr.Data(),"VZEROCmult", "", kFALSE, 300, 0.0, 20000., AliCMEVarManager::kVZEROCTotalMult);
//      man->AddHistogram(classStr.Data(),"VZEROmult", "", kFALSE, 300, 0.0, 30000., AliCMEVarManager::kVZEROTotalMult);
//      man->AddHistogram(classStr.Data(),"VZEROAemptyChannels", "", kFALSE, 64, 0.0, 64., AliCMEVarManager::kVZEROAemptyChannels);
//      man->AddHistogram(classStr.Data(),"VZEROCemptyChannels", "", kFALSE, 64, 0.0, 64., AliCMEVarManager::kVZEROCemptyChannels);
//      man->AddHistogram(classStr.Data(),"VZEROmult_VtxContributors", "", kFALSE,
//		   200, 0.0, 10000., AliCMEVarManager::kVZEROTotalMult, 200, 0.0, 1000., AliCMEVarManager::kNVtxContributors);
//      man->AddHistogram(classStr.Data(),"VZEROAmult_VtxContributors", "", kFALSE,
//		   200, 0.0, 10000., AliCMEVarManager::kVZEROATotalMult, 200, 0.0, 1000., AliCMEVarManager::kNVtxContributors);
//      man->AddHistogram(classStr.Data(),"VZEROCmult_VtxContributors", "", kFALSE,
//		   200, 0.0, 10000., AliCMEVarManager::kVZEROCTotalMult, 200, 0.0, 1000., AliCMEVarManager::kNVtxContributors);
//      for(Int_t i=0; i<10; ++i)
//	man->AddHistogram(classStr.Data(), Form("ZDCNenergy_ch%d",i), "", kFALSE, 200, 0.0, (i==0 || i==5 ? 100000.0 : 20000.), AliCMEVarManager::kZDCnEnergyCh+i);
//      for(Int_t i=0; i<10; ++i)
//	man->AddHistogram(classStr.Data(), Form("ZDCPenergy_ch%d",i), "", kFALSE, 200, 0.0, (i==0 || i==5 ? 40000.0 : 20000.), AliCMEVarManager::kZDCpEnergyCh+i);
//      for(Int_t i=0; i<26; ++i)
//	man->AddHistogram(classStr.Data(), Form("TZEROamp_ch%d",i), "", kFALSE, 200, 0.0, 100.0, AliCMEVarManager::kTZEROAmplitudeCh+i);
//      for(Int_t i=0; i<3; ++i)
//	man->AddHistogram(classStr.Data(), Form("TZERO_TOF%d",i), "", kFALSE, 200, -1000.0, 1000.0, AliCMEVarManager::kTZEROTOF+i);
//      for(Int_t i=0; i<3; ++i)
//	man->AddHistogram(classStr.Data(), Form("TZERO_TOFbest%d",i), "", kFALSE, 200, -1000.0, 1000.0, AliCMEVarManager::kTZEROTOFbest+i);
//      man->AddHistogram(classStr.Data(),"TZEROzVtx", "", kFALSE, 300, -15.0, 15., AliCMEVarManager::kTZEROzVtx);
//      man->AddHistogram(classStr.Data(),"TZEROstartTime", "", kFALSE, 1000, -1000.0, 1000., AliCMEVarManager::kTZEROstartTime);
//      man->AddHistogram(classStr.Data(),"TZEROpileup", "TZERO pileup", kFALSE, 2, -0.5, 1.5, AliCMEVarManager::kTZEROpileup);
//      man->AddHistogram(classStr.Data(),"TZEROsatellite", "TZERO satellite", kFALSE, 2, -0.5, 1.5, AliCMEVarManager::kTZEROsatellite);
//
//      man->AddHistogram(classStr.Data(),"NTracksTPCoutVsITSout","",kFALSE,100,0.,15.,AliCMEVarManager::kNTracksTPCoutVsITSout);
//      man->AddHistogram(classStr.Data(),"NTracksTRDoutVsITSout","",kFALSE,100,0.,5.,AliCMEVarManager::kNTracksTRDoutVsITSout);
//      man->AddHistogram(classStr.Data(),"NTracksTOFoutVsITSout","",kFALSE,100,0.,1.,AliCMEVarManager::kNTracksTOFoutVsITSout);
//      man->AddHistogram(classStr.Data(),"NTracksTRDoutVsTPCout","",kFALSE,100,0.,1.,AliCMEVarManager::kNTracksTRDoutVsTPCout);
//      man->AddHistogram(classStr.Data(),"NTracksTOFoutVsTPCout","",kFALSE,100,0.,0.6,AliCMEVarManager::kNTracksTOFoutVsTPCout);
//      man->AddHistogram(classStr.Data(),"NTracksTOFoutVsTRDout","",kFALSE,100,0.,2.,AliCMEVarManager::kNTracksTOFoutVsTRDout);
//      man->AddHistogram(classStr.Data(),"NTracksITSoutVsSPDtracklets","",kFALSE,100,0.,50.,AliCMEVarManager::kNTracksITSoutVsSPDtracklets);
//      man->AddHistogram(classStr.Data(),"NTracksTPCoutVsSPDtracklets","",kFALSE,100,0.,50.,AliCMEVarManager::kNTracksTPCoutVsSPDtracklets);
//      man->AddHistogram(classStr.Data(),"NTracksTRDoutVsSPDtracklets","",kFALSE,100,0.,30.,AliCMEVarManager::kNTracksTRDoutVsSPDtracklets);
//      man->AddHistogram(classStr.Data(),"NTracksTOFoutVsSPDtracklets","",kFALSE,100,0.,10.,AliCMEVarManager::kNTracksTOFoutVsSPDtracklets);
//      continue;
    }  // end if className contains "Event"
//
//    if(classStr.Contains("EvtTags")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      TString tagNames = "";
//      for(Int_t i=0; i<64; ++i) {tagNames += Form("Tag %d", i); tagNames+=";";}
//      man->AddHistogram(classStr.Data(), "EventTags", "Event tags", kFALSE,
//	           64, -0.5, 63.5, AliCMEVarManager::kEventTag, 0, 0.0, 0.0, AliCMEVarManager::kNothing, 0, 0.0, 0.0, AliCMEVarManager::kNothing, tagNames.Data());
//      continue;
//    }
//    if(classStr.Contains("L0TriggerInput")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      TString trigInputs = "";
//      for(Int_t i=0; i<32; ++i) {trigInputs += Form("L0 Input %d", i); trigInputs+=";";}
//      man->AddHistogram(classStr.Data(), "L0TriggerInputs", "L0 trigger inputs", kFALSE,
//	           32, -0.5, 31.5, AliCMEVarManager::kL0TriggerInput, 0, 0.0, 0.0, AliCMEVarManager::kNothing, 0, 0.0, 0.0, AliCMEVarManager::kNothing, trigInputs.Data());
//      continue;
//    }
//    if(classStr.Contains("L1TriggerInput")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      TString trigInputs = "";
//      for(Int_t i=0; i<32; ++i) {trigInputs += Form("L1 Input %d", i); trigInputs+=";";}
//      man->AddHistogram(classStr.Data(), "L1TriggerInputs", "L1 trigger inputs", kFALSE,
//	           32, -0.5, 31.5, AliCMEVarManager::kL1TriggerInput, 0, 0.0, 0.0, AliCMEVarManager::kNothing, 0, 0.0, 0.0, AliCMEVarManager::kNothing, trigInputs.Data());
//      continue;
//    }
//    if(classStr.Contains("L2TriggerInput")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      TString trigInputs = "";
//      for(Int_t i=0; i<16; ++i) {trigInputs += Form("L2 Input %d", i); trigInputs+=";";}
//      man->AddHistogram(classStr.Data(), "L2TriggerInputs", "L2 trigger inputs", kFALSE,
//	           16, -0.5, 15.5, AliCMEVarManager::kL2TriggerInput, 0, 0.0, 0.0, AliCMEVarManager::kNothing, 0, 0.0, 0.0, AliCMEVarManager::kNothing, trigInputs.Data());
//      continue;
//    }
//
//    // Offline trigger histograms
//    if(classStr.Contains("OnlineTriggers")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      TString triggerNames = "";
//      for(Int_t i=0; i<64; ++i) {triggerNames += AliCMEVarManager::fgkOfflineTriggerNames[i]; triggerNames+=";";}
//
//      man->AddHistogram(classStr.Data(), "Triggers", "", kFALSE,
//	           64, -0.5, 63.5, AliCMEVarManager::kOnlineTrigger, 2, -0.5, 1.5, AliCMEVarManager::kOnlineTriggerFired, 0, 0.0, 0.0, AliCMEVarManager::kNothing, triggerNames.Data(), "off;on");
//      man->AddHistogram(classStr.Data(), "Triggers2", "", kFALSE,
//	           64, -0.5, 63.5, AliCMEVarManager::kOnlineTriggerFired2, 0, 0.0, 0.0, AliCMEVarManager::kNothing, 0, 0.0, 0.0, AliCMEVarManager::kNothing, triggerNames.Data());
//      man->AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "", kFALSE,
//	           64, -0.5, 63.5, AliCMEVarManager::kOnlineTriggerFired2, 20, 0.0, 100.0, AliCMEVarManager::kCentrality, 0, 0.0, 0.0, AliCMEVarManager::kNothing, triggerNames.Data());
//      continue;
//    }
//
//    if(classStr.Contains("OnlineTriggers_vs_L0TrigInputs")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      TString triggerNames = "";
//      for(Int_t i=0; i<64; ++i) {triggerNames += AliCMEVarManager::fgkOfflineTriggerNames[i]; triggerNames+=";";}
//      TString trigInputs = "";
//      for(Int_t i=0; i<32; ++i) {trigInputs += Form("L0 Input %d", i); trigInputs+=";";}
//      man->AddHistogram(classStr.Data(), "OnlineTriggers_L0TriggerInputs", "", kFALSE,
//	           64, -0.5, 63.5, AliCMEVarManager::kOnlineTriggerFired2, 32, -0.5, 31.5, AliCMEVarManager::kL0TriggerInput, 0, 0.0, 0.0, AliCMEVarManager::kNothing,
//		   triggerNames.Data(), trigInputs);
//    }
//    if(classStr.Contains("OnlineTriggers_vs_L1TrigInputs")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      TString triggerNames = "";
//      for(Int_t i=0; i<64; ++i) {triggerNames += AliCMEVarManager::fgkOfflineTriggerNames[i]; triggerNames+=";";}
//      TString trigInputs = "";
//      for(Int_t i=0; i<32; ++i) {trigInputs += Form("L1 Input %d", i); trigInputs+=";";}
//      man->AddHistogram(classStr.Data(), "OnlineTriggers_L1TriggerInputs", "", kFALSE,
//	           64, -0.5, 63.5, AliCMEVarManager::kOnlineTriggerFired2, 32, -0.5, 31.5, AliCMEVarManager::kL1TriggerInput, 0, 0.0, 0.0, AliCMEVarManager::kNothing,
//		   triggerNames.Data(), trigInputs);
//    }
//    if(classStr.Contains("OnlineTriggers_vs_L2TrigInputs")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      TString triggerNames = "";
//      for(Int_t i=0; i<64; ++i) {triggerNames += AliCMEVarManager::fgkOfflineTriggerNames[i]; triggerNames+=";";}
//      TString trigInputs = "";
//      for(Int_t i=0; i<16; ++i) {trigInputs += Form("L2 Input %d", i); trigInputs+=";";}
//      man->AddHistogram(classStr.Data(), "OnlineTriggers_L2TriggerInputs", "", kFALSE,
//	           64, -0.5, 63.5, AliCMEVarManager::kOfflineTriggerFired2, 16, -0.5, 15.5, AliCMEVarManager::kL2TriggerInput, 0, 0.0, 0.0, AliCMEVarManager::kNothing,
//		   triggerNames.Data(), trigInputs);
//    }
//
//    if(classStr.Contains("ITSclusterMap")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      man->AddHistogram(classStr.Data(), "ITSlayerHit", "Hits in the ITS layers", kFALSE,
//                   13, -6.5, 6.5, AliCMEVarManager::kITSlayerHit);
//      man->AddHistogram(classStr.Data(), "ITSlayerHit_Phi", "Hits in the ITS layers vs #varphi", kFALSE,
//                   100, 0.0, 6.29, AliCMEVarManager::kPhi, 13, -6.5, 6.5, AliCMEVarManager::kITSlayerHit);
//      man->AddHistogram(classStr.Data(), "ITSlayerHit_Eta", "Hits in the ITS layers vs #eta", kFALSE,
//                   100, -1.0, 1.0, AliCMEVarManager::kEta, 13, -6.5, 6.5, AliCMEVarManager::kITSlayerHit);
//      man->AddHistogram(classStr.Data(), "ITSlayerHit_DCAxy", "Hits in the ITS layers vs DCA_{xy}", kFALSE,
//                   1000, -0.5, 0.5, AliCMEVarManager::kDcaXY, 13, -6.5, 6.5, AliCMEVarManager::kITSlayerHit);
//      man->AddHistogram(classStr.Data(), "ITSlayerHit_DCAz", "Hits in the ITS layers vs DCA_{z}", kFALSE,
//                   1800, -1.0, 1.0, AliCMEVarManager::kDcaZ, 13, -6.5, 6.5, AliCMEVarManager::kITSlayerHit);
//      continue;
//    }  // end of ITSclusterMap histogram definitions
//
//    if(classStr.Contains("TPCclusterMap")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      man->AddHistogram(classStr.Data(), "TPCclusterMap", "TPC cluster map", kFALSE,
//                   8, -0.5, 7.5, AliCMEVarManager::kTPCclusBitFired);
//      man->AddHistogram(classStr.Data(), "TPCclusterMap_Phi", "TPC cluster map vs #varphi", kFALSE,
//                   360, 0.0, 6.29, AliCMEVarManager::kPhi, 8, -0.5, 7.5, AliCMEVarManager::kTPCclusBitFired);
//      man->AddHistogram(classStr.Data(), "TPCclusterMap_Eta", "TPC cluster map vs #eta", kFALSE,
//                   100, -1.0, 1.0, AliCMEVarManager::kEta, 8, -0.5, 7.5, AliCMEVarManager::kTPCclusBitFired);
//      man->AddHistogram(classStr.Data(), "TPCclusterMap_Pt", "TPC cluster map vs p_{T}", kFALSE,
//                   100, 0.0, 10.0, AliCMEVarManager::kPt, 8, -0.5, 7.5, AliCMEVarManager::kTPCclusBitFired);
//      continue;
//    }  // end of TPCclusterMap histogram definitions
//
//    // Calorimeter cluster histograms
//    if(classStr.Contains("CaloClusters")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      man->AddHistogram(classStr.Data(), "Energy", "Cluster energy", kFALSE,
//	           2, 0.5, 2.5, AliCMEVarManager::kEMCALdetector, 500, 0.0, 20.0, AliCMEVarManager::kEMCALclusterEnergy, 0, 0.0, 0.0, AliCMEVarManager::kNothing, "EMCAL;PHOS");
//      man->AddHistogram(classStr.Data(), "Dx", "Cluster dx", kFALSE,
//	           2, 0.5, 2.5, AliCMEVarManager::kEMCALdetector, 600, -300.0, 300.0, AliCMEVarManager::kEMCALclusterDx, 0, 0.0, 0.0, AliCMEVarManager::kNothing, "EMCAL;PHOS");
//      man->AddHistogram(classStr.Data(), "Dz", "Cluster dz", kFALSE,
//	           2, 0.5, 2.5, AliCMEVarManager::kEMCALdetector, 600, -300.0, 300.0, AliCMEVarManager::kEMCALclusterDz, 0, 0.0, 0.0, AliCMEVarManager::kNothing, "EMCAL;PHOS");
//      man->AddHistogram(classStr.Data(), "M20", "Cluster M20", kFALSE,
//	           2, 0.5, 2.5, AliCMEVarManager::kEMCALdetector, 400, -2.0, 200.0, AliCMEVarManager::kEMCALm20, 0, 0.0, 0.0, AliCMEVarManager::kNothing, "EMCAL;PHOS");
//      man->AddHistogram(classStr.Data(), "M02", "Cluster M02", kFALSE,
//	           2, 0.5, 2.5, AliCMEVarManager::kEMCALdetector, 400, -2.0, 200.0, AliCMEVarManager::kEMCALm02, 0, 0.0, 0.0, AliCMEVarManager::kNothing, "EMCAL;PHOS");
//      man->AddHistogram(classStr.Data(), "Dispersion", "Cluster dispersion", kFALSE,
//	           2, 0.5, 2.5, AliCMEVarManager::kEMCALdetector, 400, 0.0, 50.0, AliCMEVarManager::kEMCALdispersion, 0, 0.0, 0.0, AliCMEVarManager::kNothing, "EMCAL;PHOS");
//      man->AddHistogram(classStr.Data(), "M20_M02", "", kFALSE, 400, -2.0, 200., AliCMEVarManager::kEMCALm02, 400, -2.0, 200., AliCMEVarManager::kEMCALm20);
//
//      man->AddHistogram(classStr.Data(), "CaloClusterEnergy_Run", "Cluster <energy> vs run", kTRUE,
//                   runNBins, runHistRange[0], runHistRange[1], AliCMEVarManager::kRunNo, 2, 0.5, 2.5, AliCMEVarManager::kEMCALdetector, 200, 0.0, 50.0, AliCMEVarManager::kEMCALclusterEnergy, "", "EMCAL;PHOS");
//      man->AddHistogram(classStr.Data(), "CaloClusterDx_Run", "Cluster <dx> vs run", kTRUE,
//                   runNBins, runHistRange[0], runHistRange[1], AliCMEVarManager::kRunNo, 2, 0.5, 2.5, AliCMEVarManager::kEMCALdetector, 200, -300.0, 300.0, AliCMEVarManager::kEMCALclusterDx, "", "EMCAL;PHOS");
//      man->AddHistogram(classStr.Data(), "CaloClusterDz_Run", "Cluster <dz> vs run", kTRUE,
//                   runNBins, runHistRange[0], runHistRange[1], AliCMEVarManager::kRunNo, 2, 0.5, 2.5, AliCMEVarManager::kEMCALdetector, 200, -300.0, 300.0, AliCMEVarManager::kEMCALclusterDz, "", "EMCAL;PHOS");
//    }  // end if "CaloClusters" histograms
//
//    TString trkFlagNames = "";
//    for(Int_t iflag=0; iflag<AliCMEVarManager::kNTrackingStatus; ++iflag) {
//      trkFlagNames += AliCMEVarManager::fgkTrackingStatusNames[iflag];
//      trkFlagNames += ";";
//    }
//    TString trackQualityFlagNames = "EP;#gamma;K^{0}_{S};#Lambda;#bar{#Lambda};kink0;kink1;kink2;";
//    trackQualityFlagNames += "pure #gamma;pure K^{0}_{S};pure #Lambda;pure #bar{#Lambda};-kink0;-kink1;-kink2;";
//    trackQualityFlagNames += "bayes e>0.5;bayes #pi>0.5;bayes K>0.5;bayes p>0.5;bayes>0.7;bayes>0.8;bayes>0.9";
//
//    if(classStr.Contains("TrackingFlags")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      man->AddHistogram(classStr.Data(), "TrackingFlags", "Tracking flags;;", kFALSE,
//	                AliCMEVarManager::kNTrackingFlags, -0.5, AliCMEVarManager::kNTrackingFlags-0.5, AliCMEVarManager::kTrackingFlag,
//			0, 0.0, 0.0, AliCMEVarManager::kNothing, 0, 0.0, 0.0, AliCMEVarManager::kNothing, trkFlagNames.Data());
//      continue;
//    }
//
//    if(classStr.Contains("TrackQualityFlags")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      man->AddHistogram(classStr.Data(), "TrackQualityFlags", "Track quality flags;;", kFALSE,
//	                64, -0.5, 63.5, AliCMEVarManager::kTrackQualityFlag, 0, 0.0, 0.0, AliCMEVarManager::kNothing,
//			0, 0.0, 0.0, AliCMEVarManager::kNothing, trackQualityFlagNames.Data());
//      continue;
//    }
//
//    // Track histograms
    if(classStr.Contains("TrackQA")) {
      man->AddHistClass(classStr.Data());
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution", kFALSE, 1000, 0.0, 50.0, AliCMEVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Eta", "", kFALSE, 1000, -1.5, 1.5, AliCMEVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Phi", "", kFALSE, 1000, 0.0, 6.3, AliCMEVarManager::kPhi);
//      man->AddHistogram(classStr.Data(), "PtTPC", "p_{T} (TPC) distribution", kFALSE, 1000, 0.0, 50.0, AliCMEVarManager::kPtTPC);
//      man->AddHistogram(classStr.Data(), "EtaTPC", "#eta (TPC)", kFALSE, 1000, -1.5, 1.5, AliCMEVarManager::kEtaTPC);
//      man->AddHistogram(classStr.Data(), "PhiTPC", "#varphi (TPC)", kFALSE, 1000, 0.0, 6.3, AliCMEVarManager::kPhiTPC);
      man->AddHistogram(classStr.Data(), "Pinner", "p inner param", kFALSE, 1000, 0.0, 50.0, AliCMEVarManager::kPIn);
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE, 1000, -10.0, 10.0, AliCMEVarManager::kImpactParXY);
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE, 1000, -10.0, 10.0, AliCMEVarManager::kImpactParZ);
      man->AddHistogram(classStr.Data(), "DCA2d", "DCA2d", kFALSE, 1000, 0.0, 10.0, AliCMEVarManager::kImpactPar2D);
//      man->AddHistogram(classStr.Data(), "DCAxyTPC", "DCAxyTPC", kFALSE, 1000, -10.0, 10.0, AliCMEVarManager::kDcaXYTPC);
//      man->AddHistogram(classStr.Data(), "DCAzTPC", "DCAzTPC", kFALSE, 1000, -10.0, 10.0, AliCMEVarManager::kDcaZTPC);
//      man->AddHistogram(classStr.Data(), "TrackLength", "Track length", kFALSE, 1000, 0.0, 1000.0, AliCMEVarManager::kTrackLength);
//      man->AddHistogram(classStr.Data(), "PtVsCentrality", "p_{T} distribution", kFALSE, 100, 0.0, 10.0, AliCMEVarManager::kPt, 20, 0.0, 100.0, AliCMEVarManager::kCentrality);
//      man->AddHistogram(classStr.Data(), "PtTPCVsCentrality", "p_{T} (TPC) distribution", kFALSE, 100, 0.0, 10.0, AliCMEVarManager::kPtTPC, 20, 0.0, 100.0, AliCMEVarManager::kCentrality);
//      if(classStr.Contains("ITS")) {
//        man->AddHistogram(classStr.Data(),"ITSncls", "ITS nclusters", kFALSE, 7,-0.5,6.5,AliCMEVarManager::kITSncls);
//        man->AddHistogram(classStr.Data(),"Eta_Phi_ITSncls_prof","ITS <nclusters> vs (#eta,#phi)",kTRUE,
//                     192, -1.2, 1.2, AliCMEVarManager::kEta, 126, 0.0, 6.3, AliCMEVarManager::kPhi, 7, -0.5, 6.5, AliCMEVarManager::kITSncls);
//	man->AddHistogram(classStr.Data(),"ITSsignal", "ITS dE/dx", kFALSE, 400,0.0,1000.0, AliCMEVarManager::kITSsignal);
//	man->AddHistogram(classStr.Data(),"Eta_Phi_ITSsignal_prof","ITS <dE/dx> vs (#eta,#phi)",kTRUE,
//                     192, -1.2, 1.2, AliCMEVarManager::kEta, 126, 0.0, 6.3, AliCMEVarManager::kPhi, 100, 0.0, 1000, AliCMEVarManager::kITSsignal);
//	man->AddHistogram(classStr.Data(),"ITSchi2", "ITS #chi^{2}", kFALSE, 200,0.0,20.0, AliCMEVarManager::kITSchi2);
//	man->AddHistogram(classStr.Data(),"Eta_Phi_ITSchi2_prof","ITS <#chi^{2}> vs (#eta,#phi)",kTRUE,
//                     192, -1.2, 1.2, AliCMEVarManager::kEta, 126, 0.0, 6.3, AliCMEVarManager::kPhi, 100, 0.0, 1000, AliCMEVarManager::kITSchi2);
//      }  // end if ITS histograms
     if(classStr.Contains("TPC")) {
       man->AddHistogram(classStr.Data(),"TPCncls","",kFALSE,160,-0.5,159.5,AliCMEVarManager::kNclsTPC);
//	man->AddHistogram(classStr.Data(),"TPCcrossedRows","", kFALSE, 160,-0.5,159.5,AliCMEVarManager::kTPCcrossedRows);
//	man->AddHistogram(classStr.Data(),"TPCnclsF","",kFALSE, 160,-0.5,159.5,AliCMEVarManager::kTPCnclsF);
//	man->AddHistogram(classStr.Data(),"TPCnclsShared","",kFALSE, 160,-0.5,159.5,AliCMEVarManager::kTPCnclsShared);
	man->AddHistogram(classStr.Data(),"TPCchi2","",kFALSE, 200,0.0,10.0,AliCMEVarManager::kTPCchi2Cl);
//        man->AddHistogram(classStr.Data(),"Eta_Phi_TPCncls_prof","TPC <nclusters> vs (#eta,#phi)",kTRUE,
//                     192, -1.2, 1.2, AliCMEVarManager::kEta, 126, 0.0, 6.3, AliCMEVarManager::kPhi, 160, -0.5, 159.5, AliCMEVarManager::kTPCncls);
//	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCcrossedRows_prof","TPC <n crossed rows> vs (#eta,#phi)",kTRUE,
//                     192, -1.2, 1.2, AliCMEVarManager::kEta, 126, 0.0, 6.3, AliCMEVarManager::kPhi, 160, -0.5, 159.5, AliCMEVarManager::kTPCcrossedRows);
//	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCnclsF_prof","TPC <nclusters findable> vs (#eta,#phi)",kTRUE,
//                     192, -1.2, 1.2, AliCMEVarManager::kEta, 126, 0.0, 6.3, AliCMEVarManager::kPhi, 160, -0.5, 159.5, AliCMEVarManager::kTPCnclsF);
//	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCchi2_prof","TPC <#chi^{2}> vs (#eta,#phi)",kTRUE,
//                     192, -1.2, 1.2, AliCMEVarManager::kEta, 126, 0.0, 6.3, AliCMEVarManager::kPhi, 160, -0.5, 159.5, AliCMEVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"TPCsignal_Pin","TPC dE/dx vs. inner param P",kFALSE,
                     400,0.0,4.0,AliCMEVarManager::kPIn,500,-0.5,499.5,AliCMEVarManager::kTPCsignal);
	man->AddHistogram(classStr.Data(),"TPCsignalN","",kFALSE,160,-0.5,159.5,AliCMEVarManager::kTPCsignalN);
//	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCsignalN_prof","TPC <nclusters pid> vs (#eta,#phi)",kTRUE,
//                     192, -1.2, 1.2, AliCMEVarManager::kEta, 126, 0.0, 6.3, AliCMEVarManager::kPhi, 160, -0.5, 159.5, AliCMEVarManager::kTPCsignalN);
//        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Pin","TPC N_{#sigma} electron vs. inner param P",kFALSE,
//                     200,0.0,20.0,AliCMEVarManager::kPin,100,-5.0,5.0,AliCMEVarManager::kTPCnSig+AliCMEVarManager::kElectron);
//	man->AddHistogram(classStr.Data(),"TPCnsigPion_Pin","TPC N_{#sigma} pion vs. inner param P",kFALSE,
//                     200,0.0,20.0,AliCMEVarManager::kPin,100,-5.0,5.0,AliCMEVarManager::kTPCnSig+AliCMEVarManager::kPion);
//        man->AddHistogram(classStr.Data(),"TPCnsigKaon_Pin","TPC N_{#sigma} kaon vs. inner param P",kFALSE,
//                     200,0.0,20.0,AliCMEVarManager::kPin,100,-5.0,5.0,AliCMEVarManager::kTPCnSig+AliCMEVarManager::kKaon);
//        man->AddHistogram(classStr.Data(),"TPCnsigProton_Pin","TPC N_{#sigma} proton vs. inner param P",kFALSE,
//                     200,0.0,20.0,AliCMEVarManager::kPin,100,-5.0,5.0,AliCMEVarManager::kTPCnSig+AliCMEVarManager::kProton);
//	man->AddHistogram(classStr.Data(),"TPCnsigElectron_Run","TPC N_{#sigma} electron vs. run",kTRUE,
//                     runNBins, runHistRange[0], runHistRange[1], AliCMEVarManager::kRunNo,100,-5.0,5.0,AliCMEVarManager::kTPCnSig+AliCMEVarManager::kElectron);
//	man->AddHistogram(classStr.Data(),"TPCnsigPion_Run","TPC N_{#sigma} pion vs. run",kTRUE,
//                     runNBins, runHistRange[0], runHistRange[1], AliCMEVarManager::kRunNo,100,-5.0,5.0,AliCMEVarManager::kTPCnSig+AliCMEVarManager::kPion);
//	man->AddHistogram(classStr.Data(),"TPCnsigKaon_Run","TPC N_{#sigma} kaon vs. run",kTRUE,
//                     runNBins, runHistRange[0], runHistRange[1], AliCMEVarManager::kRunNo,100,-5.0,5.0,AliCMEVarManager::kTPCnSig+AliCMEVarManager::kKaon);
//	man->AddHistogram(classStr.Data(),"TPCnsigProton_Run","TPC N_{#sigma} proton vs. run",kTRUE,
//                     runNBins, runHistRange[0], runHistRange[1], AliCMEVarManager::kRunNo,100,-5.0,5.0,AliCMEVarManager::kTPCnSig+AliCMEVarManager::kProton);
      }      // end if TPC histograms
//      if(classStr.Contains("TOF")) {
//	man->AddHistogram(classStr.Data(),"TOFdeltaBC","TOF #delta BC", kFALSE, 7000, -3500.,3500.,AliCMEVarManager::kTOFdeltaBC);
//	man->AddHistogram(classStr.Data(),"TOFtime","TOF time", kFALSE, 500, 0., 1.0e+5,AliCMEVarManager::kTOFtime);
//	man->AddHistogram(classStr.Data(),"TOFdx","TOF dx", kFALSE, 500, -5., 5., AliCMEVarManager::kTOFdx);
//	man->AddHistogram(classStr.Data(),"TOFdz","TOF dz", kFALSE, 500, -5., 5., AliCMEVarManager::kTOFdz);
//	man->AddHistogram(classStr.Data(),"TOFmismatchProbab","TOF mismatch probability", kFALSE, 120, -0.1, 1.1, AliCMEVarManager::kTOFmismatchProbability);
//	man->AddHistogram(classStr.Data(),"TOFchi2","TOF #chi^{2}", kFALSE, 500, -1., 10., AliCMEVarManager::kTOFchi2);
//        man->AddHistogram(classStr.Data(),"TOFbeta_P","TOF #beta vs P",kFALSE,
//                     200,0.0,20.0,AliCMEVarManager::kP, 220,0.0,1.1,AliCMEVarManager::kTOFbeta);
//        man->AddHistogram(classStr.Data(),"TOFnsigElectron_P","TOF N_{#sigma} electron vs. P",kFALSE,
//                     200,0.0,20.0,AliCMEVarManager::kP, 100,-5.0,5.0,AliCMEVarManager::kTOFnSig+AliCMEVarManager::kElectron);
//        man->AddHistogram(classStr.Data(),"TOFnsigPion_P","TOF N_{#sigma} pion vs. P",kFALSE,
//                     200,0.0,20.0,AliCMEVarManager::kP, 100,-5.0,5.0,AliCMEVarManager::kTOFnSig+AliCMEVarManager::kPion);
//        man->AddHistogram(classStr.Data(),"TOFnsigKaon_P","TOF N_{#sigma} kaon vs. P",kFALSE,
//                     200,0.0,20.0,AliCMEVarManager::kP, 100,-5.0,5.0,AliCMEVarManager::kTOFnSig+AliCMEVarManager::kKaon);
//        man->AddHistogram(classStr.Data(),"TOFnsigProton_P","TOF N_{#sigma} proton vs. P",kFALSE,
//                     200,0.0,20.0,AliCMEVarManager::kP, 100,-5.0,5.0,AliCMEVarManager::kTOFnSig+AliCMEVarManager::kProton);
//      }    // end if TOF histograms
//      if(classStr.Contains("TRD")) {
//        man->AddHistogram(classStr.Data(),"TRDntracklets","TRD ntracklets",kFALSE,
//                     7,-0.5,6.5,AliCMEVarManager::kTRDntracklets);
//	man->AddHistogram(classStr.Data(),"TRDntrackletsPID","TRD ntracklets PID",kFALSE,
//                     7,-0.5,6.5,AliCMEVarManager::kTRDntrackletsPID);
//	man->AddHistogram(classStr.Data(),"TRDprobabElectronLQ1D","TRD electron probability LQ1D",kFALSE,
//                     500,0.0,1.0,AliCMEVarManager::kTRDpidProbabilitiesLQ1D);
//	man->AddHistogram(classStr.Data(),"TRDprobabPionLQ1D","TRD pion probability LQ1D",kFALSE,
//                     500,0.0,1.0,AliCMEVarManager::kTRDpidProbabilitiesLQ1D+1);
//	man->AddHistogram(classStr.Data(),"TRDprobabElectronLQ2D","TRD electron probability LQ2D",kFALSE,
//                     500,0.0,1.0,AliCMEVarManager::kTRDpidProbabilitiesLQ2D);
//	man->AddHistogram(classStr.Data(),"TRDprobabPionLQ2D","TRD pion probability LQ2D",kFALSE,
//                     500,0.0,1.0,AliCMEVarManager::kTRDpidProbabilitiesLQ2D+1);
//	man->AddHistogram(classStr.Data(),"Eta_Phi_TRDntracklets_prof","TRD <ntracklets> vs (#eta,#phi)",kTRUE,
//                     192, -1.2, 1.2, AliCMEVarManager::kEta, 126, 0.0, 6.3, AliCMEVarManager::kPhi, 7, -0.5, 6.5, AliCMEVarManager::kTRDntracklets);
//        man->AddHistogram(classStr.Data(),"Eta_Phi_TRDntrackletsPID_prof","TRD <ntracklets PID> vs (#eta,#phi)",kTRUE,
//                     192, -1.2, 1.2, AliCMEVarManager::kEta, 126, 0.0, 6.3, AliCMEVarManager::kPhi, 7, -0.5, 6.5, AliCMEVarManager::kTRDntrackletsPID);
//      }   // end if TRD histograms
//      man->AddHistogram(classStr.Data(), "MatchedCaloClusterId", "EMCAL/PHOS matched cluster id", kFALSE,
//                   5000, 0.0, 5000.0, AliCMEVarManager::kEMCALmatchedClusterId);
//      if(classStr.Contains("EMCAL")) {
//	man->AddHistogram(classStr.Data(), "MatchedEnergy", "Energy from the calorimeter matched cluster", kFALSE,
//	             200, 0.0, 50.0, AliCMEVarManager::kEMCALmatchedEnergy);
//	man->AddHistogram(classStr.Data(), "MatchedEOverP", "E/P from the calorimeter matched cluster", kFALSE,
//	             300, 0.0, 1.5, AliCMEVarManager::kEMCALmatchedEOverP);
//	man->AddHistogram(classStr.Data(), "MatchedEOverP_P", "E/P from the calorimeter matched cluster vs P", kFALSE,
//	             10, 0.0, 10.0, AliCMEVarManager::kPin, 300, 0.0, 1.5, AliCMEVarManager::kEMCALmatchedEOverP);
//	man->AddHistogram(classStr.Data(),"MatchedEnergy_Run","Calorimeter matched energy vs run",kTRUE,
//                     runNBins, runHistRange[0], runHistRange[1], AliCMEVarManager::kRunNo, 100, 0.0, 50.0, AliCMEVarManager::kEMCALmatchedEnergy);
//	man->AddHistogram(classStr.Data(),"EMCALmatchedEOverP_Run","Calorimeter matched E/P vs run",kTRUE,
//                     runNBins, runHistRange[0], runHistRange[1], AliCMEVarManager::kRunNo, 100, 0.0, 2.0, AliCMEVarManager::kEMCALmatchedEOverP);
//	man->AddHistogram(classStr.Data(),"EMCALmatchedEnergy_P","Calorimeter matched energy vs momentum",kFALSE,
//                     100, 0.0, 10.0, AliCMEVarManager::kP, 100, 0.0, 10.0, AliCMEVarManager::kEMCALmatchedEnergy);
//	man->AddHistogram(classStr.Data(),"EMCALmatchedEnergy_CentVZERO","Calorimeter matched energy vs centrality VZERO",kTRUE,
//                     20, 0.0, 100.0, AliCMEVarManager::kCentrality, 100, 0.0, 10.0, AliCMEVarManager::kEMCALmatchedEnergy);
//	man->AddHistogram(classStr.Data(),"Eta_Phi_MatchedEnergy_prof","EMCAL matched cluster <energy> vs (#eta,#phi)",kTRUE,
//                     192, -1.2, 1.2, AliCMEVarManager::kEta, 126, 0.0, 6.3, AliCMEVarManager::kPhi, 100, 0.0, 10.0, AliCMEVarManager::kEMCALmatchedEnergy);
//      }  // end if EMCAL histograms
//      continue;
    }  // end if "TrackQA"
//
//    if(classStr.Contains("PairQualityFlags")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      TString pairQualityFlagNames = "";
//      man->AddHistogram(classStr.Data(), "PairQualityFlags", "Pair quality flags;;", kFALSE,
//	                32, -0.5, 31.5, AliCMEVarManager::kPairQualityFlag, 0, 0.0, 0.0, AliCMEVarManager::kNothing, 0, 0.0, 0.0, AliCMEVarManager::kNothing, pairQualityFlagNames.Data());
//      continue;
//    }
//
//    Double_t massBinWidth = 0.001;     // *GeV/c^2
//    Double_t massRange[2] = {0.0,5.0};
//    Int_t nMassBins = TMath::Nint((massRange[1]-massRange[0])/massBinWidth);
//
//    TString candidateNames = "#gamma#rightarrow e^{+}e^{-};K^{0}_{S}#rightarrow#pi^{+}#pi^{-};";
//    candidateNames += "#Lambda#rightarrow p#pi^{-};#bar{#Lambda}#rightarrow #bar{p}#pi^{+};";
//
//    // Histograms for pairs
//    if(classStr.Contains("Pair")) {
//      man->AddHistClass(classStr.Data());
//      cout << classStr.Data() << endl;
//      if(classStr.Contains("QA")) {
//	man->AddHistogram(classStr.Data(), "CandidateId", "Candidate id", kFALSE,
//                          5, -0.5, 4.5, AliCMEVarManager::kCandidateId, 0, 0.0, 0.0, AliCMEVarManager::kNothing, 0, 0.0, 0.0, AliCMEVarManager::kNothing, candidateNames.Data());
//        man->AddHistogram(classStr.Data(), "PairType", "Pair type", kFALSE, 4, -0.5, 3.5, AliCMEVarManager::kPairType);
//	man->AddHistogram(classStr.Data(), "PairChi2", "Pair #chi^{2}", kFALSE, 200, 0.0, 50, AliCMEVarManager::kPairChisquare);
//	man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE, nMassBins, massRange[0], massRange[1], AliCMEVarManager::kMass);
//        man->AddHistogram(classStr.Data(), "Mass_V0K0s", "Invariant mass, K^{0}_{s} assumption", kFALSE,
//                     nMassBins, massRange[0], massRange[1], AliCMEVarManager::kMassV0);
//        man->AddHistogram(classStr.Data(), "Mass_V0Lambda", "Invariant mass, #Lambda^{0} assumption", kFALSE,
//                     nMassBins, massRange[0], massRange[1], AliCMEVarManager::kMassV0+1);
//        man->AddHistogram(classStr.Data(), "Mass_V0ALambda", "Invariant mass, #bar{#Lambda^{0}} assumption", kFALSE,
//                     nMassBins, massRange[0], massRange[1], AliCMEVarManager::kMassV0+2);
//        man->AddHistogram(classStr.Data(), "Mass_V0Gamma", "Invariant mass, #gamma conversion assumption", kFALSE,
//                     nMassBins, massRange[0], massRange[1], AliCMEVarManager::kMassV0+3);
//        man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 1000, 0.0, 10.0, AliCMEVarManager::kPt);
//	man->AddHistogram(classStr.Data(), "Px", "", kFALSE, 1000, 0.0, 10.0, AliCMEVarManager::kPx);
//	man->AddHistogram(classStr.Data(), "Py", "", kFALSE, 1000, 0.0, 10.0, AliCMEVarManager::kPy);
//	man->AddHistogram(classStr.Data(), "Pz", "", kFALSE, 1000, 0.0, 10.0, AliCMEVarManager::kPz);
//	man->AddHistogram(classStr.Data(), "P", "", kFALSE, 1000, 0.0, 10.0, AliCMEVarManager::kP);
//	man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 240, -1.2, 1.2, AliCMEVarManager::kRap);
//	man->AddHistogram(classStr.Data(), "Eta", "Pseudo-rapidity #eta", kFALSE, 240, -2.0, 2.0, AliCMEVarManager::kEta);
//	man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 315, 0.0, 6.3, AliCMEVarManager::kPhi);
//	man->AddHistogram(classStr.Data(), "Theta", "#theta distribution", kFALSE, 1000, 0.0, 3.2, AliCMEVarManager::kTheta);
//	man->AddHistogram(classStr.Data(), "LxyOrR", "L_{xy}/Decay Radius", kFALSE, 1000, -10.0, 10.0, AliCMEVarManager::kPairLxy);
//	man->AddHistogram(classStr.Data(), "OpeningAngle", "Opening angle", kFALSE, 1000, 0.0, 3.2, AliCMEVarManager::kPairOpeningAngle);
//        man->AddHistogram(classStr.Data(), "PointingAngle", "Pointing angle", kFALSE, 1000, 0.0, 3.2, AliCMEVarManager::kPairPointingAngle);
//      }   // end if "QA"
//
//      continue;
//    }   // end if for Pair classes of histograms

  }  // end loop over histogram classes
}
