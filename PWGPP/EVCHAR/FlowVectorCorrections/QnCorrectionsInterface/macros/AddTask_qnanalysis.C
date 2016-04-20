// qn analysis example task
// author: Jaap Onderwaater, jacobus.onderwaater@cern.ch
//         2016/Apr/20


void DefineHistogramsQnAnalysis(AliQnCorrectionsHistos* histos, TString histClass);


AliAnalysisTask* AddTask_qnanalysis() {
  /* temporal flag to use multiplicity instead of centrality and to inhibit detectors for 2015 dataset */
  Bool_t bUseMultiplicity = kFALSE;
  Bool_t b2015DataSet = kFALSE;

  gROOT->LoadMacro("AddTask_ep.C");
  AliAnalysisDataContainer * qvectors = (AliAnalysisDataContainer*) AddTask_ep();

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    mgr = new AliAnalysisManager("AnalysisManagerQnAnalysis");

    //Error("AddTask_jonderw_qnanalysis", "No analysis manager found.");
    //return 0;
  }

  AliAnalysisTaskQnAnalysis* taskQn = new AliAnalysisTaskQnAnalysis("QnAnalysis");

  AliQnCorrectionsCuts *eventCuts = new AliQnCorrectionsCuts();
  eventCuts->AddCut(AliQnCorrectionsVarManager::kVtxZ,-10.0,10.0);
  if (bUseMultiplicity) {
    eventCuts->AddCut(AliQnCorrectionsVarManager::kVZEROMultPercentile,0.0,90.0);
  }
  else {
    eventCuts->AddCut(AliQnCorrectionsVarManager::kCentVZERO,0.0,80.0);
  }



  taskQn->SetEventCuts(eventCuts);
  if (!b2015DataSet) {
    taskQn->SelectCollisionCandidates(AliVEvent::kMB);  // Events passing trigger and physics selection for analysis
  }
  else
    taskQn->SelectCollisionCandidates(AliVEvent::kMB|AliVEvent::kINT7);  // Events passing trigger and physics selection for analysis


  AliQnCorrectionsHistos* hists = taskQn->GetHistograms();
  TString histClass = "";
  histClass += "Event_NoCuts;";
  histClass += "Event_Analysis;";
  DefineHistogramsQnAnalysis(hists, histClass);

  mgr->AddTask(taskQn);
  

  //create output container
  AliAnalysisDataContainer *cOutputQnAna=
    mgr->CreateContainer("QnAnalysis",
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        "QnAnalysis.root");

  mgr->ConnectInput(taskQn,  0, mgr->GetCommonInputContainer());
  mgr->ConnectInput(taskQn,  1, qvectors);
  mgr->ConnectOutput(taskQn, 1, cOutputQnAna);

  return taskQn;

}




//__________________________________________________________________
void DefineHistogramsQnAnalysis(AliQnCorrectionsHistos* histos, TString histClass) {
  //
  // define the histograms
  //
  //#define VAR AliQnCorrectionsVarManager

  const Char_t* histClasses = histClass.Data();

  cout << "Defining histograms ..." << endl;
  cout << "histogram classes: " << histClass<< endl;

  //fHistosFile=new TFile(output,"RECREATE");

  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");

  const Int_t kNRunBins = 3000;
  Double_t runHistRange[2] = {137000.,140000.};

  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    cout << "hist class: " << classStr.Data() << endl;

    // Event wise histograms
    if(classStr.Contains("Event")) {
      histos->AddHistClass(classStr.Data());
      histos->AddHistogram(classStr.Data(),"RunNo","Run numbers;Run", kFALSE, kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo);
      histos->AddHistogram(classStr.Data(),"BC","Bunch crossing;BC", kFALSE,3000,0.,3000.,AliQnCorrectionsVarManager::kBC);
      histos->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag;;", kFALSE,
          2,-0.5,1.5,AliQnCorrectionsVarManager::kIsPhysicsSelection, 0,0.0,0.0,AliQnCorrectionsVarManager::kNothing, 0,0.0,0.0,AliQnCorrectionsVarManager::kNothing, "off;on");

      histos->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)", kFALSE,300,0.0,0.0,AliQnCorrectionsVarManager::kVtxZ);
      //histos->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)", kFALSE,300,-15.,15.,AliQnCorrectionsVarManager::kVtxZ);
      histos->AddHistogram(classStr.Data(),"VtxX","Vtx X;vtx X (cm)", kFALSE,300,-1.,1.,AliQnCorrectionsVarManager::kVtxX);
      histos->AddHistogram(classStr.Data(),"VtxY","Vtx Y;vtx Y (cm)", kFALSE,300,-1.,1.,AliQnCorrectionsVarManager::kVtxY);


      histos->AddHistogram(classStr.Data(),"CentVZEROvsMultPVZERO","Multiplicity percentile (VZERO);multiplicity VZERO (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kVZEROMultPercentile, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentVZEROvsCentSPD","Centrality(VZERO);centrality VZERO (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentSPD","Centrality(TPC);centrality TPC (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentTPC, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentVZERO","Centrality(TPC);centrality TPC (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentTPC, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentZDC","Centrality(TPC);centrality TPC (percents);centrality ZDC (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentTPC, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentZDC);
      histos->AddHistogram(classStr.Data(),"CentZDCvsCentVZERO","Centrality(ZDC);centrality ZDC (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentZDC, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentZDCvsCentSPD","Centrality(ZDC);centrality ZDC (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentZDC, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentSPD);

      histos->AddHistogram(classStr.Data(),"MultVZEROvsCentVZERO","Multiplicity;multiplicity VZERO;VZERO centrality", kFALSE,
          100, 0.0, 32000.0, AliQnCorrectionsVarManager::kVZEROTotalMult, 100,0.,100., AliQnCorrectionsVarManager::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"MultSPDvsCentPSD","Multiplicity;SPD tracklets;SPD centrality", kFALSE,
          100, 0.0, 3000.0, AliQnCorrectionsVarManager::kSPDntracklets, 100,0.,100., AliQnCorrectionsVarManager::kCentSPD);
      histos->AddHistogram(classStr.Data(),"MultTPCvsCentSPD","Multiplicity;TPC selected tracks;TPC centrality", kFALSE,
          100, 0.0, 3500.0, AliQnCorrectionsVarManager::kNtracksSelected, 100,0.,100., AliQnCorrectionsVarManager::kCentTPC);
      histos->AddHistogram(classStr.Data(),"MultZDCvsCentZDC","Multiplicity;multiplicity ZDC;ZDC centrality", kFALSE,
          100, 0.0, 300000.0, AliQnCorrectionsVarManager::kZDCTotalEnergy, 100,0.,100., AliQnCorrectionsVarManager::kCentZDC);


      histos->AddHistogram(classStr.Data(),"MultTPCvsMultVZERO","Multiplicity;tracks TPC;multiplicity VZERO", kFALSE,
          100, 0.0, 3500.0, AliQnCorrectionsVarManager::kNtracksSelected, 100, 0.0, 32000.0, AliQnCorrectionsVarManager::kVZEROTotalMult);
      histos->AddHistogram(classStr.Data(),"MultTPCvsMultSPD","Multiplicity;tracklets SPD;tracks TPC", kFALSE,
          100, 0.0, 3500.0, AliQnCorrectionsVarManager::kNtracksSelected, 100, 0.0, 3000.0, AliQnCorrectionsVarManager::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"MultSPDvsMultVZERO","Multiplicity;tracklets SPD;multiplicity VZERO", kFALSE,
          100, 0.0, 32000.0, AliQnCorrectionsVarManager::kVZEROTotalMult, 100, 0.0, 3000.0, AliQnCorrectionsVarManager::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"MultTPCvsMultZDC","Multiplicity;tracks TPC;energy ZDC", kFALSE,
          100, 0.0, 3500.0, AliQnCorrectionsVarManager::kNtracksSelected, 100, 0.0, 300000.0, AliQnCorrectionsVarManager::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultVZEROvsMultZDC","Multiplicity;multiplicity VZERO;energy ZDC", kFALSE,
          100, 0.0, 32000.0, AliQnCorrectionsVarManager::kVZEROTotalMult, 100, 0.0, 300000.0, AliQnCorrectionsVarManager::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultSPDvsMultZDC","Multiplicity;tracklets SPD;energy ZDC", kFALSE,
          100, 0.0, 3000.0, AliQnCorrectionsVarManager::kSPDntracklets, 100, 0.0, 300000.0, AliQnCorrectionsVarManager::kZDCTotalEnergy);



      histos->AddHistogram(classStr.Data(),"MultVZERO","Multiplicity;multiplicity VZERO", kFALSE,
          320, 0.0, 25000.0, AliQnCorrectionsVarManager::kVZEROTotalMult);
      histos->AddHistogram(classStr.Data(),"MultVZEROA","Multiplicity;multiplicity VZEROA", kFALSE,
          250, 0.0, 0.0, AliQnCorrectionsVarManager::kVZEROATotalMult);//10000.0
      histos->AddHistogram(classStr.Data(),"MultVZEROC","Multiplicity;multiplicity VZEROC", kFALSE,
          250, 0.0, 0.0, AliQnCorrectionsVarManager::kVZEROCTotalMult);//15000.0
      histos->AddHistogram(classStr.Data(),"MultZDC","Multiplicity;multiplicity ZDC", kFALSE,
          200, 0.0, 300000.0, AliQnCorrectionsVarManager::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultZDCA","Multiplicity;multiplicity ZDCA", kFALSE,
          200, 0.0, 150000.0, AliQnCorrectionsVarManager::kZDCATotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultZDCC","Multiplicity;multiplicity ZDCC", kFALSE,
          200, 0.0, 150000.0, AliQnCorrectionsVarManager::kZDCCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultFMD1","Multiplicity;multiplicity FMD1", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManager::kFMD1TotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD2I","Multiplicity;multiplicity FMD2I", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManager::kFMD2ITotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD2O","Multiplicity;multiplicity FMD2O", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManager::kFMD2OTotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD3I","Multiplicity;multiplicity FMD3I", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManager::kFMD3ITotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD3O","Multiplicity;multiplicity FMD3O", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManager::kFMD3OTotalMult);
      histos->AddHistogram(classStr.Data(),"MultTZEROA","Multiplicity;multiplicity TZEROA", kFALSE,
          300, 0.0, 3000.0, AliQnCorrectionsVarManager::kTZEROATotalMult);
      histos->AddHistogram(classStr.Data(),"MultTZEROC","Multiplicity;multiplicity TZEROC", kFALSE,
          300, 0.0, 3000.0, AliQnCorrectionsVarManager::kTZEROCTotalMult);




      histos->AddHistogram(classStr.Data(),"MultPercentVZERO","Multiplicity percentile (VZERO);multiplicity VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kVZEROMultPercentile);
      histos->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC);centrality TPC (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentTPC);
      histos->AddHistogram(classStr.Data(),"CentZDC","Centrality(ZDC);centrality ZDC (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentZDC);

      histos->AddHistogram(classStr.Data(),"CentQuality","Centrality quality;centrality quality", kFALSE,
          100, -50.5, 49.5, AliQnCorrectionsVarManager::kCentQuality);
      histos->AddHistogram(classStr.Data(),"CentVZERO_Run_prof","<Centrality(VZERO)> vs run;Run; centrality VZERO (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentSPD_Run_prof","<Centrality(SPD)> vs run;Run; centrality SPD (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPC_Run_prof","<Centrality(TPC)> vs run;Run; centrality TPC (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentTPC);
      histos->AddHistogram(classStr.Data(),"CentZDC_Run_prof","<Centrality(ZDC)> vs run;Run; centrality ZDC (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentZDC);


      histos->AddHistogram(classStr.Data(),"NV0sTotal","Number of V0 candidates per event;# pairs", kFALSE,
          1000,0.,30000.,AliQnCorrectionsVarManager::kNV0total);
      histos->AddHistogram(classStr.Data(),"NV0sSelected","Number of selected V0 candidates per event;# pairs", kFALSE,
          1000,0.,10000.,AliQnCorrectionsVarManager::kNV0selected);
      histos->AddHistogram(classStr.Data(),"NPairs","Number of candidates per event;# pairs", kFALSE,
          5000,0.,5000.,AliQnCorrectionsVarManager::kNdielectrons);
      histos->AddHistogram(classStr.Data(),"NPairsSelected", "Number of selected pairs per event; #pairs", kFALSE,
          5000,0.,5000.,AliQnCorrectionsVarManager::kNpairsSelected);
      histos->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event;# tracks", kFALSE,
          1000,0.,30000.,AliQnCorrectionsVarManager::kNtracksTotal);
      histos->AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event;# tracks", kFALSE,
          1000,0.,30000.,AliQnCorrectionsVarManager::kNtracksSelected);
      histos->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets; tracklets", kFALSE,
          3000, -0.5, 2999.5, AliQnCorrectionsVarManager::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"SPDnSingleClusters", "SPD #single clusters; tracklets", kFALSE,
          3000, -0.5, 2999.5, AliQnCorrectionsVarManager::kSPDnSingleClusters);

      histos->AddHistogram(classStr.Data(),"NV0total_Run_prof", "<Number of total V0s> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManager::kNV0total);
      histos->AddHistogram(classStr.Data(),"NV0selected_Run_prof", "<Number of selected V0s> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManager::kNV0selected);
      histos->AddHistogram(classStr.Data(),"Ndielectrons_Run_prof", "<Number of dielectrons> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManager::kNdielectrons);
      histos->AddHistogram(classStr.Data(),"NpairsSelected_Run_prof", "<Number of selected pairs> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManager::kNpairsSelected);
      histos->AddHistogram(classStr.Data(),"NTracksTotal_Run_prof", "<Number of tracks> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManager::kNtracksTotal);
      histos->AddHistogram(classStr.Data(),"NTracksSelected_Run_prof", "<Number of selected tracks> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManager::kNtracksSelected);
      histos->AddHistogram(classStr.Data(),"SPDntracklets_Run_prof", "<SPD ntracklets> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManager::kSPDntracklets);

      histos->AddHistogram(classStr.Data(),"VtxZ_CentVZERO","Centrality(VZERO) vs vtx. Z;vtx Z (cm); centrality VZERO (%)", kFALSE,
          300,-15.,15.,AliQnCorrectionsVarManager::kVtxZ, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"VtxZ_CentSPD","Centrality(SPD) vs vtx. Z;vtx Z (cm); centrality SPD (%)", kFALSE,
          300,-15.,15.,AliQnCorrectionsVarManager::kVtxZ, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentSPD);
      histos->AddHistogram(classStr.Data(),"VtxZ_CentTPC","Centrality(TPC) vs vtx. Z;vtx Z (cm); centrality TPC (%)", kFALSE,
          300,-15.,15.,AliQnCorrectionsVarManager::kVtxZ, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentTPC);
      continue;
    }  // end if className contains "Event"    


    // Offline trigger histograms
    if(classStr.Contains("OfflineTriggers")) {
      histos->AddHistClass(classStr.Data());

      TString triggerNames = "";
      for(Int_t i=0; i<64; ++i) {triggerNames += Form("%s",AliQnCorrectionsVarManager::fOfflineTriggerNames[i]); triggerNames+=";";}

      histos->AddHistogram(classStr.Data(), "Triggers", "Offline triggers fired; ; ;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManager::kOfflineTrigger, 2, -0.5, 1.5, AliQnCorrectionsVarManager::kOfflineTriggerFired, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, triggerNames.Data(), "off;on");
      histos->AddHistogram(classStr.Data(), "Triggers2", "Offline triggers fired; ; ;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManager::kOfflineTriggerFired2, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "Offline triggers fired vs centrality VZERO; ; centrality VZERO;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManager::kOfflineTriggerFired2, 20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentTPC_Triggers2", "Offline triggers fired vs centrality TPC; ; centrality TPC;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManager::kOfflineTriggerFired2, 20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentTPC, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentSPD_Triggers2", "Offline triggers fired vs centrality SPD; ; centrality SPD;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManager::kOfflineTriggerFired2, 20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentSPD, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentZDC_Triggers2", "Offline triggers fired vs centrality ZDC; ; centrality ZDC;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManager::kOfflineTriggerFired2, 20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentZDC, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "VtxZ_Triggers2", "Offline triggers fired vs vtxZ; ; vtx Z (cm.);", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManager::kOfflineTriggerFired2, 200, -20.0, 20.0, AliQnCorrectionsVarManager::kVtxZ, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, triggerNames.Data());
      continue;
    }

    // Track histograms
    if(classStr.Contains("Tracks")) {
      histos->AddHistClass(classStr.Data());
      for(Int_t ih=0; ih<6; ++ih) {
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
            20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 500, -1., 1., AliQnCorrectionsVarManager::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
            20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 500, -1.0, 1.0, AliQnCorrectionsVarManager::kSinNPhi+ih);
      }
    }


    // Track histograms
    if(classStr.Contains("TrackQA")) {
      histos->AddHistClass(classStr.Data());

      histos->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution; p_{T} (GeV/c^{2});", kFALSE,
          1000, 0.0, 50.0, AliQnCorrectionsVarManager::kPt);
      histos->AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE,
          1000, -1.5, 1.5, AliQnCorrectionsVarManager::kEta);
      histos->AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE,
          1000, 0.0, 6.3, AliQnCorrectionsVarManager::kPhi);
      histos->AddHistogram(classStr.Data(), "DCAxy", "DCAxy; DCAxy (cm.)", kFALSE,
          1000, -10.0, 10.0, AliQnCorrectionsVarManager::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz", "DCAz; DCAz (cm.)", kFALSE,
          1000, -10.0, 10.0, AliQnCorrectionsVarManager::kDcaZ);
      histos->AddHistogram(classStr.Data(), "TPCncls", "TPCncls; TPCncls", kFALSE,
          160, 0.0, 160.0, AliQnCorrectionsVarManager::kTPCncls);

      // run dependence
      histos->AddHistogram(classStr.Data(), "Pt_Run", "<p_{T}> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 1000, 0.0, 50.0, AliQnCorrectionsVarManager::kPt);
      histos->AddHistogram(classStr.Data(), "Eta_Run", "<#eta> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 1000, -1.5, 1.5, AliQnCorrectionsVarManager::kEta);      
      histos->AddHistogram(classStr.Data(), "Phi_Run", "<#varphi> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 1000, 0.0, 6.3, AliQnCorrectionsVarManager::kPhi);      
      histos->AddHistogram(classStr.Data(), "DCAxy_Run", "<DCAxy> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 1000, -10.0, 10.0, AliQnCorrectionsVarManager::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz_Run", "<DCAz> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 1000, -10.0, 10.0, AliQnCorrectionsVarManager::kDcaZ);

      // correlations between parameters
      histos->AddHistogram(classStr.Data(), "Eta_Pt_prof", "<p_{T}> vs #eta; #eta; p_{T} (GeV/c);", kTRUE,
          300, -1.5, +1.5, AliQnCorrectionsVarManager::kEta, 100, 0.0, 10.0, AliQnCorrectionsVarManager::kPt);
      histos->AddHistogram(classStr.Data(), "Phi_Pt", "p_{T} vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kFALSE,
          300, -0.01, 6.3, AliQnCorrectionsVarManager::kPhi, 100, 0.0, 2.2, AliQnCorrectionsVarManager::kPt);
      histos->AddHistogram(classStr.Data(), "Phi_Pt_prof", "<p_{T}> vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kTRUE,
          300, 0.0, 6.3, AliQnCorrectionsVarManager::kPhi, 100, 0.0, 10.0, AliQnCorrectionsVarManager::kPt);
      histos->AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE,
          200, -1.0, +1.0, AliQnCorrectionsVarManager::kEta, 100, 0.0, 6.3, AliQnCorrectionsVarManager::kPhi);
      histos->AddHistogram(classStr.Data(), "TPCncls_Eta_Phi_prof", "<TPC ncls> vs #varphi vs #eta; #eta; #varphi (rad.);TPC ncls", kTRUE,
          200, -1.0, +1.0, AliQnCorrectionsVarManager::kEta, 100, 0.0, 6.3, AliQnCorrectionsVarManager::kPhi, 10, 0.0, 200., AliQnCorrectionsVarManager::kTPCncls);
      histos->AddHistogram(classStr.Data(), "DCAxy_Eta_Phi_prof", "<DCAxy> vs #varphi vs #eta; #eta; #varphi (rad.);DCAxy (cm)", kTRUE,
          200, -1.0, +1.0, AliQnCorrectionsVarManager::kEta, 100, 0.0, 6.3, AliQnCorrectionsVarManager::kPhi, 10, 0.0, 200., AliQnCorrectionsVarManager::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz_Eta_Phi_prof", "<DCAz> vs #varphi vs #eta; #eta; #varphi (rad.);DCAz (cm)", kTRUE,
          200, -1.0, +1.0, AliQnCorrectionsVarManager::kEta, 100, 0.0, 6.3, AliQnCorrectionsVarManager::kPhi, 10, 0.0, 200., AliQnCorrectionsVarManager::kDcaZ);
      histos->AddHistogram(classStr.Data(), "Pt_DCAxy", "DCAxy vs p_{T}; p_{T} (GeV/c); DCA_{xy} (cm)", kFALSE,
          100, 0.0, 10.0, AliQnCorrectionsVarManager::kPt, 500, -2.0, 2.0, AliQnCorrectionsVarManager::kDcaXY);
      histos->AddHistogram(classStr.Data(), "Pt_DCAz", "DCAz vs p_{T}; p_{T} (GeV/c); DCA_{z} (cm)", kFALSE,
          100, 0.0, 10.0, AliQnCorrectionsVarManager::kPt, 500, -2.0, 2.0, AliQnCorrectionsVarManager::kDcaZ);
      histos->AddHistogram(classStr.Data(), "Eta_DCAxy", "DCAxy vs #eta; #eta; DCA_{xy} (cm)", kFALSE,
          100, -1.0, 1.0, AliQnCorrectionsVarManager::kEta, 500, -2.0, 2.0, AliQnCorrectionsVarManager::kDcaXY);
      histos->AddHistogram(classStr.Data(), "Eta_DCAz", "DCAz vs #eta; #eta; DCA_{z} (cm)", kFALSE,
          100, -1.0, 1.0, AliQnCorrectionsVarManager::kEta, 500, -2.0, 2.0, AliQnCorrectionsVarManager::kDcaZ);

      for(Int_t ih=0; ih<6; ++ih) {
        //histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
        //                   20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 500, -1., 1., AliQnCorrectionsVarManager::kCosNPhi+ih);
        //histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
        //                   20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 500, -1.0, 1.0, AliQnCorrectionsVarManager::kSinNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_Pt_Eta",ih+1), Form("<cos%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <cos%d #varphi>", ih+1, ih+1), kTRUE,
            20, -1.0, 1.0, AliQnCorrectionsVarManager::kEta, 30, 0.0, 3.0, AliQnCorrectionsVarManager::kPt, 500, -1.0, 1.0, AliQnCorrectionsVarManager::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_Pt_Eta",ih+1), Form("<sin%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <sin%d #varphi>", ih+1, ih+1), kTRUE,
            20, -1.0, 1.0, AliQnCorrectionsVarManager::kEta, 30, 0.0, 3.0, AliQnCorrectionsVarManager::kPt, 500, -1.0, 1.0, AliQnCorrectionsVarManager::kSinNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO_VtxZ",ih+1), Form("<cos%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
            30, -15.0, 15.0, AliQnCorrectionsVarManager::kVtxZ, 20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 500, -1., 1., AliQnCorrectionsVarManager::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO_VtxZ",ih+1), Form("<sin%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
            30, -15.0, 15.0, AliQnCorrectionsVarManager::kVtxZ, 20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 500, -1.0, 1.0, AliQnCorrectionsVarManager::kSinNPhi+ih);
      }
    }

    // Tracklet histograms
    if(classStr.Contains("TrackletQA")) {
      histos->AddHistClass(classStr.Data());

      histos->AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE,
          1000, -3.0, 3.0, AliQnCorrectionsVarManager::kSPDtrackletEta);
      histos->AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE,
          300, -0.01, 6.3, AliQnCorrectionsVarManager::kSPDtrackletPhi);
      histos->AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE,
          200, -3.0, +3.0, AliQnCorrectionsVarManager::kSPDtrackletEta, 100, 0.0, 6.3, AliQnCorrectionsVarManager::kSPDtrackletPhi);
    }

  }

  cout << " done" << endl;
}







