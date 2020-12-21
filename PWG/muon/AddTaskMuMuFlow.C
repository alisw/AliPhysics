
///
/// Configuration example of a task to use v2 task for dimuons analysis with MuMu
///
/// \author: A. Francisco (Subatech)
///

AliAnalysisTask* AddTaskMumuFlow(const char* outputname,
                             const char* beamYear,
                             Bool_t simulations,
                             AliAnalysisDataContainer* qntask = 0x0,
                             TString mapESE= "",
                             Bool_t raa = kFALSE)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager(); // Get the manager

  if (!mgr) {
    ::Error("AddTaskMuMu", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMuMu", "This task requires an input event handler");
    return NULL;
  }

  // AccEff correction map
  //===========================================================================
  TFile *accEff = new TFile("AccEffi2D_Cent20to40_03.root");
  TH2 *map = (TH2*)accEff->Get("histoAccEffi2D");

  // q2 map
  //===========================================================================
  TList *q2histo = 0x0;
  if(mapESE!="") {
    TFile* q2Map = TFile::Open(mapESE,"READ");
    if(!q2Map) {
      cout << "ERROR: map q2 for ESE: file not found!" << endl;
    }
    else{
    std::cout << "ese :" << mapESE << std::endl;
    q2histo = new TList();
    q2histo->Add((TH1*)q2Map->Get("smallq2"));
    q2histo->Add((TH1*)q2Map->Get("largeq2"));
    q2histo->Print();
    }
  }
  else std::cout << "no ESE map "<< std::endl;
  // Triggers
  //===========================================================================
   TList* triggerClassesToConsider = new TList; // Create pointer for trigger list
   triggerClassesToConsider->SetOwner(kTRUE); // Give rights to trigger liser
  // if (!isMC) triggers->Add(new TObjString("CINT7-B-NOPF-MUFAST"));
  if(!simulations) triggerClassesToConsider->Add(new TObjString("CMUL7-B-NOPF-MUFAST||CMLL7-B-NOPF-MUFAST"));

  // Add new trigger in case of Simulation
  //===========================================================================
  if (simulations && triggerClassesToConsider )
  {
      triggerClassesToConsider->Add(new TObjString("ANY"));
  }
//
// for dpmjet simulations (at least) we have the following "triggers" :
//    C0T0A,C0T0C,MB1,MBBG1,V0L,V0R,MULow,EMPTY,MBBG3,MULL,MULU,MUHigh
  //}
  //===========================================================================

  //  Configure inputmaps (Default on is in AliMuonEventCuts)
  //===========================================================================
  TList* triggerInputsMap = new TList();
  triggerInputsMap->SetOwner(kTRUE);
  //value = logbook -1
  triggerInputsMap->Add(new TObjString("0MSL:5"));
  triggerInputsMap->Add(new TObjString("0MSH:6"));
  triggerInputsMap->Add(new TObjString("0MUL:13"));
  triggerInputsMap->Add(new TObjString("0MUH:14"));
  triggerInputsMap->Add(new TObjString("0MLL:15"));
  triggerInputsMap->Add(new TObjString("0MLH:16"));

  //===========================================================================
  //   //**** Qncorr framework *****//
  // gROOT->LoadMacro("AliFlowQnVectorCorrections.C");
  // AliAnalysisDataContainer *corrTask = AliFlowQnVectorCorrections();

  //===========================================================================
  //  Configure task
  //===========================================================================
  AliAnalysisTaskMuMu       * task = new AliAnalysisTaskMuMu; // Call the task

  // task->UseLegacyCentrality(); // run 1
  AliAnalysisMuMuEventCutter* eventCutter = new AliAnalysisMuMuEventCutter(triggerClassesToConsider,triggerInputsMap); // To handle cuts on event
  AliAnalysisMuMuCutRegistry* cr = task->CutRegistry(); // Set CutRegistry

  // Default cuts on trigger and event level
  AliAnalysisMuMuCutElement * eventTrue = cr->AddEventCut(*eventCutter,"IsTrue","const AliVEvent&","");
  AliAnalysisMuMuCutElement * triggerSelection = cr->AddTriggerClassCut(*eventCutter,"SelectTriggerClass","const TString&,TString&,UInt_t,UInt_t,UInt_t","");
  AliAnalysisMuMuCutElement * ps2 = eventTrue;

  if (!simulations)
  {
    // ps = cr->AddEventCut(*eventCutter,"IsPhysicsSelectedMUL","const AliInputEventHandler&","");
    ps2 = cr->AddEventCut(*eventCutter,"IsPhysicsSelectedMULORMLL","const AliInputEventHandler&","");
    // ps2 = cr->AddEventCut(*eventCutter,"IsPhysicsSelectedINT7","const AliInputEventHandler&","");
  }

  AliAnalysisMuMuCutElement* zvtx = cr->AddEventCut(*eventCutter,"IsSPDzVertexInRange","AliVEvent&, const Double_t&, const Double_t&","-10.,10.");

  // Apply default cut
  // cr->AddCutCombination(triggerSelection);
  // cr->AddCutCombination(eventTrue);
  // cr->AddCutCombination(ps,triggerSelection);
  // cr->AddCutCombination(ps,triggerSelection,zvtx);
  cr->AddCutCombination(ps2,triggerSelection,zvtx);
  // cr->AddCutCombination(ps1,triggerSelection);
  // cr->AddCutCombination(ps2,triggerSelection);
  // if(qntask){
    // cr->AddCutCombination(ps, zvtx);
  // }
  task->SetBeamYear(beamYear);
  // TH2* AccEffHisto;

  AliAnalysisMuMuGlobal* globalAnalysis = 0x0; //new AliAnalysisMuMuGlobal; // Basic histograms analysis;
  AliAnalysisMuMuSingle* singleAnalysis = new AliAnalysisMuMuSingle;// Analysis dealing with single muon
  AliAnalysisMuMuMinv  * minvAnalysis = new AliAnalysisMuMuMinv(map);// Analysis creating invariant mass spectrum
  AliAnalysisMuMuFlow  * flowAnalysis = new AliAnalysisMuMuFlow(map,q2histo);// Analysis creating invariant mass spectrum
  // AliAnalysisMuMuFlowSP  * flowAnalysisSP = new AliAnalysisMuMuFlowSP(map);// Analysis creating invariant mass spectrum
  // AliAnalysisMuMuFlow* flowAnalysis = 0x0;
  AliAnalysisMuMuNch   * nchAnalysis = 0x0; //new AliAnalysisMuMuNch;
  // if(qntask) minvAnalysis->FillFlowHisto();
  // Configure sub analysis
  //===========================================================================
  if ( globalAnalysis )
  {
    // Cuts on trigger level
    AliAnalysisMuMuCutElement* triggerAll = cr->AddTriggerClassCut(*globalAnalysis,"SelectAnyTriggerClass","const TString&,TString&","");
    // Adding this cut on trigger level
    cr->AddCutCombination(triggerAll);
    task->AdoptSubAnalysis(globalAnalysis);
  }
  if ( nchAnalysis )
  {
    task->AdoptSubAnalysis(nchAnalysis);
  }
  // if ( flowAnalysis )
  // {
  //   task->AdoptSubAnalysis(flowAnalysis);
  // }

  if ( singleAnalysis )
  {
    // Cuts on tracks
    AliAnalysisMuMuCutElement* trackTrue = cr->AddTrackCut(*cr,"AlwaysTrue","const AliVParticle&",""); // Apply "AlwaysTrue" cut on AliVParticle derived from AliAnalysisMuMuSingle
    AliAnalysisMuMuCutElement* rabs = cr->AddTrackCut(*singleAnalysis,"IsRabsOK","const AliVParticle&","");
    AliAnalysisMuMuCutElement* matchlow = cr->AddTrackCut(*singleAnalysis,"IsMatchingTriggerLowPt","const AliVParticle&","");
    AliAnalysisMuMuCutElement* eta = cr->AddTrackCut(*singleAnalysis,"IsEtaInRange","const AliVParticle&","");
    // AliAnalysisMuMuCutElement* pdca = cr->AddTrackCut(*singleAnalysis,"IsPDCAOK","const AliVParticle&","");


    if ( minvAnalysis )
    {
      minvAnalysis->FillMeanPtHisto();
      // Array of cut elements
      TObjArray cutElements;
      TObjArray cutElements_sq2;
      TObjArray cutElements_Lq2;

      // Cuts on track level
      AliAnalysisMuMuCutElement* pairTrue = cr->AddTrackPairCut(*cr,"AlwaysTrue","const AliVParticle&,const AliVParticle&","");// Apply "AlwaysTrue" cut on AliVParticle derived from AliAnalysisMuMuMinv
      AliAnalysisMuMuCutElement* pairy = cr->AddTrackPairCut(*minvAnalysis,"IsRapidityInRange","const AliVParticle&,const AliVParticle&","");
      AliAnalysisMuMuCutElement* ptRange = cr->AddTrackPairCut(*minvAnalysis,"IsPtInRange","const AliVParticle&,const AliVParticle&,Double_t&,Double_t&","0.,12.");

      if(raa){
        AliAnalysisMuMuCutElement* inPlane = cr->AddTrackPairCut(*flowAnalysis,"IsDPhiInPlane","const AliVParticle&,const AliVParticle&","");
        AliAnalysisMuMuCutElement* outofPlane = cr->AddTrackPairCut(*flowAnalysis,"IsDPhiOutOfPlane","const AliVParticle&,const AliVParticle&","");
        cutElements_sq2.Add(inPlane);
        cutElements_Lq2.Add(outofPlane);
      }
      if(q2histo){
        AliAnalysisMuMuCutElement* smallq2 = cr->AddEventCut(*flowAnalysis,"Isq2InSmallRange","const AliVEvent&, Double_t&,Double_t&","0.,0.05");
        AliAnalysisMuMuCutElement* bigq2 = cr->AddEventCut(*flowAnalysis,"Isq2InLargeRange","const AliVEvent&,Double_t&,Double_t&","0.1,1.");
        cutElements_sq2.Add(smallq2);
        cutElements_Lq2.Add(bigq2);
      }

      cutElements.Add(pairTrue);
      cutElements.Add(ptRange);
      cutElements.Add(rabs);
      cutElements.Add(matchlow);
      cutElements.Add(eta);

      cutElements_sq2.Add(pairTrue);
      cutElements_sq2.Add(ptRange);
      cutElements_sq2.Add(rabs);
      cutElements_sq2.Add(matchlow);
      cutElements_sq2.Add(eta);

      cutElements_Lq2.Add(pairTrue);
      cutElements_Lq2.Add(ptRange);
      cutElements_Lq2.Add(rabs);
      cutElements_Lq2.Add(matchlow);
      cutElements_Lq2.Add(eta);

      // cutElements.Add(ps);
      // add them
      cr->AddCutCombination(cutElements);
      if(q2histo || raa) {
        cr->AddCutCombination(cutElements_sq2);
        cr->AddCutCombination(cutElements_Lq2);
      }
      // Adding the sub analysis
      task->AdoptSubAnalysis(minvAnalysis);
      task->AdoptSubAnalysis(flowAnalysis);
      // task->AdoptSubAnalysis(flowAnalysisSP);
    }
  }

  /// below are the kind of configurations that can be performed :
  /// - adding cuts (at event, track or pair level)
  /// - adding bins (in pt, y, centrality, etc...) for minv (and meanpt)

  // Configure sub analysis
  //===========================================================================
  AliAnalysisMuMuBinning* binning = task->Binning(); // Create and set the "binning manager"

  if (minvAnalysis)
  {
    // Integrated
   binning->AddBin("psi","integrated");

     // pt binning
    // binning->AddBin("psi","pt", 0.0, 4.0,"L");
    // binning->AddBin("psi","pt", 4.0, 8.0,"L");
    // binning->AddBin("psi","pt", 8.0, 12.0,"L");

    binning->AddBin("psi","pt", 0.0, 2.0,"S");
    binning->AddBin("psi","pt", 2.0, 4.0,"S");
    binning->AddBin("psi","pt", 4.0, 6.0,"S");
    binning->AddBin("psi","pt", 6.0, 8.0,"S");
    binning->AddBin("psi","pt", 8.0, 12.0,"S");

    // binning->AddBin("psi","pt", 3.0, 4.0,"AF");
    // binning->AddBin("psi","pt", 4.0, 5.0,"AF");
    // binning->AddBin("psi","pt", 5.0, 6.0,"AF");
    // binning->AddBin("psi","pt", 6.0, 7.0,"AF");
    // binning->AddBin("psi","pt", 7.0, 8.0,"AF");
    // binning->AddBin("psi","pt", 8.0, 9.0,"AF");
    // binning->AddBin("psi","pt", 9.0, 10.0,"AF");
    // binning->AddBin("psi","pt", 10.0, 11.0,"AF");
    // binning->AddBin("psi","pt", 11.0, 12.0,"AF");
  }

  // binning->AddBin("centrality","V0M",0.,5.);
  // binning->AddBin("centrality","V0M",5.,20.);
  binning->AddBin("centrality","V0M",20.,40.);
  // binning->AddBin("centrality","V0M",40.,60.);
  // binning->AddBin("centrality","V0M",60.,90.);
  // binning->AddBin("centrality","V0M",0.,90.);
  //===========================================================================
  // task->SetExpectedCorrectionPass("align");  //raw, plain, rec, align, twist and rescale
  // task->SetAlternativeCorrectionPass("latest");

  // add the configured task to the analysis manager
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  AliAnalysisDataContainer *coutputHC =
  mgr->CreateContainer("OC",AliMergeableCollection::Class(),AliAnalysisManager::kOutputContainer,outputname);

  AliAnalysisDataContainer *coutputCC =
  mgr->CreateContainer("CC",AliCounterCollection::Class(),AliAnalysisManager::kOutputContainer,outputname);

  AliAnalysisDataContainer* cparam =
  mgr->CreateContainer("BIN", AliAnalysisMuMuBinning::Class(),AliAnalysisManager::kParamContainer,outputname);

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  if(qntask) mgr->ConnectInput(task, 1, qntask);
  mgr->ConnectOutput(task, 1, coutputHC);
  mgr->ConnectOutput(task, 2, coutputCC);
  mgr->ConnectOutput(task, 3, cparam);

  return task;
}

