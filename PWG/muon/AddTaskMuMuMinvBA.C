///
/// Configuration example of a task to get invariant mass spectrum of dimuons
///
/// \author: B. Audurier (Subatech) (benjamin.audurier - at - subatech.in2p3.fr)
///

AliAnalysisTask* AddTaskMuMuMinvBA(const char* outputname,
                                   const char     * triggerClassesToConsider,
                                   const char     * triggerInputsToUse,
                                   const char     * beamYear,
                                   Bool_t simulations)
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

  //  Configure task
  //===========================================================================
  AliAnalysisTaskMuMu       * task = new AliAnalysisTaskMuMu; // Call the task
  task->SetBeamYear(beamYear);

  AliAnalysisMuMuEventCutter* eventCutter = new AliAnalysisMuMuEventCutter(triggerClassesToConsider,triggerInputsToUse); // To handle cuts on event
  AliAnalysisMuMuCutRegistry* cr = task->CutRegistry(); // Set CutRegistry

  // Default cuts on trigger and event level
  AliAnalysisMuMuCutElement* eventTrue = cr->AddEventCut(*eventCutter,"IsTrue","const AliVEvent&",""); 
  AliAnalysisMuMuCutElement* triggerSelection = cr->AddTriggerClassCut(*eventCutter,"SelectTriggerClass","const TString&,TString&,UInt_t,UInt_t,UInt_t","");
  AliAnalysisMuMuCutElement* ps = eventTrue;
  
  if (!simulations)
  {
    ps = cr->AddEventCut(*eventCutter,"IsPhysicsSelected","const AliInputEventHandler&","");
  }

  // Apply default cut
  cr->AddCutCombination(eventTrue);
  cr->AddCutCombination(ps);
  cr->AddCutCombination(triggerSelection);

  
//  AliAnalysisMuMuGlobal* globalAnalysis =  new AliAnalysisMuMuGlobal; // Basic histograms analysis;
  AliAnalysisMuMuSingle* singleAnalysis = new AliAnalysisMuMuSingle;// Analysis dealing with single muon
  AliAnalysisMuMuMinv  * minvAnalysis = new AliAnalysisMuMuMinv;// Analysis creating invariant mass spectrum

  
  // Configure sub analysis
  //===========================================================================
//  if ( globalAnalysis )
//  {
//    // Cuts on trigger level
//    AliAnalysisMuMuCutElement* triggerAll = cr->AddTriggerClassCut(*globalAnalysis,"SelectAnyTriggerClass","const TString&,TString&","");
//    // Adding this cut on trigger level  
//    cr->AddCutCombination(triggerAll);
//    task->AdoptSubAnalysis(globalAnalysis);
//  }
  if ( singleAnalysis )
  {
    // Cuts on tracks
    AliAnalysisMuMuCutElement* trackTrue = cr->AddTrackCut(*cr,"AlwaysTrue","const AliVParticle&",""); // Apply "AlwaysTrue" cut on AliVParticle derived from AliAnalysisMuMuSingle
    AliAnalysisMuMuCutElement* rabs = cr->AddTrackCut(*singleAnalysis,"IsRabsOK","const AliVParticle&","");
    AliAnalysisMuMuCutElement* matchlow = cr->AddTrackCut(*singleAnalysis,"IsMatchingTriggerLowPt","const AliVParticle&","");
    AliAnalysisMuMuCutElement* eta = cr->AddTrackCut(*singleAnalysis,"IsEtaInRange","const AliVParticle&","");
    AliAnalysisMuMuCutElement* pdca = cr->AddTrackCut(*singleAnalysis,"IsPDCAOK","const AliVParticle&","");
    
    // Create combination of cuts to apply
    cr->AddCutCombination(trackTrue);
    cr->AddCutCombination(matchlow);
    cr->AddCutCombination(rabs); 
    cr->AddCutCombination(eta); 
    cr->AddCutCombination(pdca); 
    // Adding the sub analysis
    task->AdoptSubAnalysis(singleAnalysis); 

    if ( minvAnalysis )
    {
      // Array of cut elements
      TObjArray cutElements;

      // Cuts on track level
      AliAnalysisMuMuCutElement* pairTrue = cr->AddTrackPairCut(*cr,"AlwaysTrue","const AliVParticle&,const AliVParticle&","");// Apply "AlwaysTrue" cut on AliVParticle derived from AliAnalysisMuMuMinv
      AliAnalysisMuMuCutElement* pairy = cr->AddTrackPairCut(*minvAnalysis,"IsRapidityInRange","const AliVParticle&,const AliVParticle&","");
      AliAnalysisMuMuCutElement* ptRange = cr->AddTrackPairCut(*minvAnalysis,"IsPtInRange","const AliVParticle&,const AliVParticle&,Double_t&,Double_t&","0.,10.");

      cutElements.Add(pairTrue);
      cutElements.Add(pairy);
      cutElements.Add(ptRange);
      cutElements.Add(rabs);
      cutElements.Add(matchlow);
      cutElements.Add(eta);
      cutElements.Add(pdca);

      // add them
      cr->AddCutCombination(cutElements);    
      // Adding the sub analysis
      task->AdoptSubAnalysis(minvAnalysis); 
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
    minvAnalysis->DefineMinvRange(2,8,0.025);
    // Integrated
    binning->AddBin("psi","integrated");
    // pt binning for low pt exces
//    binning->AddBin("psi","pt", 0.0, 0.1,"BENJ");
//    binning->AddBin("psi","pt", 0.1, 0.2,"BENJ");
//    binning->AddBin("psi","pt", 0.2, 0.3,"BENJ");
//    binning->AddBin("psi","pt", 0.3, 0.4,"BENJ");
//    binning->AddBin("psi","pt", 0.4, 0.5,"BENJ");
//    binning->AddBin("psi","pt", 0.5, 0.6,"BENJ");
//    binning->AddBin("psi","pt", 0.6, 0.7,"BENJ");
//    binning->AddBin("psi","pt", 0.7, 0.8,"BENJ");
//    binning->AddBin("psi","pt", 0.8, 0.9,"BENJ");
//    binning->AddBin("psi","pt", 0.9, 1.0,"BENJ");
//    binning->AddBin("psi","pt", 1.0, 1.1,"BENJ");
//    binning->AddBin("psi","pt", 1.1, 1.2,"BENJ");
//    binning->AddBin("psi","pt", 1.2, 1.3,"BENJ");
//    binning->AddBin("psi","pt", 1.3, 1.4,"BENJ");
//    binning->AddBin("psi","pt", 1.4, 1.5,"BENJ");
//    binning->AddBin("psi","pt", 1.5, 1.6,"BENJ");
//    binning->AddBin("psi","pt", 1.6, 1.7,"BENJ");
//    binning->AddBin("psi","pt", 1.7, 1.8,"BENJ");
//    binning->AddBin("psi","pt", 0.0, 0.3,"BENJ");
//    binning->AddBin("psi","pt", 0.3, 1.0,"BENJ");
//    binning->AddBin("psi","pt", 1.0, 8.0,"BENJ");
    // pt binning
    binning->AddBin("psi","pt", 0.0, 1.0,"BENJ");
    binning->AddBin("psi","pt", 1.0, 2.0,"BENJ");
    binning->AddBin("psi","pt", 2.0, 3.0,"BENJ");
    binning->AddBin("psi","pt", 3.0, 4.0,"BENJ");
    binning->AddBin("psi","pt", 4.0, 5.0,"BENJ");
    binning->AddBin("psi","pt", 5.0, 6.0,"BENJ");
    binning->AddBin("psi","pt", 6.0, 8.0,"BENJ");
    // y binning
    binning->AddBin("psi","y",-4,-3.75,"BENJ");
    binning->AddBin("psi","y",-3.75,-3.5,"BENJ");
    binning->AddBin("psi","y",-3.5,-3.25,"BENJ");
    binning->AddBin("psi","y",-3.25,-3,"BENJ");
    binning->AddBin("psi","y",-3,-2.75,"BENJ");
    binning->AddBin("psi","y",-2.75,-2.5,"BENJ");

     // phi binning
   binning->AddBin("psi","phi",0,0.4,"BENJ");
   binning->AddBin("psi","phi",0.4,0.8,"BENJ");
   binning->AddBin("psi","phi",0.8,1.2,"BENJ");
   binning->AddBin("psi","phi",1.2,1.6,"BENJ");
   binning->AddBin("psi","phi",1.6,2,"BENJ");
   binning->AddBin("psi","phi",2,2.4,"BENJ");
   binning->AddBin("psi","phi",2.4,2.8,"BENJ");
   binning->AddBin("psi","phi",2.4,2.8,"BENJ");
   binning->AddBin("psi","phi",2.8,3.14,"BENJ");
  }
  // v0 centrality binning
  binning->AddBin("centrality","V0M",0.,90.);
  binning->AddBin("centrality","V0M",0.,10.);
  binning->AddBin("centrality","V0M",10.,20.);
  binning->AddBin("centrality","V0M",20.,30.);
  binning->AddBin("centrality","V0M",30.,40.);
  binning->AddBin("centrality","V0M",40.,50.);
  binning->AddBin("centrality","V0M",50.,60.);
  binning->AddBin("centrality","V0M",60.,70.);
  binning->AddBin("centrality","V0M",70.,80.);
  binning->AddBin("centrality","V0M",80.,90.);
  // v0 centrality binning for low pt exces
//  binning->AddBin("centrality","V0M",0.,10.);
//  binning->AddBin("centrality","V0M",10.,30.);
//  binning->AddBin("centrality","V0M",30.,50.);
//  binning->AddBin("centrality","V0M",50.,70.);
//  binning->AddBin("centrality","V0M",70.,90.);

  // add the configured task to the analysis manager
  mgr->AddTask(task);
  
  TString output;
  output.Form("%s:%s",AliAnalysisManager::GetCommonFileName(),outputname);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  AliAnalysisDataContainer *coutputHC =
  mgr->CreateContainer(Form("OC_%s",outputname),AliMergeableCollection::Class(),AliAnalysisManager::kOutputContainer,output.Data());
  
  AliAnalysisDataContainer *coutputCC =
  mgr->CreateContainer(Form("CC_%s",outputname),AliCounterCollection::Class(),AliAnalysisManager::kOutputContainer,output.Data());
  
  AliAnalysisDataContainer* cparam =
  mgr->CreateContainer(Form("BIN_%s",outputname), AliAnalysisMuMuBinning::Class(),AliAnalysisManager::kParamContainer,output.Data());
  
  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutputHC);
  mgr->ConnectOutput(task, 2, coutputCC);
  mgr->ConnectOutput(task, 3, cparam);
  
  return task;
}

