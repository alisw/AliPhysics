///
/// Configuration example of a task to get invariant mass spectrum of dimuons
///
/// \author: L. Aphecetche (Subatech) (laurent.aphecetche - at - subatech.in2p3.fr)
///

AliAnalysisTask* AddTaskMuMu(const char* outputname,
                             TList* triggerClassesToConsider,
                             const char* beamYear,
                             Bool_t simulations)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMuMu", "No analysis manager to connect to.");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMuMu", "This task requires an input event handler");
    return NULL;
  }
  
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  AliAnalysisTaskMuMu* task = new AliAnalysisTaskMuMu;
  task->SetBeamYear(beamYear);
  
  AliAnalysisMuMuCutRegistry* cr = task->CutRegistry();
  
  AliAnalysisMuMuEventCutter* eventCutter = new AliAnalysisMuMuEventCutter(triggerClassesToConsider);
  
  // to the very list we need a cut combination at the event level
  // (i.e. if we select no event, nothing will happen...)
  // we use by default the usual physics selection for real data
  // and a trivial selection (i.e. always true) for simulations
  //
  AliAnalysisMuMuCutElement* eventTrue = cr->AddEventCut(*eventCutter,"IsTrue","const AliVEvent&","");
  
  AliAnalysisMuMuCutElement* ps = eventTrue;
  
  if (!simulations)
  {
    ps = cr->AddEventCut(*eventCutter,"IsPhysicsSelected","const AliInputEventHandler&","");
//    cr->AddCutCombination(ps);
  }
  
  cr->AddCutCombination(eventTrue);
  
  AliAnalysisMuMuGlobal* globalAnalysis = new AliAnalysisMuMuGlobal;
  
  AliAnalysisMuMuCutElement* triggerSelection = cr->AddTriggerClassCut(*eventCutter,"SelectTriggerClass","const TString&,TString&,UInt_t,UInt_t,UInt_t","");
  
  cr->AddCutCombination(triggerSelection);
  
  AliAnalysisMuMuMinv* minvAnalysis = 0x0; // new AliAnalysisMuMuMinv;
  
  AliAnalysisMuMuSingle* singleAnalysis = new AliAnalysisMuMuSingle;
  
  AliAnalysisMuMuCutElement* trackTrue = cr->AddTrackCut(*cr,"AlwaysTrue","const AliVParticle&","");
  
  if ( singleAnalysis )
  {
    AliAnalysisMuMuCutElement* rabs = cr->AddTrackCut(*singleAnalysis,"IsRabsOK","const AliVParticle&","");
    AliAnalysisMuMuCutElement* eta = cr->AddTrackCut(*singleAnalysis,"IsEtaInRange","const AliVParticle&","");
    AliAnalysisMuMuCutElement* matchlow = cr->AddTrackCut(*singleAnalysis,"IsMatchingTriggerLowPt","const AliVParticle&","");
    
    cr->AddCutCombination(trackTrue);
    cr->AddCutCombination(rabs);
    cr->AddCutCombination(eta);
    cr->AddCutCombination(rabs,eta);
    cr->AddCutCombination(rabs,eta,matchlow);
    
    if ( minvAnalysis )
    {
      AliAnalysisMuMuCutElement* pairTrue = cr->AddTrackPairCut(*cr,"AlwaysTrue","const AliVParticle&, const AliVParticle&","");
      
      AliAnalysisMuMuCutElement* pairy = cr->AddTrackPairCut(*minvAnalysis,"IsRapidityInRange","const AliVParticle&,const AliVParticle&","");
      
      cr->AddCutCombination(rabs,eta,pairy);
      cr->AddCutCombination(rabs,eta,matchlow,pairy);
      cr->AddCutCombination(rabs,pairy);
      cr->AddCutCombination(pairy);
      
      cr->AddCutCombination(pairTrue);
      
    }
  }
  
  if ( globalAnalysis) task->AdoptSubAnalysis(globalAnalysis);
  if ( minvAnalysis ) task->AdoptSubAnalysis(minvAnalysis);
  if ( singleAnalysis ) task->AdoptSubAnalysis(singleAnalysis);
  
  /// below are the kind of configurations that can be performed :
  /// - adding cuts (at event, track or pair level)
  /// - adding bins (in pt, y, centrality, etc...) for minv (and meanpt)
  
  AliAnalysisMuMuBinning* binning = task->Binning();
  
  if ( minvAnalysis )
  {
    binning->AddBin("psi","integrated");
    
    binning->AddBin("psi","pt", 0.0, 1.0);
    binning->AddBin("psi","pt", 1.0, 2.0);
    binning->AddBin("psi","pt", 2.0, 3.0);
    binning->AddBin("psi","pt", 3.0, 4.0);
    binning->AddBin("psi","pt", 4.0, 5.0);
    binning->AddBin("psi","pt", 5.0, 6.0);
    binning->AddBin("psi","pt", 6.0, 7.0);
    binning->AddBin("psi","pt", 7.0, 9.0);
    binning->AddBin("psi","pt", 9.0,11.0);
    binning->AddBin("psi","pt",11.0,15.0);
    binning->AddBin("psi","pt",15.0,25.0);
  }
  
  binning->AddBin("centrality","pp");
  
  // disable some histograms if we don't want them
  //   task->DisableHistograms("^V02D");
  //   task->DisableHistograms("^dca");
  //   task->DisableHistograms("^Chi12");
  //   task->DisableHistograms("^Rabs12");
  // //  task->DisableHistograms("BCX");
  
  //   task->DisableHistograms("SPDZvertexResolutionNContributors");
  //   task->DisableHistograms("ZvertexMinusSPDZvertexNContributors");
  
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
  mgr->ConnectOutput(task, 1, coutputHC);
  mgr->ConnectOutput(task, 2, coutputCC);
  mgr->ConnectOutput(task, 3, cparam);
  
  return task;
}

