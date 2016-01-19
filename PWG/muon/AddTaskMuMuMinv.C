///
/// Configuration example of a task to get invariant mass spectrum of dimuons
///
/// \author: L. Aphecetche (Subatech) (laurent.aphecetche - at - subatech.in2p3.fr)
///

AliAnalysisTask* AddTaskMuMuMinv(const char* foldername,
                                 const char* triggerClassesToConsider,
                                 const char* triggerInputsToUse,
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
  
  AliAnalysisMuMuEventCutter* eventCutter = new AliAnalysisMuMuEventCutter(triggerClassesToConsider,triggerInputsToUse);
  
  // to the very least we need a cut combination at the event level
  // (i.e. if we select no event, nothing will happen...)
  // we use by default the usual physics selection for real data
  // and a trivial selection (i.e. always true) for simulations
  //
  AliAnalysisMuMuCutElement* eventTrue = cr->AddEventCut(*eventCutter,"IsTrue","const AliVEvent&","");
  
  AliAnalysisMuMuCutElement* ps = eventTrue;
  
  if (!simulations)
  {
    ps = cr->AddEventCut(*eventCutter,"IsPhysicsSelected","const AliInputEventHandler&","");
    cr->AddCutCombination(ps);
  }
  
  cr->AddCutCombination(eventTrue);
  
  AliAnalysisMuMuGlobal* globalAnalysis = new AliAnalysisMuMuGlobal;
  
  AliAnalysisMuMuCutElement* triggerSelection = cr->AddTriggerClassCut(*eventCutter,"SelectTriggerClass","const TString&,TString&,UInt_t,UInt_t,UInt_t","");
  
  cr->AddCutCombination(triggerSelection);
  
  AliAnalysisMuMuMinv* minvAnalysis = new AliAnalysisMuMuMinv;
  
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
    minvAnalysis->DefineMinvRange(0,16,0.025);

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

  if ( TString(beamYear).Contains("pbpb",TString::kIgnoreCase) )
  {
    binning->AddBin("centrality","V0M");
    binning->AddBin("centrality","V0M",0,90);
    binning->AddBin("centrality","V0M",0,10);
    binning->AddBin("centrality","V0M",10,20);
    binning->AddBin("centrality","V0M",20,30);
    binning->AddBin("centrality","V0M",30,40);
    binning->AddBin("centrality","V0M",40,50);
    binning->AddBin("centrality","V0M",50,60);
    binning->AddBin("centrality","V0M",60,70);
    binning->AddBin("centrality","V0M",70,80);
    binning->AddBin("centrality","V0M",80,90);
  }
  else
  {
    binning->AddBin("centrality","pp");
  }
  
  // add the configured task to the analysis manager
  mgr->AddTask(task);
  
  TString output;
  output.Form("%s:%s",AliAnalysisManager::GetCommonFileName(),foldername);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  AliAnalysisDataContainer *coutputHC =
  mgr->CreateContainer(Form("OC_%s",foldername),AliMergeableCollection::Class(),AliAnalysisManager::kOutputContainer,output.Data());
  
  AliAnalysisDataContainer *coutputCC =
  mgr->CreateContainer(Form("CC_%s",foldername),AliCounterCollection::Class(),AliAnalysisManager::kOutputContainer,output.Data());
  
  AliAnalysisDataContainer* cparam =
  mgr->CreateContainer(Form("BIN_%s",foldername), AliAnalysisMuMuBinning::Class(),AliAnalysisManager::kParamContainer,output.Data());
  
  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutputHC);
  mgr->ConnectOutput(task, 2, coutputCC);
  mgr->ConnectOutput(task, 3, cparam);
  
  return task;
}

