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
    cr->AddCutCombination(ps);
  }

  cr->AddCutCombination(eventTrue);

  
  AliAnalysisMuMuGlobal* globalAnalysis = new AliAnalysisMuMuGlobal;

  AliAnalysisMuMuCutElement* triggerSelection = cr->AddTriggerClassCut(*eventCutter,"SelectTriggerClass","const TString&,TString&,UInt_t,UInt_t,UInt_t","");
  
  cr->AddCutCombination(triggerSelection);
  
  task->AdoptSubAnalysis(globalAnalysis);

  /// below are the kind of configurations that can be performed :
  /// - adding cuts (at event, track or pair level)
  /// - adding bins (in pt, y, centrality, etc...) for minv (and meanpt)
  

  AliAnalysisMuMuBinning* binning = task->Binning();
  
  binning->AddBin("centrality","pp");
  
   // disable some histograms if we don't want them
  task->DisableHistograms("^V02D");
  task->DisableHistograms("^dca");
  task->DisableHistograms("^Chi12");
  task->DisableHistograms("^Rabs12");
  task->DisableHistograms("BCX");
  
  task->DisableHistograms("SPDZvertexResolutionNContributors");
  task->DisableHistograms("ZvertexMinusSPDZvertexNContributors");

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

