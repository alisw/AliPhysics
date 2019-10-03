///
/// Configuration example of a task to get basic distributions from MC generated events
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
  }
  
  cr->AddCutCombination(ps);
  
  AliAnalysisMuMuMCGene* analysis = new AliAnalysisMuMuMCGene;

  // we select all triggers
  AliAnalysisMuMuCutElement* triggerAll = cr->AddTriggerClassCut(*analysis,"SelectAnyTriggerClass","const TString&,TString&","");
  
  cr->AddCutCombination(triggerAll);
  
  task->AdoptSubAnalysis(analysis);
  
  /// below are the kind of configurations that can be performed :
  /// - adding cuts (at event, track or pair level)
  /// - adding bins (in pt, y, centrality, etc...) for minv (and meanpt)
  
  AliAnalysisMuMuBinning* binning = task->Binning();
  
  // v0 centrality binning
  
  binning->AddBin("centrality","v0a");
  //  binning->AddBin("centrality","v0a",0,5);
  
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

