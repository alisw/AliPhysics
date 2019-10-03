///
/// Configuration example of a MuMu task to get simple counting of events, without histogramming
///
/// \author: L. Aphecetche (Subatech) (laurent.aphecetche - at - subatech.in2p3.fr)
///

AliAnalysisTask* AddTaskMuMuCount(const char* foldername,
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
  AliAnalysisMuMuCutElement* psI = eventTrue;
  AliAnalysisMuMuCutElement* psIiM = eventTrue;
  AliAnalysisMuMuCutElement* psMUL = eventTrue;
  AliAnalysisMuMuCutElement* psMSL = eventTrue;

  if (!simulations)
  {
  ps = cr->AddEventCut(*eventCutter,"IsPhysicsSelectedANY","const AliInputEventHandler&","");
  cr->AddCutCombination(ps);

  psI = cr->AddEventCut(*eventCutter,"IsPhysicsSelectedINT7","const AliInputEventHandler&","");
  cr->AddCutCombination(psI);
  
  psIiM = cr->AddEventCut(*eventCutter,"IsPhysicsSelectedINT7inMUON","const AliInputEventHandler&","");
  cr->AddCutCombination(psIiM);
  
  psMUL = cr->AddEventCut(*eventCutter,"IsPhysicsSelectedMUL","const AliInputEventHandler&","");
  cr->AddCutCombination(psMUL);

  psMSL = cr->AddEventCut(*eventCutter,"IsPhysicsSelectedMSL","const AliInputEventHandler&","");
  cr->AddCutCombination(psMSL);

  }

  cr->AddCutCombination(eventTrue);

  AliAnalysisMuMuCutElement* triggerSelection = cr->AddTriggerClassCut(*eventCutter,"SelectTriggerClass","const TString&,TString&,UInt_t,UInt_t,UInt_t","");
  
  cr->AddCutCombination(triggerSelection);
  
  AliAnalysisMuMuBinning* binning = task->Binning();
  
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

