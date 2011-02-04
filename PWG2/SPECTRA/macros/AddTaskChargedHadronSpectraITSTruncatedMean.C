 AliAnalysisChargedHadronSpectraITSTruncatedMeanTask* AddTaskChargedHadronSpectraITSTruncatedMean(Float_t lowcut=-1.0,Float_t upcut=-1.0,Int_t mc=0)
 {
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) 
	{
		::Error("AddTaskITSsaTracks", "No analysis manager to connect to.");
		return NULL;
	}   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
	if (!mgr->GetInputEventHandler()) 
	{
		::Error("AddTaskITSsaTracks", "This task requires an input event handler");
		return NULL;
	}   
  
  	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	if(type.Contains("AOD"))
	{
		::Error("AddTaskITSsaTracks", "This task requires to run on ESD");
		return NULL;
	}
	TString outputFileName = AliAnalysisManager::GetCommonFileName();
	outputFileName += ":PWG2SpectraITSTPC";
	
	gROOT->LoadMacro("./config_ChargedHadronSpectraITSTruncatedMeanTask.C");
	AliAnalysisChargedHadronSpectraITSTruncatedMeanTask* task = GetAliAnalysisChargedHadronSpectraITSTruncatedMeanTask(mc);
	task->SetMultiplicityCut(lowcut,upcut);
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput = mgr->CreateContainer("output", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
	mgr->ConnectInput(task, 0, cinput);
 	mgr->ConnectOutput(task, 1, coutput);
	mgr->AddTask(task);
	return task;
  }