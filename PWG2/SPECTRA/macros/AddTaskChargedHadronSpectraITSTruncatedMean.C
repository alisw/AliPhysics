 AliAnalysisChargedHadronSpectraITSTruncatedMeanTask* AddTaskChargedHadronSpectraITSTruncatedMean(Int_t lowcut=-1,Int_t upcut=-1,Int_t mc=0)
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
	gROOT->LoadMacro("./configChargedHadronSpectraITSTruncatedMeanTask.C");
	AliAnalysisChargedHadronSpectraITSTruncatedMeanTask* task = GetAliAnalysisChargedHadronSpectraITSTruncatedMeanTask(mc);
	mgr->AddTask(task);
	Int_t upint=0;
	Int_t lowint=0;
	
	
	task->SetMultiplicityCut(lowcut,upcut);
	
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();	
	AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("outputlow%dup%dHI0",lowcut,upcut),TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
	mgr->ConnectInput(task, 0, cinput);
 	mgr->ConnectOutput(task, 1, coutput);
	
	return task;
  }

AliAnalysisChargedHadronSpectraITSTruncatedMeanTask* AddTaskChargedHadronSpectraITSTruncatedMean(Float_t lowcut=0.0,Float_t upcut=100.0,Int_t mc=0, Int_t hi=1)
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
	gROOT->LoadMacro("./configChargedHadronSpectraITSTruncatedMeanTask.C");
	AliAnalysisChargedHadronSpectraITSTruncatedMeanTask* task = GetAliAnalysisChargedHadronSpectraITSTruncatedMeanTask(mc);
	mgr->AddTask(task);
	Int_t upint=0;
	Int_t lowint=0;
	
	upint=(Int_t)upcut;
	lowint=(Int_t)lowcut;
	
	task->SetCentralityCut(lowcut,upcut);
	task->SetHImode();
	
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();	
	AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("outputlow%dup%dHI1",lowint,upint),TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
	mgr->ConnectInput(task, 0, cinput);
 	mgr->ConnectOutput(task, 1, coutput);
	
	return task;
  }