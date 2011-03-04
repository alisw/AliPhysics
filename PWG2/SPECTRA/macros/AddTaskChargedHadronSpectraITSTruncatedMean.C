AliAnalysisChargedHadronSpectraITSTruncatedMeanTask* AddTaskChargedHadronSpectraITSTruncatedMean(Int_t lowcut=-1,Int_t upcut=-1,Int_t mc=0,Int_t hi=0 ,TString filename="./configChargedHadronSpectraITSTruncatedMeanTask.C")
{
	//pp settings 	
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
	gROOT->LoadMacro(filename.Data());
	AliAnalysisChargedHadronSpectraITSTruncatedMeanTask* task = GetAliAnalysisChargedHadronSpectraITSTruncatedMeanTask(mc);
	mgr->AddTask(task);
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();	
	AliAnalysisDataContainer *coutput =0x0;
	if(hi)
	{
		Float_t lowcencut=(Float_t)lowcut;
		Float_t upcencut=(Float_t)upcut;
		if(lowcencut<0.0)
			lowcencut=0.0;
		if(upcencut<0.0||upcencut>100.0)
			upcencut=100.0;
		task->SetCentralityCut(lowcencut,upcencut);
		task->SetHImode();	
		coutput =mgr->CreateContainer(Form("outputlow%dup%dHI1",lowcut,upcut),TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
			
	}
	else
	{
		task->SetMultiplicityCut(lowcut,upcut);
		coutput =mgr->CreateContainer(Form("outputlow%dup%dHI0",lowcut,upcut),TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
	}
	mgr->ConnectInput(task, 0, cinput);
 	mgr->ConnectOutput(task, 1, coutput);
	
	return task;
  }
