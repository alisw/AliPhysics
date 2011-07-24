
AliAnalysisTaskLRC *AddTaskLRC(Bool_t RunKine=kFALSE, TString OutputRootFolder=":PWG2LRC")
{
// This macro adds AliAnalysisTaskLRC to existing AnalysisManager
// RunKine paramiter switch task to kinematics analysis 

// Author : Andrey Ivanov , St.Peterburg State University
// Email: Andrey.Ivanov@cern.ch

// Version line : 3.6
// Version 3.6.7

	gROOT->LoadMacro("configLRCAnalysis.C");
	
	//gROOT->LoadMacro("AliAnalysisTaskIA.cxx+g"); 
	

	// A. Get the pointer to the existing analysis manager via the static access method.
	//==============================================================================
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskLRC", "No analysis manager to connect to.");
		return NULL;
	}  
	
	// B. Check the analysis type using the event handlers connected to the analysis
	//    manager. The availability of MC handler cann also be checked here.
	//==============================================================================
	if (!mgr->GetInputEventHandler()) {
		Error("AddTaskLRC", "This task requires an input event handler");
		return NULL;
	}  
	TString type = mgr->GetInputEventHandler()->GetDataType(); 
	cout<<" # TaskLRC - input :"<<type<<"\n";

	// MB - Global2

	taskLRC = createLRCtaskSkeleton("Task_LRC_MB_Global2",RunKine);
	taskLRC->SetTrackCuts(createAliLRCcuts("Global2"));
	addAliLRCProcessors(taskLRC,kTRUE);
	taskLRC->SetVtxDiamond(0.4,0.4,7.0);
	taskLRC->SetMaxPtLimit(1.5);
	taskLRC->SetMinPtLimit(0.3);
	mgr->AddTask(taskLRC);
	configureLRCtaskOutput(taskLRC,":PWG2LRC");

return taskLRC;
}

