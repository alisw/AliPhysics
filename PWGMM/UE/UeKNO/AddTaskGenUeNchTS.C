// 
// AddTaskGenUeNchTS
// 
AliAnalysisTask *AddTaskGenUeNchTS(TString suffixName =""){

	// Create the task, add it to manager and configure it
	//===========================================================================

	AliAnalysisTaskGenUeNchTS* taskUeNchTSMM = new  AliAnalysisTaskGenUeNchTS("AliAnalysisTaskGenUeNchTS");
	taskUeNchTSMM -> SetGenerator("Pythia8");  
	taskUeNchTSMM -> SetYRange(0.5);
	taskUeNchTSMM -> SetMinPtLeading(5.0);

	// Get the pointer to the existing analysis manager via the static access method
	//===========================================================================

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if(!mgr){ Printf("AliAnalysisTaskSimSpectraLF: No analysis manager to connect to."); return NULL; }

	// Check the analysis type using the event handlers connected to the analysis manager
	//===========================================================================

	if(!mgr->GetMCtruthEventHandler()){ Printf("AliAnalysisTaskSimSpectraLF: This task requires an input MC event handler."); return NULL; }

	// ADD the task
	//===========================================================================
	mgr -> AddTask(taskUeNchTSMM);  

	// Create containers for input/output

	TString finDirname	= "";
	TString inname	= "cinput";
	TString outBasic	= "cList";

	finDirname	+= suffixName.Data();
	inname	+= finDirname.Data();
	outBasic	+= finDirname.Data();


	// Input and Output Slots
	//===========================================================================

	TString outputfile = AliAnalysisManager::GetCommonFileName();
	outputfile += ":PWGMM_SimUeNchTSMM";

	AliAnalysisDataContainer *coutSim = mgr->CreateContainer(outBasic,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

	mgr->ConnectInput (taskUeNchTSMM, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(taskUeNchTSMM, 1, coutSim);

	return taskUeNchTSMM;

}

