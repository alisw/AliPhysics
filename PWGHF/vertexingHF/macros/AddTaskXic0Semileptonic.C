// Macro for analisys task of preselected-central-diffractive events
//------------------------------------------------------------------
// When there is no time to wait for legotrain results,
// please use this macro to run jobs on grid.
// Recommended for pp or pA runs, not for AA runs.
//------------------------------------------------------------------
AliAnalysisTaskSEXic0Semileptonic *AddTaskXic0Semileptonic(
		const char *taskname = "Xic"
		, const char *option = "DataAOD" // when scanning AOD, add "AOD"
		, const char *coll = "PP"
		, TString outputFileName = ""
		, bool UseTrigHM = false
		)
{
	// analysis manager
	AliAnalysisManager* mgr =  AliAnalysisManager::GetAnalysisManager();
	if (!mgr)
	{
		::Error("AddTaskXic0pp13TeV", "No analysis manager to connect to.");
		return NULL;
	}  

	AliAnalysisTaskSEXic0Semileptonic *task = new AliAnalysisTaskSEXic0Semileptonic(
			taskname, Form("%s%s_task",taskname,option));
	TString foption = option;
	TString fcoll = coll;
	bool ismc = kFALSE;
	bool ispa = kFALSE;
	if(foption.Contains("MC")) ismc = true;
	if(fcoll.Contains("PA")) ispa = true;
	task->SetMC(ismc);
	task->SetPA(ispa);
	
	task->UseTrig_kINT7();
	if (UseTrigHM)
	{
		task->UseTrig_kHMV0();
		task->UseTrig_kHMSPD();
	}

	TString outputfile = AliAnalysisManager::GetCommonFileName();
	outputfile += ":PWG3_D2H_Xic02eXipp13TeV";
	outputfile += outputFileName.Data();
	if (UseTrigHM) outputfile += "_HM";
	// Create containers for input/output

	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histogram", TDirectory::Class(), AliAnalysisManager::kOutputContainer,outputfile.Data());
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("cut", TList::Class(), AliAnalysisManager::kOutputContainer,outputfile.Data());
	AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("MCXicTree", TTree::Class(), AliAnalysisManager::kOutputContainer,outputfile.Data());
	AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("PairTree", TTree::Class(), AliAnalysisManager::kOutputContainer,outputfile.Data());
	AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("MCPairTree", TTree::Class(), AliAnalysisManager::kOutputContainer,outputfile.Data());
	AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("EventTree", TTree::Class(), AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutput7 = mgr->CreateContainer("eleXiCounter",AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); //counter
	
    mgr->AddTask(task);
	mgr->ConnectInput(task, 0, cinput);
	mgr->ConnectOutput(task, 1, coutput1);
	mgr->ConnectOutput(task, 2, coutput2);
	mgr->ConnectOutput(task, 3, coutput3);
	mgr->ConnectOutput(task, 4, coutput4);
	mgr->ConnectOutput(task, 5, coutput5);
	mgr->ConnectOutput(task, 6, coutput6);
	mgr->ConnectOutput(task, 7, coutput7);
	mgr->SetDebugLevel(2);
	if (!mgr->InitAnalysis()) return 0x0;

	return task;
}
