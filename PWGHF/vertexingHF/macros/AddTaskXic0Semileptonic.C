AliAnalysisTaskSEXic0Semileptonic *AddTaskXic0Semileptonic(
		const char *taskname   = "Xic",
		const char *option     = "DataAOD", //When scanning AOD, add "AOD"
		const char *coll       = "PP",
		TString outputFileName = "",
		bool UseTrigHM         = true //Eanble if high multiplicity trigger required
		)
{
	//Analysis manager
	AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) { ::Error("AddTaskXic0pp13TeV", "No analysis manager to connect to."); return NULL; }  

	//This task
	AliAnalysisTaskSEXic0Semileptonic *task = new AliAnalysisTaskSEXic0Semileptonic(
			taskname, Form("%s%s_task", taskname, option));

	//-------------------------------------------

	//Data or MC
	bool ismc = kFALSE;
	TString foption = option;
	if (foption.Contains("MC")) ismc = true;
	task->SetMC(ismc);

	//Collision type
	TString fcoll = coll;
	bool ispa = kFALSE;
	if (fcoll.Contains("PA")) ispa = true;
	task->SetPA(ispa);

	//Trigger to use: kINT7 is default but can be disabled if needed
	task->UseTrig_kINT7();
	if (UseTrigHM)
	{
		task->UseTrig_kHMV0();
		//task->UseTrig_kHMSPD();
	}

	//-------------------------------------------

	//Output file and subdir
	TString outputfile = AliAnalysisManager::GetCommonFileName();
	outputfile += ":PWG3_D2H_Xic02eXipp13TeV";
	outputfile += outputFileName.Data();
	if (UseTrigHM) outputfile += "_HM";

	//Create containers for input/output
	AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histogram",
			TDirectory::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
	AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("cut",
			TList::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
	AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("MCXicTree",
			TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
	AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("PairTree",
			TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
	AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("MCPairTree",
			TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
	AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("EventTree",
			TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
	AliAnalysisDataContainer *coutput7 = mgr->CreateContainer("eleXiCounter",
			AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());

	//kimc, Mar. 18, additional AliNormalizationCounters
	AliAnalysisDataContainer *coutput8 = mgr->CreateContainer("OldANC_MB_0to100",
			AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
	AliAnalysisDataContainer *coutput9 = mgr->CreateContainer("OldANC_MB_0p1to30",
			AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
	AliAnalysisDataContainer *coutput10 = mgr->CreateContainer("OldANC_MB_30to100",
			AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
	AliAnalysisDataContainer *coutput11 = mgr->CreateContainer("OldANC_HMV0_0to0p1",
			AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());

	//kimc, June 22, additional AliNormalizationCounters
	AliAnalysisDataContainer *coutput12 = mgr->CreateContainer("OldANCINEL0_MB_0to100",
			AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
	AliAnalysisDataContainer *coutput13 = mgr->CreateContainer("OldANCINEL0_MB_0p1to30",
			AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
	AliAnalysisDataContainer *coutput14 = mgr->CreateContainer("OldANCINEL0_MB_30to100",
			AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());
	AliAnalysisDataContainer *coutput15 = mgr->CreateContainer("OldANCINEL0_HMV0_0to0p1",
			AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outputfile.Data());

	//-------------------------------------------

	mgr->AddTask(task);
	mgr->ConnectInput (task, 0, cinput);
	mgr->ConnectOutput(task, 1, coutput1);
	mgr->ConnectOutput(task, 2, coutput2);
	mgr->ConnectOutput(task, 3, coutput3);
	mgr->ConnectOutput(task, 4, coutput4);
	mgr->ConnectOutput(task, 5, coutput5);
	mgr->ConnectOutput(task, 6, coutput6);
	mgr->ConnectOutput(task, 7, coutput7);

	mgr->ConnectOutput(task, 8, coutput8); //Mar. 18
	mgr->ConnectOutput(task, 9, coutput9);
	mgr->ConnectOutput(task, 10, coutput10);
	mgr->ConnectOutput(task, 11, coutput11);

	mgr->ConnectOutput(task, 12, coutput12); //June 22
	mgr->ConnectOutput(task, 13, coutput13);
	mgr->ConnectOutput(task, 14, coutput14);
	mgr->ConnectOutput(task, 15, coutput15);

	//mgr->SetDebugLevel(2);
	//if (!mgr->InitAnalysis()) return 0x0;

	return task;
}
