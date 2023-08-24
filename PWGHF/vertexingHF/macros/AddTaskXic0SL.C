AliAnalysisTaskSEXic0SL* AddTaskXic0SL(const char* name, const char* option)
{
	AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr || !mgr->GetInputEventHandler()) return NULL;

	//-------------------------------------------

	AliAnalysisTaskSEXic0SL* task = new AliAnalysisTaskSEXic0SL(name, option);
	if (!task) return NULL;

	//Setup (default)
	task->SetMC(false);
	task->SetPA(false);
	task->SetTrigMB(true);
	task->SetTrigHMV0(true);
	task->SetCutsLegacy(false); //Set conditions same to the "old" code
	task->SetCutsByFile(false); //Set conditions using external ROOT file
	task->SetValidEvtOnly(false);

	//Setup (via invoke)
	TString taskOpt = option;
	if (taskOpt.Contains("MC")) task->SetMC(true);
	if (taskOpt.Contains("PA")) task->SetPA(true);
	if (taskOpt.Contains("LEGACY")) task->SetCutsLegacy(true);
	if (taskOpt.Contains("VALIDONLY")) task->SetValidEvtOnly(true);

	mgr->AddTask(task);

	//-------------------------------------------

	TString outFile = AliAnalysisManager::GetCommonFileName();
	outFile += Form(":PWG3_D2H_Xic0SL%s", option);

	AliAnalysisDataContainer* ADC[] =
	{
		mgr->GetCommonInputContainer(),
		mgr->CreateContainer(Form("Histo%s", option),
				TDirectory::Class(), AliAnalysisManager::kOutputContainer, outFile.Data()),
		mgr->CreateContainer(Form("Tree%s", option),
				TTree::Class(), AliAnalysisManager::kOutputContainer, outFile.Data()),
		mgr->CreateContainer(Form("ANC_MB_0to100%s", option),
				AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outFile.Data()),
		mgr->CreateContainer(Form("ANC_MB_30to100%s", option),
				AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outFile.Data()),
		mgr->CreateContainer(Form("ANC_MB_0p1to30%s", option),
				AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outFile.Data()),
		mgr->CreateContainer(Form("ANC_HMV0_0to0p1%s", option),
				AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outFile.Data()),
		mgr->CreateContainer(Form("ANC_INEL0_MB_0to100%s", option),
				AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outFile.Data()),
		mgr->CreateContainer(Form("ANC_INEL0_MB_30to100%s", option),
				AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outFile.Data()),
		mgr->CreateContainer(Form("ANC_INEL0_MB_0p1to30%s", option),
				AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outFile.Data()),
		mgr->CreateContainer(Form("ANC_INEL0_HMV0_0to0p1%s", option),
				AliNormalizationCounter::Class(), AliAnalysisManager::kOutputContainer, outFile.Data())
	};

	const int nADC = sizeof(ADC)/sizeof(ADC[0]);
	for (int i=0; i<nADC; i++) { (i==0) ? mgr->ConnectInput(task,i,ADC[i]) : mgr->ConnectOutput(task,i,ADC[i]); }

	return task;
}//AddTask
