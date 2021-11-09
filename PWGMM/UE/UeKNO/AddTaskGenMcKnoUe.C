///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskGenMcKnoUe Macro to run on grids            //
//            Last update: 5/11/2021                             //
//                                                               //
///////////////////////////////////////////////////////////////////

AliAnalysisTask *AddTaskGenMcKnoUe(Double_t minpT=0.5, Bool_t isPP=kTRUE, TString suffixName ="")
{
	AliAnalysisTaskGenMcKnoUe* taskUE = new AliAnalysisTaskGenMcKnoUe("taskKno");
	taskUE->SetPtMin(minpT);
	taskUE->SetIsPP(isPP); 

	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Printf("AliAnalysisTaskSimSpectraLF: No analysis manager to connect to.");
		return 0x0;
	}
	// get the input event handler this handler is part of the managing system and feeds events to your task
	if (!mgr->GetMCtruthEventHandler()) {
		Printf("AliAnalysisTaskSimSpectraLF: This task requires an input MC event handler.");
		return 0x0;
	}

	mgr->AddTask(taskUE);

	// Create containers for input/output

	TString finDirname    = "";
	TString inname    = "cinput";
	TString outBasic    = "cList";

	finDirname    += suffixName.Data();
	inname    += finDirname.Data();
	outBasic    += finDirname.Data();

	// Input and Output Slots
	//===========================================================================

	TString outputfile = AliAnalysisManager::GetCommonFileName();
	outputfile += Form(":PWGMM_GenKnoUe_%1.2f",minpT);

	AliAnalysisDataContainer *coutSim = mgr->CreateContainer(outBasic,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

	mgr->ConnectInput (taskUE, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(taskUE, 1, coutSim);

	return taskUE;
}
