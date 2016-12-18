AliAnalysisTask *AddTaskPsi3HBTsystFitrange()
{
	//get the current analysis manager
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) {
		Error("AddTaskPsi3HBTsystFitrange", "No analysis manager found.");
		return NULL;
	}
	if (!mgr->GetInputEventHandler()) {
		::Error("AddTaskPsi3HBTsystFitrange", "This task requires an input event handler");
		return NULL;
	}
	TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

	Bool_t MCthere=kFALSE;
	AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
	if(!mcH){
		MCthere=kFALSE;
	}else{
		MCthere=kTRUE;
	}

	// add your task
	AliAnalysisTaskPsi3HBTsystFitrange *hbtpsi3 = new AliAnalysisTaskPsi3HBTsystFitrange("TaskHBT");
	mgr->AddTask(hbtpsi3);

	Bool_t	AnaCentral	= kTRUE; //true : 0-10%, false : 10 - 50%
	Int_t		TrigSel			= 3; //0:MB, 1:Semi-central, 2:Central, 3:All(MB+Semi+Central)

	Char_t trigSel_name[4][10] = {"MB", "Semi", "Central", "All"};
	Char_t anaCent_name[2][10] = {"1050", "010"};

	hbtpsi3 -> SetCentralityCentral(AnaCentral)	;
	hbtpsi3 -> SetCollisionCandidates(TrigSel)	;

	// set name of output directory
	TString containerName = mgr->GetCommonFileName();
	containerName += ":PWGCF_psi3hbtFitrange";
	TString SubcontainerName = Form("Trig%s_Cent%s", trigSel_name[TrigSel], anaCent_name[AnaCentral]);
	AliAnalysisDataContainer *cinput		= mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput1	= mgr->CreateContainer(SubcontainerName, TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
	mgr->ConnectInput(hbtpsi3,	0, cinput)	;
	mgr->ConnectOutput(hbtpsi3, 1, coutput1);

	return NULL;
}
