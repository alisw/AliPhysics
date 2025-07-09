#include "AliAnalysisTaskpOmegaDibaryon.h"

AliAnalysisTaskpOmegaDibaryon *AddTaskpOmegaDibaryon(
		const TString taskname = "MSDibaryons",
		//const UInt_t trigger = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral,	
		const UInt_t trigger = AliVEvent::kCentral,
		//const UInt_t trigger = AliVEvent::kCentral,	
		const TString estimator = "V0M",
		const Float_t cenmin = 0,
		const Float_t cenmax = 90
) {

	// Connection to the analysis manager
	//================================================================
	AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
	if(!mgr){
		printf("ERROR: No analysis manager to connect to.\n");
		return NULL;
	}

	// Task creation
	//================================================================
	AliAnalysisTaskpOmegaDibaryon *task=new AliAnalysisTaskpOmegaDibaryon(taskname.Data());
	task->SetTrigger(trigger);
	//task->SelectCollisionCandidates(trigger);
	mgr->AddTask(task);

	TString outlistname = Form("outlist_msdibaryons_%02d%02d",int(cenmin),int(cenmax));
	// Connect input/output
	//================================================================
	AliAnalysisDataContainer *cinput =mgr->GetCommonInputContainer();
	AliAnalysisDataContainer *coutput=mgr->CreateContainer(outlistname,THashList::Class(),
							       AliAnalysisManager::kOutputContainer,
							       AliAnalysisManager::GetCommonFileName());
	
	AliAnalysisDataContainer *coutput2=mgr->CreateContainer("fTree_omegam",TTree::Class(),
								AliAnalysisManager::kOutputContainer,
								AliAnalysisManager::GetCommonFileName());

	AliAnalysisDataContainer *coutput3=mgr->CreateContainer("fTree_omegap",TTree::Class(),
								AliAnalysisManager::kOutputContainer,
								AliAnalysisManager::GetCommonFileName());
	
	AliAnalysisDataContainer *coutput4=mgr->CreateContainer("fTree_sidebandm",TTree::Class(),
								AliAnalysisManager::kOutputContainer,
								AliAnalysisManager::GetCommonFileName());

	AliAnalysisDataContainer *coutput5=mgr->CreateContainer("fTree_sidebandp",TTree::Class(),
								AliAnalysisManager::kOutputContainer,
								AliAnalysisManager::GetCommonFileName());
	
	AliAnalysisDataContainer *coutput6=mgr->CreateContainer("fTree_cascade",TTree::Class(),
								AliAnalysisManager::kOutputContainer,
								AliAnalysisManager::GetCommonFileName());

	AliAnalysisDataContainer *coutput7=mgr->CreateContainer("fTree_proton",TTree::Class(),
								AliAnalysisManager::kOutputContainer,
								AliAnalysisManager::GetCommonFileName());

	AliAnalysisDataContainer *coutput8=mgr->CreateContainer("fTree_antiproton",TTree::Class(),
								AliAnalysisManager::kOutputContainer,
								AliAnalysisManager::GetCommonFileName());

	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);
	mgr->ConnectOutput(task,2,coutput2);
	mgr->ConnectOutput(task,3,coutput3);
	mgr->ConnectOutput(task,4,coutput4);
	mgr->ConnectOutput(task,5,coutput5);
	mgr->ConnectOutput(task,6,coutput6);
	mgr->ConnectOutput(task,7,coutput7);
	mgr->ConnectOutput(task,8,coutput8);

	return task;
}
