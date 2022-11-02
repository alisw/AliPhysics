#include "AliAnalysisTaskMSDibaryons.h"

AliAnalysisTaskMSDibaryons *AddTaskMSDibaryons(
		const TString taskname = "MSDibaryons",
		//const UInt_t trigger = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral,	
		const UInt_t trigger = AliVEvent::kCentral,	
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
	AliAnalysisTaskMSDibaryons *task=new AliAnalysisTaskMSDibaryons(taskname.Data());
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
	mgr->ConnectInput(task,0,cinput);
	mgr->ConnectOutput(task,1,coutput);

	return task;

}
