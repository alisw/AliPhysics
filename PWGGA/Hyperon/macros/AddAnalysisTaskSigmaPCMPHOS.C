#include "AliVEvent.h"

AliAnalysisTaskSigmaPCMPHOS *AddAnalysisTaskSigmaPCMPHOS(TString name = "SigmaPCMPHOS")
{
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	if (!mgr || !mgr->GetInputEventHandler())
		return 0x0;

	TString fileName = AliAnalysisManager::GetCommonFileName();
	fileName += ":MyTask"; // creates a folder in the ROOT file

	AliAnalysisTaskSigmaPCMPHOS *task = new AliAnalysisTaskSigmaPCMPHOS(name.Data());

	if (!task)
		return 0x0;
	mgr->AddTask(task);
	task->SelectCollisionCandidates(AliVEvent::kINT7); // minimum vias trigger V0M

	mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task, 1, mgr->CreateContainer("Histos", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
	mgr->ConnectOutput(task, 2, mgr->CreateContainer("Cuts", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
	mgr->ConnectOutput(task, 3, mgr->CreateContainer("MC", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

	return task;
}