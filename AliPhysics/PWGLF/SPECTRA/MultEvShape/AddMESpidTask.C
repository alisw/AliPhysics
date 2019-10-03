#if ! defined (__CINT__) || defined (__MAKECINT__)
#include <TTree.h>
#include <TError.h>
#include <AliLog.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliMESbaseTask.h>
#include <AliMESpidTask.h>
#endif

AliMESpidTask *AddMESpidTask(Bool_t mc)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliMESpidTask *pid = new AliMESpidTask("MESpid");
  mgr->AddTask(pid);
  pid->SetPostProcess(kFALSE);
  pid->SetMCdata(mc);
  pid->SetDebugLevel(0);

  mgr->ConnectInput(pid, 0, mgr->GetCommonInputContainer()); // connect main (ESD) container

  // get exchange containers
  AliAnalysisDataContainer *ci[AliMESbaseTask::kNcontainers] = {NULL};
  ci[0] = (AliAnalysisDataContainer*)(AliAnalysisManager::GetAnalysisManager()->GetContainers()->FindObject("MESeventInfo"));
  ci[1] = (AliAnalysisDataContainer*)(AliAnalysisManager::GetAnalysisManager()->GetContainers()->FindObject("MEStracks"));
  if(mc){
	  ci[2] = (AliAnalysisDataContainer*)(AliAnalysisManager::GetAnalysisManager()->GetContainers()->FindObject("MESMCeventInfo"));
	  ci[3] = (AliAnalysisDataContainer*)(AliAnalysisManager::GetAnalysisManager()->GetContainers()->FindObject("MESMCtracks"));
  }
  Info("AddMESpidTask", Form("Inputs : [0]=\"%s\" [1]=\"%s\" [2]=\"%s\" [3]=\"%s\"",
							 ci[0]?ci[0]->GetName():"none",
							 ci[1]?ci[1]->GetName():"none",
							 ci[2]?ci[2]->GetName():"none",
							 ci[3]?ci[3]->GetName():"none"));

  // connect containers
  mgr->ConnectInput(pid, AliMESbaseTask::kEventInfo, ci[0]); // connect event info
  mgr->ConnectInput(pid, AliMESbaseTask::kTracks, ci[1]);    // connect track info container
  if(mc){
    mgr->ConnectInput(pid, AliMESbaseTask::kMCeventInfo, ci[2]); // connect MC event info
    mgr->ConnectInput(pid, AliMESbaseTask::kMCtracks, ci[3]);    // connect MC tracks container
  }
  mgr->ConnectOutput(pid, AliMESbaseTask::kQA, mgr->CreateContainer("pidQA", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:MES", mgr->GetCommonFileName())));

  return pid;
}

