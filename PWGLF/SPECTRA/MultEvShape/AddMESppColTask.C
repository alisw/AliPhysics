#if ! defined (__CINT__) || defined (__MAKECINT__)
#include <TTree.h>
#include <TError.h>
#include <AliLog.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliMESbaseTask.h>
#include <AliMESppColTask.h>
#endif

AliMESppColTask *AddMESppColTask(Bool_t mc)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliMESppColTask *ppCol = new AliMESppColTask("MESppCol");
  mgr->AddTask(ppCol);
  ppCol->SetPostProcess(kFALSE);
  ppCol->SetMCdata(mc);
  ppCol->SetDebugLevel(1);

  mgr->ConnectInput(ppCol, 0, mgr->GetCommonInputContainer()); // connect main (ESD) container

  // get exchange containers
  AliAnalysisDataContainer *ci[AliMESbaseTask::kNcontainers] = {NULL};
  ci[0] = (AliAnalysisDataContainer*)(AliAnalysisManager::GetAnalysisManager()->GetContainers()->FindObject("MESeventInfo"));
  ci[1] = (AliAnalysisDataContainer*)(AliAnalysisManager::GetAnalysisManager()->GetContainers()->FindObject("MEStracks"));
  if(mc){
	  ci[2] = (AliAnalysisDataContainer*)(AliAnalysisManager::GetAnalysisManager()->GetContainers()->FindObject("MESMCeventInfo"));
	  ci[3] = (AliAnalysisDataContainer*)(AliAnalysisManager::GetAnalysisManager()->GetContainers()->FindObject("MESMCtracks"));
  }
  Info("AddMESppColTask", Form("Inputs : [0]=\"%s\" [1]=\"%s\" [2]=\"%s\" [3]=\"%s\"",
							   ci[0]?ci[0]->GetName():"none",
							   ci[1]?ci[1]->GetName():"none",
							   ci[2]?ci[2]->GetName():"none",
							   ci[3]?ci[3]->GetName():"none"));

  // connect containers
  mgr->ConnectInput(ppCol, AliMESbaseTask::kEventInfo, ci[0]); // connect event info
  mgr->ConnectInput(ppCol, AliMESbaseTask::kTracks, ci[1]);    // connect track info container
  if(mc){
    mgr->ConnectInput(ppCol, AliMESbaseTask::kMCeventInfo, ci[2]); // connect MC event info
    mgr->ConnectInput(ppCol, AliMESbaseTask::kMCtracks, ci[3]);    // connect MC tracks container
  }
  mgr->ConnectOutput(ppCol, AliMESbaseTask::kQA, mgr->CreateContainer("ppColQA", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:MES", mgr->GetCommonFileName())));

  return ppCol;
}

