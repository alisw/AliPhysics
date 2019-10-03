#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TList.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTask.h"
#include "AliTRDonlineTrackletQA.h"
#endif

AliAnalysisTask *AddTRDonlineTrackletQA(AliAnalysisManager *mgr)
{
  if (!mgr) {
    cerr << "No Analysis manager available" << endl;
    return 0x0;
  }

  AliTRDonlineTrackletQA *task = new AliTRDonlineTrackletQA("TRDtrackletQA");
  mgr->AddTask(task);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->GetCommonOutputContainer();
  AliAnalysisDataContainer *ctracklets =
    (AliAnalysisDataContainer*) AliAnalysisManager::GetAnalysisManager()->GetContainers()->FindObject("TRDtrackletFilter");

  AliAnalysisDataContainer *ctrklqa =
    mgr->CreateContainer("TRDtrackletQA", TList::Class(),
			 AliAnalysisManager::kOutputContainer,
			 Form("%s:TRD_Performance", mgr->GetCommonFileName()));

  task->ConnectInput(0, cinput);
  task->ConnectInput(1, ctracklets);

  task->ConnectOutput(1, ctrklqa);

  return task;
}
