#if ! defined (__CINT__) || defined (__MAKECINT__)
#include <TTree.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliAnalysisTask.h>
#include <AliTRDonlineTrackletFilter.h>
#endif

AliAnalysisTask *AddTRDonlineTrackletFilter(AliAnalysisManager *mgr)
{
  if (!mgr) {
    cerr << "No Analysis manager available" << endl;
    return 0x0;
  }

  AliTRDonlineTrackletFilter *task = new AliTRDonlineTrackletFilter("TRDtrackletfilter");
  mgr->AddTask(task);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  //AliAnalysisDataContainer *coutput = mgr->GetCommonOutputContainer();

  task->ConnectInput(0, cinput);
  task->ConnectOutput(0, mgr->CreateContainer("sim.onl.tracklets", TClonesArray::Class(), AliAnalysisManager::kExchangeContainer));
/*  AliAnalysisDataContainer *ctracklets =
    mgr->CreateContainer("TRDtrackletFilter", TTree::Class(),
                         AliAnalysisManager::kOutputContainer,
                         Form("%s:TRD_Performance", mgr->GetCommonFileName()));

  task->ConnectOutput(1, ctracklets);
*/
  return task;
}
