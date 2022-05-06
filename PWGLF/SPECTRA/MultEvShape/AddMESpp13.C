#ifdef __CLING__
#include <TTree.h>
#include <TError.h>
#include <AliLog.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliMESbaseTask.h>
#include <AliMESpp13.h>
#include <AliMESeventInfo.h>
#endif

AliMESpp13 *AddMESpp13(Bool_t mc, Int_t configuration = 1)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliMESpp13 *task = new AliMESpp13((char *)"MESpp13");
  mgr->AddTask(task);

  // task set-up
  task->SetMCdata(mc);
  switch (configuration)
  {
  case 0:
    task->ConfigTask(AliMESpp13::AliMESconfigTender::k7TeV,                        // event cuts
                     AliMESpp13::AliMESconfigTender::kStandardITSTPCTrackCuts2010, // track cuts
                     AliMESpp13::AliMESconfigTender::kIterative);                  // PID priors
    break;
  case 1:
    task->ConfigTask(AliMESpp13::AliMESconfigTender::k13TeV,                       // event cuts
                     AliMESpp13::AliMESconfigTender::kStandardITSTPCTrackCuts2011, // track cuts
                     AliMESpp13::AliMESconfigTender::kNoPP);                       // PID priors
    break;
  default:
    printf("Configuration not defined\n");
    break;
  }
  task->SetPriors(); // always call this after ConfigTask !!

  // connect input
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // create output containers
  AliAnalysisDataContainer *co[AliMESpp13::kNcontainers] = {NULL};
  co[0] = mgr->CreateContainer("taskQA", TList::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  co[AliMESpp13::kEventInfo] = mgr->CreateContainer("MESeventInfo", AliMESeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliMESpp13::kTracks] = mgr->CreateContainer("MEStracks", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliMESpp13::kEventTree] = mgr->CreateContainer("MES-ev", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  co[AliMESpp13::kTracksTree] = mgr->CreateContainer("MES-trk", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  if (mc)
  {
    co[AliMESpp13::kMCeventInfo] = mgr->CreateContainer("MESMCeventInfo", AliMESeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
    co[AliMESpp13::kMCtracks] = mgr->CreateContainer("MESMCtracks", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
    co[AliMESpp13::kMCeventTree] = mgr->CreateContainer("MES-MCev", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
    co[AliMESpp13::kMCtracksTree] = mgr->CreateContainer("MES-MCtrk", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
    co[AliMESpp13::kMCGenTracksTree] = mgr->CreateContainer("MES-genTrk", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
    co[AliMESpp13::kMCMissedTracksTree] = mgr->CreateContainer("MES-missedTrk", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  }

  // connect output
  for (Int_t ios(0); ios < AliMESpp13::kNcontainers; ios++)
    if (co[ios])
      mgr->ConnectOutput(task, ios + 1, co[ios]);

  return task;
}
