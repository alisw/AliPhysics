#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TError.h>
#include <AliLog.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliMESbaseTask.h>
#include <AliMEStenderV2.h>
#endif

AliMEStenderV2 *AddMEStenderV2(Bool_t mc, Int_t configuration = 1, const AliVEvent *event)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliCDBManager *man = AliCDBManager::Instance();
  if (man)
  {
    if (!man->IsDefaultStorageSet())
    {
      man->SetDefaultStorage(Form("alien://folder=/alice/data/2015/OCDB?cacheFolder=%s/local", gSystem->ExpandPathName("$HOME")));
      man->SetRun(event->GetRunNumber());
      AliGRPObject *grp = (AliGRPObject *)man->Get("GRP/GRP/Data")->GetObject();
    }
  }

  AliMEStenderV2 *tender = new AliMEStenderV2((char *)"MEStenderV2");
  mgr->AddTask(tender);

  // task set-up
  tender->SetMCdata(mc);
  tender->SetDebugLevel(1);
  switch (configuration)
  {
  case 0:
    tender->ConfigTask(AliMEStenderV2::AliMESconfigTender::k7TeV,                        // event cuts
                       AliMEStenderV2::AliMESconfigTender::kStandardITSTPCTrackCuts2010, // track cuts
                       AliMEStenderV2::AliMESconfigTender::kIterative);                  // PID priors
    break;
  case 1:
    tender->ConfigTask(AliMEStenderV2::AliMESconfigTender::k13TeV,                       // event cuts
                       AliMEStenderV2::AliMESconfigTender::kStandardITSTPCTrackCuts2011, // track cuts
                       AliMEStenderV2::AliMESconfigTender::kNoPP);                       // PID priors
    break;
  default:
    printf("Configuration not defined\n");
    break;
  }
  tender->SetPriors(); // always call this after ConfigTask !!

  // connect input
  mgr->ConnectInput(tender, 0, mgr->GetCommonInputContainer());

  // create output containers
  AliAnalysisDataContainer *co[AliMESbaseTask::kNcontainers] = {NULL};
  co[0] = mgr->CreateContainer("tenderQA", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:MES", mgr->GetCommonFileName()));
  co[AliMESbaseTask::kEventInfo] = mgr->CreateContainer("MESeventInfo", AliMESeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliMESbaseTask::kTracks] = mgr->CreateContainer("MEStracks", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  if (mc)
  {
    co[AliMESbaseTask::kMCeventInfo] = mgr->CreateContainer("MESMCeventInfo", AliMESeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
    co[AliMESbaseTask::kMCtracks] = mgr->CreateContainer("MESMCtracks", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  }

  // connect output
  for (Int_t ios(0); ios < AliMESbaseTask::kNcontainers; ios++)
    if (co[ios])
      mgr->ConnectOutput(tender, ios + 1, co[ios]);

  return tender;
}
