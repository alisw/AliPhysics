#if ! defined (__CINT__) || defined (__MAKECINT__)
#include <TTree.h>
#include <TError.h>
#include <AliLog.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliMESbaseTask.h>
#include <AliMEStender.h>
#endif

AliMEStender *AddMEStender(Bool_t mc)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliMEStender *tender = new AliMEStender((char*)"MEStender");
  mgr->AddTask(tender);

  // task set-up
  tender->SetMCdata(mc);
  tender->SetDebugLevel(0);
  tender->ConfigTask( AliMEStender::AliMESconfigTender::kStandard,                    // event cuts
					  AliMEStender::AliMESconfigTender::kStandardITSTPCTrackCuts2010, // track cuts
					  AliMEStender::AliMESconfigTender::kIterative);                  // PID priors
  tender->SetPriors();  // always call this after ConfigTask !!

  // connect input
  mgr->ConnectInput (tender, 0, mgr->GetCommonInputContainer());

  // create output containers
  AliAnalysisDataContainer *co[AliMESbaseTask::kNcontainers] = {NULL};
  co[0]                          = mgr->CreateContainer("tenderQA", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:MES", mgr->GetCommonFileName()));
  co[AliMESbaseTask::kEventInfo] = mgr->CreateContainer("MESeventInfo", AliMESeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliMESbaseTask::kTracks]    = mgr->CreateContainer("MEStracks", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  if(mc){
	  co[AliMESbaseTask::kMCeventInfo] = mgr->CreateContainer("MESMCeventInfo", AliMESeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
	  co[AliMESbaseTask::kMCtracks]    = mgr->CreateContainer("MESMCtracks", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  }

  // connect output
  for(Int_t ios(0);ios<AliMESbaseTask::kNcontainers;ios++)
	  if(co[ios]) mgr->ConnectOutput(tender, ios+1, co[ios]);


  return tender;

}

