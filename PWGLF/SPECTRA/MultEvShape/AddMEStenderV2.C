#if ! defined (__CINT__) || defined (__MAKECINT__)
#include <TTree.h>
#include <TError.h>
#include <AliLog.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliMESbaseTask.h>
#include <AliMEStenderV2.h>
#include <AliMESeventInfo.h>
#endif

AliMEStenderV2 *AddMEStenderV2(Bool_t mc, Int_t configuration = 1)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  AliMEStenderV2 *tender = new AliMEStenderV2((char*)"MEStenderV2");
  mgr->AddTask(tender);

  // task set-up
  tender->SetMCdata(mc);
  tender->SetDebugLevel(1);
	switch (configuration) {
		case 0:
	 		tender->ConfigTask( AliMEStenderV2::AliMESconfigTender::k7TeV,      // event cuts
					  AliMEStenderV2::AliMESconfigTender::kStandardITSTPCTrackCuts2010, // track cuts
					  AliMEStenderV2::AliMESconfigTender::kIterative);                  // PID priors
  		break;
		case 1:
			tender->ConfigTask( AliMEStenderV2::AliMESconfigTender::k13TeV,       // event cuts
					AliMEStenderV2::AliMESconfigTender::kStandardITSTPCTrackCuts2011, // track cuts
					AliMEStenderV2::AliMESconfigTender::kNoPP);                  // PID priors
			break;
		default: printf("Configuration not defined\n");
			break;
	}
	tender->SetPriors();  // always call this after ConfigTask !!

  // connect input
  mgr->ConnectInput (tender, 0, mgr->GetCommonInputContainer());

  // create output containers
  AliAnalysisDataContainer *co[AliMEStenderV2::kNcontainers] = {NULL};
  co[0]                          = mgr->CreateContainer("tenderQA", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:MES", mgr->GetCommonFileName()));
  co[AliMEStenderV2::kEventInfo] = mgr->CreateContainer("MESeventInfo", AliMESeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliMEStenderV2::kTracks]    = mgr->CreateContainer("MEStracks", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliMEStenderV2::kTree] = mgr->CreateContainer("MEStree", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  if(mc){
	  co[AliMEStenderV2::kMCeventInfo] = mgr->CreateContainer("MESMCeventInfo", AliMESeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
	  co[AliMEStenderV2::kMCtracks]    = mgr->CreateContainer("MESMCtracks", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
    co[AliMEStenderV2::kTreeMC]      = mgr->CreateContainer("MESMCtree", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName() );
    co[AliMEStenderV2::kTreeGen]     = mgr->CreateContainer("MESGen", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
    co[AliMEStenderV2::kTreeMiss]    = mgr->CreateContainer("MESMiss", TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName());
  }

  // connect output
  for(Int_t ios(0);ios<AliMEStenderV2::kNcontainers;ios++)
	  if(co[ios]) mgr->ConnectOutput(tender, ios+1, co[ios]);


  return tender;

}
