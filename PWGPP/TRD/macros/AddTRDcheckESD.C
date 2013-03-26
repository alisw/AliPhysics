
#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisCuts.h"
#include "AliESDtrackCuts.h"
#include "PWGPP/TRD/AliTRDcheckESD.h"
#endif

AliESDtrackCuts* SetupESDcuts();

void AddTRDcheckESD(AliAnalysisManager *mgr)
{
  //AliLog::SetClassDebugLevel("AliTRDcheckESD", 5);
  //  AliInfo("aaaaaa6666666666");
  AliTRDcheckESD *checkESD = new AliTRDcheckESD((char*)"TRDcheckESD");
  checkESD->SetRefTrackFilter(SetupESDcuts());
  mgr->AddTask(checkESD);
  Bool_t mc = mgr->GetMCtruthEventHandler();
  checkESD->SetMC(mc);
  checkESD->SetCollision(/*kFALSE*/);
  checkESD->SetDebugLevel(0);
  checkESD->AddUserTrigger("WU");
  checkESD->AddUserTrigger("QU");
  checkESD->AddUserTrigger("SE");
  checkESD->AddUserTrigger("JT");
  
  mgr->ConnectInput(checkESD,  0, mgr->GetCommonInputContainer());  
  mgr->ConnectOutput(checkESD, 1, mgr->CreateContainer(checkESD->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance",mgr->GetCommonFileName())));
}

AliESDtrackCuts* SetupESDcuts() {
  // Setup ESD cuts for the TPC reference tracks
  AliESDtrackCuts* esdCuts = new AliESDtrackCuts;
  esdCuts->SetPtRange(0.5, 100.0);
  esdCuts->SetEtaRange(-0.9, +0.9);
  esdCuts->SetRequireTPCRefit(kTRUE);
  esdCuts->SetAcceptKinkDaughters(kFALSE);
  esdCuts->SetMaxDCAToVertexXY(1.);
  esdCuts->SetMaxDCAToVertexZ(3.);
  esdCuts->SetMinNClustersTPC(70);
  esdCuts->SetRequireITSRefit(kTRUE);
  esdCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  return esdCuts;
}

