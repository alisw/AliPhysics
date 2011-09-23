#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisCuts.h"
#include "AliESDtrackCuts.h"
#include "PWG1/TRD/AliTRDcheckESD.h"
#endif

AliESDtrackCuts* SetupESDcuts();

void AddTRDcheckESD(AliAnalysisManager *mgr)
{
  //AliLog::SetClassDebugLevel("AliTRDcheckESD", 5);
  AliTRDcheckESD *checkESD = new AliTRDcheckESD((char*)"TRDcheckESD");
  checkESD->SetRefTrackFilter(SetupESDcuts());
  mgr->AddTask(checkESD);
  Bool_t mc = mgr->GetMCtruthEventHandler();
  checkESD->SetMC(mc);
  checkESD->SetCollision(/*kFALSE*/);
  checkESD->SetDebugLevel(0);

  mgr->ConnectInput(checkESD,  0, mgr->GetCommonInputContainer());  
  mgr->ConnectOutput(checkESD, 1, mgr->CreateContainer(checkESD->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance",mgr->GetCommonFileName())));
}

AliESDtrackCuts* SetupESDcuts() {
  // Setup ESD cuts for the TPC reference tracks
  AliESDtrackCuts* esdCuts = new AliESDtrackCuts;
  esdCuts->SetPtRange(0.2, 100.0);
  esdCuts->SetEtaRange(-0.9, +0.9);
  esdCuts->SetRequireTPCRefit(kTRUE);
  esdCuts->SetAcceptKinkDaughters(kFALSE);
  esdCuts->SetMaxDCAToVertexXY(40.);
  esdCuts->SetMaxDCAToVertexZ(15.);
  esdCuts->SetMinNClustersTPC(70);
  return esdCuts;
}

