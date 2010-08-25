#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWG1/TRD/AliTRDcheckESD.h"
#endif

void AddTRDcheckESD(AliAnalysisManager *mgr)
{
  //AliLog::SetClassDebugLevel("AliTRDcheckESD", 5);
  AliTRDcheckESD *checkESD = new AliTRDcheckESD((char*)"TRDcheckESD");
  mgr->AddTask(checkESD);
  Bool_t mc = mgr->GetMCtruthEventHandler();
  checkESD->SetMC(mc);
  checkESD->SetCollision(/*kFALSE*/);
  checkESD->SetDebugLevel(0);

  mgr->ConnectInput(checkESD,  0, mgr->GetCommonInputContainer());  
  mgr->ConnectOutput(checkESD, 1, mgr->CreateContainer(checkESD->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance",mgr->GetCommonFileName())));
}

