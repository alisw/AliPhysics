#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWG1/TRD/AliTRDcheckESD.h"
#endif

void AddTRDcheckESD(AliAnalysisManager *mgr)
{
  AliTRDcheckESD *checkESD = new AliTRDcheckESD("checkESD");
  Bool_t hasMCtruth = (mgr->GetMCtruthEventHandler() != 0);
  
  checkESD->SetMC(hasMCtruth);
  mgr->AddTask(checkESD);
  mgr->ConnectInput(checkESD,  0, mgr->GetCommonInputContainer());  
  mgr->ConnectOutput(checkESD, 1, mgr->CreateContainer(checkESD->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.Performance.root"));
}

