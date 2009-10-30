#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "TRD/qaRec/AliTRDcheckESD.h"
#endif

void AddTRDcheckESD(AliAnalysisManager *mgr)
{
  AliTRDcheckESD *checkESD = new AliTRDcheckESD();
  checkESD->SetMC(mgr->GetMCtruthEventHandler());
  mgr->AddTask(checkESD);
  mgr->ConnectInput(checkESD, 0, mgr->GetCommonInputContainer());  mgr->ConnectOutput(checkESD, 0, mgr->CreateContainer(checkESD->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.Performance.root"));
}

