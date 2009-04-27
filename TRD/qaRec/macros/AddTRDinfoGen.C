#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "TRD/qaRec/AliTRDinfoGen.h"
#include "TRD/qaRec/info/AliTRDeventInfo.h"
#include "TRD/qaRec/macros/AliTRDperformanceTrain.h"
#endif

void AddTRDinfoGen(AliAnalysisManager *mgr, Char_t *trd, AliAnalysisDataContainer **/*ci*/, AliAnalysisDataContainer **co)
{
  Int_t map = ParseOptions(trd);
  if(!(TSTBIT(map, kInfoGen))) return;

  AliTRDinfoGen *info = 0x0;
  mgr->AddTask(info = new AliTRDinfoGen());
  info->SetDebugLevel(0);
  info->SetMCdata(mgr->GetMCtruthEventHandler());
  mgr->ConnectInput( info, 0, mgr->GetCommonInputContainer());
  co[0] = mgr->CreateContainer("trackInfo", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[1] = mgr->CreateContainer("eventInfo", AliTRDeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
  mgr->ConnectOutput(info, 0, co[0]);
  mgr->ConnectOutput(info, 1, co[1]);
}

