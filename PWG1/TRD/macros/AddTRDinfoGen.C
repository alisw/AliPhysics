#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWG1/TRD/AliTRDinfoGen.h"
#include "PWG1/TRD/info/AliTRDeventInfo.h"
#include "PWG1/TRD/macros/AliTRDperformanceTrain.h"
#endif

#include "PWG1/TRD/macros/helper.C"
void AddTRDinfoGen(AliAnalysisManager *mgr, Char_t *trd, AliAnalysisDataContainer **/*ci*/, AliAnalysisDataContainer **co)
{
  Int_t map = ParseOptions(trd);
  if(!(TSTBIT(map, kInfoGen))) return;
  
  //AliLog::SetClassDebugLevel("AliTRDinfoGen", 5);  
  AliTRDinfoGen *info(NULL);
  info = new AliTRDinfoGen("genInfo");
  mgr->AddTask(info);
  info->SetDebugLevel(0);
  info->SetMCdata(mgr->GetMCtruthEventHandler());
  AliAnalysisDataContainer* cin   = mgr->CreateContainer("dummy", TObjArray::Class(), AliAnalysisManager::kInputContainer);
  

  mgr->ConnectInput( info, 0, mgr->GetCommonInputContainer());
  // settings for collisions
  info->SetCollision();
  if(info->IsCollision()){
    info->SetTrigger(
      "CINT1B-ABCE-NOPF-ALL"
      " CSCO1-ABCE-NOPF-CENT" // cosmic SPD trigger
    );
    info->SetLocalEvSelection();
    info->SetLocalTrkSelection();
  }
  co[0] = mgr->CreateContainer("trackInfo", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[1] = mgr->CreateContainer("eventInfo", AliTRDeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
  co[2] = mgr->CreateContainer("v0Info",    TObjArray::Class(),       AliAnalysisManager::kExchangeContainer);

  mgr->ConnectInput (info, 0, mgr->GetCommonInputContainer());
  mgr->ConnectInput (info, 1, cin  );   // Dummy to avoid orphan
  mgr->ConnectOutput(info, 1, co[0]);
  mgr->ConnectOutput(info, 2, co[1]);
  mgr->ConnectOutput(info, 3, co[2]);
}

