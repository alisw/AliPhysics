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
  Bool_t mc(mgr->GetMCtruthEventHandler()?kTRUE:kFALSE);
  //AliLog::SetClassDebugLevel("AliTRDinfoGen", 2);
  AliTRDinfoGen *info(NULL);
  mgr->AddTask(info = new AliTRDinfoGen((char*)"genInfo"));
  info->SetDebugLevel(0);
  info->SetMCdata(mc);
  info->SetLocalTrkSelection();

  // settings for collisions
  info->SetCollision(/*kFALSE*/);
  if(info->IsCollision()){
    if(!mc) info->SetTrigger(
      "CINT1B-ABCE-NOPF-ALL"
      " CSCO1-ABCE-NOPF-CENT" // cosmic SPD trigger
    );
    info->SetLocalEvSelection();
  }
  
  // Connect IO slots
  mgr->ConnectInput (info, 0, mgr->GetCommonInputContainer());
  co[kEventInfo] = mgr->CreateContainer("eventInfo", AliTRDeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
  co[kTracksBarrel] = mgr->CreateContainer("tracksBarrel", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[kTracksSA] = mgr->CreateContainer("tracksSA", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[kTracksKink] = mgr->CreateContainer("tracksKink", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[kV0List] = mgr->CreateContainer("v0List", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  for(Int_t ios(1);ios<kNOutSlots-1;ios++) mgr->ConnectOutput(info, ios, co[ios]);
  
  // add last monitor container
  AliAnalysisDataContainer *mon=mgr->CreateContainer("infoGen", TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance",mgr->GetCommonFileName()));
  mgr->ConnectOutput(info, kNOutSlots-1, mon);
}

