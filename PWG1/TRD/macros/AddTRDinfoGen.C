#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWG1/TRD/AliTRDinfoGen.h"
#include "PWG1/TRD/AliTRDpwg1Helper.h"
#include "PWG1/TRD/info/AliTRDeventInfo.h"
#endif

void AddTRDinfoGen(AliAnalysisManager *mgr, Int_t /*map*/, AliAnalysisDataContainer **/*ci*/, AliAnalysisDataContainer **co)
{
  Bool_t mc=mgr->GetMCtruthEventHandler();
  //AliLog::SetClassDebugLevel("AliTRDinfoGen", 2);
  AliTRDinfoGen *info(NULL);
  mgr->AddTask(info = new AliTRDinfoGen((char*)"TRDinfoGen"));
  info->SetDebugLevel(0);
  info->SetMCdata(mc);
  info->SetLocalTrkSelection();
  info->SetOCDB("alien://folder=/alice/data/2010/OCDB");
  // settings for collisions
  info->SetCollision(/*kFALSE*/);
  if(info->IsCollision()){
    if(!mc) info->SetTrigger(
      "CINT1B-ABCE-NOPF-ALL"
      " CINT5-B-NOPF-ALL"
      " CINT1WU-B-NOPF-ALL"
      " CINT5WU-B-NOPF-ALL"
      " CSCO1-ABCE-NOPF-CENT" // cosmic SPD trigger
    );
    info->SetLocalEvSelection();
  }
  
  // Connect IO slots
  mgr->ConnectInput (info, 0, mgr->GetCommonInputContainer());
  co[AliTRDpwg1Helper::kEventInfo] = mgr->CreateContainer("eventInfo", AliTRDeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliTRDpwg1Helper::kTracksBarrel] = mgr->CreateContainer("tracksBarrel", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliTRDpwg1Helper::kTracksSA] = mgr->CreateContainer("tracksSA", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliTRDpwg1Helper::kTracksKink] = mgr->CreateContainer("tracksKink", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliTRDpwg1Helper::kV0List] = mgr->CreateContainer("v0List", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  for(Int_t ios(1);ios<AliTRDpwg1Helper::kNOutSlots-1;ios++) mgr->ConnectOutput(info, ios, co[ios]);
  
  // add last monitor container
  AliAnalysisDataContainer *mon=mgr->CreateContainer(info->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance",mgr->GetCommonFileName()));
  mgr->ConnectOutput(info, AliTRDpwg1Helper::kNOutSlots-1, mon);
}

