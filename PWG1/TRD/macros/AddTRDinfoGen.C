#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWG1/TRD/AliTRDinfoGen.h"
#include "PWG1/TRD/AliTRDpwg1Helper.h"
#include "PWG1/TRD/info/AliTRDeventInfo.h"
#include "PWG1/TRD/info/AliTRDeventCuts.h"
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
    AliTRDeventCuts ec;
    if(mc){
      ec.SetEventType(0);
//      ec.AddTrigger("MB1");
    } else {
/*      Int_t bunches[] = ;
      ec.SetBunchSelection();*/
/*      ec.AddTrigger("CINT1B-ABCE-NOPF-ALL");
      ec.AddTrigger("CINT5-B-NOPF-ALL");
      ec.AddTrigger("CINT1WU-B-NOPF-ALL");
      ec.AddTrigger("CINT5WU-B-NOPF-ALL");
      ec.AddTrigger("CINT7WU-I-NOPF-ALL");
      ec.AddTrigger("CSCO1-ABCE-NOPF-CENT"); // cosmic SPD trigger
*/
    }
    info->SetLocalEvSelection(ec);
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

