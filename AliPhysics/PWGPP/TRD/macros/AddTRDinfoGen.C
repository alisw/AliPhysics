void AddTRDinfoGen(AliAnalysisDataContainer **co)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) return;

  Bool_t mc = (Bool_t)mgr->GetMCtruthEventHandler();
  //AliLog::SetClassDebugLevel("AliTRDinfoGen", 2);
  AliTRDinfoGen *info(NULL);
  mgr->AddTask(info = new AliTRDinfoGen((char*)"TRDinfoGen"));
  info->SetDebugLevel(0);
  info->SetMCdata(mc);
  info->SetLocalTrkSelection();
  info->UseTrackPoints(kFALSE); // set it to true if track points for alignment are to be saved in trackInfo object
//  info->SetOCDB("alien://folder=/alice/data/2012/OCDB?cacheFolder=/home/niham/abercuci/local");
  info->SetOCDB(Form("local://%s/local/alice/data/%d/OCDB", gSystem->ExpandPathName("$HOME"), AliTRDpwgppHelper::GetRunYear()));
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
  co[AliTRDpwgppHelper::kEventInfo] = mgr->CreateContainer("eventInfo", AliTRDeventInfo::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliTRDpwgppHelper::kTracksBarrel] = mgr->CreateContainer("tracksBarrel", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliTRDpwgppHelper::kTracksITS] = mgr->CreateContainer("tracksITS", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliTRDpwgppHelper::kTracksSA] = mgr->CreateContainer("tracksSA", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliTRDpwgppHelper::kTracksKink] = mgr->CreateContainer("tracksKink", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliTRDpwgppHelper::kV0List] = mgr->CreateContainer("v0List", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliTRDpwgppHelper::kTracklets] = mgr->CreateContainer("onl.tracklets", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  co[AliTRDpwgppHelper::kClusters] = mgr->CreateContainer("clusters", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  for(Int_t ios(1);ios<AliTRDpwgppHelper::kNOutSlots-1;ios++) mgr->ConnectOutput(info, ios, co[ios]);
  
  // add last monitor container
  AliAnalysisDataContainer *mon=mgr->CreateContainer(info->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance",mgr->GetCommonFileName()));
  mgr->ConnectOutput(info, AliTRDpwgppHelper::kNOutSlots-1, mon);
}

