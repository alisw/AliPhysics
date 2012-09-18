//
// Usage :
//  Select the function that adapts to your needs, then
//  AliOADBPhysicsSelection *oadb = function();
//  physSelTask->GetPhysicsSelection()->SetCustomOADBObjects(oadb,0);
//

AliOADBPhysicsSelection *OADBSelection_CINT5_V0OR(){

  AliOADBPhysicsSelection * oadbLHC12g = new AliOADBPhysicsSelection("oadbLHC12g");
  oadbLHC12g->AddCollisionTriggerClass   ( AliVEvent::kCINT5,"+CINT5-B-NOPF-ALLNOTRD","B",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-ACE-NOPF-ALLNOTRD","ACE",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-A-NOPF-ALLNOTRD","A",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-C-NOPF-ALLNOTRD","C",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC12g->SetHardwareTrigger         ( 0,"V0A || V0C");
  oadbLHC12g->SetOfflineTrigger          ( 0,"(V0A || V0C) && !V0ABG && !V0CBG && !ZNABG && !ZNCBG && !TPCLaserWarmUp");

  return oadbLHC12g;
}

AliOADBPhysicsSelection *OADBSelection_CINT5_V0A(){

  AliOADBPhysicsSelection * oadbLHC12g = new AliOADBPhysicsSelection("oadbLHC12g");
  oadbLHC12g->AddCollisionTriggerClass   ( AliVEvent::kCINT5,"+CINT5-B-NOPF-ALLNOTRD","B",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-ACE-NOPF-ALLNOTRD","ACE",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-A-NOPF-ALLNOTRD","A",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-C-NOPF-ALLNOTRD","C",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC12g->SetHardwareTrigger         ( 0,"V0A || V0C");
  oadbLHC12g->SetOfflineTrigger          ( 0,"V0A && !V0ABG && !V0CBG && !ZNABG && !ZNCBG && !TPCLaserWarmUp");

  return oadbLHC12g;
}

AliOADBPhysicsSelection *OADBSelection_CINT5_V0AND(){

  AliOADBPhysicsSelection * oadbLHC12g = new AliOADBPhysicsSelection("oadbLHC12g");
  oadbLHC12g->AddCollisionTriggerClass   ( AliVEvent::kCINT5,"+CINT5-B-NOPF-ALLNOTRD","B",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-AC-NOPF-ALLNOTRD","AC",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-ACE-NOPF-ALLNOTRD","ACE",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-A-NOPF-ALLNOTRD","A",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-C-NOPF-ALLNOTRD","C",0);
  oadbLHC12g->AddBGTriggerClass          ( AliVEvent::kCINT5,"+CINT5-E-NOPF-ALLNOTRD","E",0);  
  oadbLHC12g->SetHardwareTrigger         ( 0,"V0A || V0C");
  oadbLHC12g->SetOfflineTrigger          ( 0,"(V0A && V0C) && !V0ABG && !V0CBG && !ZNABG && !ZNCBG && !TPCLaserWarmUp");

  return oadbLHC12g;
}
