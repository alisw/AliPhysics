AliPHOSEventCuts *CreatePHOSEventCuts(Bool_t isMC, AliPHOSTriggerHelper *obj=0x0)
{

  AliPHOSEventCuts *cuts = new AliPHOSEventCuts("PHOSEventCuts");
  cuts->SetMCFlag(isMC);
  cuts->SetMaxAbsZvtx(10.);
  cuts->SetRejectPileup(kTRUE);
  cuts->SetRejectDAQIncompleteEvent(kTRUE);
  cuts->DoPHOSTriggerAnalysis(kFALSE);

  if(obj){
    //this is necessary for only kPHI7 analysis.
    cuts->DoPHOSTriggerAnalysis(kTRUE,obj);
  }

  return cuts;

}
//________________________________________________
