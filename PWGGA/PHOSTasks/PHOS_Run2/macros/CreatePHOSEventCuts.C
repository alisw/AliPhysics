AliPHOSEventCuts *CreatePHOSEventCuts(Bool_t isMC)
{

  AliPHOSEventCuts *cuts = new AliPHOSEventCuts("PHOSEventCuts");
  cuts->SetMCFlag(isMC);
  cuts->SetMaxAbsZvtx(10.);
  cuts->SetRejectPileup(kTRUE);
  cuts->SetRejectDAQIncompleteEvent(kTRUE);
  return cuts;
}
//________________________________________________
