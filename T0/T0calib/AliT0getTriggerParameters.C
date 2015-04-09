void AliT0getTriggerParameters(Int_t run)
{ 
 // Arguments:
  AliCDBManager* man = AliCDBManager::Instance();
   man->SetDefaultStorage("raw://");
  //  man->SetDefaultStorage("local:///home/alla/alice/Jul14/OCDB/");
  man->SetRun(run);
  AliCDBEntry *entry = man->Get("T0/Calib/TriggerParam");

  AliT0TriggerParameters *calibda=(AliT0TriggerParameters*)entry->GetObject();
  calibda->Print();

}
