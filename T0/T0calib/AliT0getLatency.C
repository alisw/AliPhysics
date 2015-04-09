void AliT0getLatency(Int_t run)
{ 
 // Arguments:
  AliCDBManager* man = AliCDBManager::Instance();
   man->SetDefaultStorage("raw://");
  //  man->SetDefaultStorage("local:///home/alla/alice/Jul14/OCDB/");
  man->SetRun(run);
  AliCDBEntry *entry = man->Get("T0/Calib/Latency");

  AliT0CalibLatency *calibda=(AliT0CalibLatency*)entry->GetObject();
  calibda->Print();

}
