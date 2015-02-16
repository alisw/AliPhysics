void AliT0getTimeDelay(Int_t run)
{
  // Read Calib/TimeDelay entry with channel equalizing and mean CFD position 

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  man->SetRun(run);
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("T0/Calib/TimeDelay");
  AliT0CalibTimeEq *clb = (AliT0CalibTimeEq*)entry->GetObject();
  clb->Print();

}
