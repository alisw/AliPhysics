void AliT0setTriggerParameters (Int_t firstRun, Int_t lastRun)
{

 // Arguments:
   //Store calibration data into database
  
  AliT0TriggerParameters *calibda=new AliT0TriggerParameters();
  for (Int_t ipmt=0; ipmt<24; ipmt++) {
    calibda-> SetPMTstatus(ipmt, 1);
    calibda->SetThreshold(ipmt,1);
  }
   calibda->SetAmpCentr(30);
   calibda->SetAmpSemiCentr(5);
   calibda->SetTimeWindowLow(6090);
   calibda->SetTimeWindowHigh(6210);

 
  AliCDBMetaData md;
  md.SetResponsible("Alla");
  TString fPath="T0/Calib/TriggerParameters";
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
  if(storage) {
    AliCDBId id(fPath.Data(),firstRun, lastRun );
    storage->Put(calibda, id, &md);
    calibda->Print();
  }
  
}

