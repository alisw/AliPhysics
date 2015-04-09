void AliT0setLatency (Int_t firstRun, Int_t lastRun)
{

 // Arguments:
   //Store calibration data into database
  
  AliT0CalibLatency *calibda=new AliT0CalibLatency("T0");
  calibda->SetLatencyHPTDC(9000);

  //125097 
  calibda-> SetLatencyL1 (8.91406e+03)  ;
  calibda-> SetLatencyL1A( 8.91401e+03);
  calibda-> SetLatencyL1C (8.91412e+03) ;
  
  AliCDBMetaData md;
  md.SetResponsible("Alla");
  TString fPath="T0/Calib/Latency";
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
  if(storage) {
    AliCDBId id(fPath.Data(),firstRun, lastRun );
    storage->Put(calibda, id, &md);
    calibda->Print();
  }
  
}

