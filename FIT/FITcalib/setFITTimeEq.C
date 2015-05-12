void setFITTimeEq(Int_t firstRun, Int_t lastRun)
{
  // Writing calibration coefficients into the Calibration DB
  // Arguments:
  Int_t beamPeriod =  1;
  char*   objFormat = "FIT initial time delay";
  
  AliT0CalibTimeEq *clb = new AliT0CalibTimeEq();
  for(Int_t ipmt=0; ipmt<200; ipmt++)  {
    clb->SetTimeEq(ipmt, 0);
    clb->SetCFDvalue(ipmt,0);
  }
  clb->Print();
   //Store calibration data into database
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Alla");
  TString fPath="FIT/Calib/TimeDelay";
  cout<<fPath.Data()<<endl;
 AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
  if(storage) {
     AliCDBId id(fPath.Data(),firstRun, lastRun );
    storage->Put(clb, id, &md);
  }
}
