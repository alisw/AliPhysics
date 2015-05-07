void setlookUpTable()
{
  // Writing Lookup table into the Calibration DB
  // Arguments:

  TString DBFolder;
  Int_t firstRun   =  0;
  Int_t lastRun    = AliCDBRunRange::Infinity();
  Int_t beamPeriod =  1;
  char* objFormat = "T0 Lookup Table";

  AliT0CalibData *calibda=new AliT0CalibData("T0");

//  calibda->ReadAsciiLookup("lookUpTable.txt");
  calibda->ReadAsciiLookup("lookup_run2.txt");
  cout<<" @@@@ TRM "<<calibda->GetNumberOfTRMs()<<endl;
  //Store calibration data into database
  // AliCDBManager::Instance()->SetDefaultStorage("local:///scratch/alla/alice/Aug10/TestCDB");
   AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
   
  //  AliCDBManager::Instance()->SetSpecificStorage("T0",DBFolder.Data());

  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Alla");
  TString fPath="T0/Calib/LookUp_Table";


  // AliCDBStorage* storage = AliCDBManager::Instance()->GetSpecificStorage("T0");
  AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
  if(storage) {
    AliCDBId id(fPath.Data(),firstRun,lastRun);
    storage->Put(calibda, id, &md);
  }
}
