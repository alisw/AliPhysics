void MakeADCalibEntry(Int_t run=1,const char *inputCDB = "raw://"){

  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage(inputCDB);
  man->SetRun(run);

  AliADCalibData *calibda = new AliADCalibData();

  // Creation of the object AD Calibration as a MetaData
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Michal Broz");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("AD Calibration test");
  AliCDBId id("AD/Calib/Data",0,AliCDBRunRange::Infinity());

  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBStorage *storLoc = man->GetDefaultStorage();
  storLoc->Put(calibda, id, md);

}
