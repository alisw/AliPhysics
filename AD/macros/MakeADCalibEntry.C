void MakeADCalibEntry(const char *outputCDB = "local://$ALICE_ROOT/../AliRoot/OCDB"){

  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage(outputCDB);

  AliADCalibData *calibda = new AliADCalibData();

  // Creation of the object AD Calibration as a MetaData
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Michal Broz");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("AD Calibration test");
  AliCDBId id("AD/Calib/Data", 0, AliCDBRunRange::Infinity());

  man->Put(calibda, id, md);

}
