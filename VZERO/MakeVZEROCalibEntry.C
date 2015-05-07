void MakeVZEROCalibEntry(Int_t run, const char *inputCDB = "raw://"){

  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage(inputCDB);
  man->SetRun(run);

  AliCDBEntry *entry = man->Get("VZERO/Calib/Data");
  AliVZEROCalibData *calibdaorg = (AliVZEROCalibData*)entry->GetObject();
  AliVZEROCalibData *calibda = new AliVZEROCalibData(*calibdaorg);

  for (Int_t i = 0; i < 64; ++i) {
    calibda->SetTimeOffset(5.0,i i);
  }
  for (Int_t i = 0; i < 8; ++i) {
    calibda->SetWidthResolution(2, i);
  }

  // Creation of the object VZERO Calibration as a MetaData
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("VZERO Calibration from RAW OCDB");
  AliCDBId id("VZERO/Calib/Data", 0, AliCDBRunRange::Infinity());

  AliCDBStorage *storLoc = man->GetStorage("local://$ALICE_ROOT/OCDB");
  storLoc->Put(calibda, id, md);

}
