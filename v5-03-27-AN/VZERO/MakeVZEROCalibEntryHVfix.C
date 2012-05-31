void MakeVZEROCalibEntryHVfix(Int_t runIn, Int_t runOut, const char *inputCDB = "raw://"){

  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage(inputCDB);
  man->SetRun(runIn);

  AliCDBEntry *entry = man->Get("VZERO/Calib/Data");
  AliVZEROCalibData *calibdaorg = (AliVZEROCalibData*)entry->GetObject();
  AliVZEROCalibData *calibda = new AliVZEROCalibData(*calibdaorg);

  // Creation of the object VZERO Calibration as a MetaData
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("VZERO Calibration from RAW OCDB (repaired manually for HV ramping-down before EOR)");
  AliCDBId id("VZERO/Calib/Data",runOut,runOut);

  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBStorage *storLoc = man->GetDefaultStorage();
  storLoc->Put(calibda, id, md);

}
