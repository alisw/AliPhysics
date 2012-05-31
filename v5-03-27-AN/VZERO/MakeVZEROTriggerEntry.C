void MakeVZEROTriggerEntry(Int_t run,const char *inputCDB = "raw://"){

  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage(inputCDB);
  man->SetRun(run);

  AliCDBEntry *entry = man->Get("VZERO/Trigger/Data");
  AliVZEROTriggerData *trigdata = (AliVZEROTriggerData*)entry->GetObject();
  entry->SetObject(NULL);
  entry->SetOwner(kTRUE);

  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("VZERO Trigger Conditions Data from RAW OCDB");
  AliCDBId id("VZERO/Trigger/Data",0,AliCDBRunRange::Infinity());

  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBStorage *storLoc = man->GetDefaultStorage();
  storLoc->Put(trigdata, id, md);

}
