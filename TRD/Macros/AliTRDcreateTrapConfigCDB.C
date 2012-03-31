void AliTRDcreateTrapConfigCDB(const TString &dirname = ".") {

  TString initName("initialize.r3610");
  TRegexp initCheck("initialize.r[0-9]*$");
  TRegexp cfgCheck("cf_.*.r[0-9]*$");

  AliTRDCalTrapConfig *caltrap = new AliTRDCalTrapConfig();

  TList cfgList;

  void *dirhandle = gSystem->OpenDirectory(dirname);

  const char* filename;
  while (filename = gSystem->GetDirEntry(dirhandle)) {
    TString file(filename);
    if (file.Contains(initCheck)) {
      initName = file;
      continue;
    }

    if (file.Contains(cfgCheck)) {
      cfgList.Add(new TObjString(file));
    }
  }

  TIter cfgIter(&cfgList);

  TObjString *cfgName = 0x0;
  while (cfgName = (TObjString*) cfgIter()) {
    ::Info("createTrapConfigCDB", Form("adding config: %s", cfgName->GetString().Data()));
    AliTRDtrapConfig *cfg = new AliTRDtrapConfig(cfgName->GetString(), cfgName->GetString());
    AliTRDtrapConfigHandler cfgHandler(cfg);
    cfgHandler.ResetMCMs();
    cfgHandler.Init();
    cfgHandler.LoadConfig(dirname + "/" + initName);
    cfgHandler.LoadConfig(dirname + "/" + cfgName->GetString());
    caltrap->Add(cfg);
  }

  AliCDBMetaData *cdbMeta = new AliCDBMetaData("Jochen Klein <jochen.klein@cern.ch>");
  AliCDBId cdbId("TRD/Calib/TrapConfig", 0, 999999999);;

  AliCDBStorage *storage = AliCDBManager::Instance()->GetStorage("local://$ALICE_ROOT/OCDB");
  storage->Put(caltrap, cdbId, cdbMeta);
}
