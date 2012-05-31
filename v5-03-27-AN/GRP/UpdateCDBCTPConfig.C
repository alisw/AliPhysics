void UpdateCDBCTPConfig(const char *CTPcfg, const char* cdbUri, const char* cfgFile){
  // Produces a trigger configuration starting from a 'cfg' file in the
  // GRP/CTP folder (CTPcfg argument). Stores the configuration in the specified CDB.
  // The third argument allows to specify the config file against which
  // the trigger confiuration is checked.

  AliCDBManager* cdb = AliCDBManager::Instance();
  // we set the default storage to the repository because some dets require
  // already at the time of geometry creation to find calibration objects in the cdb
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBStorage* storage = cdb->GetStorage(cdbUri);
  cdb->SetRun(0);
  AliCDBId id("GRP/CTP/Config",0,AliCDBRunRange::Infinity());
  AliCDBMetaData *md= new AliCDBMetaData();

  // Get root and AliRoot versions
  const char* rootv = gROOT->GetVersion();
  gROOT->ProcessLine(".L $ALICE_ROOT/macros/GetARversion.C");
  TString av(GetARversion());

  md->SetAliRootVersion(av.Data());
  md->SetComment(Form("Default CTP configuration for p-p mode produced with root version %s and AliRoot version %s",rootv,av.Data()));
  
  // construct the CTP configuration starting from GRP/CTP/<CTPcfg>.cfg file
  AliTriggerConfiguration *trconfig = AliTriggerConfiguration::LoadConfiguration(CTPcfg);
  if (!trconfig) {
    Printf("Invalid cfg file! Exiting...");
    return;
  }
  if (!trconfig->CheckConfiguration(cfgFile)) {
    Printf("CTP configuration is incompatible with the specified Config.C and AliRoot version! Exiting...");
    return;
  }
  
  Printf("Storing in CDB geometry produced with root version %s and AliRoot version %s",rootv,av.Data());
  storage->Put(trconfig,id,md);
  // This is to allow macros lauched after this one in the same session to find the
  // newly produced geometry.
  storage->QueryCDB(cdb->GetRun());

}


