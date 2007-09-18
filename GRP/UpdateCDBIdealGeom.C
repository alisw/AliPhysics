void UpdateCDBIdealGeom(const char* cdbUri, const char* cfgFile){
  // Produce the ideal geometry and store it in the specified CDB
  // The second argument allows to specify the config file to be used
  // in particular for giving the choice to generate either a full or
  // a partial geometry.
  //

  AliCDBManager* cdb = AliCDBManager::Instance();
  // we set the default storage to the repository because some dets require
  // already at the time of geometry creation to find calibration objects in the cdb
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBStorage* storage = cdb->GetStorage(cdbUri);
  cdb->SetRun(0);
  AliCDBId id("GRP/Geometry/Data",0,AliCDBRunRange::Infinity());
  AliCDBMetaData *md= new AliCDBMetaData();

  // Get root and AliRoot versions
  const char* rootv = gROOT->GetVersion();
  gROOT->ProcessLine(".L $ALICE_ROOT/macros/GetARversion.C");
  TString av(GetARversion());

  md->SetAliRootVersion(av.Data());
  md->SetComment(Form("Geometry produced with root version %s and AliRoot version %s",rootv,av.Data()));
  
  gAlice->Init(cfgFile);
  
  if(!gGeoManager){
    Printf("Unable to produce a valid geometry to be put in the CDB!");
    return;
  }
  
  Printf("Storing in CDB geometry produced with root version %s and AliRoot version %s",rootv,av.Data());
  storage->Put(gGeoManager,id,md);
  // This is to allow macros lauched after this one in the same session to find the
  // newly produced geometry.
  storage->QueryCDB(cdb->GetRun());

}


