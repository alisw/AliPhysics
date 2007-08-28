void MakeFMDResMisAlignment()
{
  // Create TClonesArray of residual misalignment objects for FMD
  //
  const char* macroname = "MakeFMDResMisAlignment.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT");
  cdb->SetRun(0);
  
  TString Storage;
  
  if( gSystem->Getenv("TOCDB") == TString("kTRUE") ){
    Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBPath path("GRP","Geometry","Data");
    AliCDBEntry *entry = storage->Get(path.GetPath(),cdb->GetRun());
    if(!entry) Fatal(macroname,"Could not get the specified CDB entry!");
    entry->SetOwner(0);
    TGeoManager* geom = (TGeoManager*) entry->GetObject();
    AliGeomManager::SetGeometry(geom);
  }else{
    AliGeomManager::LoadGeometry(); //load geom from default CDB storage
  }    
  
  gSystem->Load("libFMDutil.so");
  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
    // save on file
    AliFMDAlignFaker faker(AliFMDAlignFaker::kAll, "geometry.root","FMDresidualMisalignment.root");
  }else{
    // save in CDB storage
    AliFMDAlignFaker faker(AliFMDAlignFaker::kAll, "geometry.root", Storage.Data());
  }

  // fRunMax should be changed in the constructor

  faker.SetSensorDisplacement(-0.005, -0.005, -0.005, 0.005, 0.005, 0.005);
  faker.SetSensorRotation(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  faker.SetHalfDisplacement(-0.25, -0.25, -0.25, 0.25, 0.25, 0.25);
  faker.SetHalfRotation(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  faker.Exec();


}
