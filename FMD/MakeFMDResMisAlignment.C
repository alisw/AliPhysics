void MakeFMDResMisAlignment()
{
  // Create TClonesArray of residual misalignment objects for FMD
  //
  const char* macroname = "MakeFMDResMisAlignment.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT");
  cdb->SetRun(0);

  Bool_t  toCdb   = TString(gSystem->Getenv("TOCDB")) == TString("kTRUE");
  TString storage = gSystem->Getenv("STORAGE");
  
  if(toCdb) {
    if(!storage.BeginsWith("local://") && 
       !storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE=\"%s\" is not valid. Exiting\n", storage.Data());
      return;
    }

    AliCDBStorage* store = cdb->GetStorage(storage.Data());
    if(!store){
      Error(macroname,"Unable to open storage %s\n", storage.Data());
      return;
    }

    AliCDBPath   path("GRP","Geometry","Data");
    AliCDBEntry* entry = storage->Get(path.GetPath(),cdb->GetRun());
    if(!entry) Fatal(macroname,"Could not get the specified CDB entry!");

    entry->SetOwner(0);
    TGeoManager* geom = static_cast<TGeoManager*>(entry->GetObject());
    AliGeomManager::SetGeometry(geom);
  }else
    //load geom from default CDB storage
    AliGeomManager::LoadGeometry(); 
  
  gSystem->Load("libFMDutil.so");
  AliFMDAlignFaker* faker = (toCdb ? 
			     // save on file
			     new AliFMDAlignFaker(AliFMDAlignFaker::kAll, 
						  "geometry.root",
						  "FMDfullMisalignment.root") :
			     // save in CDB storage
			     new AliFMDAlignFaker(AliFMDAlignFaker::kAll, 
						  "geometry.root", 
						  storage.Data()));

  faker->SetSensorDisplacement(-0.005, -0.005, -0.005, 0.005, 0.005, 0.005);
  faker->SetSensorRotation(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  faker->SetHalfDisplacement(-0.25, -0.25, -0.25, 0.25, 0.25, 0.25);
  faker->SetHalfRotation(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  faker->Exec();
  delete faker;

}
//
// EOF
//
