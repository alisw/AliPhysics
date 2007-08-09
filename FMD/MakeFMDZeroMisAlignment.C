void MakeFMDZeroMisAlignment()
{
  // Create TClonesArray of zero-misalignment objects for FMD
  //
  if(!AliGeomManager::GetGeometry()){
    if(!(AliCDBManager::Instance())->IsDefaultStorageSet())
      AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
      AliCDBManager::Instance()->SetRun(0);
    AliGeomManager::LoadGeometry();
  }
  
  gSystem->Load("libFMDutil.so");
  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
    // save on file
    AliFMDAlignFaker faker(AliFMDAlignFaker::kAll, "geometry.root","FMDAlignObjs.root");
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("STORAGE");
    AliFMDAlignFaker faker(AliFMDAlignFaker::kAll, "geometry.root", Storage);
  }

  // fRunMax should be changed in the constructor

  faker.SetSensorDisplacement(0., 0., 0., 0., 0., 0.);
  faker.SetSensorRotation(0., 0., 0., 0., 0., 0.);
  faker.SetHalfDisplacement(0., 0., 0., 0., 0., 0.);
  faker.SetHalfRotation(0., 0., 0., 0., 0., 0.);
  faker.Exec();


}
