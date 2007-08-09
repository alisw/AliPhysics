void MakeFMDFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for FMD
  //
  if(!AliGeomManager::GetGeometry()){
    if(!(AliCDBManager::Instance())->IsDefaultStorageSet())
      AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
      AliCDBManager::Instance()->SetRun(0);
    AliGeomManager::LoadGeometry();
  }
  // needed for the constructors with local coordinates not to fail
  
  gSystem->Load("libFMDutil.so");
  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
    // save on file
    AliFMDAlignFaker faker(AliFMDAlignFaker::kAll, "geometry.root","FMDfullMisalignment.root");
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("STORAGE");
    AliFMDAlignFaker faker(AliFMDAlignFaker::kAll, "geometry.root", Storage);
  }

  // fRunMax should be changed in the constructor

  faker.SetSensorDisplacement(-0.005, -0.005, -0.005, 0.005, 0.005, 0.005);
  faker.SetSensorRotation(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  faker.SetHalfDisplacement(-0.25, -0.25, -0.25, 0.25, 0.25, 0.25);
  faker.SetHalfRotation(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  faker.Exec();

}
