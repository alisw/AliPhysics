void MakeFMDFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for FMD
  //
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail
  
  gSystem->Load("libFMDutil.so");
  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    AliFMDAlignFaker faker(AliFMDAlignFaker::kAll, "geometry.root","FMDfullMisalignment.root");
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliFMDAlignFaker faker(AliFMDAlignFaker::kAll, "geometry.root", Storage);
  }

  // fRunMax should be changed in the constructor

  faker.SetSensorDisplacement(-0.005, -0.005, -0.005, 0.005, 0.005, 0.005);
  faker.SetSensorRotation(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  faker.SetHalfDisplacement(-0.25, -0.25, -0.25, 0.25, 0.25, 0.25);
  faker.SetHalfRotation(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  faker.Exec();

}
