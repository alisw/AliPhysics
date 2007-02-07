void MakeFMDZeroMisAlignment()
{
  // Create TClonesArray of zero-misalignment objects for FMD
  //
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  
  gSystem->Load("libFMDutil.so");
  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    AliFMDAlignFaker faker(AliFMDAlignFaker::kAll, "geometry.root","FMDAlignObjs.root");
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliFMDAlignFaker faker(AliFMDAlignFaker::kAll, "geometry.root", Storage);
  }

  // fRunMax should be changed in the constructor

  faker.SetSensorDisplacement(0., 0., 0., 0., 0., 0.);
  faker.SetSensorRotation(0., 0., 0., 0., 0., 0.);
  faker.SetHalfDisplacement(0., 0., 0., 0., 0., 0.);
  faker.SetHalfRotation(0., 0., 0., 0., 0., 0.);
  faker.Exec();


}
