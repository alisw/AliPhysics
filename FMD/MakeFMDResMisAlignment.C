void MakeFMDResMisAlignment()
{
  // Create TClonesArray of residual misalignment objects for FMD
  //
  const char* macroname = "MakeFMDResMisAlignment.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);

  Bool_t    toCdb   = TString(gSystem->Getenv("TOCDB")) == TString("kTRUE");
  TString   storage = gSystem->Getenv("STORAGE");
  TString   output  = "FMDfullMisalignment.root";
  if(toCdb) output  = storage;
  
  gSystem->Load("libFMDutil.so");
  AliFMDAlignFaker::GetGeometry(toCdb, storage);
  AliFMDAlignFaker* faker = new AliFMDAlignFaker(AliFMDAlignFaker::kAll, 
						 "geometry.root", 
						 output.Data());
  

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
