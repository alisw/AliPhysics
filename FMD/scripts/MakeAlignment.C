//____________________________________________________________________
//
// $Id$
//
// Make fake alignment data.
//
void
MakeAlignment()
{
  if (!TGeoManager::Import("geometry.root")) 
    gAlice->Init("$ALICE_ROOT/FMD/Config.C");
  AliCDBManager* cdb   = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://cdb");

  gSystem->Load("libFMDutil.so");
  AliFMDAlignFaker f;
  f.RemoveAlign(AliFMDAlignFaker::kHalves);
  f.SetSensorDisplacement(0, 0, 0, 0, 0, 0);
  f.SetSensorRotation(0, 0, 0, 3, 3, 3);
  f.Exec();
}
//____________________________________________________________________
//
// EOF
//
