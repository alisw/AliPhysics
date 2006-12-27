//____________________________________________________________________
//
// $Id$
//
// Make fake alignment data.
//
/** Make fake alignment data 
    @ingroup simple_script
 */
void
MakeAlignment()
{
  if (!TGeoManager::Import("geometry.root")) 
    gAlice->Init("$ALICE_ROOT/FMD/Config.C");
  AliCDBManager* cdb   = AliCDBManager::Instance();
  cdb->SetDefaultStorage("$ALICE_ROOT");
  AliLog::SetModuleDebugLevel("FMD", 1);
  gSystem->Load("libFMDutil.so");
  AliFMDAlignFaker f(AliFMDAlignFaker::kAll, "geometry.root", 0);
  // f.RemoveAlign(AliFMDAlignFaker::kHalves);
  f.SetSensorDisplacement(0, 0, 0, 0, 0, 0);
  f.SetSensorRotation(0, 0, 0, 0, 0, 0);
  f.SetHalfDisplacement(0, 0, 0, 0, 0, 0);
  f.SetHalfRotation(0, 0, 0, 0, 0, 0);
  f.Exec();
}
//____________________________________________________________________
//
// EOF
//
