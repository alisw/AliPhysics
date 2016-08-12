//____________________________________________________________________
//
// $Id$
//
// Make fake alignment data.
//
/** Make fake alignment data 
    @ingroup FMD_simple_script
 */
void
MakeAlignment()
{
  if (!TGeoManager::Import("geometry.root")) 
    gAlice->Init("$ALICE_ROOT/FMD/Config.C");
  AliCDBManager* cdb   = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) 
    cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  AliLog::SetModuleDebugLevel("FMD", 1);
  gSystem->Load("libFMDutil");
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
