//____________________________________________________________________
//
// $Id$
//
// Make fake residual alignment data.
//
/** Make fake residual alignment data 
    @ingroup simple_script
 */
void
MakeResidualAlignment()
{
  if (!TGeoManager::Import("geometry.root")) 
    gAlice->Init("$ALICE_ROOT/FMD/Config.C");
  AliCDBManager* cdb   = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT");
  
  gSystem->Load("libFMDutil.so");
  AliFMDAlignFaker f(AliFMDAlignFaker::kAll, "geometry.root", "residual.root");
  f.SetComment("Residual alignment for PDC06");
  // f.RemoveAlign(AliFMDAlignFaker::kHalves);
  f.SetSensorDisplacement(-0.005, -0.005, -0.005, 0.005, 0.005, 0.005);
  f.SetSensorRotation(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  f.SetHalfDisplacement(-0.25, -0.25, -0.25, 0.25, 0.25, 0.25);
  f.SetHalfRotation(-0.5, -0.5, -0.5, 0.5, 0.5, 0.5);
  f.Exec();
}
//____________________________________________________________________
//
// EOF
//
