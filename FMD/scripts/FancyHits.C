//____________________________________________________________________
//
// $Id$
//
// Draw hits in the specialised FMD event fancy 
//
/** Fancy hits 
    @ingroup FMD_script
 */
void
FancyHits()
{
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT");
  gSystem->Load("libFMDutil.so");
  AliFMDFancy* d = new AliFMDFancy;
  d->AddLoad(AliFMDInput::kHits);
  d->AddLoad(AliFMDInput::kKinematics);
  d->Run();
}

//____________________________________________________________________
//
// EOF
//
