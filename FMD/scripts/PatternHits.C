//____________________________________________________________________
//
// $Id$
//
// Draw hits in the specialised FMD event display 
//
/** Pattern hits 
    @ingroup FMD_script
 */
void
PatternHits()
{
  // AliCDBManager* cdb = AliCDBManager::Instance();
  // cdb->SetDefaultStorage("local://cdb");
  gSystem->Load("libFMDutil.so");
  AliFMDPattern* d = new AliFMDPattern;
  d->AddLoad(AliFMDInput::kHits);
  // d->AddLoad(AliFMDInput::kKinematics);
  d->Run();
}

//____________________________________________________________________
//
// EOF
//
