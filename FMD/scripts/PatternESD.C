//____________________________________________________________________
//
// $Id$
//
// Draw hits in the specialised FMD event display 
//
/** Display hits 
    @ingroup FMD_script
 */
void
PatternESD()
{
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT");
  gSystem->Load("libFMDutil.so");
  AliFMDPattern* d = new AliFMDPattern;
  d->SetMultiplicityCut(0);
  d->AddLoad(AliFMDInput::kESD);
  d->Run();
}

//____________________________________________________________________
//
// EOF
//
