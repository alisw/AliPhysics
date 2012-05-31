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
DisplayRecs()
{
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  gSystem->Load("libFMDutil.so");
  AliFMDDisplay* d = new AliFMDDisplay;
  d->AddLoad(AliFMDInput::kRecPoints);
  // d->AddLoad(AliFMDInput::kDigits);
  // d->AddLoad(AliFMDInput::kKinematics);
  d->Run();
}

//____________________________________________________________________
//
// EOF
//
