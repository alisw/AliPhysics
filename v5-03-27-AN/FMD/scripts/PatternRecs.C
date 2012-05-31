//____________________________________________________________________
//
// $Id$
//
// Draw hits in the specialised FMD event pattern 
//
/** Pattern hits 
    @ingroup FMD_script
 */
void
PatternRecs()
{
  AliLog::SetModuleDebugLevel("FMD", 1);
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  AliFMDParameters::Instance()->Init();
  gSystem->Load("libFMDutil.so");
  AliFMDPattern* d = new AliFMDPattern;
  d->AddLoad(AliFMDInput::kRecPoints);
  // d->AddLoad(AliFMDInput::kDigits);
  // d->AddLoad(AliFMDInput::kKinematics);
  d->Run();
}

//____________________________________________________________________
//
// EOF
//
