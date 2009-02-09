//____________________________________________________________________
//
// $Id: PatternDigits.C 28055 2008-08-18 00:33:20Z cholm $
//
// Draw hits in the specialised FMD event display 
//
/** Display hits 
    @ingroup FMD_script
 */
void
PatternRaw(const char* file="raw.root")
{
  // AliLog::SetModuleDebugLevel("FMD", 8);
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  AliFMDParameters::Instance()->Init();
  gSystem->Load("libFMDutil.so");
  AliFMDPattern* d = new AliFMDPattern;
  d->AddLoad(AliFMDInput::kRaw);
  d->SetRawFile(file);
  d->SetName("raw");
  d->SetTitle("Raw");
  // d->AddLoad(AliFMDInput::kKinematics);
  d->Run();
}

//____________________________________________________________________
//
// EOF
//
