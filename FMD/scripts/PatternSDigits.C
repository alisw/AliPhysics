//____________________________________________________________________
//
// $Id: PatternDigits.C 16307 2006-12-27 13:37:57Z cholm $
//
// Draw hits in the specialised FMD event display 
//
/** Display hits 
    @ingroup FMD_script
 */
void
PatternSDigits()
{
  AliLog::SetModuleDebugLevel("FMD", 1);
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  AliFMDParameters::Instance()->Init();
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libFMDanalysis.so");
  gSystem->Load("libFMDutil.so");
  AliFMDPattern* d = new AliFMDPattern;
  d->SetName("sdigit");
  d->SetTitle("Summable digits");
  d->AddLoad(AliFMDInput::kSDigits);
  d->Run();
}

//____________________________________________________________________
//
// EOF
//
