//____________________________________________________________________
//
// $Id$
//
// Read in digits, and convert them to raw data files.  This is mainly
// for testing. 
//
/** Convert digits to Raw data
    @ingroup FMD_simple_script
*/
void
Convert2Raw()
{
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetRun(0);
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliRunLoader* runLoader = 
    AliRunLoader::Open("galice.root", "Alice", "read");
  runLoader->LoadgAlice();
  AliRun* run = runLoader->GetAliRun();
  AliLoader* fmdLoader = runLoader->GetLoader("FMDLoader");
  AliFMD* fmd = static_cast<AliFMD*>(run->GetDetector("FMD"));
  AliLog::SetModuleDebugLevel("FMD", 5);

  AliFMDParameters::Instance()->Init();
  AliFMDRawWriter rw(fmd);
  rw.Exec();
}
//____________________________________________________________________
// 
// EOF
//
