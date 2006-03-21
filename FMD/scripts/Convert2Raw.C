//____________________________________________________________________
//
// $Id$
//
// Read in digits, and convert them to raw data files.  This is mainly
// for testing. 
//
void
Convert2Raw()
{
  AliRunLoader* runLoader = 
    AliRunLoader::Open("galice.root", "Alice", "read");
  runLoader->LoadgAlice();
  AliRun* run = runLoader->GetAliRun();
  AliLoader* fmdLoader = runLoader->GetLoader("FMDLoader");
  AliFMD* fmd = static_cast<AliFMD*>(run->GetDetector("FMD"));
  AliLog::SetModuleDebugLevel("FMD", 1);

  AliFMDParameters::Instance()->Init();
  AliFMDRawWriter rw(fmd);
  rw.Exec();
}
//____________________________________________________________________
// 
// EOF
//
