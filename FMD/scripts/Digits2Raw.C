void
Digits2Raw()
{
  // AliLog::SetModuleDebugLevel("FMD", 10);
  AliLog::SetModuleDebugLevel("RAW", 1);
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root", "Alice", "read");
  if (!runLoader) {
    AliError("Coulnd't read the file galice.root");
    return;
  }
  
  if  (runLoader->LoadgAlice()) return;
  AliRun* run = runLoader->GetAliRun();
  
  // Get the FMD 
  AliFMD* fmd = static_cast<AliFMD*>(run->GetDetector("FMD"));
  if (!fmd) {
    AliError("Failed to get detector FMD from loader");
    return;
  }
  
  // Get the FMD loader
  AliLoader* loader = runLoader->GetLoader("FMDLoader");
  if (!loader) {
    AliError("Failed to get detector FMD loader from loader");
    return;
  }
  if (runLoader->LoadHeader()) { 
    AliError("Failed to get event header information from loader");
    return;
  }
  TTree* treeE = runLoader->TreeE();
  
  AliCDBManager::Instance()->SetRun(0);
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  AliFMDParameters::Instance()->Init(kFALSE,
				     (AliFMDParameters::kPulseGain|
				     (AliFMDParameters::kPedestal|
				      AliFMDParameters::kDeadMap|
				      AliFMDParameters::kSampleRate|
				      AliFMDParameters::kAltroMap|
				      AliFMDParameters::kStripRange)));
  AliFMDParameters::Instance()->SetZeroSuppression(1);
  // AliFMDParameters::Instance()->SetPedestal(100);
  // AliFMDParameters::Instance()->SetPedestalWidth(0);
  // AliFMDParameters::Instance()->SetZSPreSamples(0);
  // AliFMDParameters::Instance()->SetZSPostSamples(0);
  
  fmd->Digits2Raw();
}
