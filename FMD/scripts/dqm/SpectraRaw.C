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
SpectraRaw(const char* file="raw.root", Int_t runno=0)
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  // gSystem->Load("libFMDanalysis.so");
  gSystem->Load("$(HOME)/scripts/foo.so");
  gSystem->Load("libFMDutil.so");

  AliCDBManager* cdb = AliCDBManager::Instance();
  const char* cdbUri = gSystem->Getenv("AMORE_CDB_URI");
  cdb->SetDefaultStorage(cdbUri);
  cdb->SetRun(runno);

  AliFMDParameters::Instance()->Init(kFALSE,
				     AliFMDParameters::kPulseGain|
				     AliFMDParameters::kPedestal|
				     AliFMDParameters::kDeadMap|
				     AliFMDParameters::kZeroSuppression|
				     AliFMDParameters::kAltroMap|
				     AliFMDParameters::kStripRange);
  AliFMDParameters::Instance()->SetSampleRate(2);
  
  AliFMDSpectraDisplay* d = new AliFMDSpectraDisplay;
  d->AddLoad(AliFMDInput::kRaw);
  d->SetRawFile(file);
  d->SetName("raw");
  d->SetTitle("Raw");
  d->Run();
}

//____________________________________________________________________
//
// EOF
//
