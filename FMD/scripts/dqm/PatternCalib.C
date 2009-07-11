//____________________________________________________________________
//
// $Id: PatternDigits.C 28055 2008-08-18 00:33:20Z cholm $
//
// Draw hits in the specialised FMD event display 
//
/** Display hits 
    @ingroup FMD_script
 */
Bool_t
CheckFile(const char* prefix, int number, TString& f)
{
  f = (Form("%s%d.csv", prefix, number));
  std::cout << "Checking if " << f << " exists ... " << std::flush;
  f = gSystem->Which("$(HOME)/calib/", f.Data());
  std::cout << '"' << f << '"' << std::endl;
  return !f.IsNull();
}
void
PatternCalib(const char* file="raw.root", Int_t runno=0)
{
  // AliLog::SetModuleDebugLevel("FMD", 1);
  gSystem->Load("libFMDutil.so");

  AliCDBManager* cdb = AliCDBManager::Instance();
  const char* cdbUri = gSystem->Getenv("AMORE_CDB_URI");
  cdb->SetDefaultStorage(cdbUri);
  cdb->SetRun(runno);

  AliFMDCalibStripRange* range = new AliFMDCalibStripRange;
  AliFMDCalibSampleRate* rate  = new AliFMDCalibSampleRate;
  AliFMDCalibPedestal*   peds  = new AliFMDCalibPedestal;
  AliFMDCalibGain*       gains = new AliFMDCalibGain;
  Bool_t gotConds = kFALSE;
  Bool_t gotPeds  = kFALSE;
  Bool_t gotGains = kFALSE;
  for (Int_t i = 1; i <= 3; i++) { 
    TString f;
    if (CheckFile("conditions", i, f)) {
      gotConds = kTRUE;
      std::cout << "Reading conditions for FMD" <<i<< " from " <<f<< std::endl;
      std::ifstream in(f.Data());
      range->ReadFromFile(in);
      rate->ReadFromFile(in);
    }
    if (CheckFile("peds", i, f)) {
      gotPeds = kTRUE;
      std::cout << "Reading pedestals for FMD" <<i<< " from " <<f<< std::endl;
      std::ifstream in(f.Data());
      peds->ReadFromFile(in);
    }
    if (CheckFile("gains", i, f)) {
      gotGains = kTRUE;
      std::cout << "Reading gains for FMD" <<i<< " from " <<f<< std::endl;
      std::ifstream in(f.Data());
      gains->ReadFromFile(in);
    }
  }

  Int_t mask = (AliFMDParameters::kDeadMap|
		AliFMDParameters::kZeroSuppression|
		AliFMDParameters::kAltroMap);


  if (!gotConds) mask |= AliFMDParameters::kStripRange;
  if (!gotConds) mask |= AliFMDParameters::kSampleRate;
  if (!gotPeds)  mask |= AliFMDParameters::kPedestal;
  if (!gotGains) mask |= AliFMDParameters::kPulseGain;

  AliFMDParameters* pars = AliFMDParameters::Instance();
  pars->Init(kFALSE, mask);

  if (gotConds) pars->SetStripRange(range);
  if (gotConds) pars->SetSampleRate(rate);
  if (gotPeds)  pars->SetPedestal(peds);
  if (gotGains) pars->SetGain(gains);
  
  // pars->Print("pedestal");

  AliFMDPattern* d = new AliFMDPattern;
  d->AddLoad(AliFMDInput::kRawCalib);
  d->SetRawFile(file);
  d->SetName("rawCalib");
  d->SetTitle("Calibrated Raw");
  d->Run();
}

//____________________________________________________________________
//
// EOF
//
