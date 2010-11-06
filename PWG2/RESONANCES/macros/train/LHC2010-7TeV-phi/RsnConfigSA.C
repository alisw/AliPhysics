//
// Configuration: (see loaded macro for details)
//
// - PID: realistic (full)
// - ITS: included
// - dip: not included
//
Bool_t RsnConfigSA(const char *taskName, const char *options)
{
  gROOT->LoadMacro("$(ALICE_ROOT)/PWG2/RESONANCES/macros/train/LHC2010-7TeV-phi/RsnConfig.C");
  return RsnConfig(taskName, options, "realistic+its");
}
