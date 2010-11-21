//
// Configuration: (see loaded macro for details)
//
// - PID: realistic (full)
// - ITS: not included
// - dip: included
//
Bool_t RsnConfigDipNoSA(const char *taskName, const char *options, const char *path)
{
  gROOT->LoadMacro("$(ALICE_ROOT)/PWG2/RESONANCES/macros/train/LHC2010-7TeV-phi/RsnConfig.C");
  return RsnConfig(taskName, options, "realistic+dip", path);
}
