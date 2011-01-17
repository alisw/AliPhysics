//
// Configuration: (see loaded macro for details)
//
// - PID: realistic (full)
// - ITS: included
// - dip: included
//
Bool_t RsnConfigDipSA(const char *taskName, const char *options, const char *path, Int_t multMin = 0, Int_t multMax = 0)
{
  gROOT->LoadMacro("$(ALICE_ROOT)/PWG2/RESONANCES/macros/train/LHC2010-7TeV-phi/RsnConfig.C");
  return RsnConfig(taskName, options, "pid+its+dip", path, multMin, multMax);
}
