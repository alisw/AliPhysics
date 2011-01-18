//
// Configuration: (see loaded macro for details)
//
// - PID: realistic (full)
// - ITS: not included
//
Bool_t RsnConfigNoSA(const char *taskName, const char *options, const char *path)
{
  gROOT->LoadMacro(Form("%s/RsnConfig.C", path));
  return RsnConfig(taskName, options, "pid", path);
}
