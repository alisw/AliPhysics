//
// Configuration: (see loaded macro for details)
//
// - PID: realistic (full)
// - ITS: included
//
Bool_t RsnConfigSA(const char *taskName, const char *options, const char *path)
{
  gROOT->LoadMacro(Form("%s/RsnConfig.C", path));
  return RsnConfig(taskName, options, "pid+its", path);
}
