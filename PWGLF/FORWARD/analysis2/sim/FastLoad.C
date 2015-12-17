void
LoadOne(const char* simdir, const char* script, const char* opt="+g")
{
  gROOT->LoadMacro(Form("%s/%s.C+%s", simdir, script, opt));
}

void
FastLoad()
{
  const char* simdir = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/sim";
  gSystem->AddIncludePath(Form("-I%s", simdir));
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");

  const char*  scripts[] = { "FastShortHeader",
			     "FastMonitor",
			     "FastCentEstimators",
			     "FastSim",
			     "FastAnalysis",
			     "FastCentHelper",
			     "dNdetaAnalysis",
			     "dNdyAnalysis",
			     "spectraAnalysis",
			     0 };
  const char** ptr = scripts;
  while (*ptr) {
    LoadOne(simdir, *ptr);
    ptr++;
  }
}
