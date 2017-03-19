void
LoadOne(const char* simdir, const char* script, const char* opt="+g")
{
  gROOT->LoadMacro(Form("%s/%s.C+%s", simdir, script, opt));
}

void
FastLoad(const char* opt="+g")
{
  TString simdir = "$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/sim";
  if (gSystem->Getenv("ANA_SRC"))  
    simdir.Form("%s/sim", gSystem->Getenv("ANA_SRC"));
  gSystem->AddIncludePath(Form("-I%s", simdir.Data()));
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
    LoadOne(simdir, *ptr, opt);
    ptr++;
  }
}
