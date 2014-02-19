makePlots(const char *inputFile)
{
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  //gROOT->LoadMacro("$ALICE_ROOT/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx++");
  gROOT->LoadMacro("AliHighPtTreeAnalysis.cxx++");

  AliHighPtTreeAnalysis *a = new AliHighPtTreeAnalysis( inputFile );
  a->Loop();
}
