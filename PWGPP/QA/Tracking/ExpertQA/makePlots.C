makePlots(const char *inputFile)
{
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx++");
  gROOT->LoadMacro("AliHighPtTreeAnalysis.cxx++");

  AliHighPtTreeAnalysis *a = new AliHighPtTreeAnalysis( inputFile );
  a->Loop();
}
