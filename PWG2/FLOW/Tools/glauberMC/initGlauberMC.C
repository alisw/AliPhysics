{
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$PWD");
  gROOT->LoadMacro("AliGlauberNucleon.cxx++");
  gROOT->LoadMacro("AliGlauberNucleus.cxx++");
  gROOT->LoadMacro("AliGlauberMC.cxx++");
}
