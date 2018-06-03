{
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gROOT->ProcessLine(".L TInfo.cxx+");
  gROOT->ProcessLine(".L LInfo.cxx+");
}
