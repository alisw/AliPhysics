{
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gROOT->ProcessLine(".L TInfo.cxx+");
  gROOT->ProcessLine(".L LInfo.cxx+");
  gROOT->ProcessLine(".L readOCDB_Temperature.C+");
  gROOT->ProcessLine(".L readOCDB_LED.C+");
}
