{
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gROOT->ProcessLine(".L TInfo.cxx+");
  gROOT->ProcessLine(".L LInfo.cxx+");
  gROOT->ProcessLine(".L readOCDB_Temperature.C+");
  gROOT->ProcessLine(".L readOCDB_LED.C+");
  gROOT->ProcessLine(".L readOCDB.C+");
  gROOT->ProcessLine(".L plotOCDB_Temperature.C+");
  gROOT->ProcessLine(".L plotOCDB_LED.C+");
  gROOT->ProcessLine(".L createTree.C+");
}
