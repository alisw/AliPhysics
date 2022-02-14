{
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->AddIncludePath("-INonSeparationAnalysis");

  gROOT->LoadMacro("AliVdMTree.cxx+");
  gROOT->LoadMacro("AliVdMScanData.cxx+");
  gROOT->LoadMacro("proc.C+");
  gROOT->LoadMacro("proc_pileup.C+");
  gROOT->LoadMacro("proc_bkgd.C+");
}
