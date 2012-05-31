void runAnalysis() {
  TStopwatch timer;
  timer.Start();

  gSystem->AddIncludePath("-I\"$ALICE_ROOT/include\"");
  gSystem->Load("libANALYSIS.so");

  gROOT->LoadMacro("testEvent.C+");
//  generate();
  filter_reco();

  timer.Stop();
  timer.Print();
}
