void runProcess() {
  TStopwatch timer;
  timer.Start();
  gSystem->AddIncludePath("-I\"$ALICE_ROOT/include\"");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISRL.so");

  gROOT->LoadMacro("AliAnalysisTaskRLPt.cxx+");
  gROOT->LoadMacro("demoLocal.C");
  demoLocal();

  timer.Stop();
  timer.Print();
}


