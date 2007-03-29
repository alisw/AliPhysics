void runBatchProcess() {
  TStopwatch timer;
  timer.Start();

  printf("*** Connect to AliEn ***\n");
  TGrid::Connect("alien://");
  gSystem->Load("libProofPlayer.so");
  gSystem->AddIncludePath("-I\"$ALICE_ROOT/include\"");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISRL.so");

  gROOT->LoadMacro("AliAnalysisTaskRLPt.cxx+");
  gROOT->LoadMacro("demoBatch.C");
  demoBatch();

  //gROOT->LoadMacro("CreateXML.C");
  //CreateXML();

  timer.Stop();
  timer.Print();
}

