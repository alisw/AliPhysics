void runAnalysis() {
  TStopwatch timer;
  timer.Start();
  
  gSystem->AddIncludePath("-I\"$ALICE_ROOT/include\"");
  printf("****** Connect to PROOF *******\n");
  TProof::Open("proof://<username>@lxb6046.cern.ch"); 

  // Enable the Analysis Package
  gProof->UploadPackage("ESD.par");
  gProof->EnablePackage("ESD");
  gProof->UploadPackage("AOD.par");
  gProof->EnablePackage("AOD");
  gProof->UploadPackage("ANALYSIS.par");
  gProof->EnablePackage("ANALYSIS");

  gProof->GetManager()->ShowROOTVersions();
  gProof->ShowEnabledPackages();

  // You should get this macro and the txt file from:
  // http://aliceinfo.cern.ch/Offline/Analysis/CAF/
  gROOT->LoadMacro("CreateESDChain.C");
  TChain* chain = 0x0;
  chain = CreateESDChain("ESD1.txt",100);

  gROOT->LoadMacro("AliAnalysisTaskPt.cxx+");
  gProof->Load("AliAnalysisTaskPt.cxx+");
  gROOT->LoadMacro("demoCAF.C");
  demoCAF(chain,"proof");
 
  gSystem->Exec("rm -rf ESD ANALYSIS");

  timer.Stop();
  timer.Print();
}
