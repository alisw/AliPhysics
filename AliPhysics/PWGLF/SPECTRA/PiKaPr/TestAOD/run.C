void run()
{


  LoadLibs();

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gStyle->SetPalette(1);
  gStyle->SetFillColor(kWhite);
  
  gROOT->LoadMacro("AliSpectraAODTrackCuts.cxx+g");
  gROOT->LoadMacro("AliSpectraAODEventCuts.cxx+g");
  gROOT->LoadMacro("AliSpectraAODHistoManager.cxx+g");
  gROOT->LoadMacro("AliSpectraAODPID.cxx+g");
  gROOT->LoadMacro("AliAnalysisTaskSpectraAOD.cxx+g");

}
