void run()
{

   gSystem->Load("libTree.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libPhysics.so");
   gSystem->Load("libSTEERBase.so");
   gSystem->Load("libESD.so");
   gSystem->Load("libAOD.so");
   gSystem->Load("libANALYSIS.so");
   gSystem->Load("libANALYSISalice.so");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");
   gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
   gSystem->AddIncludePath("-I$ALICE_ROOT/include");
   gStyle->SetPalette(1);
   gStyle->SetFillColor(kWhite);

   gROOT->LoadMacro("AliSpectraAODTrackCuts.cxx+g");
   gROOT->LoadMacro("AliSpectraAODEventCuts.cxx+g");
   gROOT->LoadMacro("AliSpectraAODHistoManager.cxx+g");
   gROOT->LoadMacro("AliAnalysisTaskSpectraAOD.cxx+g");








}
